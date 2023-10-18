#include <array>
#include <cassert>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

#include <aio.h>
#include <fcntl.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef __BMI2__
#  include <x86intrin.h>
#endif

#ifndef _LARGEFILE64_SOURCE
#  define _LARGEFILE64_SOURCE
#endif
#define _FILE_OFFSET_BITS 64

#ifndef UNROLL
#  define UNROLL 1
#endif

namespace fastdna {

using namespace std::chrono;

using fsize_t = off_t;  // File sizes and offsets beyond 32 bits
using binary_buffer_t = std::vector<uint8_t>;
using word_t = uint64_t; // The largest size that can be bit-processed in parallel
using bit_mapping_t = std::array<uint8_t, 256>; // Mapping from DNA character to 2 bits
using char_mapping_t = std::array<uint32_t, 256>; // Mapping from 8 bits to 4 DNA characters


// Constants:

constexpr auto HEADER_SZ = sizeof(uint64_t); // How many bytes to store input size
constexpr auto BPW = sizeof(word_t); // bytes per word, for use with pext
constexpr auto BPPB = 4; // Base-pairs per byte (two bits each)

constexpr auto UNROLL_FACTOR = UNROLL;
constexpr fsize_t MIN_CHUNK_SZ = 1 << 10;
constexpr fsize_t DEFAULT_CHUNK_SZ = 1 << 17;
constexpr char char_of[4] = { 'A', 'C', 'T', 'G' };
constexpr fsize_t MAX_BUFFER_SZ = (1 << 30); // Buffers >= 2GB are problematic for I/O

constexpr char_mapping_t
letters_of_code()
{
  char_mapping_t ret;
  for (int i = 0; i < 256; i++) {
    const uint8_t bits = static_cast<uint8_t>(i);
    const uint32_t decoded =
      char_of[uint8_t(bits << 0) >> 6] << 0  |
      char_of[uint8_t(bits << 2) >> 6] << 8  |
      char_of[uint8_t(bits << 4) >> 6] << 16 |
      char_of[uint8_t(bits << 6) >> 6] << 24;
    ret[i] = decoded;
  }
  return ret;
}

constexpr bit_mapping_t
codes_of_letters()
{
  bit_mapping_t ret;
  ret['A'] = 0; // 0b00
  ret['C'] = 1; // 0b01
  ret['T'] = 2; // 0b10
  ret['G'] = 3; // 0b11
  return ret;
}

constexpr auto bit_map = codes_of_letters();
constexpr auto char_map = letters_of_code();

// Global configuration
struct Configuration {
  fsize_t chunk_sz = DEFAULT_CHUNK_SZ; // Size of I/O operations (in bytes)
  fsize_t task_sz = 0;     // Size of input handled by one thread (in bytes)
  unsigned nthreads = 1;   // No. of threads to execute with
  bool async = false;      // Should file I/O be asynchronous?
  bool skip_io = false;    // Should we skip all read/write operations?
  bool ordered = true;     // Should we encode bits in order of base pairs in a word?
  bool use_pext = false;   // Should we use bit-parallel operations?
  bool encode = true;      // Whether operation is encoding or decoding
} g_config;

/*****************************************************************************
    Utilities
******************************************************************************/

//////////////////////////////////////////////////////////////////////////////
// Round up or down to the nearest multiple:
constexpr auto round_up = [](auto num, auto multiple) {
  assert(multiple > 1);
  return multiple? multiple * ((num + multiple - 1) / multiple) : num;
};

constexpr auto round_down = [](auto num, auto multiple) {
  return (num / multiple) * multiple;
};

constexpr auto encoded_byte_count = [](auto unencoded_count) {
  return round_up(unencoded_count, BPPB) / BPPB;
};

//////////////////////////////////////////////////////////////////////////////
// Where to start writing the output for an input that starts at 'input_pos'
inline constexpr fsize_t
output_pos(fsize_t input_pos)
{
  return HEADER_SZ + input_pos / BPPB;
}

//////////////////////////////////////////////////////////////////////////////
// Figure out what is the amount of bytes each thread should work on.
// The returned task size must be a multiple of BPW (unless the file is shorter).
void
compute_thread_count(fsize_t filesize)
{
  // If there's too little work to be done, no point in multithreading:
  if (filesize < MIN_CHUNK_SZ) {
    g_config.nthreads = 1;
    g_config.task_sz = filesize;
    return;
  }

  const auto equal_work = filesize / g_config.nthreads;
  g_config.task_sz = round_down(equal_work, BPW * UNROLL_FACTOR);
  assert(g_config.task_sz > 0);

  // If chunk size is larger than a task, shrink it to avoid waste:
  if (g_config.chunk_sz > g_config.task_sz) {
    g_config.chunk_sz = std::max(round_up(g_config.task_sz, 8), MIN_CHUNK_SZ);
  }
}


//////////////////////////////////////////////////////////////////////////////
// Open input and output files, return their descriptors and input size
// If encoding, also writes the input size to the first 8 bytes of the output.
std::tuple<int, int, fsize_t>
open_files(std::string fn)
{
  const fsize_t insize = std::filesystem::file_size(fn);
  int infile = open(fn.c_str(), O_RDONLY | O_LARGEFILE);
  assert(infile > 0);

  const auto outfn = fn + (g_config.encode? ".2b" : ".dec");
  // Initialize and resize output file, write header (input size)
  int outfile = open(outfn.c_str(), O_CREAT | O_RDWR | O_LARGEFILE, S_IRUSR | S_IWUSR);
  assert(outfile > 0);

  if (g_config.encode) {
    [[maybe_unused]] int err = ftruncate(outfile, encoded_byte_count(insize) + HEADER_SZ);
    assert(err == 0);

    [[maybe_unused]] auto written = write(outfile, &insize, HEADER_SZ);
    assert(written == HEADER_SZ);
  }

  return { infile, outfile, insize };
}


//////////////////////////////////////////////////////////////////////////////
// Initialize an aio struct
struct aiocb
create_aio(int fd, char* buf)
{
  struct aiocb ret;
  ret.aio_fildes = fd;
  ret.aio_buf = buf;
  ret.aio_reqprio = 0;
  ret.aio_sigevent.sigev_notify = SIGEV_NONE;
  return ret;
}

//////////////////////////////////////////////////////////////////////////////
// Wait for an aio operation to complete (blocking)
void
inline aio_wait(struct aiocb* io)
{
  const struct aiocb* aio_list[1] = { io };
  [[maybe_unused]] auto status = aio_suspend(aio_list, 1, NULL);
  assert(status == 0);
}

//////////////////////////////////////////////////////////////////////////////
// Given a current byte position in input, return either a chunk size, or if
// that exceed the last byte to be read, the difference between the two.
fsize_t
next_read_size(fsize_t cur_byte, fsize_t last_byte)
{
  return std::min(g_config.chunk_sz, last_byte - cur_byte);
}

//////////////////////////////////////////////////////////////////////////////
void
usage(char**argv)
{
  std::cerr << argv[0] <<" usage: [-t threads] [-c chunk_size] [-a] encode|decode filename\n";
  std::cerr << "\t-t: Maximum number of threads to use (may be lower)\n";
  std::cerr << "\t-c: I/O chunk size (0 for entire file, i.e., no chunks; default chunk size=" << DEFAULT_CHUNK_SZ << ")\n";
  std::cerr << "\t-a: Perform asynchronous I/O\n";
  std::cerr << "\t-i: Skip I/O operations, synchronous only (using uninitialized data)\n";
  std::cerr << "\t-p: Use bit-parallel operations instead of one byte at a time.\n";
  std::cerr << "\t-u: Encode bits unordered in word (faster)\n";
  std::cerr << "\te[ncode] or [d]ecode: operation to perform\n";
  std::cerr << "\tfilename of the file to encode/decode\n";
  std::exit(1);
}

//////////////////////////////////////////////////////////////////////////////
// Parse command-line options. Sets configuration options in g_config.
// Decides how many threads to invoke and how much work (bytes) each thread
// should handle (which is further divided by chunk_sz in each thread).
// Returns an open file handle for the input and output file, as well as input size.
std::tuple<int, int, fsize_t>
parse_cmdline(int argc, char** argv)
{
  if (argc < 3) usage(argv);

  int opt;
  while ((opt = getopt(argc, argv, "aiupt:c:")) != -1) {
    switch (opt) {
      case 'a':  g_config.async = true; break;
      case 'i':  g_config.skip_io = true; break;
      case 'u':  g_config.ordered = false; break;
      case 'p':  g_config.use_pext = true; break;
      case 't':  g_config.nthreads = atoi(optarg); break;
      case 'c':  g_config.chunk_sz = atoll(optarg); break;
      default:   usage(argv);
    }
  }

  if (optind != argc - 2) usage(argv);
  switch(argv[optind++][0]) {
    case 'e':  g_config.encode = true;  break;
    case 'd':  g_config.encode = false; break;
    default:   usage(argv);
  }

  assert(g_config.ordered || g_config.use_pext); // Must have pext for unordered
  assert(!g_config.async || !g_config.skip_io); // Can't have async I/O with skipping
  assert(!g_config.async || g_config.nthreads == 1); // Can't have async I/O with multithreading

  const auto fn = std::string(argv[optind]);
  const auto [infile, outfile, size] = open_files(fn);

  // Too-small chunk sizes are unsupported due to division errros
  if (g_config.chunk_sz < MIN_CHUNK_SZ) g_config.chunk_sz = MIN_CHUNK_SZ;

  compute_thread_count(size);
  assert(g_config.chunk_sz % 8 == 0); // Chunks have to be multiples of 8 for perfromance

  return { infile, outfile, size };
}

/*****************************************************************************
    Encoding algorithms
******************************************************************************/

//////////////////////////////////////////////////////////////////////////////
// Encode a range of 64-bit strings (words) into two bytes using PEXT
template <bool ORDERED>
inline char*
encode_pext(const word_t* words, fsize_t limit, char* out)
{
#ifndef __BMI2__
  assert(false && "BMI2 instruction set not supported!");
  return nullptr;
#else

  constexpr word_t MASK = 0b00000110'00000110'00000110'00000110'00000110'00000110'00000110'00000110;
  auto wide_out = reinterpret_cast<uint16_t*>(out);

#pragma GCC unroll UNROLL_FACTOR
  for (fsize_t i = 0; i < limit; i += BPW) {
    if constexpr (ORDERED) {
      const auto encoded = _pext_u64(__builtin_bswap64(*words++), MASK);
      *out++ = (encoded >> 8) & 0xFF;
      *out++ = encoded & 0xFF;
    } else {
      *wide_out++ = _pext_u64(*words++, MASK);
    }
  }

  return ORDERED? out : reinterpret_cast<char*>(wide_out);
#endif
}

//////////////////////////////////////////////////////////////////////////////
// Encode any number of bytes into two-bit representation in memory buffers.
// At first values are encoded in 8-bytes batches by encode_words.
// Then, bytes are encoded one-by-one into bytes.
// If the remaining bytes is not a multiple of BPPB, the last few values are
// condensed to the left of a byte.
char*
encode(const char* in, fsize_t insize, char* out)
{
  fsize_t current_base = 0;
  if (g_config.use_pext) {
    const auto words = reinterpret_cast<const word_t*>(in);
    current_base = round_down(insize, UNROLL_FACTOR * BPW);
    out = g_config.ordered ?
          encode_pext<true>(words, current_base, out) :
          encode_pext<false>(words, current_base, out);
  }

  while (insize - current_base >= BPPB) {
    *out = bit_map[in[current_base++]];
    *out = (*out << 2) | bit_map[in[current_base++]];
    *out = (*out << 2) | bit_map[in[current_base++]];
    *out = (*out << 2) | bit_map[in[current_base++]];
    out++;
  }

  assert(insize - current_base < BPPB);
  if (current_base == insize) {
    return out;
  }

  *out = 0;
  switch (insize - current_base) {
    case 3: *out |= bit_map[in[current_base++]] << 2; [[fallthrough]];
    case 2: *out |= bit_map[in[current_base++]] << 4; [[fallthrough]];
    case 1: *out |= bit_map[in[current_base++]] << 6; [[fallthrough]];
    default: break;
  }

  return out;
}


//////////////////////////////////////////////////////////////////////////////
// Read, encode, and write a subset of the input using synchronous I/O
void
encode_file_chunk_sync(int infile, int outfile, fsize_t from, fsize_t to)
{
  auto input = new char[g_config.chunk_sz];
  auto output = new char[encoded_byte_count(g_config.chunk_sz)];
  fsize_t total_io = 0;

  while (from < to) {
    const auto count = next_read_size(from, to);

    // Read next `count` bytes, possibly in blocks of size MAX_BUFFER_SZ
    if (!g_config.skip_io) {
      for (total_io = 0; total_io < count; ) {
        total_io += pread64(infile, input + total_io, std::min(count - total_io, MAX_BUFFER_SZ), from + total_io);
      }
    }

    encode(input, count, output);

    if (!g_config.skip_io) {
      for (total_io = 0; total_io < encoded_byte_count(count); ) {
        total_io += pwrite64(outfile, output + total_io, std::min(encoded_byte_count(count) - total_io, MAX_BUFFER_SZ), output_pos(from) + total_io);
      }
    }

    from += count;
  }

  delete[] output;
  delete[] input;
}

//////////////////////////////////////////////////////////////////////////////
// Read, encode, and write a subset of the input using asynchronous I/O
// This is a concurrent algorithm that works with two sets of input and output
// buffers, so that a previous write and next read can take place concurrently
// with the current encoding.
void
encode_file_chunk_async(int infile, int outfile, fsize_t from, fsize_t to)
{
  // Two sets of input and output data buffers, for simultaneous operations:
  struct aiocb read1 = create_aio(infile, new char[g_config.chunk_sz]);
  struct aiocb read2 = create_aio(infile, new char[g_config.chunk_sz]);
  struct aiocb write1 = create_aio(outfile, new char[encoded_byte_count(g_config.chunk_sz)]);
  struct aiocb write2 = create_aio(outfile, new char[encoded_byte_count(g_config.chunk_sz)]);
  struct aiocb *cur_read = &read1, *nxt_read = &read2, *cur_write = &write1, *nxt_write = &write2;

  // Can't handle too-large chunks in aio:
  if (g_config.chunk_sz > MAX_BUFFER_SZ) {
    g_config.chunk_sz = MAX_BUFFER_SZ;
  }

  // Before main loop, read in first chunk of data:
  cur_read->aio_offset = from;
  cur_read->aio_nbytes = next_read_size(from, to);
  assert(cur_read->aio_nbytes <= MAX_BUFFER_SZ);   // Overly-large buffers not supported in aio
  aio_read(cur_read);

  // main loop:
  while (from < to) {

    // Start reading the next input chunk:
    nxt_read->aio_offset = cur_read->aio_offset + cur_read->aio_nbytes;
    nxt_read->aio_nbytes = next_read_size(nxt_read->aio_offset, to);
    aio_read(nxt_read); // Could be zero bytes

    // Make sure previous write is completed before encoding into that buffer:
    aio_wait(cur_write);
    std::swap(nxt_write, cur_write);

    // Ensure current input has been read, then encode it and start writing it:
    aio_wait(cur_read);
    encode(static_cast<char*>(const_cast<void*>(cur_read->aio_buf)),
           cur_read->aio_nbytes,
           static_cast<char*>(const_cast<void*>(nxt_write->aio_buf)));
    nxt_write->aio_offset = output_pos(cur_read->aio_offset);
    nxt_write->aio_nbytes = encoded_byte_count(cur_read->aio_nbytes);
    aio_write(nxt_write);

    // Prepare for the next read at the next iteration:
    from += cur_read->aio_nbytes;
    std::swap(cur_read, nxt_read);
  }

  aio_wait(nxt_write);   // All done, wait for last write

  delete[] static_cast<char*>(const_cast<void*>(cur_read->aio_buf));
  delete[] static_cast<char*>(const_cast<void*>(nxt_read->aio_buf));
  delete[] static_cast<char*>(const_cast<void*>(cur_write->aio_buf));
  delete[] static_cast<char*>(const_cast<void*>(nxt_write->aio_buf));
}

//////////////////////////////////////////////////////////////////////////////
// Multithreaded implementation of encoding a file.
// Takes file handles for the input and output and the input file size.
void
encode_file(int infile, int outfile, fsize_t insize)
{
  std::vector<std::thread> threads;

  for (unsigned i = 0; i < g_config.nthreads; ++i) {
    threads.push_back(std::thread(
          g_config.async? encode_file_chunk_async : encode_file_chunk_sync,
          infile,
          outfile,
          i * g_config.task_sz,
          i == (g_config.nthreads - 1)? insize : (i + 1) * g_config.task_sz)
        );
  }

  for (auto& thr : threads) {
    thr.join();
  }

  close(infile);
  close(outfile);
}


/*****************************************************************************
    Decoding algorithms
******************************************************************************/

//////////////////////////////////////////////////////////////////////////////
// Decode two bytes of input at a time using PDEP and reordering
template <bool ORDERED> // Does encoding maintain original basepair order?
inline char*
decode_pdep(const char* in, fsize_t remaining_bases, char* out)
{
#ifndef __BMI2__
  assert(false && "BMI2 instruction set not supported!");
  return nullptr;
#else

  constexpr word_t OR =
      0b01000001'01000001'01000001'01000001'01000001'01000001'01000001'01000001;
  constexpr word_t MASK1 =
      0b00000110'00000110'00000110'00000110'00000110'00000110'00000110'00000110;
  constexpr word_t MASK2 =
      0b00000001'00000001'00000001'00000001'00000001'00000001'00000001'00000001;
  auto wide_in = reinterpret_cast<const uint16_t*>(in);
  auto wide_out = reinterpret_cast<uint64_t*>(out);

#pragma GCC unroll UNROLL_FACTOR
  for (fsize_t i = 0; i < remaining_bases; i += sizeof(*wide_in) * BPPB, ++wide_out) {
    if constexpr (ORDERED) {
      *wide_out = __builtin_bswap64(OR | _pdep_u64(__builtin_bswap16(*wide_in++), MASK1));
    } else {
      *wide_out = OR | _pdep_u64(*wide_in++, MASK1);
    }

    // Translate 'E's to 'T's using bit manipulation
    const auto is_e = (*wide_out >> 2) & ((~*wide_out) >> 1) & MASK2;
    *wide_out ^= is_e | (is_e << 4);
  }

  return reinterpret_cast<char*>(wide_out);
#endif
}

//////////////////////////////////////////////////////////////////////////////
// Decode num_bases bytes of an an input buffer to an output buffer
char*
decode(const char* in,  char* out, fsize_t remaining_bases)
{
  if (g_config.use_pext && remaining_bases >= fsize_t(BPW)) {
    const fsize_t bases = round_down(remaining_bases, UNROLL_FACTOR * BPW);
    out = g_config.ordered? decode_pdep<true>(in, bases, out) : decode_pdep<false>(in, bases, out);
    remaining_bases -= bases;
    in += bases / BPPB;
  }

  // write any whole bytes left in the input
  auto wide_out = reinterpret_cast<uint32_t*>(out);
  out += round_down(remaining_bases, BPPB);

  for (; remaining_bases >= BPPB; remaining_bases -= BPPB) {
    *wide_out++ = char_map[static_cast<uint8_t>(*in++)];
  }

  assert(remaining_bases < BPPB);
  // if there is an input byte left, it must not be full (<4 encoded chars in the byte)
  switch (remaining_bases)
  {
    case 3:  *out++ = char_of[(*in >> 2) & 0x3];  [[fallthrough]];
    case 2:  *out++ = char_of[(*in >> 4) & 0x3];  [[fallthrough]];
    case 1:  *out++ = char_of[(*in >> 6) & 0x3];  [[fallthrough]];
    default: break;
 }

  return out;
}


//////////////////////////////////////////////////////////////////////////////
// Read, decode, and write a subset of the input using synchronous I/O. Params:
// infile, outfile: file handles for input and output.
// from, to: indices of first byte to read from input, one after last byte.
// remaining_bases: how many base pairs (bytes) remain to be decoded.
void
decode_file_chunk_sync(int infile, int outfile, fsize_t from, fsize_t to, fsize_t remaining_bases)
{
  auto input = new char[g_config.chunk_sz];
  auto output = new char[g_config.chunk_sz * BPPB];
  fsize_t total_io = 0;

  while (from < to) {
    const auto nbytes_read = next_read_size(from, to);

    // Read next `nbytes_read` bytes, possibly in blocks of size MAX_BUFFER_SZ
    if (!g_config.skip_io) {
      for (total_io = 0; total_io < nbytes_read; ) {
        total_io += pread64(infile, input + total_io, std::min(nbytes_read - total_io, MAX_BUFFER_SZ), from + total_io);
      }
    }

    const auto bases = std::min(remaining_bases, nbytes_read * BPPB);
    decode(input, output, bases);

    if (!g_config.skip_io) {
      for (total_io =- 0; total_io < bases; ) {
        total_io += pwrite64(outfile, output + total_io, std::min(bases - total_io, MAX_BUFFER_SZ), BPPB * (from-HEADER_SZ) + total_io);
      }
    }

    remaining_bases -= bases;
    from += nbytes_read;
  }

  delete[] output;
  delete[] input;
}

//////////////////////////////////////////////////////////////////////////////
// Read, decode, and write a subset of the input using asynchronous I/O
// This is a concurrent algorithm that works with two sets of input and output
// buffers, so that a previous write and next read can take place concurrently
// with the current encoding.
void
decode_file_chunk_async(int infile, int outfile, fsize_t from, fsize_t to, fsize_t remaining_bases)
{
  // Two sets of input and output data variables, for simultaneous operations:
  struct aiocb read1 = create_aio(infile, new char[g_config.chunk_sz]);
  struct aiocb read2 = create_aio(infile, new char[g_config.chunk_sz]);
  struct aiocb write1 = create_aio(outfile, new char[g_config.chunk_sz * BPPB]);
  struct aiocb write2 = create_aio(outfile, new char[g_config.chunk_sz * BPPB]);
  struct aiocb *cur_read = &read1, *nxt_read = &read2, *cur_write = &write1, *nxt_write = &write2;

  // Can't handle too-large chunks in aio:
  if (g_config.chunk_sz > MAX_BUFFER_SZ / 4) {
    g_config.chunk_sz = MAX_BUFFER_SZ / 4;
  }

  // Before main loop, read in first chunk of data:
  cur_read->aio_offset = from;
  cur_read->aio_nbytes = next_read_size(from, to);
  aio_read(cur_read);

  // main loop:
  while (from < to) {

    // Start reading the next input chunk:
    nxt_read->aio_offset = cur_read->aio_offset + cur_read->aio_nbytes;
    nxt_read->aio_nbytes = next_read_size(nxt_read->aio_offset, to);
    aio_read(nxt_read); // Could be zero bytes

    // Make sure previous write is completed before encoding into that buffer:
    aio_wait(cur_write);
    std::swap(nxt_write, cur_write);

    // Ensure current input has been read, then encode it and start writing it:
    const auto bases = std::min(remaining_bases, fsize_t(cur_read->aio_nbytes * BPPB));
    aio_wait(cur_read);
    decode(static_cast<char*>(const_cast<void*>(cur_read->aio_buf)),
           static_cast<char*>(const_cast<void*>(nxt_write->aio_buf)), bases);
    nxt_write->aio_offset = (cur_read->aio_offset - HEADER_SZ) * BPPB;
    nxt_write->aio_nbytes = std::min(fsize_t(BPPB * cur_read->aio_nbytes), bases);
    assert(nxt_write->aio_nbytes <= MAX_BUFFER_SZ);   // Overly-large buffers not supported in aio
    aio_write(nxt_write);

    // Prepare for the next read at the next iteration:
    from += cur_read->aio_nbytes;
    remaining_bases -= bases;
    std::swap(cur_read, nxt_read);
  }

  aio_wait(nxt_write);   // All done, wait for last write:

  delete[] static_cast<char*>(const_cast<void*>(cur_read->aio_buf));
  delete[] static_cast<char*>(const_cast<void*>(nxt_read->aio_buf));
  delete[] static_cast<char*>(const_cast<void*>(cur_write->aio_buf));
  delete[] static_cast<char*>(const_cast<void*>(nxt_write->aio_buf));
}

//////////////////////////////////////////////////////////////////////////////
// Multithreaded implementation of decoding a file.
// Takes file handles for the input and output and the input file size.
void decode_file(int infile, int outfile, fsize_t insize)
{
  std::vector<std::thread> threads;

  int64_t outsize = 0;
  if (read(infile, reinterpret_cast<uint64_t*>(&outsize), sizeof(outsize)) !=
      sizeof(outsize)) {
    std::cerr << "Failed to read size from input file" << std::endl;
    std::exit(3);
  }

  // We need to track how many encoded bases (bytes) are left in the input
  // to account for the last byte, which may not be full
  fsize_t remaining_bases = outsize;

  for (unsigned i = 0; i < g_config.nthreads; ++i) {
    const fsize_t from = i * g_config.task_sz + HEADER_SZ;
    const fsize_t to = (i == (g_config.nthreads - 1)) ?
                        insize :
                       (i + 1) * g_config.task_sz + HEADER_SZ;
    const auto bases = std::min(remaining_bases, (to - from) * BPPB);

    threads.push_back(std::thread(
         g_config.async? decode_file_chunk_async : decode_file_chunk_sync,
         infile,
         outfile,
         from,
         to,
         bases));
    remaining_bases -= bases;
  }

  for (auto &thr : threads) {
    thr.join();
  }

  close(infile);
  close(outfile);
}

/*****************************************************************************
    main
******************************************************************************/
} // namespace fastdna

int
main(int argc, char** argv)
{
  using namespace fastdna;

  // Initialize configuration:
  auto [infile, outfile, size] = parse_cmdline(argc, argv);

  // Start the clock, do the op
  std::cout << (g_config.encode? "Encoding" : "Decoding") <<
    " with configuration:" <<
    " asynchronous I/O: " << g_config.async <<
    " Skipping I/O: " << g_config.skip_io <<
    " Use pext: " << g_config.use_pext <<
    " Ordered bits: " << g_config.ordered <<
    " I/O chunk size: " << g_config.chunk_sz <<
    " threads: " << g_config.nthreads <<
    " task size: " << g_config.task_sz << std::endl;

  const auto t0 = high_resolution_clock::now();
  if (g_config.encode) {
    encode_file(infile, outfile, size);
  } else {
    decode_file(infile, outfile, size);
  }

  const auto run_time = duration_cast<microseconds>(high_resolution_clock::now() - t0).count();
  assert(run_time > 0);
  std::cout << "Time: " << run_time << " μs. ";
  std::cout << "Input size: " << size << " Bytes. ";
  std::cout << std::fixed << "Bandwidth: " << (1e6 * double(size) / run_time) << " Bytes/s\n";

  return 0;
}
