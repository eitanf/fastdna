The `cal.csv` file contains all of the experimental results for the CAL paper, one row per obersvation (there can be multiples copies or observation per experiment).

Each column represents another input or output variable, as detailed next:

 * `experiment` (string): A name for the experiment.
 * `filesystem` (string): The path to the filesystem for input and output files.
 * `model` (string): Model name for filesystem's physical device.
 * `input_size` (integer): Size (in bytes) of the raw input file.
 * `output_size` (integer): Size (in bytes) of the encoded bases file.
 * `threads` (integer): Number of threads used in encoding and decoding.
 * `chunk_sz` (integer): Size (in bytes) of I/O units; 0 signifies whole file.
 * `async` (Boolean): Whether encoder's I/O was asynchronous.
 * `skip_io` (Boolean): Whether encoder and decoder used only RAM buffers.
 * `pext` (Boolean): Whether `pext` was used in encoding and `pdep` in decoding.
 * `unordered` (Boolean): Whether `pext`/`pdep` operations reordered bases.
 * `copy` (integer): repetition or copy number of the same repeated experiment.
 * `eoncoding_time` (float): Time to encode input (seconds).
 * `enoding_bw` (float): Encoder throughput (Bases/sec).
 * `decoding_time` (float): Time to encode input (seconds).
 * `decoding_bw` (float): Decoder throughput (Bases/sec).
