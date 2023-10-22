#!/usr/bin/env -S python3 -u
#
# This script runs all the experiments for the paper in sequence.
# It generates its own inputs.

import csv
import os
import random
import re
import sys
import subprocess
import tempfile

inputs = {}
repeats = 100
input_size = 3*(1<<30)
chunk_size = 1*(1<<17)
max_threads = 16
exprfn = "cal.csv"

if len(sys.argv) > 1:
    repeats = int(sys.argv[1])

############
# Generate a random input file of specified size and location
# This is quite slow at large sizes, so it's cached.
def generate_random_input(size, path):
    if size in inputs:
        return inputs[size]

    fn = f"{path}/input.{size}"
    if os.path.isfile(fn) and os.path.getsize(fn) == size:
        inputs[size] = fn
        return inputs[size]

    '''
    faster alternative?
    os.system(f"tr -dc 'ACGT' < /dev/urandom | head -c {size} > {fn}")
    inputs[size] = fn
    '''

    alpha = ['A', 'C', 'G', 'T']
    with open(fn, 'w') as tmp:
        #    with tempfile.NamedTemporaryFile(mode = "w", dir=path, delete=False) as tmp:
        for _ in range(size):
            tmp.write(random.choice(alpha))

    inputs[size] = tmp.name
    return inputs[size]

############
# Figure out filesystem for a path
def get_mount(path):
    path = os.path.realpath(os.path.abspath(path))
    while path != os.path.sep:
        if os.path.ismount(path):
            return path
        path = os.path.abspath(os.path.join(path, os.pardir))
    return path


############
# Calculate the expected compressed file size for a given raw input size
def encoded_size(input_size):
    return 8 + int(((3 + input_size) / 4))


############
# Reset OS file cache
def reset_cache():
    cmd = [ "sudo", "sh", "-c", "/usr/bin/sync; /sbin/sysctl vm.drop_caches=3" ]
    try:
        subprocess.run(cmd, stdout=subprocess.DEVNULL)
    except PermissionError:
        print("Can't reset caches. Consider adding this line to /etc/sudoers:")
        print(f"{os.getlogin()}\tALL=NOPASSWD:\t/sbin/sysctl vm.drop_caches=3")
        raise


############
# Encode a file, then decode it, report back the performance of both.
# Checks that the encoded file is of the correct size and that the decoded file
# matches the raw input file.
# Deletes the encoded and decoded files.
def encode_decode(fn,
                  chunk_sz,
                  threads = 1,
                  asynch = False,
                  skip_io = False,
                  unordered = False,
                  pext = False):

    opts = []
    ret = {}
    if threads > 0: opts += ["-t", str(threads)]
    if chunk_sz >= 0: opts += ["-c", str(chunk_sz)]
    if asynch: opts += ["-a"]
    if skip_io: opts += ["-i"]
    if pext: opts += ["-p"]
    if unordered: opts += ["-u"]

    # run encoder:
    if not skip_io: reset_cache()
    cmd = ["./fastdna"] + opts + ["e"] + [fn]
    out = subprocess.run(cmd, stdout = subprocess.PIPE).stdout.decode('utf-8')
#    print(" ".join(cmd), out)
    m = re.match(".*chunk size: (\d+) threads: (\d+).*\nTime: (\d+).*Bandwidth: (\d*\.?\d+).*", out)
    if not m:   print(f"\n{cmd}\n{out}")
    assert m, "Couldn't parse encoder output"
    enc_ch = int(m.group(1))
    enc_thr = int(m.group(2))
    enc_t = float(m.group(3)) / 1.e6
    enc_bw = float(m.group(4))

    # Check that encoded output is of correct size:
    inp_sz = os.path.getsize(fn)
    enc_sz = os.path.getsize(fn + ".2b")
    assert enc_sz == encoded_size(inp_sz), "Encoded file doesn't match expected size"

    # run decoder:
    if not skip_io: reset_cache()
    cmd = ["./fastdna"] + opts + ["d"] + [fn+".2b"]
    out = subprocess.run(cmd, stdout = subprocess.PIPE).stdout.decode('utf-8')

    m = re.match(".*chunk size: (\d+) threads: (\d+).*\nTime: (\d+).*Bandwidth: (\d*\.?\d+).*", out)
    if not m:   print(cmd, out)
    assert m, "Couldn't find time or bandwidth in decoder output"
    dec_t = int(m.group(3)) / 1.e6
    dec_bw = float(m.group(4))

    # Check that decoded file is same as input, then delete it:
    if not skip_io:
        code = subprocess.run(["cmp", fn, fn + ".2b.dec"]).returncode
        assert code == 0, "Decompressed output is different than input!"
        os.remove(fn + ".2b")
        os.remove(fn + ".2b.dec")

    return enc_t, enc_bw, dec_t, dec_bw, enc_sz, enc_ch, enc_thr


############
def get_filesystem_partition(path):
    cmd = f"/usr/bin/df {path} | tail -1 | cut -f1 -d' '"
    return subprocess.check_output(cmd, shell=True, text=True).strip()


############
def get_disk_model(partition):
    disk = re.sub("p\d+$", "", partition)
    cmd = f"lsblk -o MODEL {disk} | head -2 | tail -1"
    return subprocess.check_output(cmd, shell=True, text=True).strip()


############
# Save a set of runs to file
def write_runs(fn, runs):
    needs_header = not os.path.isfile(fn)
    with open(fn, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames = runs[0].keys())
        if needs_header:
            writer.writeheader()

        for row in runs:
            writer.writerow(row)



experiment_num = 0
############
# Run a fully-specified experiment for a number of copies (repetitions)
def run_experiment(experiment,
                   size = input_size,
                   copies = repeats,
                   chunk_sz = input_size,
                   threads = 1,
                   asynch = False,
                   skip_io = False,
                   unordered = False,
                   pext = False,
                   path = "/tmp"):

    global experiment_num
    if chunk_sz == -1: chunk_sz = size
    runs = []
    tfile = generate_random_input(size, path)

    experiment_num += 1
    print(f"{experiment_num}. {experiment}", end="")
    for copy in range(copies):
        enc_t, enc_bw, dec_t, dec_bw, enc_size, enc_ch, enc_thr = \
            encode_decode(tfile, chunk_sz, threads, asynch, skip_io, unordered, pext)
        tmp = {'experiment': experiment,
               'filesystem': get_filesystem_partition(path),
               'model': get_disk_model(get_filesystem_partition(path)),
               'input_size': size,
               'output_size': enc_size,
               'threads': enc_thr,
               'chunk_sz': enc_ch,
               'async': asynch,
               'skip_io': skip_io,
               'pext': pext,
               'unordered': unordered,
               'copy': copy + 1,
               'encoding_time': enc_t,
               'encoding_bw': enc_bw,
               'decoding_time': dec_t,
               'decoding_bw': dec_bw
               }
        runs.append(tmp)
        print(".", end="")

    print()
    # os.remove(tfile)
    write_runs(exprfn, runs)
    return runs


############
# Main:

run = run_experiment
# Arc of experiments: First, using naive encoding, maximize File I/O
# Then, maximize RAM BW and then cache BW and finally ILP + unrolling


os.system("rm fastdna")
os.system('make fastdna OPTFLAGS="-O3 -DNDEBUG -DUNROLL=1"')

short=repeats
min_ch=10
max_ch=27
step=2

run("Baseline", copies=short)
for pow2 in range(10, 32):
    run("Baseline", size = (1 << pow2), copies=short)

for pow2 in range(10, 28):
    run("Chunked I/O", threads=1, chunk_sz=(1<<pow2), copies=short)

for pow2 in range(10, 28):
    run("Asynchronous I/O", asynch=True, chunk_sz=(1<<pow2), copies=short)

for thr in range(1, max_threads+1):
    run("Multithreading", threads=thr, chunk_sz=chunk_size, copies=short)

run("Other device", threads=8, path="/home/eitan/", copies=short)

for thr in range(1, max_threads+1):
    run("No I/O", threads=thr, skip_io=True)

for pow2 in range(min_ch, max_ch, step):
    run("memory chunked", threads=1, chunk_sz=(1<<pow2), skip_io=True)
    run("memory chunked", threads=2, chunk_sz=(1<<pow2), skip_io=True)
    run("memory chunked", threads=4, chunk_sz=(1<<pow2), skip_io=True)
    run("memory chunked", threads=8, chunk_sz=(1<<pow2), skip_io=True)
    run("memory chunked", threads=16, chunk_sz=(1<<pow2), skip_io=True)

for pow2 in range(min_ch, max_ch, step):
    run("Bit-level parallelism", chunk_sz=(1<<pow2), threads=1, skip_io=True, pext=True)
    run("Bit-level parallelism", chunk_sz=(1<<pow2), threads=2, skip_io=True, pext=True)
    run("Bit-level parallelism", chunk_sz=(1<<pow2), threads=4, skip_io=True, pext=True)
    run("Bit-level parallelism", chunk_sz=(1<<pow2), threads=8, skip_io=True, pext=True)
    run("Bit-level parallelism", chunk_sz=(1<<pow2), threads=16, skip_io=True, pext=True)

for pow2 in range(min_ch, max_ch, step):
    run("Unordered", threads=1, skip_io=True, pext=True, unordered=True, chunk_sz=(1<<pow2))
    run("Unordered", threads=2, skip_io=True, pext=True, unordered=True, chunk_sz=(1<<pow2))
    run("Unordered", threads=4, skip_io=True, pext=True, unordered=True, chunk_sz=(1<<pow2))
    run("Unordered", threads=8, skip_io=True, pext=True, unordered=True, chunk_sz=(1<<pow2))
    run("Unordered", threads=16, skip_io=True, pext=True, unordered=True, chunk_sz=(1<<pow2))

os.system("rm fastdna")
os.system('make fastdna OPTFLAGS="-O3 -DNDEBUG -DUNROLL=2"')

for pow2 in range(min_ch, max_ch, step):
    run("Unroll2", threads=1, skip_io=True, pext=True, unordered=True, chunk_sz=(1<<pow2))
    run("Unroll2", threads=2, skip_io=True, pext=True, unordered=True, chunk_sz=(1<<pow2))
    run("Unroll2", threads=4, skip_io=True, pext=True, unordered=True, chunk_sz=(1<<pow2))
    run("Unroll2", threads=8, skip_io=True, pext=True, unordered=True, chunk_sz=(1<<pow2))
    run("Unroll2", threads=16, skip_io=True, pext=True, unordered=True, chunk_sz=(1<<pow2))

