# FastDNA 4:1 encoder

The fastdna program can encode DNA ASCII files into binary files that are one quarter of the original file, or decode them back.
The raw inputs are ASCII files that contain only the bases (characters) 'A', 'C', 'G', or 'T'.
The encoded files (`.2b` extension) are binary files, which each base encoded with a unique two-bit combination.

The program has various options to control its speed. Use `./fastdna -h` to list all these options.

## Building it

The program is designed for x86-64 architectures, ideally with the BMI2 instruction set.
To build it on such an architecture, simply type `make`.
You may want to override the compiler choice using CXX, e.g., `make CXX=clang++-12`.

## Testing it

The script `runtests.sh` will go through a set of simple combinations of inputs and parameters.
If all the tests pass, you won't see any error messages.

## Measuring it

Evaluating the full performance potential of fastdna requires large inputs sizes (in the order of GBs.).
But because the contents of the input file don't affect performance, you can generate a random input for benchmarking. One easy way to generate a random file (say, of size 1TB), might be:

```
tr -dc 'ACGT' < /dev/urandom | head -c 1T > inputfile
```

## Rebuilding the document

The document `cal.pdf` can be built from the source markdown file `cal.Rmd`. It requires R with a few packages, as well as the `cal.csv` data file with the results of the experiment.
Both `cal.csv` and `cal.pdf` can be built from scratch using `make` once the prerequisites have been met.
