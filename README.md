## EINSim: Statistically Inferring ECC Characteristics Using Only Post-Correction Errors 

The Error Inference Simulator (EINSim) uses Monte-Carlo simulation to model how different error-correction code (ECC) schemes transform the spatial distribution of errors in a memory device (e.g., DRAM data-retention errors). EINSim simulates an ECC dataword's lifetime throughout ECC encoding, error injection, and ECC decoding. By simulating enough datawords (e.g., > 1e6), the statistical distribution of post-correction errors approximates that of a real device. 

EINSim is a key part of the EIN statistical inference methodology originally proposed in our DSN 2019 academic paper [1]. EIN uses maximum-a-posteriori (MAP) estimation over different ECC models (computed using EINSim) to infer a device's pre-correction error characteristics from post-correction error distributions. This information can be used to:

- Reverse-engineer characteristics of an unknown ECC scheme (e.g., ECC data-/code-word length, ECC correction capability)
- Infer pre-correction error rates that are obfuscated by ECC

EIN is especially useful for devices in which the pre-correction error characteristics are hidden (e.g., DRAM with on-die ECC). Our paper [1] explains EIN and its potential applications in much greater detail.

*Please send questions to Minesh Patel at minesh.patelh@gmail.com*

## Codebase Overview

Current version: 2.0.0

EINSim comprises two parts:
- A C++ command-line tool that simulates the raw statistical distributions
- Python scripts that apply the simulation results to real experimental data

We use Doxygen to document the source code and provide a Doxyfile for building HTML and LaTeX documentation. To build the documentation, simply issue:

```
$ doxygen
```
when in the project directory, or point ```doxygen``` to the provided Doxyfile. The HTML documentation will be built under ```doxygen/html/index.html```. 

## Dependencies

Building and running EINSim requires the following working installations:

- A C++11 toolchain (e.g., GCC, Clang, MSVC)
	- pthreads

- Python 2/3
	- Matplotlib 2.1.1+ (for graph generation only)
	- SciPy
	- NumPy

EINSim has been tested with:

	- GCC 9.1.0
	- GCC 8.3.0
	- GCC 7.4.0
	- GCC 6.3.0
	- Apple LLVM 10.0.0 (clang-1000.11.45.5)
	- Apple LLVM 9.1.0 (clang-902.0.39.2)
	- Python 3.7.3 (Anaconda 3)
	- Python 3.7.2
	- Python 3.7.1
	- Python 3.5.3
	- Python 2.7.13

## Building in Linux

```
$ make [-j <# threads>] [other make options] <target>
```

The makefile has various targets, described as follows:

- ```release``` builds ```einsim``` with full optimizations
- ```debug``` builds ```einsim.d``` with no optimization and debug symbols
- ```doc``` builds Doxygen documentation in the directory ```doxygen```
- ```all``` builds both ```release``` and ```debug```
- ```clean``` cleans build and binary files for both ```release``` and ```debug```

Omitting the ```target``` argument defaults to the ```release``` configuration. 

## Usage

EINSim runs as a command-line tool with several CLI options that are shown when running EINSim without options:

```
$ ./path/to/einsim
```

Quickly testing whether EINSim works correctly can be done by issuing:
```
$ ./path/to/einsim -m t -T fast [-t <# threads>]
```
This will run a brief test of all implemented ECC schemes across various code parameters and will report success or failure. 

EINSim has three modes of operation:

- ```Simulation``` (```./path/to/einsim -m s```) is the primary mode of operation and simulates how the spatial distribution of errors is transformed by different ECC schemes.

- ```Test``` (```./path/to/einsim -m t```) is used primarily for sanity-checking and regression testing. It sweeps through a set of pre-baked configurations for each ECC scheme that are designed to roughly cover the range of important parameters for the scheme. This helps check for obvious errors when modifying, updating, or adding ECC schemes.  

- ```Debug``` (```./path/to/einsim -m d```) is used for debugging various components of EINSim. It is similar to the *simulation* mode, but provides a simpler core simulation loop with more controls over an ECC word. While *debug* mode is designed to make working with EINSim easier, it need not be used if the user is comfortable with a different debugging approach.

## EINSim high-level overview

Our paper [1] provides a detailed overview of EINSim's purpose and design, and we highly recommend the interested user to refer to it when delving deeply into EINSim.

At a high level, EINSim simulates the error characteristics of configurable-length *words*, which are referred to as *bursts* in the context of DRAM. Bursts are typically several bytes long (e.g., 32B) and may comprise multiple ECC data words. Therefore, EINSim subdivides each burst into one or more ECC datawords as per the particular ECC schemes that have been chosen for simulation. 

In the core simulation loop, EINSim repeatedly simulates individual bursts through the following steps:

1) Writing the burst with a configurable data pattern (e.g., 0xFF, RANDOM, custom data)
2) ECC encoding
3) Error injection according to error-mechanism parameters (e.g., error distribution, DRAM true-/anti-cell layout, data pattern)
4) ECC decoding
5) Measuring some property of the output (e.g., counting the number of errors in the decoded burst with respect to the originally written data pattern)

EINSim then outputs the final measured property using the output format described in the following section.

## JSON representation of ECC schemes

EINSim represents ECC schemes using JSON with the following members:

- **s**: type of ECC scheme (e.g., HSC for Hamming single-error correcting)
- **k**: dataword length (bits)
- **p**: random seed used to initially generate the ECC code
- **uid**: 64-bit unsigned identifier for this ECC scheme
- **G**: generator matrix
- **H**: parity-check matrix
- **R**: degenerator matrix
- **R**: degenerator matrix
- **miscorrection_profile**: miscorrection profile as used by [BEER](https://github.com/CMU-SAFARI/BEER) (ignored by EINSim)

As an example, a (7, 4, 3) Hamming code:

```
{"s":"HSC","k":4,"p":0,"uid":3285674560482760615,"G":[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[1,1,0,1],[0,1,1,1],[1,0,1,1]],"H":[[1,1,0,1,1,0,0],[0,1,1,1,0,1,0],[1,0,1,1,0,0,1]],"R":[[1,0,0,0,0,0,0],[0,1,0,0,0,0,0],[0,0,1,0,0,0,0],[0,0,0,1,0,0,0]]}
```

Additional examples are provided under the ```examples``` directory.

## EINSim output format

The EINSim output format is verbose by design. Each line of output contains the pre- and post-correction error distribution when simulating N words for a single system configuration. This keeps it relatively easy to search for data concerning a particular configuration using common tools such as ```awk``` or ```grep```.

Since a large N (e.g., N > 1e6) is required to achieve a statistically representative sample of the true error distributions (discussed further in Sections 4 and 5 of Patel et al. [1]), EINSim is designed for high parallelization. EINSim typically divides N into smaller units, and each unit outputs its results as soon as it finishes. The outputs across all units can later be combined into a single distribution using the provided Python script ```script/utils/aggregate_output.py```.

The output format of a single output line is described as follows and provides reasonable human-readability:

```
[DATA] <ECC scheme>: p:<int> t:<int> k:<int> n:<int> m:<int> rber:<float> bl:<int> bcl:<int> ps:<int> ed:<string> cd:<string> dp:<string> [ <error count>:<# words pre-correction>:<# words post-correction> <>:<>:<> ... repeat ]
```

We briefly explain each of the fields:

- **uid**: 64-bit unsigned UID computed for the ECC scheme based on its generator, parity-check, and degenerator matrices - used to identify the ECC scheme
- **nw (number of words)**: number of ECC words simulated within this line of output
- **bl (burst length)**: size of the simulated burst (data bits only; will be subdivided into ECC datawords as required)
- **bcl (burst codeword length)**: size of the simulated burst, including ECC check bits
- **ps (pad size)**: additional zero-padding used for when a burst unevenly divides into ECC words
- **em (error model)**: error model used for error-injection, including parameter values
- **cd (true-/anti-cell distribution)**: the spatial distribution of true-/anti-cells to assume when simulating bursts
- **dp (data pattern)**: the data pattern to assume for all datawords (e.g., 0xFF, RANDOM) 
- **obs (observable)**: the metric to compute across all simulated ECC words (e.g., per-bit error counts)
- **[data]**: data values corresponding to the requested observable

Example output for a BCH 3EC code with a 511-bit codeword, 484-bit dataword, and a 256-bit burst length with a UNIFORM_RANDOM pre-correction error distribution, an equal number of true-/anti-cells, and a RANDOM data pattern:

```
$ ./path/to/einsim -m s -s HSC -k 4 -b 4 -p 0 -w BLOCKS -c ALL_TRUE -n 100 -d ALL_ONES -o PER_BIT_ERROR_COUNT -e DATA_RETENTION,0.5
...
[DATA] uid:3285674560482760615 nw:100 bl:4 bcl:7 ps:0 em:DATA_RETENTION(p:0.500000) cd:ALL_TRUE dp:ALL_ONES obs:PER_BIT_ERROR_COUNT [ 52 47 42 52 : 56 45 42 53 52 53 45 ]
```

The script ```script/aggregate_output.py``` can be used to aggregate multiple outputs across multiple runs (or simply within just one run) as follows:

```
$ ./path/to/script/aggregate_output.py -o <filename: destination> [<filename: input EINSim results 0> <filename: input EINSim results 1> ...]
```

Note: it is *not* always necessary to aggregate simulation outputs, but aggregation saves considerable computation time and disk space when applying many simulation results to real experimental data.

## Applying real experimental data to EINSim output

We provide a Python script ```script/ein/ein.py``` to apply EINSim results to experimental data. ```ein.py``` can be run as follows:

```
$ ./path/to/script/ein/ein.py <simulation data filename> [<simulation data filename> ...] -e <experimental data filename>
```

Experimental data files use same format as simulation outputs for convenience. However, most of the parameters are ignored and may be set to **-1**.

## Using the provided examples

Under ```example```, we provide several example configuration files for both ECC schemes and error models. These can be used as inputs to EINSim as alternatives to passing command-line arguments and are intended to simplify the use of complicated configurations and uses that include automation.

## Running the examples given in the original publication [1]

EINSim tag "full-size-dataset-v1.0" provides full working examples of running EINSim end-to-end, including the *simulation* data used in our paper [1]. Because this EINSim release is not fully backwards-compatible with the initial release, please refer to the README within tag "full-size-dataset-v1.0" for instructions on how to run the original version.

## Licensing

The current version of the simulator is provided as-is under the MIT license.

The following libraries are used and are located under ```lib``` with their own license:
- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
- [cxxopts](https://github.com/jarro2783/cxxopts)
- libtp - Inspired by [CTPL](https://github.com/vit-vit/CTPL)
- [rapidjson](https://rapidjson.org/)
- crc64 - CRC64 routines written by Salvatore Sanfilippo

The BCH encoding/decoding helper functions are found under src/codes/bch_helpers* and are adapted from bch3.c found at http://eccpage.com/. The original license provided by the author, Robert Morelos-Zaragoza, is retained along with a brief description of the modifications we have made.

## Attribution

Please cite the following paper when using EINSim:

[\[1\] Minesh Patel, Jeremie S. Kim, Hasan Hassan, and Onur Mutlu, "Understanding and Modeling On-Die Error Correction in Modern DRAM: An Experimental Study Using Real Devices", in the Proceedings of the 49th Annual IEEE/IFIP International Conference on Dependable Systems and Networks (DSN 2019), Portland, OR, USA, June 2019.](https://people.inf.ethz.ch/omutlu/pub/understanding-and-modeling-in-DRAM-ECC_dsn19.pdf)