## EINSim: Working Examples

Using EINSim requires two types of datasets:

- **Simulation** results are those generated using EINSim
- **Experimental** results are measured from real devices (*simulation* results can be substituted for *experimental* results since they share the same format. This provides a demonstrative example of how EINSim may be used with real experimental data, and we refer to the substitute data as *proxy-experimental* results when doing so).

Due to the large number of simulation results required for high accuracy (explained in detail in our DSN 2019 paper [1]), we provide two different example datasets:

- **Full-size dataset:** contains all of the simulation data that we use in our paper [1], amounting to several gigabytes of output. These are available for separate download under GitHub "Releases" at https://github.com/CMU-SAFARI/EINSim.
- **Reduced-size dataset:** is a subset of the full-size dataset that is provided directly in this repository and contains enough simulation data to cleanly demonstrate how EINSim output may be applied to experimental data.

## Experimental data

Unfortunately, we are unable to release the original experimental data that we present in our paper [1] for confidentiality reasons. Instead, we provide *proxy-experiment* examples that are themselves generated using EINSim simulation. While these do not constitute real experimental data, they show how EINSim would be used if such data were available.

The *proxy experimental* results that we provide are located under ```examples/proxy_experimental```:

- **example/proxy_experiment/proxy_experiment**: unmodified EINSim output
- **example/proxy_experiment/proxy_experiment_noised_{0.01 - 0.10}**: EINSim output with normally-distributed noise added in proportion to each histogram bin value with a {0.01 - 0.10} multiplier

These results are generated using EINSim with the following configuration:

```
./path/to/einsim -o <outfile> -m s -b 256 -d RANDOM -e UNIFORM_RANDOM -c ALL_TRUE_OR_ALL_ANTI -r 0.1 -k 128 -s HSC -p 0 -n 1000000
```

The noised examples demonstrate how the inference is fairly tolerant to added noise up to a certain point (\~0.03), after which the log-likelihoods begin to decrease and the models begin to mispredict the proxy-experimental configuration.

## Full-size datasets (separate download)

The example simulation datasets that we use in our paper [1] are provided as GitHub "Releases" at https://github.com/CMU-SAFARI/EINSim. The distributed files are xz-compressed for space efficiency using ```xz -9e``` (note that the analysis scripts can directly consume xz/gz/bz2 compressed files).  

We provide simulation datasets:

- example/simulation/full_size_dataset/BCH_0xFF_T1.xz (11.50 MiB)        
- example/simulation/full_size_dataset/BCH_0xFF_T2.xz (13.03 MiB)        
- example/simulation/full_size_dataset/BCH_0xFF_T3.xz (13.68 MiB)        
- example/simulation/full_size_dataset/BCH_RANDOM_T1.xz (9.34 MiB)          
- example/simulation/full_size_dataset/BCH_RANDOM_T2.xz (9.93 MiB)          
- example/simulation/full_size_dataset/BCH_RANDOM_T3.xz (10.09 MiB)          
- example/simulation/full_size_dataset/HSC_0xFF_T1.xz (12.35 MiB)        
- example/simulation/full_size_dataset/HSC_RANDOM_T1.xz (9.36 MiB)          
- example/simulation/full_size_dataset/REP_0xFF_T1.xz (12.02 MiB)        
- example/simulation/full_size_dataset/REP_0xFF_T2.xz (13.85 MiB)        
- example/simulation/full_size_dataset/REP_0xFF_T3.xz (16.09 MiB)        
- example/simulation/full_size_dataset/REP_RANDOM_T1.xz (10.62 MiB)          
- example/simulation/full_size_dataset/REP_RANDOM_T2.xz (13.00 MiB)          
- example/simulation/full_size_dataset/REP_RANDOM_T3.xz (15.36 MiB)          

These datasets comprise 980000 simulated models are generated as a cross-product of the following parameters:

- RBERs from 0.0 to 0.2 in 10000 steps
- RANDOM and 0xFF data patterns
- Four different 'p' {0, 1, 2, 3} (aggregated)
- ECC schemes: Hamming (d = 3), BCH (d = 3, 5, 7), Repetition (d = 3, 5, 7)
- Data-word lengths: {4, 8, 16, 32, 64, 128, 256} (note: BCH simulations are output with non-truncated lengths)
- Burst length: 256
- Uniform-random error distribution
- Bursts that are entirely either true-cells or anti-cells with 50% probability

## Reduced-size example datasets (included)

The reduced-size example dataset comprises 120000 simulated models taken as a subset of the full-size example dataset. We include results for the following restrictions:

- ECC schemes: Hamming (d = 3), BCH (d = 5)
- Data-word lengths: {64, 128, 256}

They are found at:

- example/simulation/reduced_size_dataset/BCH_0xFF_T2_reduced.xz (4.77 MiB)
- example/simulation/reduced_size_dataset/BCH_RANDOM_T2_reduced.xz (3.74 MiB)
- example/simulation/reduced_size_dataset/HSC_0xFF_T1_reduced.xz (4.58 MiB)
- example/simulation/reduced_size_dataset/HSC_RANDOM_T1_reduced.xz (3.62 MiB)
	
## Running the examples

To run EIN on the examples, one simply specifies the relevant *simulation* and *experiment* files to the ```ein.py``` script:

```
$ ./script/ein.py <simulation file> [<simulation file> ...] -e <experimental file> -b <# bootstrapping samples>
```

```ein.py``` will compute a confidence interval (min/max) using bootstrapping if the ```-b``` argument is provided. Note that if ```-b``` is omitted, ```ein.py``` will default to no bootstrapping, and the output min/max values will be identical. Keep in mind that the script runtime will be roughly proportional to the number of bootstrap samples requested.

The provided proxy-experimental data can be run as follows:

```
$ ./script/ein.py example/simulation/reduced_size_dataset/* -e example/proxy_experiment/proxy_experiment
```

Adding the ```-g``` flag will plot the 20 best-fit models in the same style as presented in our paper [1].

## Working with and extending the example data

The example data is split among several files simply for convenience; if desired, the data may be combined into one file using a tool such as ```cat``` or subdivided further using ```grep``` (or the equivalent ```xzutils```). 

New simulation or experimental results can be incorporated simply by passing them into the ```ein.py``` script for analysis.

## Bibliography

[\[1\] Minesh Patel, Jeremie S. Kim, Hasan Hassan, and Onur Mutlu, "Understanding and Modeling On-Die Error Correction in Modern DRAM: An Experimental Study Using Real Devices", in the Proceedings of the 49th Annual IEEE/IFIP International Conference on Dependable Systems and Networks (DSN 2019), Portland, OR, USA, June 2019.](https://people.inf.ethz.ch/omutlu/pub/understanding-and-modeling-in-DRAM-ECC_dsn19.pdf)