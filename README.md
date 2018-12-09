# DIVERS: Decomposition of Variance Using Replicate Sampling
code for DIVERS (Decomposition of Variance Using Replicate Sampling), including absolute abundance calculation for spike-in sequencing and variance decompostion of absolute abundance.

## dependencies

* R (we have tested this code for R version 3.5.0)
	- argparse
	- matrixStats
	- progress

* python 3.6, jupyter 4.3.0 (for data visualization)
	- pandas
	- numpy
	- matplotlib
	- _NB: Above libraries are bundled together in the [Anaconda distribution](https://www.continuum.io/downloads)_

## Absolute abundance calculation for spike-in sequencing

### Description
```
usage: ./script/calculate_absolute_abundance.R [-h] [-s sample_list]
                                               [-i otu_count]
                                               [-w weight_table]
                                               [-p spikein_count]
                                               [-o output_prefix] [-r]

Calculate absolute abundance for spike-in sequencing

optional arguments:
  -h, --help            show this help message and exit
  -s sample_list, --samples sample_list
                        list of all samples [required]
  -i otu_count, --input otu_count
                        reads count matrix of OTUs(.csv) [required]
  -w weight_table, --weight weight_table
                        table of sample weights(mg) [required]
  -p spikein_count, --spikein spikein_count
                        numbers of reads mapped to spike-in strain [required]
  -o output_prefix, --output output_prefix
                        prefix of output files [required]
  -r, --renormalize     renormalize bacterial densities to mean of 1
```
### Input format

****[Important] avoid any delimiter (tab or blackspace) in OTU IDs and sample IDs****

**sample_list:** list of sample IDs

[example: ./test_data/test.sample_list.txt]

```
d16s1r1
d16s2r1
d16s2r2
d17s1r1
...
```

**otu_count:** matrix of reads counts for each OTU and sample

[example: ./test_data/test.OTU_readsCount.csv]

* each row is a OTU and each column is a sample

* the matrix should be provided in CSV format

* reads counts for spike-in strain are excluded

```
,d16s1r1,d16s2r1,d16s2r2,d17s1r1,...
otu_1,3906,4034,3111,7183,...
otu_2,2221,2418,2256,2232,...
otu_3,2763,3294,3067,2481,...
otu_4,553,503,378,622,...
...
```

**weight_table:** table of sample weights (mg)

[example: ./test_data/test.sample_weight_tsv]

* the first column is sample ID and second column is sample weight

* tab-delimited

* first row should be header (sample[tab]weight)

```
sample  weight
d16s1r1 41.5
d16s2r1 85.5
d16s2r2 61.6
d17s1r1 49.5
...
```

**spikein_count:** table of reads counts of spike-in strain

[example: ./test_data/test.spikein_readsCount.tsv]

* the first column is sample ID and second column is number of reads mapped to spike-in strain

* tab-delimited

* first row should be header (sample[tab]spikein)

```
sample  spikein
d16s1r1 8215
d16s2r1 8263
d16s2r2 11423
d17s1r1 4107
...
```

### Output for absolute abundance calculation

* [output_prefix].total_density.tsv
	- if -r, --renormalize is specified, the bacterial densities will be normalized to mean of 1

* [output_prefix].relative_abundance.csv
	- matrix of relative abundance for each OTU and sample

* [output_prefix].absolute_abundance.csv
	- matrix of absolute abundance for each OTU and sample
	- DIVERS will take this file as input for variance decompostion

### Example
```
mkdir ./test_output
chmod +x ./script/calculate_absolute_abundance.R

./script/calculate_absolute_abundance.R -s ./test_data/test.sample_list.txt \
					-i ./test_data/test.OTU_readsCount.csv \
					-w ./test_data/test.sample_weight_tsv \
					-p ./test_data/test.spikein_readsCount.tsv \
					-o ./test_output/test
					-r
```

## Variance decompostion of absolute abundance by DIVERS

### Description
```
usage: ./script/DIVERS.R [-h] [-i abundance_matrix] [-c configure]
                         [-o output_prefix] [-v number_variance]
                         [-n number_iteration]

DIVERS: decomposition of variance using replicate sampling

optional arguments:
  -h, --help            show this help message and exit
  -i abundance_matrix, --input abundance_matrix
                        matrix of absolute abundance(.csv) [required]
  -c configure, --config configure
                        configure file of sample hierarchy [required]
  -o output_prefix, --output output_prefix
                        prefix of output files [required]
  -v number_variance, --variance number_variance
                        depth of variance hierarchy to be decomposed, should
                        be either 3 (temporal, spatial and technical) or 2
                        (temporal_spatial and technical) [default: 3]
  -n number_iteration, --iteration number_iteration
                        number of iterations [default: 1000]
```

### Input format

****[Important] avoid any delimiter (tab or blackspace) in OTU IDs and sample IDs****

**abundance_matrix:** matrix of absolute abundance for each OTU and example

[example: ./test_output/test.absolute_abundance.csv]

* each row is a OTU and each column is a sample

* the matrix should be provided in CSV format

```
,d16s1r1,d16s2r1,d16s2r2,d17s1r1,...
otu_1,0.0629932321456365,0.0313941944754388,0.0243084036859114,0.194264427469352,...
otu_2,0.0358187323593084,0.018817838929502,0.0176276948619145,0.0603644998067095,...
otu_3,0.0445597287297474,0.0256352197823737,0.0239646011265478,0.0670987114786945,...
otu_4,0.00891839666578007,0.00391454631163751,0.0029535765327144,0.0168220066665651,...
...
```

**configure:** configure file of sample hierarchy

[example: ./test_output/test.sample_info.config]

* column 1 is sample ID and column 5 is the label of variable.

* column 2-4 are indices for different types of variance, and DIVERS will only take the information of sample hierarchy from column 1,2,5.

* if [number_variance] is specified as 3, DIVERS will expect to get exactly one sample labelled as X for each **temporal index**, one sample labelled as Y for each **temporal index** and one sample labelled as Z for each **temporal index**.

* if [number_variance] is specified as 2, DIVERS will expect to get exactly one sample labelled as X for each **temporal index** and one sample labelled as Y for each **temporal index**. Spatial index can be assigned as arbitrary value.

* tab-delimited

* first row should be header (sample[tab]temporal[tab]spatial[tab]technical[tab]variable)


**if [number_variance] is specified as 3**
```
sample  temporal    spatial technical   variable
d16s1r1 16  1   1   Z   
d16s2r1 16  2   1   X   
d16s2r2 16  2   2   Y   
d17s1r1 17  1   1   Z   
d17s2r1 17  2   1   X   
d17s2r2 17  2   2   Y 
...
```
**if [number_variance] is specified as 2**
```
sample  temporal    spatial technical   variable 
d16s2r1 16  2   1   X   
d16s2r2 16  2   2   Y  
d17s2r1 17  2   1   X
d17s2r2 17  2   2   Y
...
```

### Output for variance decompostion by DIVERS

* [output_prefix].variance_decomposition.tsv
	- result of variance decompostion, including the **mean**, **total variance** and **decomposed variances** of absolute abundance for each OTU
	- if [number_variance] is specified as 3, decomposed variances will be temporal (vars_T), spatial (vars_S) and technical(vars_N) variance.
	- if [number_variance] is specified as 2, decomposed variances will be temporal_spatial (vars_TS) and technical (vars_N) variance.

* [output_prefix].taylor_law_exponents.tsv
	- taylor's law exponents for total variance and 3 or 2 decomposed variances

* [output_prefix].covariance_total.csv
	- covariance of absolute abundance for each pair of OTU
	
* [output_prefix].covariance_decomposition_[\*].csv
	- decomposed covariance of absolute abundance for each pair of OTU
	- if [number_variance] is specified as 3, decomposed covariances will be temporal (vars_T), spatial (vars_S) and technical(vars_N) covariance.
	- if [number_variance] is specified as 2, decomposed covariances will be temporal_spatial (vars_TS) and technical (vars_N) covariance.

* [output_prefix].correlation_total.csv
	- correlation of absolute abundance for each pair of OTU
	
* [output_prefix].correlation_decomposition_[\*].csv
	- decomposed correlation of absolute abundance for each pair of OTU
	- if [number_variance] is specified as 3, decomposed correlations will be temporal (vars_T), spatial (vars_S) and technical(vars_N) correlation.
	- if [number_variance] is specified as 2, decomposed correlations will be temporal_spatial (vars_TS) and technical (vars_N) correlation.


### Example
```
chmod +x ./script/DIVERS.R

./script/DIVERS.R -i ./test_output/test.absolute_abundance.csv \
		  -c ./test_output/test.sample_info.config \
		  -o ./test_output/test \
		  -v 3
		  -n 1000
```

## data analysis

to get you started tinkering with the TRACE system, we have provided some example analysis code to investigate the resulting data. check out the demo notebook, where we analyze the 4 day temporal recording experiment from Fig. 2 and 3 in the Science manuscript: [_demo/trace_4day_analysis.ipynb_](demo/trace_4day_analysis.ipynb)
```
$ cd demo
$ jupyter notebook
```
