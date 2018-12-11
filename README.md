# DIVERS: Decomposition of Variance Using Replicate Sampling
Code for DIVERS (Decomposition of Variance Using Replicate Sampling), including absolute abundance estimation from spike-in sequencing and the variance/covariance decompostion of absolute bacterial abundances.

## Dependencies

* R (we have tested this code for R version 3.5.0)
	- argparse
	- matrixStats
	- progress

* Python 3.6, jupyter 4.3.0 (for results analysis)
	- pandas
	- numpy
	- matplotlib
	- _NB: Above libraries are bundled together in the [Anaconda distribution](https://www.continuum.io/downloads)_

## Absolute abundance estimation from spike-in sequencing

### Description
```
usage: ./script/calculate_absolute_abundance.R [-h] [-s sample_list]
                                               [-i otu_count]
                                               [-w weight_table]
                                               [-p spikein_count]
                                               [-o output_prefix] [-r]

Calculate absolute abundances from spike-in sequencing

optional arguments:
  -h, --help            show this help message and exit
  -s sample_list, --samples sample_list
                        list of all samples [required]
  -i otu_count, --input otu_count
                        reads count matrix of OTUs (spike-in OTU removed) (.csv) [required]
  -w weight_table, --weight weight_table
                        table of sample weights(mg) [required]
  -p spikein_count, --spikein spikein_count
                        numbers of reads mapped to spike-in OTU [required]
  -o output_prefix, --output output_prefix
                        prefix of output files [required]
  -r, --renormalize     renormalize total bacterial densities to mean of 1
```
### Input format

****[Important] Avoid any delimiter (tab or blackspace) in OTU IDs and sample IDs****

**sample_list:** List of sample IDs

[example: ./test_data/test.sample_list.txt]

```
d16s1r1
d16s2r1
d16s2r2
d17s1r1
...
```

**otu_count:** Matrix of reads counts for each OTU and sample. Spike-in OTU read counts should be removed.

[example: ./test_data/test.OTU_readsCount.csv]

* each row is a OTU and each column is a sample

* matrix should be provided in CSV format

* reads counts from spike-in OTU are excluded

```
,d16s1r1,d16s2r1,d16s2r2,d17s1r1,...
otu_1,3906,4034,3111,7183,...
otu_2,2221,2418,2256,2232,...
otu_3,2763,3294,3067,2481,...
otu_4,553,503,378,622,...
...
```

**weight_table:** Table of sample weights (mg)

[example: ./test_data/test.sample_weight.tsv]

* first column is sample ID and second column is sample weight

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

**spikein_count:** Table of read counts of the spike-in OTU

[example: ./test_data/test.spikein_readsCount.tsv]

* first column is sample ID and second column is number of spike-in OTU reads

* tab-delimited

* First row should be header (sample[tab]spikein)

```
sample  spikein
d16s1r1 8215
d16s2r1 8263
d16s2r2 11423
d17s1r1 4107
...
```

### Output of absolute abundance estimation

* [output_prefix].total_density.tsv
	- if -r, --renormalize is specified, the bacterial densities will be normalized to mean of 1

* [output_prefix].relative_abundance.csv
	- matrix of relative abundances for each OTU and sample

* [output_prefix].absolute_abundance.csv
	- matrix of absolute abundances for each OTU and sample
	- DIVERS will take this file as input for the variance/covariance decomposition

### Example
```
mkdir ./test_output
chmod +x ./script/calculate_absolute_abundance.R

./script/calculate_absolute_abundance.R -s ./test_data/test.sample_list.txt \
					-i ./test_data/test.OTU_readsCount.csv \
					-w ./test_data/test.sample_weight.tsv \
					-p ./test_data/test.spikein_readsCount.tsv \
					-o ./test_output/test \
					-r
```

## DIVERS variance and covariance decomposition of absolute bacterial abundances

### Description
```
usage: ./script/DIVERS.R [-h] [-i abundance_matrix] [-c configure]
                         [-o output_prefix] [-v number_variance]
                         [-n number_iteration] [-t abundance_threshold] [-cv]

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
                        (biological and technical) [default: 3]
  -n number_iteration, --iteration number_iteration
                        number of iterations [default: 500]
  -t abundance_threshold, --threshold abundance_threshold
                        abundance threshold for covariance decomposition,
                        DIVERS will only perform covariance decomposition for
                        OTUs with average abundance greater than
                        abundance_threshold [default: 1e-4]
  -cv, --covariance     write total and decomposed covariance matrices to
                        output. By default, DIVERS will only write total and
                        decomposed correlation matrices to output
```

### Input format

****[Important] Avoid any delimiter (tab or blackspace) in OTU IDs and sample IDs****

**abundance_matrix:** Matrix of absolute abundances for each OTU and example

[example: ./test_output/test.absolute_abundance.csv]

* each row is a OTU and each column is a sample

* matrix should be provided in CSV format

```
,d16s1r1,d16s2r1,d16s2r2,d17s1r1,...
otu_1,0.0629932321456365,0.0313941944754388,0.0243084036859114,0.194264427469352,...
otu_2,0.0358187323593084,0.018817838929502,0.0176276948619145,0.0603644998067095,...
otu_3,0.0445597287297474,0.0256352197823737,0.0239646011265478,0.0670987114786945,...
otu_4,0.00891839666578007,0.00391454631163751,0.0029535765327144,0.0168220066665651,...
...
```

**configure:** Configure file of sample hierarchy

 * **if [number_variance] is specified as 3**
 
   [example: ./test_output/test.sample_info.variance_3.config]
   
	- column 1 is sample ID and **column 5** is the variable label (X, Y or Z).
	
	- columns 2-4 label the time, spatial replicate, and technical replicate number of each sample. DIVERS will only take information from columns 1,2,5. Spatial (column 3) and technical (column 4) index can be assigned as arbitrary value.

	- DIVERS will expect exactly one sample labelled as X for each **temporal index**, one sample labelled as Y for each **temporal index** and one sample labelled as Z for each **temporal index**.

	- tab-delimited
	
	- first row should be header (sample[tab]temporal[tab]spatial[tab]technical[tab]variable)

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

 * **if [number_variance] is specified as 2**
 
    [example: ./test_output/test.sample_info.variance_2.config]
    
 	- column 1 is sample ID and **column 4** is the variable label (X or Y).
	
	- columns 2-3 label the biological replicate and technical replicate number of each sample. DIVERS will only take information from columns 1,2,4. Technical (column 3) index can be assigned as arbitrary value.
	
	- DIVERS will expect exactly one sample labelled as X for each **biological index** and one sample labelled as Y for each **biological index**.
	
	- tab-delimited

	- first row should be header (sample[tab]biological[tab]technical[tab]variable)

```
sample  biological  technical   variable
d16s2r1 16  1   X   
d16s2r2 16  2   Y   
d17s2r1 17  1   X   
d17s2r2 17  2   Y 
...
```

### Output of the DIVERS variance and covariance decomposition model

* [output_prefix].variance_decomposition.tsv
	- result of the variance decomposition, including the **average abundance**, **total abundance variance** and **decomposed variances** of absolute abundances for each OTU
	- if [number_variance] is specified as 3, decomposed variances reflect temporal (Temporal_variances), spatial (Spatial_variances) and technical (Technical_variances) variances.
	- if [number_variance] is specified as 2, decomposed variances reflect biological (Biological_variances) and technical (Technical_variances) variances.

* [output_prefix].correlation_total.csv
	- absolute abundance correlations for every pair of OTUs **(average abundance > [abundance_threshold])**
	
* [output_prefix].correlation_decomposition_[\*].csv
	- Decomposed abundance correlations for every pair of OTUs **(average abundance > [abundance_threshold])**
	- if [number_variance] is specified as 3, decomposed correlations will reflect temporal, spatial and technical correlations.
	- if [number_variance] is specified as 2, decomposed correlations will  reflect biological and technical correlations.

**if -cv, --covariance is specified**

* [output_prefix].covariance_total.csv
	- absolute abundance covariances for every pair of OTUs **(average abundance > [abundance_threshold])**
	
* [output_prefix].covariance_decomposition_[\*].csv
	- decomposed abundance covariances for every pair of OTUs **(average abundance > [abundance_threshold])**
	- if [number_variance] is specified as 3, decomposed covariances will reflect temporal, spatial, and technical covariances.
	- if [number_variance] is specified as 2, decomposed covariances will reflect biological and technical covariances.



### Example
```
chmod +x ./script/DIVERS.R

# For temporal, spatial and technical variance decomposition
./script/DIVERS.R -i ./test_output/test.absolute_abundance.csv \
		  -c ./test_data/test.sample_info.variance_3.config \
		  -o ./test_output/test.variance_3 \
		  -v 3 \
		  -t 1e-4
		  -n 500
		  
# For biological and technical variance decomposition
./script/DIVERS.R -i ./test_output/test.absolute_abundance.csv \
		  -c ./test_data/test.sample_info.variance_2.config \
		  -o ./test_output/test.variance_2 \
		  -v 2 \
		  -t 1e-4
		  -n 500
```

## Data analysis

We have provided some example analysis code to investigate the resulting data. Check out the demo notebook, where we analyze the DIVERS result of fecal samples spike-in sequencing from Fig. 1d and Fig. 3a/b in the manuscript [_demo/DIVERS_analysis.ipynb_](demo/DIVERS_analysis.ipynb)
```
$ cd demo
$ jupyter notebook
```
