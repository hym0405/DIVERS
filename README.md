# DIVERS: Decomposition of Variance Using Replicate Sampling
code for DIVERS (Decomposition of Variance Using Replicate Sampling), including absolute abundance calculation for spike-in sequencing and variance decompostion of absolute abundance

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
$ usage: ./script/calculate_absolute_abundance.R [-h] [-s sample_list]
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
sample_list: list of samples (./test_data/test.sample_list.txt)
```
d16s1r1
d16s2r1
d16s2r2
d17s1r1
...
```

create BLAST database for determining spacer origins. this only needs to be done once. note that if the ncbi-BLAST bin is already on your path, the script can be executed without the path argument. *we only provide the reference for the main pRec/pTrig recording strain in the _ref_ folder to save space in the repo, but references for the other recording strains can be easily recreated using plasmid sequences from the _plasmid_maps_ folder.*
```
$ ./extraction/build_blast_db.sh [bin_path (optional)]
```

search the spacers against the BLAST database. again, note that if the ncbi-BLAST bin is already on your path, the script can be executed without the path argument.
```
$ ./extraction/blast_search.sh [spacer_dir] [reference_fasta] [bin_path (optional)]
```

determine unique spacers from the BLAST search results
```
$ ./extraction/unique_spacers.py [working_directory]
```

## data analysis

to get you started tinkering with the TRACE system, we have provided some example analysis code to investigate the resulting data. check out the demo notebook, where we analyze the 4 day temporal recording experiment from Fig. 2 and 3 in the Science manuscript: [_demo/trace_4day_analysis.ipynb_](demo/trace_4day_analysis.ipynb)
```
$ cd demo
$ jupyter notebook
```
