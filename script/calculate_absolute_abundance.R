#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description="Calculate absolute abundance for spike-in sequencing")
parser$add_argument("-s", "--samples", metavar = "sample_list", 
					help = "list of all samples [required]")
parser$add_argument("-i", "--input", metavar = "otu_count", 
					help = "reads count matrix of OTUs(.csv) [required]")
parser$add_argument("-w", "--weight", metavar = "weight_table", 
					help = "table of sample weights(mg) [required]")
parser$add_argument("-p", "--spikein", metavar = "spikein_count", 
					help = "numbers of reads mapped to spike-in strain [required]")
parser$add_argument("-o", "--output", metavar = "output_prefix", 
					help = "prefix of output files [required]")
parser$add_argument("-r", "--renormalize", action = "store_true", default = FALSE,
					help = "renormalize bacterial densities to mean of 1")

argv <- parser$parse_args()
path_sampleList = argv$samples
path_OTUcount = argv$input
path_weight = argv$weight
path_spikein = argv$spikein
if_renormalize = argv$renormalize
prefix_out = argv$output

sample.list <- read.table(path_sampleList, stringsAsFactors = F, header = F)$V1
spikein.count <- read.table(path_spikein, stringsAsFactors = F, header = T, row.names = 1)[sample.list,1]
weight.list <- read.table(path_weight, stringsAsFactors = F, header = T, row.names = 1)[sample.list,1]
OTU.count <- read.csv(path_OTUcount, stringsAsFactors = F, header = T, row.names = 1)[,sample.list]

spikein_abunds <- spikein.count / (colSums(OTU.count) + spikein.count)
abs_abunds_raw <- (1 - spikein_abunds) / (spikein_abunds * weight.list)

if (if_renormalize == TRUE){
	abs_abunds <- abs_abunds_raw / mean(abs_abunds_raw)
}else{
	abs_abunds <- abs_abunds_raw
} 

data.relativeAb <- sweep(OTU.count, 2, colSums(OTU.count), "/")
data.absoluteAb <- sweep(data.relativeAb, 2, abs_abunds, "*")

message()
message(paste("Write total bacterial density to ", prefix_out, ".total_density.tsv", sep = ""))
message(paste("Write relative abundance to ", prefix_out, ".relative_abundance.csv", sep = ""))
message(paste("Write absolute abundance to ", prefix_out, ".absolute_abundance.csv", sep = ""))
write.table(data.frame(sample = sample.list, total_bacterial_density = abs_abunds), 
			paste(prefix_out,"total_density.tsv", sep = "."), quote = F, col.names = T, 
			row.names = F, sep = "\t")
write.csv(data.absoluteAb, paste(prefix_out,"absolute_abundance.csv", sep = "."), quote = F, 
			row.names = T)
write.csv(data.relativeAb, paste(prefix_out,"relative_abundance.csv", sep = "."), quote = F, 
            row.names = T)
message("Complete!")

