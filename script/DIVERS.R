#!/usr/bin/env Rscript
library(matrixStats)
suppressPackageStartupMessages(library("argparse"))

argv <- commandArgs(trailingOnly = FALSE)
base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
source(paste(base_dir, "DIVERS_functions.R", sep="/"))

parser <- ArgumentParser(description='DIVERS: decomposition of variance using replicate sampling')
parser$add_argument("-i", "--input", metavar = "abundance_matrix", 
                    help = "matrix of absolute abundance(.csv) [required]")
parser$add_argument("-c", "--config", metavar = "configure", 
                    help = "configure file of sample hierarchy [required]")
parser$add_argument("-o", "--output", metavar = "output_prefix", 
                    help = "prefix of output files [required]")
parser$add_argument("-v", "--variance",  metavar = "number_variance",
					type = "integer", default = 3,
                    help = "depth of variance hierarchy to be decomposed, should be either 3 (temporal, spatial and technical) or 2 (temporal_spatial and technical) [default: 3]")
parser$add_argument("-n", "--iteration", metavar = "number_iteration", 
					type = "integer", default = 1000,
                    help = "number of iterations [default: 1000]")


argv <- parser$parse_args()
path_abs = argv$input
path_config = argv$config
prefix_out = argv$output
types_variance = argv$variance
number_iteration = argv$iteration
message()
if (types_variance == 2){
	message("DIVERS will perform variance decomposition for sample X and Y in configure file")
	message()
}else{
	if(types_variance == 3){
		message("DIVERS will perform variance decomposition for sample X, Y and Z in configure file")
		message()
	}else{
		message("Error: -v, --variance should be either 3 or 2")
		message("Exit!")
		message()
		quit()
	}
}

data.abs <- read.csv(path_abs, stringsAsFactors = F, header = T, row.names = 1)
data.config <- read.table(path_config, stringsAsFactors = F, header = T)
colnames(data.config) <- c("sample","temporal", "spatial", "technical", "variable")

message("0.Check configure file...")
flag.configCheck <- checkConfig(data.config, path_config, types_variance)
if (flag.configCheck == 1){
	message("Exit!")
	message()
	quit()
}
message("Complete!")
message()

temporal.list <- unique(data.config$temporal)
number_temporal <- length(temporal.list)
number_OTU <- nrow(data.abs)

XYZsamples <- getXYZsamples(data.config, temporal.list, types_variance)
data.abs.X <- data.abs[,XYZsamples[1][[1]]]
data.abs.Y <- data.abs[,XYZsamples[2][[1]]]
if (types_variance == 3){
	data.abs.Z <- data.abs[,XYZsamples[3][[1]]]
}

message("1.Calculate marginal means and variances of each OTU...") 
if (types_variance == 3){
	df.marginalInfo <- calMarginal(data.abs.X, data.abs.Y, data.abs.Z, number_OTU, number_temporal, number_iteration)
}else{
	
	df.marginalInfo <- calMarginal_dual(data.abs.X, data.abs.Y, number_OTU, number_temporal, number_iteration)
}
message("Complete!")
message()

message("2.Variance decomposition of OTU abundances...")
if (types_variance == 3){
	df.varianceResult <- varianceDecomposition(data.abs.X, data.abs.Y, data.abs.Z, number_OTU, number_temporal, number_iteration)
	df.varianceResult$vars_total <- df.marginalInfo$marg_var
	df.varianceResult$means_total <- df.marginalInfo$marg_mean
	df.varianceResult <- df.varianceResult[,c("OTU", "means_total", "vars_total", "vars_T", "vars_S", "vars_N")]
}else{
	df.varianceResult <- varianceDecomposition_dual(data.abs.X, data.abs.Y, number_OTU, number_temporal, number_iteration)
	df.varianceResult$vars_total <- df.marginalInfo$marg_var
	df.varianceResult$means_total <- df.marginalInfo$marg_mean
	df.varianceResult <- df.varianceResult[,c("OTU", "means_total", "vars_total", "vars_TS", "vars_N")]
}
message(paste("  write variance decomposition result to ", prefix_out, ".variance_decomposition.tsv", sep = ""))

write.table(df.varianceResult, paste(prefix_out, "variance_decomposition.tsv", sep = "."),
									quote = F, row.names = F, col.names = T, sep = "\t")

if (types_variance == 3){
	beta.total <- calTaylorExponents(df.varianceResult$means_total, df.varianceResult$vars_total)
	beta.T <- calTaylorExponents(df.varianceResult$means_total, df.varianceResult$vars_T)
	beta.S <- calTaylorExponents(df.varianceResult$means_total, df.varianceResult$vars_S)
	beta.N <- calTaylorExponents(df.varianceResult$means_total, df.varianceResult$vars_N)
	df.betaResult <- data.frame(variance_type = c("total", "temporal(T)", "spatial(S)", "technical(N)"),
								rbind(beta.total, beta.T, beta.S, beta.N))
}else{
	beta.total <- calTaylorExponents(df.varianceResult$means_total, df.varianceResult$vars_total)
	beta.TS <-  calTaylorExponents(df.varianceResult$means_total, df.varianceResult$vars_TS)
	beta.N <-  calTaylorExponents(df.varianceResult$means_total, df.varianceResult$vars_N)
	df.betaResult <- data.frame(variance_type = c("total", "temporal_spatial(TS)", "technical(N)"),
								rbind(beta.total, beta.TS, beta.N))
}

message(paste("  write Taylor's law exponents to ", prefix_out, ".taylor_law_exponents.tsv", sep = ""))
write.table(df.betaResult, paste(prefix_out, "taylor_law_exponents.tsv", sep = "."),
								quote = F, row.names = F, col.names = T, sep = "\t")
message("Complete!")
message()

message("3.Calculate total covariances...")
if (types_variance == 3){
	covs.total <- calCovariance(data.abs.X, data.abs.Y, data.abs.Z, number_OTU, number_temporal, number_iteration)
}else{
	covs.total <- calCovariance_dual(data.abs.X, data.abs.Y, number_OTU, number_temporal, number_iteration)
}
message("Complete!")
message()

message("4.Covariance decomposition of OTU abundances...")
if (types_variance == 3){
	list.covarianceResult <- covarianceDecomposition(data.abs.X, data.abs.Y, data.abs.Z, number_OTU, number_temporal, number_iteration)
	covs.T <- list.covarianceResult[1][[1]]
	covs.S <- list.covarianceResult[2][[1]]
	covs.N <- list.covarianceResult[3][[1]]
	colnames(covs.total) <- rownames(data.abs)
	rownames(covs.total) <- rownames(data.abs)
	colnames(covs.T) <- rownames(data.abs)
	rownames(covs.T) <- rownames(data.abs)
	colnames(covs.S) <- rownames(data.abs)
	rownames(covs.S) <- rownames(data.abs)
	colnames(covs.N) <- rownames(data.abs)
	rownames(covs.N) <- rownames(data.abs)
}else{
	list.covarianceResult <- covarianceDecomposition_dual(data.abs.X, data.abs.Y, number_OTU, number_temporal, number_iteration)
	covs.TS <- list.covarianceResult[1][[1]]
	covs.N <- list.covarianceResult[2][[1]]
	colnames(covs.total) <- rownames(data.abs)
	rownames(covs.total) <- rownames(data.abs)
	colnames(covs.TS) <- rownames(data.abs)
	rownames(covs.TS) <- rownames(data.abs)
	colnames(covs.N) <- rownames(data.abs)
	rownames(covs.N) <- rownames(data.abs)
}

message(paste("  write covariance decomposition results to ", prefix_out, ".covariance_decomposition_*.csv", sep = ""))

write.csv(covs.total, paste(prefix_out, "covariance_total.csv", sep = "."),
						quote = F, row.names = T)
if (types_variance == 3){
	write.csv(covs.T, paste(prefix_out, "covariance_decomposition_temporal.csv", sep = "."),
							quote = F, row.names = T)
	write.csv(covs.S, paste(prefix_out, "covariance_decomposition_spatial.csv", sep = "."),
							quote = F, row.names = T)
	write.csv(covs.N, paste(prefix_out, "covariance_decomposition_technical.csv", sep = "."),
							quote = F, row.names = T)
}else{	
	write.csv(covs.TS, paste(prefix_out, "covariance_decomposition_temporal_spatial.csv", sep = "."),
							quote = F, row.names = T)
	write.csv(covs.N, paste(prefix_out, "covariance_decomposition_technical.csv", sep = "."),
							quote = F, row.names = T)
}

variance.sigxsigy <- sqrt(df.varianceResult$vars_total) %*% sqrt(t(df.varianceResult$vars_total))
cors.total <- covs.total / variance.sigxsigy

message(paste("  write correlation decomposition results to ", prefix_out, ".correlation_decomposition_*.csv", sep = ""))
if (types_variance == 3){
	cors.T <- covs.T / variance.sigxsigy
	cors.S <- covs.S / variance.sigxsigy
	cors.N <- covs.N / variance.sigxsigy
	colnames(cors.total) <- rownames(data.abs)
	rownames(cors.total) <- rownames(data.abs)
	colnames(cors.T) <- rownames(data.abs)
	rownames(cors.T) <- rownames(data.abs)
	colnames(cors.S) <- rownames(data.abs)
	rownames(cors.S) <- rownames(data.abs)
	colnames(cors.N) <- rownames(data.abs)
	rownames(cors.N) <- rownames(data.abs)
	write.csv(cors.total, paste(prefix_out, "correlation_total.csv", sep = "."),
							quote = F, row.names = T)
	write.csv(cors.T, paste(prefix_out, "correlation_decomposition_temporal.csv", sep = "."),
							quote = F, row.names = T)
	write.csv(cors.S, paste(prefix_out, "correlation_decomposition_spatial.csv", sep = "."),
							quote = F, row.names = T)
	write.csv(cors.N, paste(prefix_out, "correlation_decomposition_technical.csv", sep = "."),
							quote = F, row.names = T)
}else{
	cors.TS <- covs.TS / variance.sigxsigy
	cors.N <- covs.N / variance.sigxsigy
	colnames(cors.total) <- rownames(data.abs)
	rownames(cors.total) <- rownames(data.abs)
	colnames(cors.TS) <- rownames(data.abs)
	rownames(cors.TS) <- rownames(data.abs)
	colnames(cors.N) <- rownames(data.abs)
	rownames(cors.N) <- rownames(data.abs)
	write.csv(cors.total, paste(prefix_out, "correlation_total.csv", sep = "."),
							quote = F, row.names = T)
	write.csv(cors.TS, paste(prefix_out, "correlation_decomposition_temporal_sptial.csv", sep = "."),
							quote = F, row.names = T)
	write.csv(cors.N, paste(prefix_out, "correlation_decomposition_technical.csv", sep = "."),
							quote = F, row.names = T)
}
message("Complete!")
message()
message("DIVERS complete!")
message()
