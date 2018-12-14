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
                    help = "depth of variance hierarchy to be decomposed, should be either 3 (temporal, spatial and technical) or 2 (biological and technical) [default: 3]")
parser$add_argument("-n", "--iteration", metavar = "number_iteration", 
					type = "integer", default = 500,
                    help = "number of iterations [default: 500]")

parser$add_argument("-t", "--threshold", metavar = "abundance_threshold", 
					type = "double", default = 1e-4,
                    help = "abundance threshold for covariance decomposition, DIVERS will only perform covariance decomposition for OTUs with average abundance greater than abundance_threshold [default: 1e-4]")
parser$add_argument("-cv", "--covariance", action = "store_true", default = FALSE,
                                        help = "write total and decomposed covariance matrices to output. By default, DIVERS will only write total and decomposed correlation matrices to output")


argv <- parser$parse_args()
path_abs <- argv$input
path_config <- argv$config
prefix_out <- argv$output
types_variance <- argv$variance
number_iteration <- argv$iteration
abundance_threshold <- argv$threshold
if_covariance <- argv$covariance

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

if (types_variance == 3){
	if (ncol(data.config) != 5){
		message("Error: format of configure file is wrong")
		message("Exit!")
		message()
		quit()
	}else{
		colnames(data.config) <- c("sample","temporal", "spatial", "technical", "variable")
	}
}else{
	if (ncol(data.config) != 4){
		message("Error: format of configure file is wrong")
		message("Exit!")
		message()
		quit()
	}else{
		colnames(data.config) <- c("sample","biological", "technical", "variable")
	}
}

message("0.Check configure file...")
flag.configCheck <- checkConfig(data.config, path_config, types_variance)
if (flag.configCheck == 1){
	message("Exit!")
	message()
	quit()
}
message("Complete!")
message()

if (types_variance == 3){
	temporal.list <- as.character(unique(data.config$temporal))
}else{
	temporal.list <- as.character(unique(data.config$biological))
}

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
	df.varianceResult <- df.varianceResult[,c("OTU", "means_total", "vars_total", "vars_B", "vars_N")]
}
message(paste("  write variance decomposition result to ", prefix_out, ".variance_decomposition.tsv", sep = ""))

out.df.varianceResult <- df.varianceResult
if (types_variance == 3){
	colnames(out.df.varianceResult) <- c("OTU_ID", "Average_abundances", "Total_variances","Temporal_variances", "Spatial_variances", "Technical_variances")
	write.table(out.df.varianceResult, paste(prefix_out, "variance_decomposition.tsv", sep = "."),
									quote = F, row.names = F, col.names = T, sep = "\t")
}else{
	colnames(out.df.varianceResult) <- c("OTU_ID", "Average_abundances", "Total_variances","Biological_variances", "Technical_variances")
	write.table(out.df.varianceResult, paste(prefix_out, "variance_decomposition.tsv", sep = "."),
									quote = F, row.names = F, col.names = T, sep = "\t")
}

#if (types_variance == 3){
#	beta.total <- calTaylorExponents(df.varianceResult$means_total, df.varianceResult$vars_total)
#	beta.T <- calTaylorExponents(df.varianceResult$means_total, df.varianceResult$vars_T)
#	beta.S <- calTaylorExponents(df.varianceResult$means_total, df.varianceResult$vars_S)
#	beta.N <- calTaylorExponents(df.varianceResult$means_total, df.varianceResult$vars_N)
#	df.betaResult <- data.frame(variance_type = c("total", "temporal(T)", "spatial(S)", "technical(N)"),
#								rbind(beta.total, beta.T, beta.S, beta.N))
#}else{
#	beta.total <- calTaylorExponents(df.varianceResult$means_total, df.varianceResult$vars_total)
#	beta.B <-  calTaylorExponents(df.varianceResult$means_total, df.varianceResult$vars_B)
#	beta.N <-  calTaylorExponents(df.varianceResult$means_total, df.varianceResult$vars_N)
#	df.betaResult <- data.frame(variance_type = c("total", "biological(B)", "technical(N)"),
#								rbind(beta.total, beta.B, beta.N))
#}
#
#message(paste("  write Taylor's law exponents to ", prefix_out, ".taylor_law_exponents.tsv", sep = ""))
#write.table(df.betaResult, paste(prefix_out, "taylor_law_exponents.tsv", sep = "."),
#								quote = F, row.names = F, col.names = T, sep = "\t")

message("Complete!")
message()

OTUs.pass <- as.character(df.varianceResult[which(df.varianceResult$means_total > abundance_threshold), "OTU"])

if (length(OTUs.pass) == 0){
	message(paste("There is no OTU with average abundance greater than ", abundance_threshold, sep = ""))
	message("Skip covariance decompostion...")
	message()
	message("DIVERS complete!")
	message()
	quit()
}

message("3.Calculate total covariances...")
if (types_variance == 3){
	covs.total <- calCovariance(data.abs.X[OTUs.pass,], data.abs.Y[OTUs.pass,], data.abs.Z[OTUs.pass,], length(OTUs.pass), number_temporal, number_iteration)
}else{
	covs.total <- calCovariance_dual(data.abs.X[OTUs.pass,], data.abs.Y[OTUs.pass,], length(OTUs.pass), number_temporal, number_iteration)
}
message("Complete!")
message()

message("4.Covariance decomposition of OTU abundances...")
if (types_variance == 3){
	list.covarianceResult <- covarianceDecomposition(data.abs.X[OTUs.pass,], data.abs.Y[OTUs.pass,], data.abs.Z[OTUs.pass,], length(OTUs.pass), number_temporal, number_iteration)
	covs.T <- list.covarianceResult[1][[1]]
	covs.S <- list.covarianceResult[2][[1]]
	covs.N <- list.covarianceResult[3][[1]]
	colnames(covs.total) <- OTUs.pass
	rownames(covs.total) <- OTUs.pass
	colnames(covs.T) <- OTUs.pass
	rownames(covs.T) <- OTUs.pass
	colnames(covs.S) <- OTUs.pass
	rownames(covs.S) <- OTUs.pass
	colnames(covs.N) <- OTUs.pass
	rownames(covs.N) <- OTUs.pass
}else{
	list.covarianceResult <- covarianceDecomposition_dual(data.abs.X[OTUs.pass,], data.abs.Y[OTUs.pass,], length(OTUs.pass), number_temporal, number_iteration)
	covs.B <- list.covarianceResult[1][[1]]
	covs.N <- list.covarianceResult[2][[1]]
	colnames(covs.total) <- OTUs.pass
	rownames(covs.total) <- OTUs.pass
	colnames(covs.B) <- OTUs.pass
	rownames(covs.B) <- OTUs.pass
	colnames(covs.N) <- OTUs.pass
	rownames(covs.N) <- OTUs.pass
}

if (if_covariance == TRUE){
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
		write.csv(covs.B, paste(prefix_out, "covariance_decomposition_biological.csv", sep = "."),
							quote = F, row.names = T)
		write.csv(covs.N, paste(prefix_out, "covariance_decomposition_technical.csv", sep = "."),
							quote = F, row.names = T)
	}
}
variance.sigxsigy <- sqrt(df.varianceResult[OTUs.pass,]$vars_total) %*% sqrt(t(df.varianceResult[OTUs.pass,]$vars_total))
cors.total <- covs.total / variance.sigxsigy

message(paste("  write correlation decomposition results to ", prefix_out, ".correlation_decomposition_*.csv", sep = ""))
if (types_variance == 3){
	cors.T <- covs.T / variance.sigxsigy
	cors.S <- covs.S / variance.sigxsigy
	cors.N <- covs.N / variance.sigxsigy
	colnames(cors.total) <- OTUs.pass
	rownames(cors.total) <- OTUs.pass
	colnames(cors.T) <- OTUs.pass
	rownames(cors.T) <- OTUs.pass
	colnames(cors.S) <- OTUs.pass
	rownames(cors.S) <- OTUs.pass
	colnames(cors.N) <- OTUs.pass
	rownames(cors.N) <- OTUs.pass
	write.csv(cors.total, paste(prefix_out, "correlation_total.csv", sep = "."),
							quote = F, row.names = T)
	write.csv(cors.T, paste(prefix_out, "correlation_decomposition_temporal.csv", sep = "."),
							quote = F, row.names = T)
	write.csv(cors.S, paste(prefix_out, "correlation_decomposition_spatial.csv", sep = "."),
							quote = F, row.names = T)
	write.csv(cors.N, paste(prefix_out, "correlation_decomposition_technical.csv", sep = "."),
							quote = F, row.names = T)
}else{
	cors.B <- covs.B / variance.sigxsigy
	cors.N <- covs.N / variance.sigxsigy
	colnames(cors.total) <- OTUs.pass
	rownames(cors.total) <- OTUs.pass
	colnames(cors.B) <- OTUs.pass
	rownames(cors.B) <- OTUs.pass
	colnames(cors.N) <- OTUs.pass
	rownames(cors.N) <- OTUs.pass
	write.csv(cors.total, paste(prefix_out, "correlation_total.csv", sep = "."),
							quote = F, row.names = T)
	write.csv(cors.B, paste(prefix_out, "correlation_decomposition_biological.csv", sep = "."),
							quote = F, row.names = T)
	write.csv(cors.N, paste(prefix_out, "correlation_decomposition_technical.csv", sep = "."),
							quote = F, row.names = T)
}
message("Complete!")
message()
message("DIVERS complete!")
message()
