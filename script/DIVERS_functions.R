#!/usr/bin/env Rscript
library(matrixStats)
library(progress)

checkConfig <- function(df.config, path_config, number_variance){
	if (number_variance == 3){
	    temporal.list <- unique(df.config$temporal)
		labelV <- "temporal"
	}else{
		temporal.list <- unique(df.config$biological)
		labelV <- "biological"
	}
    flag = 0
    for (i in 1:length(temporal.list)){
		if (number_variance == 3){
        	tmp.df <- df.config[which(df.config$temporal == temporal.list[i]),]
		}else{
			tmp.df <- df.config[which(df.config$biological == temporal.list[i]),]
		}
        if (length(tmp.df[which(tmp.df$variable == "X"),"variable"]) > 1){
            message(paste("Error: multiple variable X are found for ", labelV, " ",temporal.list[i], " in ", path_config, sep = ""))
            flag = 1
        }
        if (length(tmp.df[which(tmp.df$variable == "X"),"variable"]) == 0){
            message(paste("Error: variable X is not found for ", labelV, " ", temporal.list[i], " in ", path_config, sep = ""))
            flag = 1
        }
        if (length(tmp.df[which(tmp.df$variable == "Y"),"variable"]) > 1){
            message(paste("Error: multiple variable Y are found for ", labelV, " ", temporal.list[i], " in ", path_config, sep = ""))
            flag = 1
        }
        if (length(tmp.df[which(tmp.df$variable == "Y"),"variable"]) == 0){
            message(paste("Error: variable Y is not found for ", labelV, " ", temporal.list[i], " in ", path_config, sep = ""))
            flag = 1
        }
		if (number_variance == 3){
	        if (length(tmp.df[which(tmp.df$variable == "Z"),"variable"]) > 1){
            	message(paste("Error: multiple variable Z are found for temporal ", temporal.list[i], " in ", path_config, sep = ""))
            	flag = 1
        	}
        	if (length(tmp.df[which(tmp.df$variable == "Z"),"variable"]) == 0){
            	message(paste("Error: variable Z is not found for temporal ", temporal.list[i], " in ", path_config, sep = ""))
            	flag = 1
        	}
		}
    }
    return(flag)
}

getXYZsamples <- function(data.config, temporal.list, types_variance){
	if (types_variance == 3){
		tmp.df.X <- data.config[which(data.config$variable == "X"),]
		rownames(tmp.df.X) <- tmp.df.X$temporal
		Xsamples <- tmp.df.X[temporal.list,"sample"]

		tmp.df.Y <- data.config[which(data.config$variable == "Y"),]
		rownames(tmp.df.Y) <- tmp.df.Y$temporal
		Ysamples <- tmp.df.Y[temporal.list,"sample"]

		tmp.df.Z <- data.config[which(data.config$variable == "Z"),]
		rownames(tmp.df.Z) <- tmp.df.Z$temporal
		Zsamples <- tmp.df.Z[temporal.list,"sample"]
		return(list(Xsamples, Ysamples, Zsamples))
	}else{
		tmp.df.X <- data.config[which(data.config$variable == "X"),]
		rownames(tmp.df.X) <- tmp.df.X$biological
		Xsamples <- tmp.df.X[temporal.list,"sample"]

		tmp.df.Y <- data.config[which(data.config$variable == "Y"),]
		rownames(tmp.df.Y) <- tmp.df.Y$biological
		Ysamples <- tmp.df.Y[temporal.list,"sample"]
		return(list(Xsamples, Ysamples))
	}
}	

calMarginal <- function(data.abs.X, data.abs.Y, data.abs.Z, number_OTU, number_temporal, number_iteration){
	marg_means <- matrix(0.0, nrow = number_OTU, ncol = 1)
	marg_vars <- matrix(0.0, nrow = number_OTU, ncol = 1)
	pb <- progress_bar$new(format = "  processing [:bar] :current/:total eta: :eta", total = number_iteration, clear = FALSE, width = 60)
	for (i in 1:number_iteration){
		pb$tick()
		tmpChoose <- sample(1:3, number_temporal, replace = T)
		data_perm <- as.matrix(cbind(data.abs.X[,which(tmpChoose == 1)], 
							data.abs.Y[,which(tmpChoose == 2)],
							data.abs.Z[,which(tmpChoose == 3)]))
		marg_means <- marg_means + rowMeans(data_perm)
		marg_vars <- marg_vars + rowVars(data_perm)
	}
	average_means <- marg_means / number_iteration
	average_vars <- marg_vars / number_iteration
	out.df <- data.frame(row.names = rownames(data.abs.X), marg_mean = average_means, marg_var = average_vars)
	return(out.df)
}


calMarginal_dual <- function(data.abs.X, data.abs.Y, number_OTU, number_temporal, number_iteration){
	marg_means <- matrix(0.0, nrow = number_OTU, ncol = 1)
	marg_vars <- matrix(0.0, nrow = number_OTU, ncol = 1)
	pb <- progress_bar$new(format = "  processing [:bar] :current/:total eta: :eta", total = number_iteration, clear = FALSE, width = 60)
	for (i in 1:number_iteration){
		pb$tick()
		tmpChoose <- sample(1:2, number_temporal, replace = T)
		data_perm <- as.matrix(cbind(data.abs.X[,which(tmpChoose == 1)], 
							data.abs.Y[,which(tmpChoose == 2)]))
		marg_means <- marg_means + rowMeans(data_perm)
		marg_vars <- marg_vars + rowVars(data_perm)
	}
	average_means <- marg_means / number_iteration
	average_vars <- marg_vars / number_iteration
	out.df <- data.frame(row.names = rownames(data.abs.X), marg_mean = average_means, marg_var = average_vars)
	return(out.df)
}

varianceDecomposition <- function(data.abs.X, data.abs.Y, data.abs.Z, number_OTU, number_temporal, number_iteration){
	covs_XZ <- matrix(0.0, nrow = number_OTU, ncol = 1)
	covs_XmZY <- matrix(0.0, nrow = number_OTU, ncol = 1)
	vars_XmY <- matrix(0.0, nrow = number_OTU, ncol = 1)
	data.abs.XY <- cbind(data.abs.X, data.abs.Y)
	pb <- progress_bar$new(format = "  processing [:bar] :current/:total eta: :eta", total = number_iteration, clear = FALSE, width = 60)
	for (i in 1:number_iteration){
		pb$tick()
		tmpChoose <- sample(1:2, number_temporal, replace = T)
		tmpChoose.X <- 1:number_temporal + (tmpChoose - 1) * number_temporal
		tmpChoose.Y <- 1:number_temporal + (2 - tmpChoose) * number_temporal
		data_X_perm <- data.abs.XY[,tmpChoose.X]
		data_Y_perm <- data.abs.XY[,tmpChoose.Y]
		covs_XZ <- covs_XZ + diag(cov(t(data_X_perm), t(data.abs.Z)))
		covs_XmZY <- covs_XmZY + diag(cov(t(data_X_perm - data.abs.Z), t(data_Y_perm)))
		vars_XmY <- vars_XmY + diag(cov(t(data_X_perm - data_Y_perm), t(data_X_perm - data_Y_perm))) / 2
	}
	vars_T <- covs_XZ / number_iteration
	vars_T[which(vars_T < 0)] = 0
	vars_S <- covs_XmZY / number_iteration
	vars_S[which(vars_S < 0)] = 0
	vars_N <- vars_XmY / number_iteration
	out.df <- data.frame(row.names = rownames(data.abs.X), OTU =  rownames(data.abs.X), vars_T = vars_T, vars_S = vars_S, vars_N = vars_N)
	return(out.df)
}

varianceDecomposition_dual <- function(data.abs.X, data.abs.Y, number_OTU, number_temporal, number_iteration){
	covs_XY <- matrix(0.0, nrow = number_OTU, ncol = 1)
	vars_XmY <- matrix(0.0, nrow = number_OTU, ncol = 1)
	data.abs.XY <- cbind(data.abs.X, data.abs.Y)
	pb <- progress_bar$new(format = "  processing [:bar] :current/:total eta: :eta", total = number_iteration, clear = FALSE, width = 60)
	for (i in 1:number_iteration){
		pb$tick()
		tmpChoose <- sample(1:2, number_temporal, replace = T)
		tmpChoose.X <- 1:number_temporal + (tmpChoose - 1) * number_temporal
		tmpChoose.Y <- 1:number_temporal + (2 - tmpChoose) * number_temporal
		data_X_perm <- data.abs.XY[,tmpChoose.X]
		data_Y_perm <- data.abs.XY[,tmpChoose.Y]

		covs_XY <- covs_XY + diag(cov(t(data_X_perm), t(data_Y_perm)))
		vars_XmY <- vars_XmY + diag(cov(t(data_X_perm - data_Y_perm), t(data_X_perm - data_Y_perm))) / 2
	}
	vars_B <- covs_XY / number_iteration
	vars_B[which(vars_B < 0)] = 0
	vars_N <- vars_XmY / number_iteration
	out.df <- data.frame(row.names = rownames(data.abs.X), OTU =  rownames(data.abs.X), vars_B = vars_B, vars_N = vars_N)
	return(out.df)
}

calTaylorExponents <- function(value_mean, value_var){
	sample.pass <- which(value_mean > 0 & value_var > 0)
	tmpModel <- lm(log10(value_var[sample.pass]) ~ log10(value_mean[sample.pass]))
	out.beta <- c(tmpModel$coefficients[2], tmpModel$coefficients[1])
	names(out.beta) <- c("exponent","log10_constant")
	return(out.beta)	
}

calCovariance <- function(data.abs.X, data.abs.Y, data.abs.Z, number_OTU, number_temporal, number_iteration){
	covs_total <- matrix(0.0, nrow = number_OTU * number_OTU, ncol = 1)
	pb <- progress_bar$new(format = "  processing [:bar] :current/:total eta: :eta", total = number_iteration, clear = FALSE, width = 60)
	for (i in 1:number_iteration){
		pb$tick()
		tmpChoose <- sample(1:3, number_temporal, replace = T)
		data_perm <- as.matrix(cbind(data.abs.X[,which(tmpChoose == 1)], 
							data.abs.Y[,which(tmpChoose == 2)],
							data.abs.Z[,which(tmpChoose == 3)]))
		covmat_perm <- cov(t(data_perm))
		covs_total <- covs_total + as.vector(covmat_perm)
	}
	out.cov <- matrix(covs_total / number_iteration, nrow = number_OTU, ncol = number_OTU)
}

calCovariance_dual <- function(data.abs.X, data.abs.Y, number_OTU, number_temporal, number_iteration){
	covs_total <- matrix(0.0, nrow = number_OTU * number_OTU, ncol = 1)
	pb <- progress_bar$new(format = "  processing [:bar] :current/:total eta: :eta", total = number_iteration, clear = FALSE, width = 60)
	for (i in 1:number_iteration){
		pb$tick()
		tmpChoose <- sample(1:2, number_temporal, replace = T)
		data_perm <- as.matrix(cbind(data.abs.X[,which(tmpChoose == 1)], 
							data.abs.Y[,which(tmpChoose == 2)]))
		covmat_perm <- cov(t(data_perm))
		covs_total <- covs_total + as.vector(covmat_perm)
	}
	out.cov <- matrix(covs_total / number_iteration, nrow = number_OTU, ncol = number_OTU)
}

covarianceDecomposition <- function(data.abs.X, data.abs.Y, data.abs.Z, number_OTU, number_temporal, number_iteration){
	crosscovs_T <- matrix(0.0, nrow = number_OTU * number_OTU, ncol = 1)
	crosscovs_S <- matrix(0.0, nrow = number_OTU * number_OTU, ncol = 1)
	crosscovs_N <- matrix(0.0, nrow = number_OTU * number_OTU, ncol = 1)
	data.abs.XY <- cbind(data.abs.X, data.abs.Y)
	pb <- progress_bar$new(format = "  processing [:bar] :current/:total eta: :eta", total = number_iteration, clear = FALSE, width = 60)
	for (i in 1:number_iteration){
		pb$tick()
		tmpChoose <- sample(1:2, number_temporal, replace = T)
		tmpChoose.X <- 1:number_temporal + (tmpChoose - 1) * number_temporal
		tmpChoose.Y <- 1:number_temporal + (2 - tmpChoose) * number_temporal
		data_X_perm <- data.abs.XY[,tmpChoose.X]
		data_Y_perm <- data.abs.XY[,tmpChoose.Y]

		covmat_XZ <- cov(t(data_X_perm), t(data.abs.Z))
		covmat_ZX <- cov(t(data.abs.Z), t(data_X_perm))
		crosscovs_T <- crosscovs_T + as.vector((covmat_XZ + covmat_ZX) / 2)
	
		covmat_XZY <- cov(t(data_X_perm - data.abs.Z), t(data_Y_perm))
		covmat_YXZ <- cov(t(data_Y_perm), t(data_X_perm - data.abs.Z))
		crosscovs_S <- crosscovs_S + as.vector((covmat_XZY + covmat_YXZ) / 2)

		covmat_XmY <- cov(t(data_X_perm - data_Y_perm), t(data_X_perm - data_Y_perm))
		crosscovs_N <- crosscovs_N + as.vector(covmat_XmY / 2)

	}
	covs_T <- matrix(crosscovs_T / number_iteration, nrow = number_OTU, ncol = number_OTU)
	covs_S <- matrix(crosscovs_S / number_iteration, nrow = number_OTU, ncol = number_OTU)
	covs_N <- matrix(crosscovs_N / number_iteration, nrow = number_OTU, ncol = number_OTU)
	return(list(covs_T, covs_S, covs_N))
}

covarianceDecomposition_dual <- function(data.abs.X, data.abs.Y, number_OTU, number_temporal, number_iteration){
	crosscovs_B <- matrix(0.0, nrow = number_OTU * number_OTU, ncol = 1)
	crosscovs_N <- matrix(0.0, nrow = number_OTU * number_OTU, ncol = 1)
	data.abs.XY <- cbind(data.abs.X, data.abs.Y)
	pb <- progress_bar$new(format = "  processing [:bar] :current/:total eta: :eta", total = number_iteration, clear = FALSE, width = 60)
	for (i in 1:number_iteration){
		pb$tick()
		tmpChoose <- sample(1:2, number_temporal, replace = T)
		tmpChoose.X <- 1:number_temporal + (tmpChoose - 1) * number_temporal
		tmpChoose.Y <- 1:number_temporal + (2 - tmpChoose) * number_temporal
		data_X_perm <- data.abs.XY[,tmpChoose.X]
		data_Y_perm <- data.abs.XY[,tmpChoose.Y]

		covmat_XY <- cov(t(data_X_perm), t(data_Y_perm))
		covmat_YX <- cov(t(data_Y_perm), t(data_X_perm))
		crosscovs_B <- crosscovs_B + as.vector((covmat_XY + covmat_YX) / 2)
	
		covmat_XmY <- cov(t(data_X_perm - data_Y_perm), t(data_X_perm - data_Y_perm))
		crosscovs_N <- crosscovs_N + as.vector(covmat_XmY / 2)
	}
	covs_B <- matrix(crosscovs_B / number_iteration, nrow = number_OTU, ncol = number_OTU)
	covs_N <- matrix(crosscovs_N / number_iteration, nrow = number_OTU, ncol = number_OTU)
	return(list(covs_B, covs_N))
}
