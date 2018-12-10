%% Load OTU data from the human gut microbiome
load gut_table.mat;
 

%% Convert to absolute abundances;
[M,N] = size(data);
data_rel = data(2:end,:) ./ repmat(sum(data(2:end,:),1),M-1,1); %relative abundances ignoring spike-in OTU
data_abs = data_rel .* abund_mat(2:end,:) * 1; %average absolute abundance (AU) across samples is normalized to 1

%%Update taxonomic info of OTUs
tax = tax(2:end); 
ptax = ptax(2:end); ptax_conf = ptax_conf(2:end);
ctax = ctax(2:end); ctax_conf = ctax_conf(2:end);
otax = otax(2:end); otax_conf = otax_conf(2:end);
ftax = ftax(2:end); ftax_conf = ftax_conf(2:end);
gtax = gtax(2:end); gtax_conf = gtax_conf(2:end);
otu_ids = otu_ids(2:end);

% Split data into two technical replicates (X and Y) 
data_X = data_abs(:,X_inds);
data_Y = data_abs(:,Y_inds);
[Mx,Nx] = size(data_X);
[My,Ny] = size(data_Y);


%% Calculate marginal means and variances of each OTU 

%Perform 1,000 different re-sampling iterations
marg_means = [];
marg_vars = [];

for i = 1:1e3
    
    data_X_perm = [];
    
    %Randomly draw a sample from data_X, data_Y or data_Z
    for j = 1:Nx
        flip = rand;
        if flip > 1/2
            data_X_perm(:,j) = data_X(:,j);
        else
            data_X_perm(:,j) = data_Y(:,j);
        end
    end
    
    %Estimate mean and variance from this iteration
    marg_means(:,i) = mean(data_X_perm,2);
    marg_vars(:,i) = var(data_X_perm')';
end

%Average over all re-sampling iterations
means = mean(marg_means,2);
vars_total = mean(marg_vars,2);


%% Variance decomposition of OTU abundances %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

covs_XY = [];
vars_XmY = [];

%Perfor 1,000 re-sampling iterations
for i = 1:1e3

    data_X_perm = [];
    data_Y_perm = [];
    
    %Randomly permute X and Y, Z stays the same
    for j = 1:Nx
        flip = rand;
        if flip > .5
            data_X_perm(:,j) = data_X(:,j);
            data_Y_perm(:,j) = data_Y(:,j);
        else
            data_X_perm(:,j) = data_Y(:,j);
            data_Y_perm(:,j) = data_X(:,j);
        end
    end
    
    %Temporal + spatial variance
    cov_XY = [];
    for j = 1:Mx
        covmat_xy = cov(data_X_perm(j,:),data_Y_perm(j,:));
        cov_xY = covmat_xy(1,2);
        cov_XY(j) = cov_xY;
    end
    covs_XY(:,i) = cov_XY';
     
    %Technical
    var_XmY = .5 * var((data_X_perm - data_Y_perm)')';
    vars_XmY(:,i) = var_XmY;
    
    i
end


%% Average over re-sampling iterations

%Temporal variance
vars_ST = mean(covs_XY,2);
vars_ST = max(vars_ST,0); %Ignore negative estimates of variability

%Technical variance
vars_N = mean(vars_XmY,2);

%Fraction of each component to total variance
vf_N = vars_N ./ (vars_N + vars_ST);
vf_ST = vars_ST ./ (vars_N + vars_ST);

%Fraction of each component to total variance
vf_N2 = vars_N ./ vars_total;
vf_ST2 = vars_ST ./ vars_total;


%% Technical variability from multiple sequencing replicates

%Multiple technical replicates
data_Nr = data_abs(:,tech_inds);
means_Nr = mean(data_Nr,2);
vars_Nr = var(data_Nr')';

%Noise floor
abunds = abund_mat(1,:);
noise_floor = var(abunds(tech_inds)) / mean(abunds(tech_inds))^2;


%% Variance decomposition of total bacterial densities %%%%%%%%%%%%%%%%%%%%

% Total absolute abundance densities
abunds = abund_mat(1,:);

%Include days 27 and 48 in this calculation
abunds_X = abunds([X_inds 84 87]);
abunds_Y = abunds([Y_inds 85 88]);
abunds_all = [abunds_X; abunds_Y];


%% Calculate marginal variance of total bacterial densities

%Perform 1,000 re-sampling iterations
for i = 1:1e3
    
    abunds_all_perm = [];
    
    %Randomly shuffle X, Y and Z
    for j = 1:length(abunds_X)
        ord = randperm(2);
        abunds_all_perm(:,j) = abunds_all(ord,j);
    end
    
    var_perm = mean(var(abunds_all_perm')');
    vars_perm(i) = var_perm;
    
end

%Average over 1,000 iterations
abunds_vars_total = mean(vars_perm);


%% Perform variance decomposition of total bacterial densities

abunds_covs_XY = [];
abunds_vars_XmY = [];

%Perform 1,000 re-sampling iterations
for i = 1:1e3
    
    abunds_X_perm = [];
    abunds_Y_perm = [];
    
    for j = 1:length(abunds_X)
        flip = rand;
        if flip > .5
            abunds_X_perm(j) = abunds_X(j);
            abunds_Y_perm(j) = abunds_Y(j);
        else
            abunds_X_perm(j) = abunds_Y(j);
            abunds_Y_perm(j) = abunds_X(j);
        end
    end 
    
    %Temporal + spatial variance
    abunds_covmat_XY = cov(abunds_X_perm,abunds_Y_perm);
    abunds_cov_XY = abunds_covmat_XY(1,2);
    abunds_covs_XY(i) = abunds_cov_XY;
    
    %Technical variance
    abunds_var_XmY = .5 * var([abunds_X_perm - abunds_Y_perm])';
    abunds_vars_XmY(i) = abunds_var_XmY;
    
    i
end


%% Average over re-sampling iterations

%Temporal + spatial variance
abunds_vars_ST = mean(abunds_covs_XY,2);
abunds_vars_ST = max(abunds_vars_ST,0);  %Ignore negative estimates of variability

%Technical variance
abunds_vars_N = mean(abunds_vars_XmY,2);

%Fraction of each component to total variance
abunds_vf_ST = abunds_vars_ST / (abunds_vars_ST + abunds_vars_N);
abunds_vf_N = abunds_vars_N / (abunds_vars_ST + abunds_vars_N);

%Standard deviation over re-sampling iterations
abunds_vf_ST_std = std(abunds_covs_XY / (abunds_vars_ST + abunds_vars_N));
abunds_vf_N_std = std(abunds_vars_XmY / (abunds_vars_ST + abunds_vars_N));

%Fraction of each component to total variance
abunds_vf_ST2 = abunds_vars_ST / (abunds_vars_total);
abunds_vf_N2 = abunds_vars_N / (abunds_vars_total);




%% Covariance decomposition of OTU abundances %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Calculate total covariances
covs_total = [];

%Perform 1,000 re-sampling iterations
for i = 1:1e3

    data_X_perm = [];

    %Randomly choose a sample from of data_X, data_Y, data_Z
    for j = 1:Nx
        flip = rand;
        if flip > 1/2
            data_X_perm(:,j) = data_X(:,j);
        else
            data_X_perm(:,j) = data_Y(:,j);
        end
    end
    
    %Mean center the data 
    data_X_mc = data_X_perm - repmat(mean(data_X_perm,2),1,Nx);
    
    %Calculate covariance matrix
    covmat_X_perm = 1/(Nx-1) * data_X_mc * data_X_mc';
    
    %Store covariance matrix for this iteration
    covs_total(:,:,i) = covmat_X_perm;  
    
    i
end
covmat_total = mean(covs_total,3);
clear covs_total;


%% Temporal + spatial covariance decomposition
crosscovs_ST = [];

%Perform 1,000 re-sampling iterations
for i = 1:1e3

    data_X_perm = [];
    data_Y_perm = [];
    
    %One realization of the data
    for j = 1:Nx
        flip = rand;
        if flip > .5
            data_X_perm(:,j) = data_X(:,j);
            data_Y_perm(:,j) = data_Y(:,j);
        else
            data_X_perm(:,j) = data_Y(:,j);
            data_Y_perm(:,j) = data_X(:,j);
        end
    end
    
    %Mean center matrices;
    data_X_mc = data_X_perm - repmat(mean(data_X_perm,2),1,Nx);
    data_Y_mc = data_Y_perm - repmat(mean(data_Y_perm,2),1,Ny);
    
    %Calculate all pairwise covariances for this realization (2
    %permutations)
    covmat_XY = 1/(Nx-1) * data_X_mc*data_Y_mc'; %(X_i,Y_j)
    covmat_YX = 1/(Nx-1) * data_Y_mc*data_X_mc'; %(Y_i,X_j)
    covmat = 1/2 * (covmat_XY + covmat_YX);
    
    %Store covmat for each realization
    crosscovs_ST(:,:,i) = covmat;    
    i
end

%Average over all iterations
covmat_ST = mean(crosscovs_ST,3);
clear crosscovs_ST;



 %% Technical covariance decomposition
crosscovs_N = [];

%Peform 1,000 re-sampling iterations
for i = 1:1e3

    data_X_perm = [];
    data_Y_perm = [];
    
    %One realization of the data
    for j = 1:Nx
    flip = rand;
        if flip > .5
            data_X_perm(:,j) = data_X(:,j);
            data_Y_perm(:,j) = data_Y(:,j);
        else
            data_X_perm(:,j) = data_Y(:,j);
            data_Y_perm(:,j) = data_X(:,j);
        end
    end
   
    data_XmY = data_X_perm - data_Y_perm;
    
    %Mean center matrices;
    data_XmY_mc = data_XmY - repmat(mean(data_XmY,2),1,Nx);
    
    %Calculate all pairwise covariances for this realization (1
    %permutation)
    covmat = 1/(Nx-1) * data_XmY_mc*data_XmY_mc';
    
    %Store covmat for each realization
    crosscovs_N(:,:,i) = .5 * covmat;   
    i
end

%Average over all iterations
covmat_N = mean(crosscovs_N,3);  
clear crosscovs_N;


%% Re-scale covariances to obtain correlations
L = length(vars_total);

%Calculate product of marginal standard deviations for each pair of OTUs
sigxsigy = [];
for i = 1:L
   for j = 1:i
       if vars_total(i) > 0 && vars_total(j) > 0
            sigxsigy(i,j) = sqrt(vars_total(i))*sqrt(vars_total(j));
            sigxsigy(j,i) = sqrt(vars_total(i))*sqrt(vars_total(j));
       else
            sigxsigy(i,j) = 0;
            sigxsigy(j,i) = 0;
       end
   end
end

cormat_total = covmat_total ./ sigxsigy;
cormat_ST = covmat_ST ./ sigxsigy;
cormat_N = covmat_N ./ sigxsigy;





    
clear covs_XY vars_XmY sigxsigy data_X_perm data_Y_perm  data_X_mc data_Y_mc data_XmY_mc data_XmY 


 save('/Users/brianji/Documents/dv_lab/Manuscripts/noise/revisions/round2/matData/gut_data_dual.mat');

