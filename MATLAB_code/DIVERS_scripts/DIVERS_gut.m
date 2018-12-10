%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script information

% This script loads the fecal OTU times series data
% (in the MAT file gut_table.mat) and performs the DIVERS variance and
% covariance decompositions of individual OTU abundances
% Variance decomposition results are saved in DIVERS_gut.txt
% The resulting data is saved into a MAT file DIVERS_gut.mat

% User input required here:

%Directory containing DIVERS files
file_dir = ['/Path/To/.../DIVERS/DIVERS_files/'];

%Directory containing DIVERS scripts
save_dir = ['/Path/To/.../DIVERS/DIVERS_scripts/'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Load OTU data from the human gut microbiome
addpath(file_dir);
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

% Split data into technical replicates (X and Y) and second spatial
% replicate (Z)
data_X = data_abs(:,X_inds);
data_Y = data_abs(:,Y_inds);
data_Z = data_abs(:,Z_inds);
[Mx,Nx] = size(data_X);
[My,Ny] = size(data_Y);
[Mz,Nz] = size(data_Z);


%% Calculate marginal means and variances of each OTU 

%Perform 1,000 different re-sampling iterations
marg_means = [];
marg_vars = [];

for i = 1:1e3
    
    data_X_perm = [];
    
    %Randomly draw a sample from data_X, data_Y or data_Z
    for j = 1:Nx
        flip = rand;
        if flip > (2/3)
            data_X_perm(:,j) = data_X(:,j);
        elseif flip > (1/3)
            data_X_perm(:,j) = data_Y(:,j);
        else
            data_X_perm(:,j) = data_Z(:,j);
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
disp(['Performing variance decomposition...'])


covs_XZ = [];
covs_XmZY = [];
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
    
    %Temporal variance
    cov_XZ = [];
    for j = 1:Mx
        covmat_xz = cov(data_X_perm(j,:),data_Z(j,:));
        cov_xz = covmat_xz(1,2);
        cov_XZ(j) = cov_xz;
    end
    covs_XZ(:,i) = cov_XZ';
    
    %Spatial sampling variance
    cov_XmZY = [];
    for j = 1:Mx
        covmat_xmzy = cov(data_X_perm(j,:)-data_Z(j,:),data_Y_perm(j,:));
        cov_xmzy = covmat_xmzy(1,2);
        cov_XmZY(j) = cov_xmzy;
    end
    covs_XmZY(:,i) = cov_XmZY';
    
    %Technical
    var_XmY = .5 * var((data_X_perm - data_Y_perm)')';
    vars_XmY(:,i) = var_XmY;
    
end
disp(['Variance decomposition complete!'])


%% Average over re-sampling iterations

%Temporal variance
vars_T = mean(covs_XZ,2);
vars_T = max(vars_T,0); %Ignore negative estimates of variability

%Spatial sampling variance
vars_S = mean(covs_XmZY,2);
vars_S = max(vars_S,0); %Ignore negative estimates of variability

%Technical variance
vars_N = mean(vars_XmY,2);

%Fraction of each component to total variance
vf_N = vars_N ./ (vars_N + vars_S + vars_T);
vf_S = vars_S ./ (vars_N + vars_S + vars_T);
vf_T = vars_T ./ (vars_N + vars_S + vars_T);

%Fraction of each component to total variance
vf_N2 = vars_N ./ vars_total;
vf_S2 = vars_S ./ vars_total;
vf_T2 = vars_T ./ vars_total;


%% Technical variability from multiple sequencing replicates

%Multiple technical replicates
data_Nr = data_abs(:,tech_inds);
means_Nr = mean(data_Nr,2);
vars_Nr = var(data_Nr')';

%Noise floor
abunds = abund_mat(1,:);
noise_floor = var(abunds(tech_inds)) / mean(abunds(tech_inds))^2;


%% Taylor's Law

inds = find(means > 0 & vars_total > 0);
inds_T = find(means > 0 & vars_T > 0);
inds_S = find(means > 0 & vars_S > 0);
inds_N = find(means > 0 & vars_N > 0);

%Taylor's law exponents
[beta,s] = polyfit(log10(means(inds)),log10(vars_total(inds)),1);
[beta_T,s_T] = polyfit(log10(means(inds_T)),log10(vars_T(inds_T)),1);
[beta_S,s_S] = polyfit(log10(means(inds_S)),log10(vars_S(inds_S)),1);
[beta_N,s_P] = polyfit(log10(means(inds_N)),log10(vars_N(inds_N)),1);


%% Variance decomposition of total bacterial densities %%%%%%%%%%%%%%%%%%%%

% Total absolute abundance densities
abunds = abund_mat(1,:);

%Include days 27 and 48 in this calculation
abunds_X = abunds([X_inds 84 87]);
abunds_Y = abunds([Y_inds 85 88]);
abunds_Z = abunds([Z_inds 83 86]);
abunds_all = [abunds_X; abunds_Y; abunds_Z];


%% Calculate marginal variance of total bacterial densities

%Perform 1,000 re-sampling iterations
for i = 1:1e3
    
    abunds_all_perm = [];
    
    %Randomly shuffle X, Y and Z
    for j = 1:length(abunds_X)
        ord = randperm(3);
        abunds_all_perm(:,j) = abunds_all(ord,j);
    end
    
    var_perm = mean(var(abunds_all_perm')');
    vars_perm(i) = var_perm;
    
end

%Average over 1,000 iterations
abunds_vars_total = mean(vars_perm);


%% Perform variance decomposition of total bacterial densities

abunds_covs_XZ = [];
abunds_covs_XmZY = [];
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
    
    %Temporal variance
    abunds_covmat_XZ = cov(abunds_X_perm,abunds_Z);
    abunds_cov_XZ = abunds_covmat_XZ(1,2);
    abunds_covs_XZ(i) = abunds_cov_XZ;
    
    %Spatial sampling variance
    abunds_covmat_XmZY = cov(abunds_X_perm-abunds_Z,abunds_Y_perm);
    abunds_cov_XmZY = abunds_covmat_XmZY(1,2);
    abunds_covs_XmZY(i) = abunds_cov_XmZY;
    
    %Technical variance
    abunds_var_XmY = .5 * var([abunds_X_perm - abunds_Y_perm])';
    abunds_vars_XmY(i) = abunds_var_XmY;
    
end

%% Average over re-sampling iterations

%Temporal variance
abunds_vars_T = mean(abunds_covs_XZ,2);
abunds_vars_T = max(abunds_vars_T,0);  %Ignore negative estimates of variability

%Spatial sampling variance
abunds_vars_S = mean(abunds_covs_XmZY,2);
abunds_vars_S = max(abunds_vars_S,0);  %Ignore negative estimates of variability

%Technical variance
abunds_vars_N = mean(abunds_vars_XmY,2);

%Fraction of each component to total variance
abunds_vf_T = abunds_vars_T / (abunds_vars_T + abunds_vars_S + abunds_vars_N);
abunds_vf_S = abunds_vars_S / (abunds_vars_T + abunds_vars_S + abunds_vars_N);
abunds_vf_N = abunds_vars_N / (abunds_vars_T + abunds_vars_S + abunds_vars_N);

%Standard deviation over re-sampling iterations
abunds_vf_T_std = std(abunds_covs_XZ / (abunds_vars_T + abunds_vars_S + abunds_vars_N));
abunds_vf_S_std = std(abunds_covs_XmZY / (abunds_vars_T + abunds_vars_S + abunds_vars_N));
abunds_vf_N_std = std(abunds_vars_XmY / (abunds_vars_T + abunds_vars_S + abunds_vars_N));

%Fraction of each component to total variance
abunds_vf_T2 = abunds_vars_T / (abunds_vars_total);
abunds_vf_S2 = abunds_vars_S / (abunds_vars_total);
abunds_vf_N2 = abunds_vars_N / (abunds_vars_total);




%% Covariance decomposition of OTU abundances %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Performing covariance decomposition...'])


%% Calculate total covariances
covs_total = [];

%Perform 1,000 re-sampling iterations
for i = 1:1e3

    data_X_perm = [];

    %Randomly choose a sample from of data_X, data_Y, data_Z
    for j = 1:Nx
        flip = rand;
        if flip > (2/3)
            data_X_perm(:,j) = data_X(:,j);
        elseif flip > (1/3)
            data_X_perm(:,j) = data_Y(:,j);
        else
            data_X_perm(:,j) = data_Z(:,j);
        end
    end
    
    %Mean center the data 
    data_X_mc = data_X_perm - repmat(mean(data_X_perm,2),1,Nx);
    
    %Calculate covariance matrix
    covmat_X_perm = 1/(Nx-1) * data_X_mc * data_X_mc';
    
    %Store covariance matrix for this iteration
    covs_total(:,:,i) = covmat_X_perm;  
    
end
covmat_total = mean(covs_total,3);
clear covs_total;


%% Temporal covariance decomposition
crosscovs_T = [];

%Perform 1,000 re-sampling iterations
for i = 1:1e3

    data_X_perm = [];
    data_Y_perm = [];
    data_Z_perm = data_Z;
    
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
    data_Z_mc = data_Z_perm - repmat(mean(data_Z_perm,2),1,Nz);
    
    %Calculate all pairwise covariances for this realization (2
    %permutations)
    covmat_XZ = 1/(Nx-1) * data_X_mc*data_Z_mc'; %(X_i,Z_j)
    covmat_ZX = 1/(Nx-1) * data_Z_mc*data_X_mc'; %(Z_i,X_j)
    covmat = 1/2 * (covmat_XZ + covmat_ZX);
    
    %Store covmat for each realization
    crosscovs_T(:,:,i) = covmat;    
    
end

%Average over all iterations
covmat_T = mean(crosscovs_T,3);
clear crosscovs_T;



%% Spatial sampling covariance decomposition
crosscovs_S = [];

%Peform 1,000 re-sampling iterations
for i = 1:1e3

    data_X_perm = [];
    data_Y_perm = [];
    data_Z_perm = data_Z;
    
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
    
    data_XmZ = data_X_perm - data_Z_perm;
    
    %Mean center matrices;
    data_XmZ_mc = data_XmZ - repmat(mean(data_XmZ,2),1,Nx);
    data_Y_mc = data_Y_perm - repmat(mean(data_Y_perm,2),1,Ny);
    
    %Calculate all pairwise covariances for this realization (2
    %permutations)
    covmat_XZY = 1/(Nx-1) * data_XmZ_mc*data_Y_mc';
    covmat_YXZ = 1/(Nx-1) * data_Y_mc*data_XmZ_mc';
    covmat = 1/2 * (covmat_XZY + covmat_YXZ);
    
    %Store covmat for each realization
    crosscovs_S(:,:,i) = covmat;  
    
end

%Average over all iterations
covmat_S = mean(crosscovs_S,3);
clear crosscovs_S;



 %% Technical covariance decomposition
crosscovs_N = [];

%Peform 1,000 re-sampling iterations
for i = 1:1e3

    data_X_perm = [];
    data_Y_perm = [];
    data_Z_perm = data_Z;
    
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
   
end

%Average over all iterations
covmat_N = mean(crosscovs_N,3);  
clear crosscovs_N;

disp(['Covariance decomposition complete!'])


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
cormat_T = covmat_T ./ sigxsigy;
cormat_S = covmat_S ./ sigxsigy;
cormat_N = covmat_N ./ sigxsigy;





    
clear covs_XZ covs_XmZY vars_XmY sigxsigy data_X_perm data_Y_perm data_Z_perm data_X_mc data_Y_mc data_Z_mc data_XmY_mc data_XmZ_mc data_XmY data_XmZ

%% Saving
 save([save_dir 'DIVERS_gut.mat']);
 
%% Write variance decomposition to output

%Variances
vars_table = table(otu_ids, means, [vars_T + vars_S + vars_N], vars_T, vars_S, vars_N, tax);
vars_table.Properties.VariableNames = {'OTU_ID','Average_abundance','Total_variances','Temporal_variances','Spatial_variances','Technical_variances','Taxonomy'};
writetable(vars_table,[save_dir 'DIVERS_gut.txt'],'Delimiter','\t');


 
 

