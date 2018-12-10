%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Script information

%   This script reads in the fecal OTU time series data
%   It then reformats this OTU data for downstream analysis in MATLAB.
%   The data is saved in the MAT file gut_table.mat.
%   This MAT file is the input for the script DIVERS_gut.m, which
%   performs the DIVERS decomposition 

% User input required here:

%Directory containing DIVERS files
file_dir = ['/Path/To/.../DIVERS/DIVERS_files/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Read in sample meta data
MD = readtable([file_dir 'Gut_Metadata.xlsx']);

%Rows 1-88 are samples we are interested in (others are control wells)
md_weights = table2array(MD(1:88,10)); %sample weights in mg
md_names = table2array(MD(1:88,16)); 
md_days = table2array(MD(1:88,5));


%% Rearrange sample orders

%Sample names
target_names = cat(1,md_names(1:60),md_names(67:88));
target_names = cat(1,target_names,md_names(61:66));

%Sample weights
target_ws = cat(1,md_weights(1:60),md_weights(67:88));
target_ws = cat(1,target_ws,md_weights(61:66));
target_weights = [];
for i = 1:length(target_ws)
    target_weights(i) = str2double(target_ws(i));
end

%Time series days
target_days = cat(1,md_days(1:60),md_days(67:88));
target_days = cat(1,target_days,md_days(61:66));

%These will be final sample indices - OTU table will be arranged like this
time_inds = 1:60;
space_inds = [40 41 61:72];
tech_inds = [41 42 73:82];


%% Read in OTU table
T = readtable([file_dir 'gut_table.txt'],'Delimiter','\t');
samples = T.Properties.VariableNames; samples = samples(2:end-1);
otu_ids = table2array(T(:,1));
tax = table2array(T(:,end));
data = table2array(T(:,2:end-1));
[M,N] = size(data);


%% Arrange OTU table in preferred order

% Find mapping between meta data and OTU samples
target_inds = [];
for i = 1:length(target_names)
   target_inds = [target_inds; find(strcmp(target_names(i),samples))]; 
end

%Rearrange OTU table in preferred order
data = data(:,target_inds);


%% Indices of spatial and technical replicates

%Indices of technical replicates in OTU table
X_inds = time_inds(2:3:60);
Y_inds = time_inds(3:3:60);

%Indices of second spatial replicate in OTU table
Z_inds = time_inds(1:3:60);


%% Calculate total bacterial densities in each sample

%Relative abundance of spike-in strain in each sample (OTU 1)
spike_otu_abunds = data(1,:) ./ sum(data,1);

%Calculate total bacterial density per sample (up to scaling constant)
abs_abunds = (1-spike_otu_abunds) ./ (spike_otu_abunds .* target_weights);

%Renormalize total bacterial densities to mean of 1
abs_abunds_norm = abs_abunds ./ mean(abs_abunds);
abund_mat = repmat(abs_abunds_norm,M,1);


%% Rename stuff before saving

samples = target_names;
weights = target_weights;
days = target_days;


%% Parse OTU taxonomies
ptax = []; %phylum
ctax = []; %class
otax = []; %order
ftax = []; %family
gtax = []; %genus

for i = 1:M
    taxrank = tax{i};
    quote_inds = strfind(taxrank,'"');
    if length(quote_inds) > 0
        taxrank(quote_inds) = '';
        tax{i} = taxrank;
    end
end

for i = 1:M
    taxrank = strsplit(tax{i},'; ');
    
    %Contains subclass and suborder
    if length(taxrank) == 18
        ptax{i} = taxrank{5};
        ctax{i} = taxrank{7};
        otax{i} = taxrank{11};
        ftax{i} = taxrank{15};
        gtax{i} = taxrank{17};
        
        %Contains subclass
    elseif length(taxrank) == 16
        ptax{i} = taxrank{5};
        ctax{i} = taxrank{7};
        otax{i} = taxrank{11};
        ftax{i} = taxrank{13};
        gtax{i} = taxrank{15};
        
    elseif length(taxrank) == 14
        ptax{i} = taxrank{5};
        ctax{i} = taxrank{7};
        otax{i} = taxrank{9};
        ftax{i} = taxrank{11};
        gtax{i} = taxrank{13};
    else
        ptax{i} = '';
        ctax{i} = '';
        otax{i} = '';
        ftax{i} = '';
        gtax{i} = '';
    end
end

ptax = ptax';
ctax = ctax';
otax = otax';
ftax = ftax';
gtax = gtax';



%% Taxonomies with confidences
ptax_conf = []; %phlum
ctax_conf = []; %class
otax_conf = []; %order
ftax_conf = []; %family
gtax_conf = []; %genus

for i = 1:M
    taxrank = strsplit(tax{i},'; ');
    
    %Contains subclass and suborder
    if length(taxrank) == 18
        ptax_conf{i} = [taxrank{5} ' (' taxrank{6} ')'];
        ctax_conf{i} = [taxrank{7} ' (' taxrank{8} ')'];
        otax_conf{i} = [taxrank{11} ' (' taxrank{12} ')'];
        ftax_conf{i} = [taxrank{15} ' (' taxrank{16} ')'];
        gtax_conf{i} = [taxrank{17} ' (' taxrank{18} ')'];
        
        %Contains subclass
    elseif length(taxrank) == 16
        ptax_conf{i} = [taxrank{5} ' (' taxrank{6} ')'];
        ctax_conf{i} = [taxrank{7} ' (' taxrank{8} ')'];
        otax_conf{i} = [taxrank{11} ' (' taxrank{12} ')'];
        ftax_conf{i} = [taxrank{13} ' (' taxrank{14} ')'];
        gtax_conf{i} = [taxrank{15} ' (' taxrank{16} ')'];
        
    elseif length(taxrank) == 14
        ptax_conf{i} = [taxrank{5} ' (' taxrank{6} ')'];
        ctax_conf{i} = [taxrank{7} ' (' taxrank{8} ')'];
        otax_conf{i} = [taxrank{9} ' (' taxrank{10} ')'];
        ftax_conf{i} = [taxrank{11} ' (' taxrank{12} ')'];
        gtax_conf{i} = [taxrank{13} ' (' taxrank{14} ')'];
        
    else
        ptax_conf{i} = '';
        ctax_conf{i} = '';
        otax_conf{i} = '';
        ftax_conf{i} = '';
        gtax_conf{i} = '';
    end
end
ptax_conf = ptax_conf';
ctax_conf = ctax_conf';
otax_conf = otax_conf';
ftax_conf = ftax_conf';
gtax_conf = gtax_conf';




%% Save data
save([file_dir 'gut_table.mat'],'data','abund_mat','samples','weights','days','X_inds','Y_inds','Z_inds','time_inds','space_inds','tech_inds','tax','gtax','ftax','otax','ctax','ptax','otu_ids','ptax_conf','ctax_conf','otax_conf','ftax_conf','gtax_conf');
