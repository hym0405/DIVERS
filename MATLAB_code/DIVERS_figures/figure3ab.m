%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script information

% This script generates Figures 3a,b of the main text

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load DIVERS gut data
load DIVERS_gut.mat;
c = @cmu.colors;

 
%% Find OTUs with high abundance;

%Cutoff based on  variance decomposition
high_inds = find(log10(means) > -4);

%Update means and variances
means = means(high_inds);
vars_total = vars_total(high_inds);
vars_T = vars_T(high_inds);
vars_S = vars_S(high_inds);
vars_N = vars_N(high_inds);
vf_T2 = vf_T2(high_inds);
vf_S2 = vf_S2(high_inds);
vf_N2 = vf_N2(high_inds);
vf_T = vf_T(high_inds);
vf_S = vf_S(high_inds);
vf_N = vf_N(high_inds);

%Update matrices
cormat_total = cormat_total(high_inds,high_inds);
cormat_T = cormat_T(high_inds,high_inds);
cormat_S = cormat_S(high_inds,high_inds);
cormat_N = cormat_N(high_inds,high_inds);

%Update taxa
tax = tax(high_inds);
ptax = ptax(high_inds);
ctax = ctax(high_inds);
otax = otax(high_inds);
ftax = ftax(high_inds);
gtax = gtax(high_inds);
otu_ids = otu_ids(high_inds);

%Update data
data_X = data_X(high_inds,:);
data_Y = data_Y(high_inds,:);
data_Z = data_Z(high_inds,:);



%% Covariance decomposition by phyla


%% Get phyla inds
bacter_inds = find(strcmp('Bacteroidetes',ptax));
proteo_inds = find(strcmp('Proteobacteria',ptax));
actino_inds = find(strcmp('Actinobacteria',ptax));
firmi_inds = find(strcmp('Firmicutes',ptax));

phyla_inds{1} = actino_inds;
phyla_inds{2} = bacter_inds;
phyla_inds{3} = firmi_inds;
phyla_inds{4} = proteo_inds;


%% Phyla variances
vars_total_actino = vars_total(actino_inds);
vars_total_bacter = vars_total(bacter_inds);
vars_total_firmi = vars_total(firmi_inds);
vars_total_proteo = vars_total(proteo_inds);


%% Get correlations within phyla

%Actinobacteria
cormat_total_actino = cormat_total(actino_inds,actino_inds);
cormat_T_actino = cormat_T(actino_inds,actino_inds);
cormat_S_actino = cormat_S(actino_inds,actino_inds);
cormat_N_actino = cormat_N(actino_inds,actino_inds);
cormat_total_actino_triu = triu(cormat_total_actino,1);
cormat_actino_inds = find(cormat_total_actino_triu);

%Bacteroidetes
cormat_total_bacter = cormat_total(bacter_inds,bacter_inds);
cormat_T_bacter = cormat_T(bacter_inds,bacter_inds);
cormat_S_bacter = cormat_S(bacter_inds,bacter_inds);
cormat_N_bacter = cormat_N(bacter_inds,bacter_inds);
cormat_total_bacter_triu = triu(cormat_total_bacter,1);
cormat_bacter_inds = find(cormat_total_bacter_triu);

%Firmicutes
cormat_total_firmi = cormat_total(firmi_inds,firmi_inds);
cormat_T_firmi = cormat_T(firmi_inds,firmi_inds);
cormat_S_firmi = cormat_S(firmi_inds,firmi_inds);
cormat_N_firmi = cormat_N(firmi_inds,firmi_inds);
cormat_total_firmi_triu = triu(cormat_total_firmi,1);
cormat_firmi_inds = find(cormat_total_firmi_triu);

%Proteobacteria
cormat_total_proteo = cormat_total(proteo_inds,proteo_inds);
cormat_T_proteo = cormat_T(proteo_inds,proteo_inds);
cormat_S_proteo = cormat_S(proteo_inds,proteo_inds);
cormat_N_proteo = cormat_N(proteo_inds,proteo_inds);
cormat_total_proteo_triu = triu(cormat_total_proteo,1);
cormat_proteo_inds = find(cormat_total_proteo_triu);


%% Get correlations within and between phyla

cormat_phyla_total = [];
cormat_phyla_T = [];
cormat_phyla_S = [];
cormat_phyla_N = [];

for i = 1:size(phyla_inds,2);
    for j = 1:i
        
        %Get phyla indices
        phyla_inds_1 = phyla_inds{i};
        phyla_inds_2 = phyla_inds{j};
        
        %If i and j belong to same phyla
        if i == j
            
            %Get the intersection of indices in correlation matrix
            cross_cormat_total = cormat_total(phyla_inds_1,phyla_inds_2);
            cross_cormat_T = cormat_T(phyla_inds_1,phyla_inds_2);
            cross_cormat_S = cormat_S(phyla_inds_1,phyla_inds_2);
            cross_cormat_N = cormat_N(phyla_inds_1,phyla_inds_2);
            
            %Get upper triangle inds
            triu_inds = find(triu(cross_cormat_total,1));
            
            %Get average cross correlation and store
            cormat_phyla_total(i,j) = mean(cross_cormat_total(triu_inds));
            cormat_phyla_total(j,i) = mean(cross_cormat_total(triu_inds));
            cormat_phyla_T(i,j) = mean(cross_cormat_T(triu_inds));
            cormat_phyla_T(j,i) = mean(cross_cormat_T(triu_inds));
            cormat_phyla_S(i,j) = mean(cross_cormat_S(triu_inds));
            cormat_phyla_S(j,i) = mean(cross_cormat_S(triu_inds));
            cormat_phyla_N(i,j) = mean(cross_cormat_N(triu_inds));
            cormat_phyla_N(j,i) = mean(cross_cormat_N(triu_inds));
            
        else
            
            %Get the intersection of indices in correlation matrix
            cross_cormat_total = cormat_total(phyla_inds_1,phyla_inds_2);
            cross_cormat_T = cormat_T(phyla_inds_1,phyla_inds_2);
            cross_cormat_S = cormat_S(phyla_inds_1,phyla_inds_2);
            cross_cormat_N = cormat_N(phyla_inds_1,phyla_inds_2);
            
            %Get average cross correlation and store
            cormat_phyla_total(i,j) = mean(cross_cormat_total(:));
            cormat_phyla_total(j,i) = mean(cross_cormat_total(:));
            cormat_phyla_T(i,j) = mean(cross_cormat_T(:));
            cormat_phyla_T(j,i) = mean(cross_cormat_T(:));
            cormat_phyla_S(i,j) = mean(cross_cormat_S(:));
            cormat_phyla_S(j,i) = mean(cross_cormat_S(:));
            cormat_phyla_N(i,j) = mean(cross_cormat_N(:));
            cormat_phyla_N(j,i) = mean(cross_cormat_N(:));

        end
    end
end



%%  Plotting %%%%%%%%%%%%%%%%%%%%


%% Total covariance decomposition
f1 = figure;

%% Boxplot of pairwise correlations (total, temporal, spsatial, technical)
subplot(2,2,1);

cormat_total_triu = triu(cormat_total,1);
cormat_total_inds = find(cormat_total_triu); %Get inds of upper triangle
all_cors = [cormat_total(cormat_total_inds) cormat_T(cormat_total_inds) cormat_S(cormat_total_inds) cormat_N(cormat_total_inds)];

%Positive total correlations
positions = [1 1.25 1.5 1.75];
b1a = boxplot(all_cors,'Notch','off','BoxStyle','outline','PlotStyle','traditional',...
    'Widths',.09,'Symbol','k.','OutlierSize',2,'Whisker',3,'positions',positions,...
    'Color',c('black'));
set(b1a,'LineWidth',1);
h = findobj(gca,'Tag','Box');
patch(get(h(4),'XData'),get(h(4),'YData'),c('gray'),'FaceAlpha',1,'EdgeColor',c('gray'));
patch(get(h(3),'XData'),get(h(3),'YData'),c('red'),'FaceAlpha',1,'EdgeColor',c('red'));
patch(get(h(2),'XData'),get(h(2),'YData'),c('deep sky blue'),'FaceAlpha',1,'EdgeColor',c('deep sky blue'));
patch(get(h(1),'XData'),get(h(1),'YData'),c('dark violet'),'FaceAlpha',1,'EdgeColor',c('dark violet'));
hold on
p1 = plot(0:.1:2,zeros(size(0:.1:2)),'--');
p1.LineWidth = 0.5;
p1.Color = c('pastel gray');
hold on
b1a = boxplot(all_cors,'Notch','off','BoxStyle','outline','PlotStyle','traditional',...
    'Widths',.09,'Symbol','k.','OutlierSize',2,'Whisker',3,'positions',positions,...
    'Color',c('black'));
set(b1a,'LineWidth',1);
hold off

set(gca,'xlim',[.85,1.85])
set(gca,'ylim',[-1,1.1]);
set(gca,'XTick',positions);
set(gca,'XTickLabel',{'Total','Temporal','Spatial sampling','Technical'});
set(gca,'XTickLabelRotation',45);
set(gca,'YTick',[-1:.25:1]);
set(gca,'YTickLabel',{'-1','','-0.5','','0','','0.5','','1'});
ylabel('Pairwise correlation');
set(gca,'FontSize',11);
set(gca,'LineWidth',1);
set(gca,'FontName','Arial');
box off




%% Average correlations within and between phyla
f2 = figure;

%Total correlations
subplot(2,2,1);

cmap1 = colorGradient(c('deep sky blue'),c('white'),500);
cmap2 = colorGradient(c('white'),c('tangerine'),500);
cmap = [cmap1;cmap2];
colormap(cmap);

p1 = imagesc(cormat_phyla_total);
p1.XData = 1:4;
p1.YData = 1:4;
set(gca, 'XTick', 1:4)
set(gca,'XTickLabel',{'Actinobacteria','Bacteroidetes','Firmicutes','Proteobacteria'});
set(gca,'XTickLabelRotation',45);
set(gca,'YTick', 1:4);
set(gca,'YTickLabel',{'Actinobacteria','Bacteroidetes','Firmicutes','Proteobacteria'});
title('Total correlations','FontWeight','Normal','FontSize',12);
set(gca,'FontSize',10);
set(gca,'LineWidth',1);
c1 = colorbar;
c1.LineWidth = 1;
caxis([-1,1]);
c1.Limits = [-1 1];

%Temporal correlations
subplot(2,2,2);

p2 = imagesc(cormat_phyla_T);
p2.XData = 1:4;
p2.YData = 1:4;
set(gca, 'XTick', 1:4)
set(gca,'XTickLabel',{'Actinobacteria','Bacteroidetes','Firmicutes','Proteobacteria'});
set(gca,'XTickLabelRotation',45);
set(gca,'YTick', 1:4);
set(gca,'YTickLabel',{'Actinobacteria','Bacteroidetes','Firmicutes','Proteobacteria'});
title('Temporal correlations','FontWeight','Normal','FontSize',12);
set(gca,'FontSize',10);
set(gca,'LineWidth',1);
c2 = colorbar;
c2.LineWidth = 1;
caxis([-1,1]);
c2.Limits = [-1 1];

%Spatial correlations
subplot(2,2,3);

p3 = imagesc(cormat_phyla_S);
p3.XData = 1:4;
p3.YData = 1:4;
set(gca, 'XTick', 1:4)
set(gca,'XTickLabel',{'Actinobacteria','Bacteroidetes','Firmicutes','Proteobacteria'});
set(gca,'XTickLabelRotation',45);
set(gca,'YTick', 1:4);
set(gca,'YTickLabel',{'Actinobacteria','Bacteroidetes','Firmicutes','Proteobacteria'});
title('Spatial correlations','FontWeight','Normal','FontSize',12);
set(gca,'FontSize',10);
set(gca,'LineWidth',1);
c3 = colorbar;
c3.LineWidth = 1;
caxis([-1,1]);
c3.Limits = [-1 1];

%Technical correlations
subplot(2,2,4);

p4 = imagesc(cormat_phyla_N);
p4.XData = 1:4;
p4.YData = 1:4;
set(gca, 'XTick', 1:4)
set(gca,'XTickLabel',{'Actinobacteria','Bacteroidetes','Firmicutes','Proteobacteria'});
set(gca,'XTickLabelRotation',45);
set(gca,'YTick', 1:4);
set(gca,'YTickLabel',{'Actinobacteria','Bacteroidetes','Firmicutes','Proteobacteria'});
title('Technical correlations','FontWeight','Normal','FontSize',12);
set(gca,'FontSize',10);
set(gca,'LineWidth',1);
c4 = colorbar;
c4.LineWidth = 1;
caxis([-1,1]);
c4.Limits = [-1 1];


%% Tests

%Wilcoxan Rank Sum of one phyla versus all others combined
cormat_total_triu = triu(cormat_total,1);
cormat_total_inds = find(cormat_total_triu); %Get inds of upper triangle
cors_total = cormat_total(cormat_total_inds);
cors_T = cormat_T(cormat_total_inds);
cors_S = cormat_S(cormat_total_inds);


%Bacteroidetes
bacter_cors_total = cormat_total_bacter(cormat_bacter_inds);
bacter_cors_T = cormat_T_bacter(cormat_bacter_inds);
bacter_cors_S = cormat_S_bacter(cormat_bacter_inds);

%Actinobacteria
actino_cors_total = cormat_total_actino(cormat_actino_inds);
actino_cors_T = cormat_T_actino(cormat_actino_inds);
actino_cors_S = cormat_S_actino(cormat_actino_inds);

%Firmicutes
firmi_cors_total = cormat_total_firmi(cormat_firmi_inds);
firmi_cors_T = cormat_T_firmi(cormat_firmi_inds);
firmi_cors_S = cormat_S_firmi(cormat_firmi_inds);

%Proteobacteria
proteo_cors_total = cormat_total_proteo(cormat_proteo_inds);
proteo_cors_T = cormat_T_proteo(cormat_proteo_inds);
proteo_cors_S = cormat_S_proteo(cormat_proteo_inds);

%Within phyla correlations versus whole community
[p, h] = ranksum(actino_cors_total, cors_total);
[p, h] = ranksum(actino_cors_T, cors_T);
[p, h] = ranksum(actino_cors_S, cors_S);

[p, h] = ranksum(bacter_cors_total,  cors_total);
[p, h] = ranksum(bacter_cors_T, cors_T);
[p, h] = ranksum(bacter_cors_S,  cors_S);

[p, h] = ranksum(firmi_cors_total, cors_total);
[p, h] = ranksum(firmi_cors_T, cors_T);
[p, h] = ranksum(firmi_cors_S, cors_S);

[p, h] = ranksum(proteo_cors_total, cors_total);
[p, h] = ranksum(proteo_cors_T,cors_T);
[p, h] = ranksum(proteo_cors_S, cors_S);







