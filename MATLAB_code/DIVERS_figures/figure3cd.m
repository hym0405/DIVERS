%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script information

% This script generates Figure 3c,d of the main text

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

covmat_total = covmat_total(high_inds,high_inds);
covmat_T = covmat_T(high_inds,high_inds);
covmat_S = covmat_S(high_inds,high_inds);
covmat_N = covmat_N(high_inds,high_inds);

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
data_XY = .5 * (data_X + data_Y);


%% Examples of positive and negative temporal and spatial correlations

%Positive temporal (inds 45 and 54)
%Negative temporal (inds 4 and 11)
%Positive spatial (inds 12 and 29)
%Negative spatial (inds 33 and 52)


%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1= figure;

triu_inds = find(triu(cormat_total,1));

%% Plot temporal versus spatial correlations for all pairs of OTUs
subplot(2,2,1);

cc = colorGradient(c('light gray'),c('black'),100);
cc(end,:) = [1 1 1];
colormap(cc);

%All pairs
p1 = scatter_kde(cormat_T(triu_inds),cormat_S(triu_inds),'filled','o');
p1.MarkerEdgeColor = c('black');
p1.SizeData = 25;
p1.MarkerFaceAlpha = .1;
p1.MarkerEdgeAlpha = .01;
hold on

%Positive temporal correlation
p1a = scatter(cormat_T(45,54),cormat_S(45,54),'o');
p1a.MarkerEdgeColor = c('tangerine');
p1a.MarkerFaceColor = c('tangerine');
p1a.SizeData = 40;
p1a.MarkerFaceAlpha = 1;
p1a.MarkerEdgeAlpha = 1;
hold on

%Negative temporal correlation
p1b = scatter(cormat_T(4,11),cormat_S(4,11),'o');
p1b.MarkerEdgeColor = c('deep sky blue');
p1b.MarkerFaceColor = c('deep sky blue');
p1b.SizeData = 40;
p1b.MarkerFaceAlpha = 1;
p1b.MarkerEdgeAlpha = 1;
hold on

%Positive spatial correlation
p1c = scatter(cormat_T(12,29),cormat_S(12,29),'o');
p1c.MarkerEdgeColor = c('neon red');
p1c.MarkerFaceColor = c('neon red');
p1c.SizeData = 40;
p1c.MarkerFaceAlpha = 1;
p1c.MarkerEdgeAlpha = 1;
hold on

set(gca,'xlim',[-1,1.3]);
set(gca,'XTick',[-1,-0.5,0,0.5,1]);
set(gca,'ylim',[-1,1.2])
set(gca,'YTick',[-1,-0.5,0,0.5,1]);
xlabel('Temporal pairwise correlation');
ylabel('Spatial pairwise correlation');
set(gca,'LineWidth',1);
set(gca,'FontSize',10);


%% Plot time series of examples
f2 = figure;

%Positive temporal
subplot(2,2,1);

p1a = plot(data_XY(45,:),':k');
p1a.LineWidth = 2.1;
p1a.Color = c('deep sky blue');
hold on
p1b = plot(data_Z(45,:),'-k');
p1b.LineWidth = 1.5;
p1b.Color = c('deep sky blue');
hold on
set(gca,'ylim',[0,.015]);
set(gca,'YTick',[0,.005,.01,.015]);
ylabel('Abundance');
set(gca,'ycolor','k')

yyaxis right
p1c = plot(data_XY(54,:),':k');
p1c.LineWidth = 2.1;
p1c.Color = c('neon red');
hold on
p1d = plot(data_Z(54,:),'-k');
p1d.LineWidth = 1.5;
p1d.Color = c('neon red');
hold off
set(gca,'ylim',[0,.012]);
set(gca,'YTick',[0,.004,.008,.012]);
ylabel('Abundance');
set(gca,'ycolor','k')

xlabel('Day');
set(gca,'xlim',[0,21]);
set(gca,'XTick',[0,5,10,15,20]);
set(gca,'XTickLabel',{'0','5','10','15','20'});
set(gca,'FontSize',10);
set(gca,'LineWidth',1);
box off
ll1 = legend([p1b,p1d],{['OTU ' num2str(otu_ids(45))],['OTU ' num2str(otu_ids(54))]});
ll1.FontSize = 10;

%Negative temporal
subplot(2,2,2);

p2a = plot(data_XY(4,:),':k');
p2a.LineWidth = 2.1;
p2a.Color = c('deep sky blue');
hold on
p2b = plot(data_Z(4,:),'-k');
p2b.LineWidth = 1.5;
p2b.Color = c('deep sky blue');
hold off
set(gca,'ylim',[0,.12]);
set(gca,'YTick',[0,.04,.08,.12]);
ylabel('Abundance');
set(gca,'ycolor','k')

yyaxis right
p2c = plot(data_XY(11,:),':k');
p2c.LineWidth = 2.1;
p2c.Color = c('neon red');
hold on
p2d = plot(data_Z(11,:),'-k');
p2d.LineWidth = 1.5;
p2d.Color = c('neon red');
hold on
set(gca,'ylim',[0,.24]);
set(gca,'YTick',[0,.08,.16,.24]);
ylabel('Abundance');
set(gca,'ycolor','k')

set(gca,'xlim',[0,21]);
set(gca,'XTick',[0,5,10,15,20]);
set(gca,'XTickLabel',{'0','5','10','15','20'});
xlabel('Day');
set(gca,'FontSize',10);
set(gca,'LineWidth',1);
box off
ll2 = legend([p2b,p2d],{['OTU ' num2str(otu_ids(4))],['OTU ' num2str(otu_ids(11))]});
ll2.FontSize = 10;

%Positive spatial
subplot(2,2,3);

p3a = plot(data_XY(12,:),':k');
p3a.LineWidth = 2.1;
p3a.Color = c('deep sky blue');
hold on
p3b = plot(data_Z(12,:),'-k');
p3b.LineWidth = 1.5;
p3b.Color = c('deep sky blue');
hold on
set(gca,'ylim',[-.01,.18]);
set(gca,'YTick',[0,.06,.12,.18]);
ylabel('Abundance');
set(gca,'ycolor','k')

yyaxis right
p3c = plot(data_XY(29,:),':k');
p3c.LineWidth = 2.1;
p3c.Color = c('neon red');
hold on
p3d = plot(data_Z(29,:),'-k');
p3d.LineWidth = 1.5;
p3d.Color = c('neon red');
hold off
set(gca,'ylim',[-.001,.03]);
set(gca,'YTick',[0,.01,.02,.03]);
ylabel('Abundance');
set(gca,'ycolor','k')

set(gca,'xlim',[0,21]);
set(gca,'XTick',[0,5,10,15,20]);
set(gca,'XTickLabel',{'0','5','10','15','20'});
xlabel('Day');
set(gca,'FontSize',10);
set(gca,'LineWidth',1);
box off
ll3 = legend([p3b,p3d],{['OTU ' num2str(otu_ids(12))],['OTU ' num2str(otu_ids(29))]});
ll3.FontSize = 10;

