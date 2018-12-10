%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script information

% This script generates Figures 1b-d of the main text

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load DIVERS gut data
load DIVERS_gut.mat;
c = @cmu.colors;
 

%% Add day 27 (42 overall) and day 48 (63 overall)
time = [1:20 24 29];


%%  Plotting %%%%%%%%%%%%%%%%%%%%
f1 = figure;

% Plot absolute abundances as a function of time
subplot(2,3,1:2);

p1a = plot(time,abunds_X,'o');
p1a.Color = c('dark violet');
p1a.MarkerFaceColor = c('dark violet');
p1a.MarkerSize = 7;
p1a.LineWidth = .01;
hold on
p1b = plot(time,abunds_Y,'^');
p1b.Color = c('dark violet');
p1b.MarkerFaceColor = c('dark violet');
p1b.MarkerSize = 7;
p1b.LineWidth = .01;
hold on
p1c = plot(time,abunds_Z,'o');
p1c.Color = c('deep sky blue');
p1c.MarkerFaceColor = c('deep sky blue');
p1c.MarkerSize = 7;
p1c.LineWidth = .01;
hold on
p1d = plot(time,mean([0.5*([abunds_X + abunds_Y]); abunds_Z],1),'-');
p1d.LineWidth = 1.5;
p1d.Color = c('pastel gray');
hold off
set(gca,'xlim',[0,32]);
set(gca,'XTick',[0:2:20, 24 29]);
set(gca,'XTickLabel',{'0','','4','','8','','12','','16','','20','27','48'});
set(gca,'ylim',[0,2.5]);
set(gca,'YTick',[0,.5,1,1.5,2,2.5]);
set(gca,'FontSize',12);
set(gca,'LineWidth',1); 
set(gca,'FontName','Arial');
xlabel('Day');
ylabel('Total bacterial density (a.u.)');
box off
ll1 = legend([p1a,p1b,p1c],{'X','Y','Z'},'orientation','vertical');
ll1.FontSize = 12;
ll1.Location = 'SouthEast';

% Total absolute density decomposition
subplot(2,3,3);
e1 = errorbar([1,2,3],[abunds_vf_N abunds_vf_S abunds_vf_T],[abunds_vf_N_std abunds_vf_S_std abunds_vf_T_std],'.');
e1.LineWidth = 2;
e1.Color = c('black');
hold on
b1a = bar(1,abunds_vf_N);
b1a.FaceColor = c('dark violet');
b1a.LineWidth = 1;
b1a.EdgeColor = c('black');
b1a.FaceAlpha = 1;
b1a.BarWidth = .6;
hold on
b1b = bar(2,abunds_vf_S);
b1b.FaceColor = c('deep sky blue');
b1b.LineWidth = 1;
b1b.EdgeColor = c('black');
b1b.FaceAlpha = 1;
b1b.BarWidth = .6;
hold on
b1c = bar(3,abunds_vf_T);
b1c.FaceColor = c('red');
b1c.LineWidth = 1;
b1c.EdgeColor = c('black');
b1c.FaceAlpha = 1;
b1c.BarWidth = .6;
hold off
set(gca,'XTick',[1,2,3]);
set(gca,'XTickLabel',{'N','S','T'});
set(gca,'xlim',[.5,3.5]);
set(gca,'YTick',[0:.25:1]);
set(gca,'YTickLabel',{'0','','0.5','','1'});
set(gca,'ylim',[0,1]);
set(gca,'FontSize',12)
set(gca,'LineWidth',1);
ylabel('Variance fraction');
box off


%% Plot variance fraction decomposition as a function of average abundance
f2 = figure; 

subplot(2,2,1);

%Bin OTUs according to average OTU abundance
b_w = 1;
bins = -7:b_w:-2;
means_vf_T = [];
means_vf_S = [];
means_vf_N = [];
sem_vf_T = [];
sem_vf_S = [];
sem_vf_N = [];
for i = 1:length(bins)
    if i == length(bins)
        inds = find(log10(means) > bins(i));
    else
        inds = find(log10(means) > bins(i) & log10(means) < bins(i+1)); 
    end
    means_vf_T(i) = nanmean(vf_T(inds));
    means_vf_S(i) = nanmean(vf_S(inds));
    means_vf_N(i) = nanmean(vf_N(inds));
    sem_vf_T(i) = nanstd(vf_T(inds))/sqrt(length(inds));
    sem_vf_S(i) = nanstd(vf_S(inds))/sqrt(length(inds));
    sem_vf_N(i) = nanstd(vf_N(inds))/sqrt(length(inds));    
end

%Create bar graphs
vf_bar = [means_vf_N; means_vf_S; means_vf_T];

b2 = bar(bins+b_w/2,vf_bar','stacked');
b2(1).FaceColor = c('matlab purple');
b2(1).EdgeColor = c('black');
b2(1).LineWidth = 1;
b2(2).FaceColor = c('deep sky blue');
b2(2).EdgeColor = c('black');
b2(2).LineWidth = 1;
b2(3).FaceColor = c('red');
b2(3).EdgeColor = c('black');
b2(3).LineWidth = 1;
b2(3).BarWidth = .5;
hold on

%Create error bars
vf_error_N = means_vf_N;
vf_error_S = means_vf_N + means_vf_S;
vf_error_T = means_vf_N + means_vf_S + means_vf_T;

e2a = errorbar(bins+b_w/2,vf_error_N,sem_vf_N,'.');
e2a.Color = c('black');
e2a.MarkerSize = .1;
e2a.LineWidth = 1;
hold on
e2b = errorbar(bins+b_w/2,vf_error_S,sem_vf_S,'.');
e2b.Color = c('black');
e2b.MarkerSize = .1;
e2b.LineWidth = 1;
hold on
e2c = errorbar(bins+b_w/2,vf_error_T,sem_vf_T,'.');
e2c.Color = c('black');
e2c.MarkerSize = .1;
e2c.LineWidth = 1;

hold off

xlabel('Mean absolute abundance (a.u.)');
ylabel('Variance fraction');
set(gca,'xlim',[-7.2,-.8])
set(gca,'XTick',[-7:-1]);
set(gca,'XTickLabel',{'10^{-7}','','10^{-5}','','10^{-3}','','10^{-1}'});
set(gca,'ylim',[0,1.1])
set(gca,'YTick',[0,.25,.5,.75,1]);
set(gca,'YTick',[0,.2,.4,.6,.8,1]);
set(gca,'FontSize',12);
set(gca,'LineWidth',1); 
set(gca,'FontName','Arial');
box off
ll2 = legend([b2(1),b2(2),b2(3)],'Technical','Spatial sampling','Temporal');
ll2.FontSize = 12;
