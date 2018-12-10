%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script information

% This script generates Figure 2 of the main text

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load DIVERS gut data
load DIVERS_gut.mat;
c = @cmu.colors;

 
%% Get phyla inds
bacter_inds = find(strcmp('Bacteroidetes',ptax));
actino_inds = find(strcmp('Actinobacteria',ptax));
firmi_inds = find(strcmp('Firmicutes',ptax));
proteo_inds = find(strcmp('Proteobacteria',ptax));


%% Find abundant, low technical noise strains
high_inds = find(log10(means) > -4); 

%Phyla
bacter_high_inds = intersect(bacter_inds,high_inds);
actino_high_inds = intersect(actino_inds,high_inds);
proteo_high_inds = intersect(proteo_inds,high_inds);
firmi_high_inds = intersect(firmi_inds,high_inds);


%% Find OTUs with high temporal and spatial variance
hsv_inds = intersect(find(vf_S > .6),high_inds);
htv_inds = intersect(find(vf_T > .8),high_inds);

%Spatial OTUs for plotting
hsv1 = find(hsv_inds == 12); %OTU 13
hsv2 = find(hsv_inds == 121); %OTU 48

%Temporal OTUs for plotting
htv1 = find(htv_inds == 11); %OTU 12
htv2 = find(htv_inds == 24); %OTU 25





%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Taylor's law (temporal variance)
f1 = figure;

subplot(2,2,1);

x = -7:.1:0;
p1a = plot(x,polyval(beta_T,x,'-'));
p1a.LineWidth = 1;
p1a.Color = c('red');
hold on
p1b = plot(log10(means), log10(vars_T),'.');
p1b.MarkerSize = 10;
p1b.Color = c('black');
hold off

set(gca,'xlim',[-7,-.5])
set(gca,'XTick',[-7:0]);
set(gca,'XTickLabel',{'10^{-7}','','10^{-5}','','10^{-3}','','10^{-1}',''});
set(gca,'ylim',[-14,-1])
set(gca,'YTick',[-13:2:-1]);
set(gca,'YTickLabel',{'10^{-13}','','10^{-9}','','10^{-5}','','10^{-1}'});
set(gca,'xscale','linear');
set(gca,'yscale','linear');
set(gca,'FontSize',12); 
set(gca,'LineWidth',1); 
set(gca,'FontName','Arial');
xlabel('Mean absolute abundance (a.u.)');
ylabel('Variance (temporal)');
box off

%% Taylor's law (spatial variance)
subplot(2,2,2);

x = -7:.1:0;
p2a = plot(x,polyval(beta_S,x,'-'));
p2a.LineWidth = 1;
p2a.Color = c('deep sky blue');
hold on
p2b = plot(log10(means), log10(vars_S),'.');
p2b.MarkerSize = 10;
p2b.Color = c('black');
hold off

set(gca,'xlim',[-7,-.5])
set(gca,'XTick',[-7:0]);
set(gca,'XTickLabel',{'10^{-7}','','10^{-5}','','10^{-3}','','10^{-1}',''});
set(gca,'ylim',[-14,-1])
set(gca,'YTick',[-13:2:-1]);
set(gca,'YTickLabel',{'10^{-13}','','10^{-9}','','10^{-5}','','10^{-1}'});
set(gca,'xscale','linear');
set(gca,'yscale','linear');
set(gca,'FontSize',12); 
set(gca,'LineWidth',1); 
set(gca,'FontName','Arial');
xlabel('Mean absolute abundance (a.u.)');
ylabel('Variance (spatial)');
box off




%% Temporal versus spatial variance fractions
f2 = figure;

%Variance fraction
subplot(2,2,1);

%Non-variable OTUs
nonv_inds = setdiff(high_inds,[htv_inds;hsv_inds]);

p3a = scatter(vf_T(nonv_inds),vf_S(nonv_inds),'o');
p3a.SizeData = 36;
p3a.MarkerEdgeColor = 'black';
p3a.MarkerEdgeAlpha = 1;
hold on
p3b = scatter(vf_T(htv_inds),vf_S(htv_inds),'o');
p3b.SizeData = 36;
p3b.MarkerEdgeColor = 'red';
p3b.MarkerEdgeAlpha = .3;
hold on
p3c = scatter(vf_T(htv_inds([htv1 htv2])),vf_S(htv_inds([htv1 htv2])),'o');
p3c.SizeData = 36;
p3c.MarkerEdgeColor = 'red';
p3c.MarkerFaceColor = 'red';
hold on
p3d = scatter(vf_T(hsv_inds),vf_S(hsv_inds),'o');
p3d.SizeData = 36;
p3d.MarkerEdgeColor = c('deep sky blue');
p3d.MarkerEdgeAlpha = .6;
hold on
p3e = scatter(vf_T(hsv_inds([hsv1 hsv2])),vf_S(hsv_inds([hsv1 hsv2])),'o');
p3e.SizeData = 36;
p3e.MarkerEdgeColor = c('deep sky blue');
p3e.MarkerFaceColor = c('deep sky blue');
hold off
set(gca,'xlim',[0,1])
set(gca,'XTick',[0:.2:1]);
set(gca,'YTick',[0:.2:1]);
set(gca,'xscale','linear');
set(gca,'yscale','linear');
set(gca,'FontSize',12); 
set(gca,'LineWidth',1); 
set(gca,'FontName','Arial');
xlabel('Variance fraction (temporal)');
ylabel('Variance fraction (spatial)');
box off



%% High temporal variance OTUs
f3 = figure;

%High temporal variance OTU 1
subplot(4,3,4);
p4x = plot(data_X(htv_inds(htv1),:),'-');
p4x.LineWidth = 1.5; 
p4x.Color = c('gray');
p4x.MarkerFaceColor = c('gray');
hold on
p4y = plot(data_Y(htv_inds(htv1),:),'-');
p4y.LineWidth = 1.5; 
p4y.Color = c('gray');
p4y.MarkerFaceColor = c('gray');
hold on
p4z = plot(data_Z(htv_inds(htv1),:),'-');
p4z.LineWidth = 1.5; 
p4z.Color = c('red');
p4z.MarkerFaceColor = c('red');
hold off
set(gca,'xlim',[0,21]);
set(gca,'XTick',[0:5:20]);
set(gca,'XTickLabel',{'0','5','10','15','20'});
set(gca,'ylim',[0,.22]);
set(gca,'YTick',[0:.05:.2]);
set(gca,'YTickLabel',{'0','','0.1','','0.2'});
set(gca,'FontSize',9);
set(gca,'LineWidth',1);
xlabel('Day');
ylabel('Abundance');
title(['OTU ' num2str(otu_ids(htv_inds(htv1))) ' (' gtax{htv_inds(htv1)} ')'],'FontSize',9,'FontWeight','Normal');
box off


%High temporal variance OTU 2
subplot(4,3,5);
p5x = plot(data_X(htv_inds(htv2),:),'-');
p5x.LineWidth = 1.5; 
p5x.Color = c('gray');
p5x.MarkerFaceColor = c('gray');
hold on
p5y = plot(data_Y(htv_inds(htv2),:),'-');
p5y.LineWidth = 1.5; 
p5y.Color = c('gray');
p5y.MarkerFaceColor = c('gray');
hold on
p5z = plot(data_Z(htv_inds(htv2),:),'-');
p5z.LineWidth = 1.5; 
p5z.Color = c('red');
p5z.MarkerFaceColor = c('red');
hold off
set(gca,'xlim',[0,21]);
set(gca,'XTick',[0:5:20]);
set(gca,'XTickLabel',{'0','5','10','15','20'});
set(gca,'ylim',[0,.045]);
set(gca,'YTick',[0:.01:.04]);
set(gca,'YTickLabel',{'0','','0.02','','0.04'});
xlabel('Day');
set(gca,'FontSize',9);
set(gca,'LineWidth',1);
%title(['OTU ' num2str(otu_ids(htv_inds(htv2))) ' (' gtax{htv_inds(htv2)} ')'],'FontSize',9,'FontWeight','Normal');
title(['OTU ' num2str(otu_ids(htv_inds(htv2))) ' (' 'Lachnospiracea incertae sedis' ')'],'FontSize',9,'FontWeight','Normal');
box off
[ll5,icons,plots,txt] = legend([p5x,p5z],{'X,Y','Z'});
ll5.FontSize = 10;
ll5.Location = 'NorthEast';
icons(3).XData = [0.35 0.6182];
icons(5).XData = [0.35 0.6182];
legend boxoff;


%% High spatial variance OTUs

%%High spatial variance OTU 1
subplot(4,3,1);
p6x = plot(data_X(hsv_inds(hsv1),:),'-');
p6x.LineWidth = 1.5; 
p6x.Color = c('gray');
p6x.MarkerFaceColor = c('gray');
hold on
p6y = plot(data_Y(hsv_inds(hsv1),:),'-');
p6y.LineWidth = 1.5; 
p6y.Color = c('gray');
p6y.MarkerFaceColor = c('gray');
hold on
p6z = plot(data_Z(hsv_inds(hsv1),:),'-');
p6z.LineWidth = 1.5; 
p6z.Color = c('deep sky blue');
p6z.MarkerFaceColor = c('deep sky blue');
hold off
set(gca,'xlim',[0,21]);
set(gca,'XTick',[0:5:20]);
set(gca,'XTickLabel',{'0','5','10','15','20'});
set(gca,'ylim',[0,.11]);
set(gca,'YTick',[0:.025:.1]);
set(gca,'YTickLabel',{'0','','0.05','','0.1'});
set(gca,'FontSize',9);
set(gca,'LineWidth',1);
xlabel('Day');
ylabel('Abundance');
title(['OTU ' num2str(otu_ids(hsv_inds(hsv1))) ' (' gtax{hsv_inds(hsv1)} ')'],'FontSize',9,'FontWeight','Normal');
box off

%High spatial variance OTU 2
subplot(4,3,2);
p7x = plot(data_X(hsv_inds(hsv2),:),'-');
p7x.LineWidth = 1.5; 
p7x.Color = c('gray');
p7x.MarkerFaceColor = c('gray');
hold on
p7y = plot(data_Y(hsv_inds(hsv2),:),'-');
p7y.LineWidth = 1.5; 
p7y.Color = c('gray');
p7y.MarkerFaceColor = c('gray');
hold on
p7z = plot(data_Z(hsv_inds(hsv2),:),'-');
p7z.LineWidth = 1.5; 
p7z.Color = c('deep sky blue');
p7z.MarkerFaceColor = c('deep sky blue');
hold off
set(gca,'xlim',[0,21]);
set(gca,'XTick',[0:5:20]);
set(gca,'XTickLabel',{'0','5','10','15','20'});
set(gca,'ylim',[0,2.2e-3]);
set(gca,'YTick',[0:.5e-3:2e-3]);
set(gca,'YTickLabel',{'0','','1e-3','','2e-3'});
xlabel('Day');
set(gca,'FontSize',9);
set(gca,'LineWidth',1);
title(['OTU ' num2str(otu_ids(hsv_inds(hsv2))) ' (' gtax{hsv_inds(hsv2)} ')'],'FontSize',9,'FontWeight','Normal');
box off
[ll2,icons,plots,txt] = legend([p7x,p7z],{'X,Y','Z'});
ll2.FontSize = 10;
icons(3).XData = [0.35 0.6182];
icons(5).XData = [0.35 0.6182];
legend boxoff;

