%% sceptre_analysis3

%% Notes:
% This analysis script uses violinplot obtained from:
% https://github.com/bastibe/Violinplot-Matlab
% This analysis script is based on the example of only having a single
% immunofluorescence channel and three FISH channels.

%% Font parameters for figure

la_fs = 20; %Font for figure labels
ax_fs = 20; %Font for figures axes
ax_lw = 3; %line width for axes 
sz = 30; %size of markers for plots

%%  figure and analysis parameters
chck = 0; % leave > 0 if random boxes have already been generated
hmark1mM = [0 4500]; % IF channel 1 axis limit for violinplots and scatters
bn = 300; %number of boxes to be randomly generated throughout each nucleus

list = 1:length(r);
gene1 = '{\it MYL6}';       %name of targeted gene in FISH channel 1
gene2 = '{\it HOXC}';       %name of targeted gene in FISH channel 1
gene3 = '{\it LINC-PINT}';  %name of targeted gene in FISH channel 1
hmark1 = 'K4me3';           %name of Immunofluorescence channel 1

% Order of channels in violinplot 
h1o = 1;    
rd1o = 2;
f1o = 3;
f2o = 4;
f3o = 5;

%% FISH cluster selection parameters

% Below is the conditions used for selecting FISH clusters from segmented
% clusters in FISH channels. For the most part, only minimum volume is
% necessary for removing dim and non-specific clusters, but other 
% parameters such as maxfluorescence may remove clusters that are 
% fluorescent contaminants.

% Parameters for FISH channel 1 

minv1 = 50; %minimum volume for selected FISH clusters
maxv1 = 1200; %maximum volume for selected FISH clusters
minfmean1 = 00; %minimum FISH mean fluorescence for selected FISH clusters
maxfmean1 = 15000; %max FISH mean fluorescence for selected FISH clusters
minfmax1 = 00;    %minimum FISH max fluorescence for selected FISH clusters
maxfmax1 = 15000; %max FISH max fluorescence for selected FISH clusters

% Parameters for FISH channel 2 

minv2 = 50; %minimum volume for selected FISH clusters
maxv2 = 1200; %maximum volume for selected FISH clusters
minfmean2 = 00; %minimum FISH mean fluorescence for selected FISH clusters
maxfmean2 = 15000; %max FISH mean fluorescence for selected FISH clusters
minfmax2 = 00;    %minimum FISH max fluorescence for selected FISH clusters
maxfmax2 = 15000; %max FISH max fluorescence for selected FISH clusters

% Parameters for FISH channel 2 

minv3 = 50; %minimum volume for selected FISH clusters
maxv3 = 1200; %maximum volume for selected FISH clusters
minfmean3 = 00; %minimum FISH mean fluorescence for selected FISH clusters
maxfmean3 = 15000; %max FISH mean fluorescence for selected FISH clusters
minfmax3 = 00;    %minimum FISH max fluorescence for selected FISH clusters
maxfmax3 = 15000; %max FISH max fluorescence for selected FISH clusters


%% extracting measurement values for analysis

hm1mean = extractfield(r,'hm1mean');

fsh1v = extractfield(r,'fsh1v');
fsh1mean = extractfield(r,'fsh1mean');
fsh1max = extractfield(r,'fsh1max');

cond1 = maxv1>fsh1v&fsh1v>minv1&maxfmax1>fsh1max&fsh1max>minfmax1&maxfmean1>fsh1mean&fsh1mean>minfmean1;

hm1infsh1mean = extractfield(r,'hm1infsh1mean');
hm1infsh1mean = hm1infsh1mean(cond1); 
hm1tofsh1edge = extractfield(r,'hm1tofsh1edge');
hm1tofsh1edge = hm1tofsh1edge(cond1);
hm1tofsh1c = extractfield(r,'hm1tofsh1c');
hm1tofsh1c = hm1tofsh1c(cond1);
hm1infsh1frov = extractfield(r,'hm1infsh1frov');
hm1infsh1frov = hm1infsh1frov(cond1);

fsh2v = extractfield(r,'fsh2v');
fsh2mean = extractfield(r,'fsh2mean');
fsh2max = extractfield(r,'fsh2max');

cond2 = maxv2>fsh2v&fsh2v>minv2&maxfmax2>fsh2max&fsh2max>minfmax2&maxfmean2>fsh2mean&fsh2mean>minfmean2;

hm1infsh2mean = extractfield(r,'hm1infsh2mean');
hm1infsh2mean = hm1infsh2mean(cond2); 
hm1tofsh2edge = extractfield(r,'hm1tofsh2edge');
hm1tofsh2edge = hm1tofsh2edge(cond2);
hm1tofsh2c = extractfield(r,'hm1tofsh2c');
hm1tofsh2c = hm1tofsh2c(cond2);
hm1infsh2frov = extractfield(r,'hm1infsh2frov');
hm1infsh2frov = hm1infsh2frov(cond2);

fsh3v = extractfield(r,'fsh3v');
fsh3mean = extractfield(r,'fsh3mean');
fsh3max = extractfield(r,'fsh3max');

cond3 = maxv3>fsh3v&fsh3v>minv3&maxfmax3>fsh3max&fsh3max>minfmax3&maxfmean3>fsh3mean&fsh3mean>minfmean3;

hm1infsh3mean = extractfield(r,'hm1infsh3mean');
hm1infsh3mean = hm1infsh3mean(cond3); 
hm1tofsh3edge = extractfield(r,'hm1tofsh3edge');
hm1tofsh3edge = hm1tofsh3edge(cond3);
hm1tofsh3c = extractfield(r,'hm1tofsh3c');
hm1tofsh3c = hm1tofsh3c(cond3);
hm1infsh3frov = extractfield(r,'hm1infsh3frov');
hm1infsh3frov = hm1infsh3frov(cond3);


boxv = mean([fsh1v(cond1),fsh2v(cond2),fsh3v(cond3)]); 
 % volume of randomly generated cubic boxes

if chck > 0
else
h = waitbar(0,strcat('OBTAINING RANDOM BOXES: ',' 1',' of ',{' '},num2str(lnm)));
setappdata(h,'canceling',0);

for k = 1:lnm
    if getappdata(h,'canceling')
        break
    end
    fname = namelist{k};
    rboxes(k) = sceptre_rbgenerator(fname,maxlevel,nuclevel,choosefilter,r(k).xROI,r(k).yROI,r(k).zROI,ch,nuch,h1ch,h2ch,imorder,fstart,smooth,boxv,bn,bnd_nuc,imdimx,imdimy);
    step = strcat('OBTAINING RANDOM BOXES: ',{' '},num2str(k+1),' of',{' '},num2str(length(namelist)));
    waitbar(k/lnm,h,step)
    pause(0.1)
end
close(h)
end

rdv = extractfield(rboxes,'rdhm1v');
rdhm1mean = extractfield(rboxes,'rdhm1mean');

%% Violin plots of the distributions of Fluorescence signal in segmented clusters


hmf1 = NaN(5000000,5);
hmf2 = hmf1;

hmf1(1:length(hm1mean),h1o)= hm1mean(:)';
hmf1(1:length(hm1infsh1mean),f1o) = hm1infsh1mean(:)';
hmf1(1:length(hm1infsh2mean),f2o) = hm1infsh2mean(:)';
hmf1(1:length(hm1infsh3mean),f3o) = hm1infsh3mean(:)';
hmf1(1:length(rdhm1mean),rd1o)= rdhm1mean(:)';

m1 = strcat('n = ',{' '},num2str(length(hm1mean)));
m2 = strcat('(n = ',{' '},num2str(length(hm1infsh1mean)),')');
m3 = strcat('(n = ',{' '},num2str(length(hm1infsh2mean)),')');
m4 = strcat('(n = ',{' '},num2str(length(hm1infsh3mean)),')');
m5 = strcat('n = ',{' '},num2str(length(rdhm1mean)));


figure (1);

v1 = violinplot(hmf1);
v1(h1o).EdgeColor = [0 0 0];
v1(h1o).BoxColor = [0 0 0];
v1(h1o).MedianColor = [1 1 1];
v1(h1o).ViolinColor = [1 0 1];
v1(h1o).BoxWidth = 0.1;
v1(h1o).WhiskerPlot.LineWidth = 1.5;
v1(h1o).ShowData = 0;
v1(h1o).ViolinAlpha = 1;

v1(f1o).EdgeColor = [0 0 0];
v1(f1o).BoxColor = [0 0 0];
v1(f1o).MedianColor = [1 1 1];
v1(f1o).ViolinColor = [0 1 0];
v1(f1o).ScatterPlot.Marker = '.';
v1(f1o).ScatterPlot.MarkerEdgeColor = [0 0 0];
v1(f1o).BoxWidth = 0.1;
v1(f1o).WhiskerPlot.LineWidth = 1.5;
v1(f1o).ViolinAlpha = 0.25;

v1(f2o).EdgeColor = [0 0 0];
v1(f2o).BoxColor = [0 0 0];
v1(f2o).MedianColor = [1 1 1];
v1(f2o).ViolinColor = [0 1 0];
v1(f2o).ScatterPlot.Marker = '.';
v1(f2o).ScatterPlot.MarkerEdgeColor = [0 0 0];
v1(f2o).BoxWidth = 0.1;
v1(f2o).WhiskerPlot.LineWidth = 1.5;
v1(f2o).ViolinAlpha = 0.25;

v1(f3o).EdgeColor = [0 0 0];
v1(f3o).BoxColor = [0 0 0];
v1(f3o).MedianColor = [1 1 1];
v1(f3o).ViolinColor = [0 1 0];
v1(f3o).ScatterPlot.Marker = '.';
v1(f3o).ScatterPlot.MarkerEdgeColor = [0 0 0];
v1(f3o).BoxWidth = 0.1;
v1(f3o).WhiskerPlot.LineWidth = 1.5;
v1(f3o).ViolinAlpha = 0.25;

v1(rd1o).EdgeColor = [0 0 0];
v1(rd1o).BoxColor = [0 0 0];
v1(rd1o).MedianColor = [1 1 1];
v1(rd1o).ViolinColor = [0.5 0.5 0.5];
v1(rd1o).ScatterPlot.Marker = '.';
v1(rd1o).ScatterPlot.MarkerEdgeColor = [0.5 0.5 0.5];
v1(rd1o).BoxWidth = 0.1;
v1(rd1o).WhiskerPlot.LineWidth = 1.5;
v1(rd1o).ShowData = 0;
v1(rd1o).ViolinAlpha = 1;


xlim([0.5 5.5])
xticks([1 2 3 4 5])
xticklabels({hmark1,'random',gene1,gene2,gene3})
yticks([0 1000  2000  3000  4000])
yticklabels({'','1000','','3000',''})
ytickangle(90)
xtickangle(20)
yL = get(gca,'YLim');
    set(gca,'fontsize',ax_fs)
    set(gca,'fontname','Arial')
    set(gca,'linewidth',ax_lw)
    set(gca,'tickdir','out')
    ylabel(strcat(hmark1,' signal',[newline '(arb.)']),'fontsize',la_fs,'fontname','Arial')
    xlabel('clusters','fontsize',la_fs,'fontname','Arial')
    set(gcf,'Position',[100,100,425,350])
    set(gcf,'Renderer','painters')    

box on 

hold off

%% Figure for comaparing alleles of the same gene

%% Gene 1 allele comparisons
for k = 1:length(r)
    fsh(k).fsh1p = ones(length(r(k).fsh1mean),1).*k;
end

fsh1p = extractfield(fsh,'fsh1p'); %The cell number each loci corresponds to
fsh1p = fsh1p(cond1);

s = 0;
for k = 1:length(r)
    if isempty(hm1infsh1mean(fsh1p == k))>0
    else
        s = s+1;
        loc1(s).hm1mean = hm1infsh1mean(fsh1p == k);
        loc1(s).hm1frov = hm1infsh1frov(fsh1p == k);
        loc1(s).pos = k;
    end
end

s = 0;
for k = 1:length(loc1)
    if length(loc1(k).hm1mean)>1 %excludes cells with less than 2 FISH clusters
        if length(loc1(k).hm1mean)>4 %excludes cells with more than 4 FISH clusters
        else
        s = s+1;
        fdm = randperm(length(loc1(k).hm1mean));
         %Random index order for FISH clusters in each cell
        loc1A(s).hm1mean = loc1(k).hm1mean(fdm(1));
        loc1A(s).hm1frov = loc1(k).hm1frov(fdm(1));
        loc1B(s).hm1mean = loc1(k).hm1mean(fdm(2));
        loc1B(s).hm1frov = loc1(k).hm1frov(fdm(2));
        end
    end
end


x1 = extractfield(loc1A,'hm1mean')';
y1 = extractfield(loc1B,'hm1mean')';
mc1 = corrcoef(x1,y1);

figure(2); 
hs2 = scatter(extractfield(loc1A,'hm1mean'),extractfield(loc1B,'hm1mean'),sz,'o','g');
set(hs2,'Marker','o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0]);
set(hs2,'LineWidth',1.5);
ylim(hmark1mM)
xlim(hmark1mM)
yL = get(gca,'Ylim');
xL = get(gca,'XLim');
hold on
line(xL,[prctile(hm1mean,5) prctile(hm1mean,5)],'Color','k','LineWidth',1);
line([prctile(hm1mean,5) prctile(hm1mean,5)],yL,'Color','k','LineWidth',1);

str = strcat('{\it r} =',{' '},num2str(mc1(1,2),'%0.2f'));
dim = [.65 .8 .2 .1];
an1 = annotation('textbox',dim,'String',str,'FitBoxToText','on');
    set(an1,'fontsize',ax_fs)
    set(an1,'fontname','arial')
    set(an1,'LineStyle','none')
    set(an1,'FaceAlpha', 0.5)
    xlabel(strcat('Allele A',{' '},hmark1,' signal',[newline '(arb.)']),'fontsize',la_fs)
    ylabel(strcat('Allele B',{' '},hmark1,' signal',[newline '(arb.)']),'fontsize',la_fs)
    set(gca,'fontsize',ax_fs)
    set(gca,'linewidth',ax_lw)
    set(gca,'tickdir','out')
    set(gca)
    set(gcf,'Position',[100,100,500,500])
box on

hold off

clear fsh

%% Gene 2 allele comparisons

for k = 1:length(r)
    fsh(k).fsh2p = ones(length(r(k).fsh2mean),1).*k;
end

fsh2p = extractfield(fsh,'fsh2p'); %The cell number each loci corresponds to
fsh2p = fsh2p(cond2);

s = 0;
for k = 1:length(r)
    if isempty(hm1infsh2mean(fsh2p == k))>0
    else
        s = s+1;
        loc2(s).hm1mean = hm1infsh2mean(fsh2p == k);
        loc2(s).hm1frov = hm1infsh2frov(fsh2p == k);
        loc2(s).pos = k;
    end
end

s = 0;
for k = 1:length(loc2)
    if length(loc2(k).hm1mean)>1 %excludes cells with less than 2 FISH clusters
        if length(loc2(k).hm1mean)>4 %excludes cells with more than 4 FISH clusters
        else
        s = s+1;
        fdm = randperm(length(loc2(k).hm1mean));
         %Random index order for FISH clusters in each cell
        loc2A(s).hm1mean = loc2(k).hm1mean(fdm(1));
        loc2A(s).hm1frov = loc2(k).hm1frov(fdm(1));
        loc2B(s).hm1mean = loc2(k).hm1mean(fdm(2));
        loc2B(s).hm1frov = loc2(k).hm1frov(fdm(2));
        end
    end
end

x2 = extractfield(loc2A,'hm1mean')';
y2 = extractfield(loc2B,'hm1mean')';

mc2 = corrcoef(x2,y2);

figure(3); 
hs2 = scatter(extractfield(loc2A,'hm1mean'),extractfield(loc2B,'hm1mean'),sz,'o','g');
set(hs2,'Marker','o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0]);
set(hs2,'LineWidth',1.5);
ylim(hmark1mM)
xlim(hmark1mM)
yL = get(gca,'Ylim');
xL = get(gca,'XLim');
hold on
line(xL,[prctile(hm1mean,5) prctile(hm1mean,5)],'Color','k','LineWidth',1);
line([prctile(hm1mean,5) prctile(hm1mean,5)],yL,'Color','k','LineWidth',1);

str = strcat('{\it r} =',{' '},num2str(mc2(1,2),'%0.2f'));
dim = [.65 .8 .2 .1];
an1 = annotation('textbox',dim,'String',str,'FitBoxToText','on');
    set(an1,'fontsize',ax_fs)
    set(an1,'fontname','arial')
    set(an1,'LineStyle','none')
    set(an1,'FaceAlpha', 0.5)
    xlabel(strcat('Allele A',{' '},hmark1,' signal',[newline '(arb.)']),'fontsize',la_fs)
    ylabel(strcat('Allele B',{' '},hmark1,' signal',[newline '(arb.)']),'fontsize',la_fs)
    set(gca,'fontsize',ax_fs)
    set(gca,'linewidth',ax_lw)
    set(gca,'tickdir','out')
    set(gca)
    set(gcf,'Position',[100,100,500,500])
box on

hold off

clear fsh

%% Gene 3 allele comparisons

for k = 1:length(r)
    fsh(k).fsh3p = ones(length(r(k).fsh3mean),1).*k;
end

fsh3p = extractfield(fsh,'fsh3p'); 
fsh3p = fsh3p(cond3);

s = 0;
for k = 1:length(r)
    if isempty(hm1infsh3mean(fsh3p == k))>0
    else
        s = s+1;
        loc3(s).hm1mean = hm1infsh3mean(fsh3p == k);
        loc3(s).hm1frov = hm1infsh3frov(fsh3p == k);
        loc3(s).pos = k;
    end
end

s = 0;
for k = 1:length(loc3)
    if length(loc3(k).hm1mean)>1 %excludes cells with less than 2 FISH clusters
        if length(loc3(k).hm1mean)>4 %excludes cells with more than 4 FISH clusters
        else
        s = s+1;
        fdm = randperm(length(loc3(k).hm1mean));
         %Random index order for FISH clusters in each cell
        loc3A(s).hm1mean = loc3(k).hm1mean(fdm(1));
        loc3A(s).hm1frov = loc3(k).hm1frov(fdm(1));
        loc3B(s).hm1mean = loc3(k).hm1mean(fdm(2));
        loc3B(s).hm1frov = loc3(k).hm1frov(fdm(2));
        end
    end
end

x3 = extractfield(loc3A,'hm1mean')';
y3 = extractfield(loc3B,'hm1mean')';
mc3 = corrcoef(x3,y3);

figure(4); 
hs2 = scatter(extractfield(loc3A,'hm1mean'),extractfield(loc3B,'hm1mean'),sz,'o','g');
set(hs2,'Marker','o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0]);
set(hs2,'LineWidth',1.5);
ylim(hmark1mM)
xlim(hmark1mM)
yL = get(gca,'Ylim');
xL = get(gca,'XLim');
hold on
line(xL,[prctile(hm1mean,5) prctile(hm1mean,5)],'Color','k','LineWidth',1);
line([prctile(hm1mean,5) prctile(hm1mean,5)],yL,'Color','k','LineWidth',1);


str = strcat('{\it r} =',{' '},num2str(mc3(1,2),'%0.2f'));
dim = [.65 .8 .2 .1];
an1 = annotation('textbox',dim,'String',str,'FitBoxToText','on');
    set(an1,'fontsize',ax_fs)
    set(an1,'fontname','arial')
    set(an1,'LineStyle','none')
    set(an1,'FaceAlpha', 0.5)
    xlabel(strcat('Allele A',{' '},hmark1,' signal',[newline '(arb.)']),'fontsize',la_fs)
    ylabel(strcat('Allele B',{' '},hmark1,' signal',[newline '(arb.)']),'fontsize',la_fs)
    set(gca,'fontsize',ax_fs)
    set(gca,'linewidth',ax_lw)
    set(gca,'tickdir','out')
    set(gca)
    set(gcf,'Position',[100,100,500,500])
box on

hold off

%% Figures for comparing Fluorescence signal between alleles of different genes

%% Comparing gene 1 vs. gene 2;
s = 0;
for k = 1:length(loc1)
    for q = 1:length(loc2)
    if loc1(k).pos == loc2(q).pos %Select cell that has loci from both genes
        s = s + 1;
        fdm = randperm(length(loc1(k).hm1mean));
        %Random permutation to select allele for first gene
        fdm2 = randperm(length(loc2(q).hm1mean));
        %Random permutation to select allele for second gene
        loc1c2(s).loc1hm1mean = loc1(k).hm1mean(fdm(1));
        loc1c2(s).loc2hm1mean = loc2(q).hm1mean(fdm2(1));
    end
    end
end

x4 = extractfield(loc1c2,'loc1hm1mean')';
y4 = extractfield(loc1c2,'loc2hm1mean')';
mc4 = corrcoef(x4,y4);


figure(5); 
hs2 = scatter(x4,y4,sz,'o','g');
set(hs2,'Marker','o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0]);
set(hs2,'LineWidth',1.5);
ylim(hmark1mM)
xlim(hmark1mM)
yL = get(gca,'Ylim');
xL = get(gca,'XLim');
hold on
line(xL,[prctile(hm1mean,5) prctile(hm1mean,5)],'Color','k','LineWidth',1);
line([prctile(hm1mean,5) prctile(hm1mean,5)],yL,'Color','k','LineWidth',1);

str = strcat('{\it r} =',{' '},num2str(mc4(1,2),'%0.2f'));
dim = [.65 .8 .2 .1];
an1 = annotation('textbox',dim,'String',str,'FitBoxToText','on');
    set(an1,'fontsize',ax_fs)
    set(an1,'fontname','arial')
    set(an1,'LineStyle','none')
    set(an1,'FaceAlpha', 0.5)
    xlabel(strcat(gene1,' alleles',[newline hmark1 ' signal (arb.)']),'fontsize',la_fs)
    ylabel(strcat(gene2,' alleles',[newline hmark1 ' signal (arb.)']),'fontsize',la_fs)
    set(gca,'fontsize',ax_fs)
    set(gca,'linewidth',ax_lw)
    set(gca,'tickdir','out')
    set(gca)
    set(gcf,'Position',[100,100,500,500])
box on

hold off

%% Comparing gene 1 vs. gene 3

s = 0;
for k = 1:length(loc1)
    for q = 1:length(loc3)
    if loc1(k).pos == loc3(q).pos %Select cell that has loci from both genes
        s = s + 1;
        fdm = randperm(length(loc1(k).fsh1v));
        %Random permutation to select allele for first gene
        fdm2 = randperm(length(loc3(q).fsh3v));
        %Random permutation to select allele for second gene gene
        loc1c3(s).loc1hm1mean = loc1(k).hm1mean(fdm(1));
        loc1c3(s).loc3hm1mean = loc3(q).hm1mean(fdm2(1));
    end
    end
end

x5 = extractfield(loc1c3,'loc1hm1mean')';
y5 = extractfield(loc1c3,'loc3hm1mean')';
mc5 = corrcoef(x5,y5);

figure(6); 
hs2 = scatter(x5,y5,sz,'o','g');
set(hs2,'Marker','o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0]);
set(hs2,'LineWidth',1.5);
ylim(hmark1mM)
xlim(hmark1mM)
yL = get(gca,'Ylim');
xL = get(gca,'XLim');
hold on
line(xL,[prctile(hm1mean,5) prctile(hm1mean,5)],'Color','k','LineWidth',1);
line([prctile(hm1mean,5) prctile(hm1mean,5)],yL,'Color','k','LineWidth',1);

str = strcat('{\it r} =',{' '},num2str(mc5(1,2),'%0.2f'));
dim = [.65 .8 .2 .1];
an1 = annotation('textbox',dim,'String',str,'FitBoxToText','on');
    set(an1,'fontsize',ax_fs)
    set(an1,'fontname','arial')
    set(an1,'LineStyle','none')
    set(an1,'FaceAlpha', 0.5)
    xlabel(strcat(gene1,' alleles',[newline hmark1 ' signal (arb.)']),'fontsize',la_fs)
    ylabel(strcat(gene3,' alleles',[newline hmark1 ' signal (arb.)']),'fontsize',la_fs)
    set(gca,'fontsize',ax_fs)
    set(gca,'linewidth',ax_lw)
    set(gca,'tickdir','out')
    set(gca)
    set(gcf,'Position',[100,100,500,500])
box on

hold off

%% Comparing gene 2 vs. gene 3

s = 0;
for k = 1:length(loc2)
    for q = 1:length(loc3)
    if loc2(k).pos == loc3(q).pos %Select cell that has loci from both genes
        s = s + 1;
        fdm = randperm(length(loc2(k).fsh2v));
        %Random permutation to select allele for first gene gene
        fdm2 = randperm(length(loc3(q).fsh3v));
        %Random permutation to select allele for second gene gene
        loc2c3(s).loc2hm1mean = loc2(k).hm1mean(fdm(1));
        loc2c3(s).loc3hm1mean = loc3(q).hm1mean(fdm2(1));
    end
    end
end

x6 = extractfield(loc2c3,'loc2hm1mean')';
y6 = extractfield(loc2c3,'loc3hm1mean')';

mc6 = corrcoef(x6,y6);


figure(7); 
hs2 = scatter(x6,y6,sz,'o','g');
set(hs2,'Marker','o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0]);
set(hs2,'LineWidth',1.5);
ylim(hmark1mM)
xlim(hmark1mM)
yL = get(gca,'Ylim');
xL = get(gca,'XLim');
hold on
line(xL,[prctile(hm1mean,5) prctile(hm1mean,5)],'Color','k','LineWidth',1);
line([prctile(hm1mean,5) prctile(hm1mean,5)],yL,'Color','k','LineWidth',1);

str = strcat('{\it r} =',{' '},num2str(mc6(1,2),'%0.2f'));
dim = [.65 .8 .2 .1];
an1 = annotation('textbox',dim,'String',str,'FitBoxToText','on');
    set(an1,'fontsize',ax_fs)
    set(an1,'fontname','arial')
    set(an1,'LineStyle','none')
    set(an1,'FaceAlpha', 0.5)
    xlabel(strcat(gene2,' alleles',[newline hmark1 ' signal (arb.)']),'fontsize',la_fs)
    ylabel(strcat(gene3,' alleles',[newline hmark1 ' signal (arb.)']),'fontsize',la_fs)
    set(gca,'fontsize',ax_fs)
    set(gca,'linewidth',ax_lw)
    set(gca,'tickdir','out')
    set(gca)
    set(gcf,'Position',[100,100,500,500])
box on

hold off


