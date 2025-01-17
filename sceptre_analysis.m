%% sceptre_analysis 

%% Notes:
% This analysis script uses violinplot obtained from:
% https://github.com/bastibe/Violinplot-Matlab
% This analysis script is based on the example of only having a single FISH
% channel and up to 2 immunofluorescence channels. Adaptations must be made
% if multiple FISH channels will be used and analysed together.
%% Font parameters for figure

la_fs = 20; %Font for figure labels
ax_fs = 20; %Font for figures axes
ax_lw = 3; %line width for axes 
sz = 30; %size of markers for plots

%%  figure and analysis parameters

contlist = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]; %contour step list
chck = 0; % leave > 0 if random boxes have already been generated
hmark1mM = [0 1600]; % IF channel 1 axis limit for violinplots and scatters
hmark2mM = [0 1800]; % IF channel 2 axis limit for violinplots and scatters
hmark1corrmM = [0 1600]; %IF channel 1 axis limit for contours
hmark2corrmM = [0 1800]; %IF channel 2 axis limit for contours
bn = 300; %number of boxes to be randomly generated throughout each nucleus

gene1 = '{\it GAPDH}'; %name of targeted gene in FISH channel 1
hmark1 = 'K4me3'; 	%name of Immunofluorescence channel 1
hmark2 = 'K27me3';  %name of Immunofluorescence channel 2

%% FISH cluster selection parameters

% Below is the conditions used for selecting FISH clusters from segmented
% clusters in FISH channels. For the most part, only minimum volume is
% necessary for removing dim and non-specific clusters, but other 
% parameters such as maxfluorescence may remove clusters that are 
% fluorescent contaminants.

minv = 50; %minimum volume for selected FISH clusters
maxv = 1200; %maximum volume for selected FISH clusters
minfmean = 00; %minimum FISH mean fluorescence for selected FISH clusters
maxfmean = 15000; %max FISH mean fluorescence for selected FISH clusters
minfmax = 00;    %minimum FISH max fluorescence for selected FISH clusters
maxfmax = 15000; %max FISH max fluorescence for selected FISH clusters



%% extracting measurement values for analysis

hm1mean = extractfield(r,'hm1mean');

if h2ch > 0
   hm2mean = extractfield(r,'hm2mean');
   hm2inhm1mean = extractfield(r,'hm2inhm1mean'); 
   hm1inhm2mean = extractfield(r,'hm1inhm2mean');
   hm1inhm2frov =extractfield(r,'hm1inhm2frov');
   hm2inhm1frov =extractfield(r,'hm2inhm1frov');
end

fsh1v = extractfield(r,'fsh1v');
fsh1mean = extractfield(r,'fsh1mean');
fsh1max = extractfield(r,'fsh1max');
cond = maxv>fsh1v&fsh1v>minv&maxfmax>fsh1max&fsh1max>minfmax&maxfmean>fsh1mean&fsh1mean>minfmean;

hm1infsh1mean = extractfield(r,'hm1infsh1mean');
hm1infsh1mean = hm1infsh1mean(cond);  
hm1infsh1frov = extractfield(r,'hm1infsh1frov');
hm1infsh1frov = hm1infsh1frov(cond);

if h2ch > 0
    hm2infsh1mean = extractfield(r,'hm2infsh1mean');
    hm2infsh1mean = hm2infsh1mean(cond);
    hm2infsh1frov = extractfield(r,'hm2infsh1frov');
    hm2infsh1frov = hm2infsh1frov(cond);
end

boxv = mean(fsh1v(cond)); % volume of randomly generated cubic boxes

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
rdhm2mean = extractfield(rboxes,'rdhm2mean');



%% Violin plots of the distributions of Fluorescence signal in segmented clusters

hmf1 = NaN(5000000,3);


hmf1(1:length(hm1mean),1)= hm1mean(:)';
hmf1(1:length(hm1infsh1mean),3) = hm1infsh1mean(:)';
hmf1(1:length(rdhm1mean),2)= rdhm1mean(:)';

figure (1);
v1 = violinplot(hmf1);
v1(1).EdgeColor = [0 0 0];
v1(1).BoxColor = [0 0 0];
v1(1).MedianColor = [1 1 1];
v1(1).ViolinColor = [1 0 0];
v1(1).BoxWidth = 0.1;
v1(1).WhiskerPlot.LineWidth = 1.5;
v1(1).ShowData = 0;
v1(1).ViolinAlpha = 1;

v1(3).EdgeColor = [0 0 0];
v1(3).BoxColor = [0 0 0];
v1(3).MedianColor = [1 1 1];
v1(3).ViolinColor = [0 1 0];
v1(3).ScatterPlot.Marker = '.';
v1(3).ScatterPlot.MarkerEdgeColor = [0 0 0];
v1(3).BoxWidth = 0.1;
v1(3).WhiskerPlot.LineWidth = 1.5;
v1(3).ViolinAlpha = 0.25;


v1(2).EdgeColor = [0 0 0];
v1(2).BoxColor = [0 0 0];
v1(2).MedianColor = [1 1 1];
v1(2).ViolinColor = [0.5 0.5 0.5];
v1(2).BoxWidth = 0.1;
v1(2).WhiskerPlot.LineWidth = 1.5;
v1(2).ShowData = 0;
v1(2).ViolinAlpha = 1;

xlim([0.5 3.5])
xticks([1 2 3])


ylim(hmark1mM)
xticklabels({hmark1,'random',gene1})
    set(gca,'fontsize',ax_fs)
    set(gca,'fontname','arial')
    set(gca,'linewidth',ax_lw)
    set(gca,'tickdir','out')
    ylabel(strcat(hmark1,' signal',[newline '(arb.)']),'fontsize',la_fs)
    xlabel('clusters','fontsize',la_fs)
set(gcf,'Position',[100,100,500,320])
set(gcf,'Renderer','painters')
box on  

if h2ch > 0
    
    hmf2 = NaN(5000000,3);
    hmf3 = hmf2;
    hmf4 = hmf2;
    
    hmf2(1:length(hm2infsh1mean),3) = hm2infsh1mean(:)';
    hmf2(1:length(hm2mean),1)= hm2mean(:)';
    hmf2(1:length(rdhm1mean),2)= rdhm2mean(:)';

    hmf3(1:length(hm1mean),1)= hm1mean(:)';
    hmf3(1:length(hm2mean),3)= hm1inhm2mean(:)';
    hmf3(1:length(rdhm1mean),2)= rdhm1mean(:)';


    hmf4(1:length(hm1mean),3)= hm2inhm1mean(:)';
    hmf4(1:length(hm2mean),1)= hm2mean(:)';
    hmf4(1:length(rdhm1mean),2)= rdhm2mean(:)';
    
    figure (2);
    v1 = violinplot(hmf2);
    v1(1).EdgeColor = [0 0 0];
    v1(1).BoxColor = [0 0 0];
    v1(1).MedianColor = [1 1 1];
    v1(1).ViolinColor = [0 0 1];
    v1(1).BoxWidth = 0.1;
    v1(1).WhiskerPlot.LineWidth = 1.5;
    v1(1).ShowData = 0;
    v1(1).ViolinAlpha = 1;

    v1(3).EdgeColor = [0 0 0];
    v1(3).BoxColor = [0 0 0];
    v1(3).MedianColor = [1 1 1];
    v1(3).ViolinColor = [0 1 0];
    v1(3).ScatterPlot.Marker = '.';
    v1(3).ScatterPlot.MarkerEdgeColor = [0 0 0];
    v1(3).BoxWidth = 0.1;
    v1(3).WhiskerPlot.LineWidth = 1.5;
    v1(3).ViolinAlpha = 0.25;


    v1(2).EdgeColor = [0 0 0];
    v1(2).BoxColor = [0 0 0];
    v1(2).MedianColor = [1 1 1];
    v1(2).ViolinColor = [0.5 0.5 0.5];
    v1(2).BoxWidth = 0.1;
    v1(2).WhiskerPlot.LineWidth = 1.5;
    v1(2).ShowData = 0;
    v1(2).ViolinAlpha = 1;

    xlim([0.5 3.5])
    xticks([1 2 3])
    
    ylim(hmark2mM)
    xticklabels({hmark2,'random',gene1})
        set(gca,'fontsize',ax_fs)
        set(gca,'fontname','arial')
        set(gca,'linewidth',ax_lw)
        set(gca,'tickdir','out')
        ylabel(strcat(hmark2,' signal',[newline '(arb.)']),'fontsize',la_fs)
        xlabel('clusters','fontsize',la_fs)
    set(gcf,'Position',[100,100,500,320])
    set(gcf,'Renderer','painters')
    box on  

    figure (3);
    v1 = violinplot(hmf3);
    v1(1).EdgeColor = [0 0 0];
    v1(1).BoxColor = [0 0 0];
    v1(1).MedianColor = [1 1 1];
    v1(1).ViolinColor = [1 0 0];
    v1(1).BoxWidth = 0.1;
    v1(1).WhiskerPlot.LineWidth = 1.5;
    v1(1).ShowData = 0;
    v1(1).ViolinAlpha = 1;

    v1(3).EdgeColor = [0 0 0];
    v1(3).BoxColor = [0 0 0];
    v1(3).MedianColor = [1 1 1];
    v1(3).ViolinColor = [0 0 1];
    v1(3).BoxWidth = 0.1;
    v1(3).WhiskerPlot.LineWidth = 1.5;
    v1(3).ShowData = 0;
    v1(3).ViolinAlpha = 1;


    v1(2).EdgeColor = [0 0 0];
    v1(2).BoxColor = [0 0 0];
    v1(2).MedianColor = [1 1 1];
    v1(2).ViolinColor = [0.5 0.5 0.5];
    v1(2).BoxWidth = 0.1;
    v1(2).WhiskerPlot.LineWidth = 1.5;
    v1(2).ShowData = 0;
    v1(2).ViolinAlpha = 1;

    xlim([0.5 3.5])
    xticks([1 2 3])

    ylim([hmark1mM])
    xticklabels({hmark1,'random',hmark2})
        set(gca,'fontsize',ax_fs)
        set(gca,'fontname','arial')
        set(gca,'linewidth',ax_lw)
        set(gca,'tickdir','out')
        ylabel(strcat(hmark1,' signal',[newline '(arb.)']),'fontsize',la_fs)
        xlabel('clusters','fontsize',la_fs)
    set(gcf,'Position',[100,100,500,320])
    set(gcf,'Renderer','painters')
    box on  


    figure (4);
    v1 = violinplot(hmf4);
    v1(1).EdgeColor = [0 0 0];
    v1(1).BoxColor = [0 0 0];
    v1(1).MedianColor = [1 1 1];
    v1(1).ViolinColor = [0 0 1];
    v1(1).BoxWidth = 0.1;
    v1(1).WhiskerPlot.LineWidth = 1.5;
    v1(1).ShowData = 0;
    v1(1).ViolinAlpha = 1;

    v1(3).EdgeColor = [0 0 0];
    v1(3).BoxColor = [0 0 0];
    v1(3).MedianColor = [1 1 1];
    v1(3).ViolinColor = [1 0 0];
    v1(3).BoxWidth = 0.1;
    v1(3).WhiskerPlot.LineWidth = 1.5;
    v1(3).ShowData = 0;
    v1(3).ViolinAlpha = 1;


    v1(2).EdgeColor = [0 0 0];
    v1(2).BoxColor = [0 0 0];
    v1(2).MedianColor = [1 1 1];
    v1(2).ViolinColor = [0.5 0.5 0.5];
    v1(2).BoxWidth = 0.1;
    v1(2).WhiskerPlot.LineWidth = 1.5;
    v1(2).ShowData = 0;
    v1(2).ViolinAlpha = 1;

    xlim([0.5 3.5])
    xticks([1 2 3])

    ylim(hmark2mM)
    xticklabels({hmark2,'random',hmark1})
        set(gca,'fontsize',ax_fs)
        set(gca,'fontname','arial')
        set(gca,'linewidth',ax_lw)
        set(gca,'tickdir','out')
        ylabel(strcat(hmark2,' signal',[newline '(arb.)']),'fontsize',la_fs)
        xlabel('clusters','fontsize',la_fs)
    set(gcf,'Position',[100,100,500,320])
    set(gcf,'Renderer','painters')
    box on  

end


%% Figures for Allele comparisons

for k = 1:length(r)

    fsh(k).fsh1p = ones(length(r(k).fsh1mean),1).*k;
end

fsh1p = extractfield(fsh,'fsh1p'); %The cell number each loci corresponds to
fsh1p = fsh1p(cond);

fsh1vol = fsh1v(cond);
fsh1m = fsh1mean(cond);
fsh1Mm = fsh1max(cond);
s = 0;
for k = 1:length(r)
    if isempty(hm1infsh1mean(fsh1p == k))>0
    else
        s = s+1;
        loc1(s).hm1mean = hm1infsh1mean(fsh1p == k);
        loc1(s).hm1frov = hm1infsh1frov(fsh1p == k);
        if h2ch > 0
            loc1(s).hm2mean = hm2infsh1mean(fsh1p == k);
            loc1(s).hm2frov = hm2infsh1frov(fsh1p == k);
        end
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
        loc1A(s).hm1frov = loc1(k).hm1frov(fdm(1));            ;
        loc1B(s).hm1mean = loc1(k).hm1mean(fdm(2));  
        loc1B(s).hm1frov = loc1(k).hm1frov(fdm(2));  
            if h2ch > 0
                loc1A(s).hm2mean = loc1(k).hm2mean(fdm(1));
                loc1A(s).hm2frov = loc1(k).hm2frov(fdm(1));
                loc1B(s).hm2mean = loc1(k).hm2mean(fdm(2));
                loc1B(s).hm2frov = loc1(k).hm2frov(fdm(2));
            end
        end
    end
end

x1 = extractfield(loc1A,'hm1mean')';
y1 = extractfield(loc1B,'hm1mean')';

mc1 = corrcoef(x1,y1);

figure(5); 

hs2 = scatter(extractfield(loc1A,'hm1mean'),extractfield(loc1B,'hm1mean'),sz,'o','g');
set(hs2,'Marker','o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0]);
set(hs2,'LineWidth',1.5);
ylim(hmark1mM)
xlim(hmark1mM)
hold on
line(hmark1mM,[prctile(hm1mean,5) prctile(hm1mean,5)],'Color','k','LineWidth',1);
line([prctile(hm1mean,5) prctile(hm1mean,5)],hmark1mM,'Color','k','LineWidth',1);

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
    set(gcf,'Position',[100,100,500,500])
box on

hold off


if h2ch > 0

    x2 = extractfield(loc1A,'hm2mean')';
    y2 = extractfield(loc1B,'hm2mean')';

    mc2 = corrcoef(x2,y2);

    figure(6); 

    hs2 = scatter(extractfield(loc1A,'hm2mean'),extractfield(loc1B,'hm2mean'),sz,'o','filled','r');
    set(hs2,'Marker','o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0]);
    set(hs2,'LineWidth',1.5);
    ylim(hmark2mM)
    xlim(hmark2mM)
    xL = get(gca,'XLim');
    hold on

    line(hmark2mM,[prctile(hm2mean,5) prctile(hm2mean,5)],'Color','k','LineWidth',1);
    line([prctile(hm2mean,5) prctile(hm2mean,5)],hmark2mM,'Color','k','LineWidth',1);

    str = strcat('{\it r} =',{' '},num2str(mc2(1,2),'%0.2f'));
    dim = [.65 .8 .2 .1];
    an1 = annotation('textbox',dim,'String',str,'FitBoxToText','on');
        set(an1,'fontsize',ax_fs)
        set(an1,'fontname','arial')
        set(an1,'LineStyle','none')
    set(gca,'fontsize',ax_fs)
        xlabel(strcat('Allele A',{' '},hmark2,' signal',[newline '(arb.)']),'fontsize',la_fs)
        ylabel(strcat('Allele B',{' '},hmark2,' signal',[newline '(arb.)']),'fontsize',la_fs)
        set(gca,'linewidth',ax_lw)
        set(gca,'tickdir','out')

        set(gcf,'Position',[100,100,500,500])
    box on

    hold off
end


%% correlation plot for gene 1
if h2ch > 0
    x1 = double(hm1infsh1mean)';
    y1 = double(hm2infsh1mean)'; 
    mc1 = corrcoef(x1,y1);
    
    figure(7); 
    hs2 = scatter(x1,y1,sz,'.','k');
    ylim(hmark2corrmM);
    xlim(hmark1corrmM);
    
    set(hs2,'Marker','o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0]);
    set(hs2,'LineWidth',1.5);

    hold on
    str = strcat('{\it r} =',{' '},num2str(mc1(1,2),'%0.2f'));
    dim = [.65 .8 .2 .1];
    an1 = annotation('textbox',dim,'String',str,'FitBoxToText','on');
        set(an1,'fontsize',ax_fs)
        set(an1,'LineStyle','none')
        xlabel(strcat(hmark1,' signal',[newline '(arb.)']),'fontsize',la_fs)
        ylabel(strcat(hmark2,' signal',[newline '(arb.)']),'fontsize',la_fs)
        set(gca,'fontsize',ax_fs)
        set(gca,'FontName','arial')
        set(gca,'linewidth',ax_lw)
        set(gca,'tickdir','out')

        set(gcf,'Position',[100,100,500,500])
    xL = get(gca,'XLim');
    yL = get(gca,'YLim');

    line(xL,[prctile(hm2mean,5) prctile(hm2mean,5)],'Color','k','LineWidth',1);
    line([prctile(hm1mean,5) prctile(hm1mean,5)],yL,'Color','k','LineWidth',1);
    box on
    hold off
end

%% Contour plots for comparisons of fluorescence signals within clusters 

% Making 2D histograms 

hmno = length(hm1mean);
hmstd = std(hm1mean);
bins1 = 3.49*hmstd/(hmno^(1/3)); % Scott's rule
xbar1 = (bins1/2):bins1:max(hm1mean)+(bins1/2);

hmno = length(hm2mean);
hmstd = std(hm2mean);
bins2 = 3.49*hmstd/(hmno^(1/3)); % Scott's rule
xbar2 = (bins2/2):bins2:max(hm2mean)+(bins2/2);

smoother = 1; 
hg = fspecial('gaussian',[5 5],smoother);
    % gaussian filter matrix for smoothing contours
clear edges

edges{2} = xbar1;
edges{1} = xbar2;

%Contours for immunofluorescence channel 1 clusters
[h1n,h1c] = hist3([hm2inhm1mean',hm1mean'],'Ctrs',edges);
h1n2 = conv2(h1n,hg,'same');

%Contours for immunofluorescence channel 2 clusters
[h2n,h2c] = hist3([hm2mean',hm1inhm2mean'],'Ctrs',edges);
h2n2 = conv2(h2n,hg,'same');

%Contours for reandomly selected region clusters
[r1n,r1c] = hist3([rdhm2mean',rdhm1mean'],'Ctrs',edges);
r1n2 = conv2(r1n,hg,'same');

h1n2m = max(h1n2(:));
h2n2m = max(h2n2(:));
r1n2m = max(r1n2(:));


%% Contours for randomly selected regions
if h2ch > 0
figure(9);

hs2 = scatter(x1,y1,sz,'.','k');
set(hs2,'Marker','.','MarkerEdgeColor',[0.25 0.25 0.25])

ylim(hmark2corrmM);
xlim(hmark1corrmM);
xL = get(gca,'XLim');
yL = get(gca,'YLim');
hold on
[n3,h3] = contourf(r1c{2},r1c{1},r1n2./r1n2m,'color','flat','LineWidth',2,'LevelList',contlist);
h3.LineWidth = 0.5;
h3.LineColor = 'k';

map = [ 3.25 3.25 3.25
        4  4  4
        4.75 4.75 4.75 
        5.5 5.5 5.5
        6.25 6.25 6.25
        7  7  7
        7.75 7.75 7.75
        8.5 8.5 8.5
        9.25 9.25 9.25
        10 10 10].*0.1;
    %Matrix map for decreasing shades of gray
colormap(map)
line(xL,[prctile(hm2mean,5) prctile(hm2mean,5)],'Color','k','LineWidth',1);
line([prctile(hm1mean,5) prctile(hm1mean,5)],yL,'Color','k','LineWidth',1);
    x1 = double(rdhm1mean)';
    y1 = double(rdhm2mean)'; 
    mc1 = corrcoef(x1,y1);

str = strcat('{\it r} =',{' '},num2str(mc1(1,2),'%0.2f'));
dim = [.65 .8 .2 .1];
an1 = annotation('textbox',dim,'String',str,'FitBoxToText','on');
    set(an1,'fontsize',ax_fs)
    set(an1,'LineStyle','none')
    xlabel(strcat(hmark1,' signal',[newline '(arb.)']),'fontsize',la_fs)
    ylabel(strcat(hmark2,' signal',[newline '(arb.)']),'fontsize',la_fs)
    set(gca,'fontsize',ax_fs)
    set(gca,'FontName','arial')
    set(gca,'linewidth',ax_lw)
    set(gca,'tickdir','out')
    
    set(gcf,'Position',[100,100,500,500])
box on

hold off

%% hmark1 plot
x1 = double(hm1mean)';
y1 = double(hm2inhm1mean)';
mc1 = corrcoef(x1,y1);

rp1 = randperm(length(x1));
x1 = x1(rp1);
y1 = y1(rp1);
x1 = x1(1:ceil(length(x1)/100));
y1 = y1(1:ceil(length(y1)/100));

map = [ 10  1  1
        10  2  2
        10  3  3
        10  4  4
        10  5  5
        10  6  6
        10  7  7
        10  8  8
        10  9  9
        10 10 10].*0.1;
    %Matrix map for decreasing shades of red
figure(11); 
hs2 = scatter(x1,y1,sz,'.','r');
ylim(hmark2corrmM);
xlim(hmark1corrmM);
hold on
str = strcat('{\it r} =',{' '},num2str(mc1(1,2),'%0.2f'));
dim = [.65 .8 .2 .1];
an1 = annotation('textbox',dim,'String',str,'FitBoxToText','on');
    set(an1,'fontsize',ax_fs)
    set(an1,'LineStyle','none')
xL = get(gca,'XLim');
yL = get(gca,'YLim');
[n1,h1] = contourf(h1c{2},h1c{1},h1n2./h1n2m,'color','flat','LineWidth',2,'LevelList',contlist);
h1.LineWidth = 0.5;
h1.LineColor = 'k';

colormap(map)

line(xL,[prctile(hm2mean,5) prctile(hm2mean,5)],'Color','k','LineWidth',1);
line([prctile(hm1mean,5) prctile(hm1mean,5)],yL,'Color','k','LineWidth',1);
    xlabel(strcat(hmark1,' signal',[newline '(arb.)']),'fontsize',la_fs)
    ylabel(strcat(hmark2,' signal',[newline '(arb.)']),'fontsize',la_fs)
    set(gca,'fontsize',ax_fs)
    set(gca,'FontName','arial')
    set(gca,'linewidth',ax_lw)
    set(gca,'tickdir','out')
    
    set(gcf,'Position',[100,100,500,500])
box on

hold off

%% hmark2 plot
x1 = double(hm1inhm2mean)';
y1 = double(hm2mean)';
mc1 = corrcoef(x1,y1);

rp1 = randperm(length(x1));
x1 = x1(rp1);
y1 = y1(rp1);
x1 = x1(1:ceil(length(x1)/100));
y1 = y1(1:ceil(length(y1)/100));

map = [  1  1 10
         2  2 10
         3  3 10
         4  4 10
         5  5 10
         6  6 10
         7  7 10
         8  8 10
         9  9 10
        10 10 10].*0.1;
    % matrix map for decreasing shades of blue   
figure(12); 
hs2 = scatter(x1,y1,sz,'.','b');
ylim(hmark2corrmM);
xlim(hmark1corrmM);
hold on
str = strcat('{\it r} =',{' '},num2str(mc1(1,2),'%0.2f'));
dim = [.65 .8 .2 .1];
an1 = annotation('textbox',dim,'String',str,'FitBoxToText','on');
    set(an1,'fontsize',ax_fs)
    set(an1,'LineStyle','none')
xL = get(gca,'XLim');
yL = get(gca,'YLim');
[n2,h2] = contourf(h2c{2},h2c{1},h2n2./h2n2m,'color','flat','LevelList',contlist);
h2.LineWidth = 0.5;
h2.LineColor = 'k';
colormap(map)
line(xL,[prctile(hm2mean,5) prctile(hm2mean,5)],'Color','k','LineWidth',1);
line([prctile(hm1mean,5) prctile(hm1mean,5)],yL,'Color','k','LineWidth',1);
    xlabel(strcat(hmark1,' signal',[newline '(arb.)']),'fontsize',la_fs)
    ylabel(strcat(hmark2,' signal',[newline '(arb.)']),'fontsize',la_fs)
    set(gca,'fontsize',ax_fs)
    set(gca,'FontName','arial')
    set(gca,'linewidth',ax_lw)
    set(gca,'tickdir','out')
    
    set(gcf,'Position',[100,100,500,500])
box on

hold off
end