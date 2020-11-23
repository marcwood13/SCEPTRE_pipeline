function data = sceptre_processing(fname,ch,imorder,fstart,smooth,bnd_nuc,edgedel,maxlevel,choosefilter,nuch,h1ch,h2ch,f1ch,f2ch,f3ch,nuclevel,hm1level,hm2level,fsh1level,fsh2level,fsh3level,minvol,imdimx,imdimy)

%% Step 1: Obtain image info and XYZ boundaries

info = imfinfo(fname);
L = length(info);
imwidth = info(1).Width;
imheight = info(1).Height;

xrange = imdimx;
xROI = ((imwidth - imdimx)/2)+1:((imwidth + imdimx)/2); 

yrange = imdimy;
yROI =  ((imheight - imdimy)/2)+1:((imheight + imdimy)/2);

zrange = (L/ch)+1-fstart;

zROI = 1:zrange;
zmat = zeros(xrange,yrange,zrange); %zero matrix for current ROI
laplace = fspecial('laplacian',0.2); %laplace filter

data.name = fname;

%% Selecting image Z indexes

if imorder>0
    
    zhm1 = h1ch+ch*(fstart-1):ch:L-(ch-h1ch);  %z index for h1ch
    
    if h2ch>0
        zhm2 = h2ch+ch*(fstart-1):ch:L-(ch-h2ch);  %z index for h2ch
    end
    
    if nuch>0
        znuc = nuch+ch*(fstart-1):ch:L-(ch-nuch); %z index for  nuch
    end
    
    if f1ch>0
        zfsh1 = f1ch+ch*(fstart-1):ch:L-(ch-f1ch); %z index for  f1ch
    end
    
    if f2ch>0
        zfsh2 = f2ch+ch*(fstart-1):ch:L-(ch-f2ch); %z index for  f2ch
    end
    
    if f3ch>0
        zfsh3 = f3ch+ch*(fstart-1):ch:L-(ch-f3ch); %z index for  f3ch
    end
    
else

    zhm1 = (L*(h1ch-1)/ch)+fstart:L*h1ch/ch; %z index for  h1ch
    
    if h2ch>0
        zhm2 = (L*(h2ch-1)/ch)+fstart:L*h2ch/ch; %z index for  h2ch
    end
    
    if nuch>0
        znuc = (L*(nuch-1)/ch)+fstart:L*nuch/ch; %z index for  nuch
    end
    
    if f1ch>0    
        zfsh1 = (L*(f1ch-1)/ch)+fstart:L*f1ch/ch; %z index for  f1ch
    end
    
    if f2ch>0    
        zfsh2 = (L*(f2ch-1)/ch)+fstart:L*f2ch/ch; %z index for  f2ch
    end
    
    if f3ch>0    
        zfsh3 = (L*(f3ch-1)/ch)+fstart:L*f3ch/ch; %z index for  f3ch
    end
end

%% Step 2: generate nuclear mask

nuc = zmat;
if nuch>0
    for k = 1:length(znuc)
         im0 = imread(fname,'tif','Index',znuc(k)); %load image
         im1 = single(im0(xROI,yROI)); %crop image
         nuc(:,:,k) = im1; 
    end 
else
    for k = 1:length(zhm1)
         im0 = imread(fname,'tif','Index',zhm1(k)); %load image
         im1 = single(im0(xROI,yROI)); %crop image
         nuc(:,:,k) = im1; 
    end
end

se1 = strel('sphere',bnd_nuc); %structure for dilation/erosion
nucs = imgaussfilt(nuc,smooth); %3D image stack smoothing
m = quantile(nucs(nucs<maxlevel),0.5); %Median of image stack
d = quantile(nucs(nucs<maxlevel),0.75)-m; %Third quartile of image stack
high = double(max(nucs(nucs<maxlevel))); %high for image contrast
low = m + nuclevel*d; %low and clipping threshold for image contrast
     data.nuclow = low;
     data.nuchigh = high;

nucadj = imadjustn(nucs./high, [low/high 1], [0 1]); 
nucbin = zmat;
if choosefilter > 0
    for k = 1:size(nucadj,3)    
        im2 = imfilter(nucadj(:,:,k), laplace);  % laplacian filter
        im3 = (im2 < 0);   % negative values of the laplacian
        nucbin(:,:,k) = im3; % first binary image of nucleus
    end
else
    nucbin = imbinarize(nucadj);
end

nucl = imdilate(nucbin,se1); %dilation to unite clusters across nucleus
nucbw = bwlabeln(nucl,18);
nucp = regionprops3(nucbw,'Volume');
mnuc = max(nucp.Volume); %largest object (i.e. nucleus) size
nucsel = bwareaopen(nucl,mnuc-1,18); %Select only nucleus
nucf0 = zmat;
for k = 1:size(nucsel,3)
    nucf0(:,:,k) = bwconvhull(nucsel(:,:,k),'union');
end
nucf0 = imerode(nucf0,se1); %erode nucleus to "reverse" dilation

clear nucl nucbw nucp mnuc nucsel

%% Limiting ROI to nuclear mask region

zpick = zeros(1,zrange);
for k = 1:size(nucf0,3)
    if sum(sum(nucf0(:,:,k)))>0
    zpick(k) = 1;
    else
        zpick(k) = 0; 
    end
end
zROId = zROI(logical(zpick)); %All planes with nuclear mask in them
mip = max(nucf0,[],3); %max projection of nuclear mask image stack
mx = max(mip'); 
my = max(mip); 
mx2 = mx.*xROI;
my2 =my.*yROI;
xROId = mx2(mx2>0); %new x dimension for image cropping
yROId = my2(my2>0); %new y dimension for y cropping 
nucf = nucf0(xROId-min(xROI)+1,yROId-min(yROI)+1,zROId); 
    %limiting nuclear mask region
nucedge = double(bwperim(nucf)); %Perimeter for removing perimiter clusters
nucbw = bwlabeln(nucf);
nucv = regionprops3(nucbw,'Volume');
data.nucv = max(nucv.Volume);
clear nucs nucadj nucbw nucbw nucp mnuc volnuc nucsel m d low high nucl thrsh nucf0 mip mx mx2 my my2 nucbin
zmat = zeros(length(xROId),length(yROId),length(zROId));

%% Segementation of immunofluorescence channel 1
hm1 = zmat; 
zhm1 = zhm1(logical(zpick));

for k = 1:length(zhm1)
     im0 = imread(fname,'tif','Index',zhm1(k)); %load  IF channel image
     im1 = single(im0(xROId,yROId)); %crop image
     hm1(:,:,k) = im1;
end

hms = imgaussfilt(hm1,smooth);  %3D image stack smoothing
m = quantile(hms(hms<maxlevel),0.5); %Median of image stack
d = quantile(hms(hms<maxlevel),0.75)-m; %Third quartile of image stack
high = double(max(hms(hms<maxlevel))); %high for image contrast
low = m + hm1level*d; %low and clipping threshold for image contrast
data.hm1low = low;
data.hm1high = high;
hmadj = imadjustn(hms./high,[low/high 1], [0 1]);   % contrast adjustment
hml = zmat;

if choosefilter > 0
    for q = 1:size(hmadj,3)    
        im2 = imfilter(hmadj(:,:,q), laplace);  % laplacian filter
        im3 = (im2 < 0);   % negative values of the laplacian
        hml(:,:,q) = im3;
    end
else
    hml = imbinarize(hmadj);
end
% Watershedding for immunofluorescence channel channel

D = bwdist(~hml);
D = -D;
D(~hml) = Inf;
L = watershed(D);
L(~hml) = 0;
hmw = hml.*(L>0);

hmm = immultiply(hmw,nucf); %Apply nuclear mask to channel
hm1f = bwareaopen(hmm,minvol,18); %Removing small morphological components
clear hmf hmw hmm hms im0 im1 im2 im3 D L m d low high hmadj thrsh

%% Segmentation of HM Channel 2

if h2ch>0
    hm2 = zmat;
    zhm = zhm2(logical(zpick));

    for k = 1:length(zhm)
         im0 = imread(fname,'tif','Index',zhm(k)); %load IF channel image
         im1 = single(im0(xROId,yROId)); %crop image
         hm2(:,:,k) = im1;
    end


    hms = imgaussfilt(hm2,smooth);
    m = quantile(hms(hms<maxlevel),0.5); %Median of image stack
    d = quantile(hms(hms<maxlevel),0.75)-m; %Third quartile of image stack
    high = double(max(hms(hms<maxlevel))); %high for image contrast
    low = m + hm2level*d; %low and clipping threshold for image contrast
    data.hm2low = low;
    data.hm2high = high;
    hmadj = imadjustn(hms./high,[low/high 1], [0 1]);   % contrast adjustment
    hml = zmat;
    
    if choosefilter > 0
        for q = 1:size(hmadj,3)    
            im2 = imfilter(hmadj(:,:,q), laplace);  % laplacian filter
            im3 = (im2 < 0);   % negative values of the laplacian
            hml(:,:,q) = im3;
        end
    else
        hml = imbinarize(hmadj);
    end
    % Watershedding for HM channel
    D = bwdist(~hml);
    D = -D;
    D(~hml) = Inf;
    L = watershed(D);
    L(~hml) = 0;
    hmw = hml.*(L>0);

    hmm = immultiply(hmw,nucf); %Apply nuclear mask to channel
    hm2f = bwareaopen(hmm,minvol,18); %Removing small morphological components
    clear hmf hmw hmm hms im0 im1 im2 im3 D L m d low high hmadj Lhm Rhm rdhma hma btlevel btper thrsh
end 

%% Segmenting FISH channel 1
if f1ch>0
    fshlevel = fsh1level;
    zfsh = zfsh1(logical(zpick));
    fsh1 = zmat;
    
    for k = 1:length(zfsh)
        im0 = imread(fname,'tif','Index',zfsh(k)); %load FISH channel 1 image
        im1 = single(im0(xROId,yROId)); %crop image
        fsh1(:,:,k) = im1;
    end

    fshs = imgaussfilt(fsh1,smooth);  %3D image stack smoothing
    fshm = immultiply(fshs,nucf); %Apply nuclear mask to channel 
    m = quantile(fshm(fshm>0&fshm<maxlevel),0.5); %median of image stack
    d = quantile(fshm(fshm>0&fshm<maxlevel),0.75)-m; %Third quartile of image stack
	low = m + fshlevel*d; %low and clipping threshold for image contrast
    high = max(fshm(fshm>0&fshm<maxlevel)); %high for image contrast
    data.fsh1low = low;
    data.fsh1high = high;
    if high>low
        fshadj = imadjustn(fshm./high,[low/high 1], [0 1]);   % contrast adjustment
        fshbin= zmat;
        if choosefilter > 0
            for k = 1:size(fshadj,3)    
                im2 = imfilter(fshadj(:,:,k), laplace);  % laplacian filter
                im3 = (im2 < 0);   % negative values of the laplacian
                fshbin(:,:,k) = im3;
            end
        else
            fshbin(:,:,k) = imbinarize(fshadj);
        end
        
        if edgedel > 0
            fshedge = logical(imadd(fshbin,nucedge)); %Adds nuclear periphery
            fshsel = fshedge - bwareaopen(fshedge,2000); 
                %removes periphery with attaching clusters
            fsh1f = bwareaopen(fshsel,minvol,18);
        else
            fsh1f = bwareaopen(fshbin,minvol,18);
        end
    else
        fsh1f = zmat; %When the channel is empty, this avoids the script from crashing
    end
    clear fshsel fshedge fshbin im3 im2 im1 fshadj high low m d fshm fshs chb
    Lfsh = bwlabeln(fsh1f,18);
    Rfsh = regionprops3(Lfsh,fsh1,'Volume','MeanIntensity','MaxIntensity','Centroid');
       data.fsh1v = single([Rfsh.Volume]); 
       data.fsh1mean = single([Rfsh.MeanIntensity]);
       data.fsh1max = single([Rfsh.MaxIntensity]);
       data.fsh1c = single([Rfsh.Centroid]);
    Rhminfsh = regionprops3(Lfsh,hm1,'MeanIntensity');
       data.hm1infsh1mean = single([Rhminfsh.MeanIntensity]); 
        % IF channel 1 fluorescence in FISH 1 channel clusters
    Rhminfsh = regionprops3(Lfsh,hm1f,'MeanIntensity');
       data.hm1infsh1frov = single([Rhminfsh.MeanIntensity]);
       % IF cluster fraction of overlap with FISH 1 channel clusters
    clear Rhminfsh hfsh Lhfsh Rhfsh
    
    if h2ch>0
        Rhminfsh = regionprops3(Lfsh,hm2,'MeanIntensity');
            data.hm2infsh1mean = single([Rhminfsh.MeanIntensity]);
            % IF channel 1 fluorescence in FISH 1 channel clusters
        Rhminfsh = regionprops3(Lfsh,hm2f,'MeanIntensity');
            data.hm2infsh1frov = single([Rhminfsh.MeanIntensity]);
            % IF cluster fraction of overlap with FISH 1 channel clusters
    end   
    clear Lfsh Rfsh fsh1 fsh1f
end
%% Segmenting FISH channel 2
if f2ch>0
    fshlevel = fsh2level;
    zfsh = zfsh2(logical(zpick));
    fsh2 = zmat;

    for k = 1:length(zfsh)
        im0 = imread(fname,'tif','Index',zfsh(k)); %load image
        im1 = single(im0(xROId,yROId)); %crop image
        fsh2(:,:,k) = im1; 
    end

    fshs = imgaussfilt(fsh2,smooth); %3D image stack smoothing
    fshm = immultiply(fshs,nucf); %Apply nuclear mask to channel 
    m = quantile(fshm(fshm>0&fshm<maxlevel),0.5); %median of image stack
    d = quantile(fshm(fshm>0&fshm<maxlevel),0.75)-m; %Third quartile of image stack
	low = m + fshlevel*d; %low and clipping threshold for image contrast
    high = max(fshm(fshm>0&fshm<maxlevel)); %high for image contrast
    data.fsh2low = low;
    data.fsh2high = high;
    if high>low
        fshadj = imadjustn(fshm./high,[low/high 1], [0 1]);   % contrast adjustment
        fshbin= zmat;
        if choosefilter > 0
            for k = 1:size(fshadj,3)    
                im2 = imfilter(fshadj(:,:,k), laplace);  % laplacian filter
                im3 = (im2 < 0);   % negative values of the laplacian
                fshbin(:,:,k) = im3;
            end
        else
            fshbin = imbinarize(fshadj);
        end

        if edgedel > 0
            fshedge = logical(imadd(fshbin,nucedge));  %Adds nuclear periphery
            fshsel = fshedge - bwareaopen(fshedge,2000);
                %removes periphery with attaching clusters
            fsh2f = bwareaopen(fshsel,minvol,18);
        else
            fsh2f = bwareaopen(fshbin,minvol,18);
        end
    else
        fsh2f = zmat; %When the channel is empty, this avoids the script from crashing   
    end
    clear fshsel fshedge fshbin im3 im2 im1 fshadj high low m d fshm fshs chb
    Lfsh = bwlabeln(fsh2f,18);
    Rfsh = regionprops3(Lfsh,fsh2,'Volume','MeanIntensity','MaxIntensity','Centroid');
       data.fsh2v = single([Rfsh.Volume]);
       data.fsh2mean = single([Rfsh.MeanIntensity]);
       data.fsh2max = single([Rfsh.MaxIntensity]);
       data.fsh2c = single([Rfsh.Centroid]);
    Rhminfsh = regionprops3(Lfsh,hm1,'MeanIntensity');
       data.hm1infsh2mean = single([Rhminfsh.MeanIntensity]);
       % IF channel 1 fluorescence in FISH 2 channel clusters
    Rhminfsh = regionprops3(Lfsh,hm1f,'MeanIntensity');
       data.hm1infsh2frov = single([Rhminfsh.MeanIntensity]);
       % IF channel 1 cluster fraction of overlap with FISH 2 channel clusters
    clear Rhminfsh hfsh Lhfsh Rhfsh
    if h2ch>0
        Rhminfsh = regionprops3(Lfsh,hm2,'MeanIntensity');
            data.hm2infsh2a = single([Rhminfsh.MeanIntensity]);
            % IF channel 2 fluorescence in FISH 2 channel clusters
        Rhminfsh = regionprops3(Lfsh,hm2f,'MeanIntensity');
            data.hm2infsh2frov = single([Rhminfsh.MeanIntensity]);
            % IF channel 2 cluster fraction of overlap with FISH 2 channel clusters
        clear Rhminfsh hfsh Lhfsh Rhfsh
    end
    clear Lfsh Rfsh fsh2 fsh2f
end
%% Segmenting FISH channel 3
if f3ch>0
    fshlevel = fsh3level;
    zfsh = zfsh3(logical(zpick));
    fsh3 = zmat;

    for k = 1:length(zfsh)
        im0 = imread(fname,'tif','Index',zfsh(k)); %load image
        im1 = single(im0(xROId,yROId)); %crop image
        fsh3(:,:,k) = im1; 
    end

    fshs = imgaussfilt(fsh3,smooth); %3D image stack smoothing
    fshm = immultiply(fshs,nucf); %Apply nuclear mask to channel 
    m = quantile(fshm(fshm>0&fshm<maxlevel),0.5); %median of image stack
    d = quantile(fshm(fshm>0&fshm<maxlevel),0.75)-m; %Third quartile of image stack
	low = m + fshlevel*d; %low and clipping threshold for image contrast
    high = max(fshm(fshm>0&fshm<maxlevel)); %high for image contrast
    data.fsh3low = low;
    data.fsh3high = high;

    if high>low
        fshadj = imadjustn(fshm./high,[low/high 1], [0 1]);   % contrast adjustment
        fshbin= zmat;
        
        for k = 1:size(fshadj,3)    
            im2 = imfilter(fshadj(:,:,k), laplace);  % laplacian filter
            im3 = (im2 < 0);   % negative values of the laplacian
            fshbin(:,:,k) = im3;
        end

        if edgedel > 0
            fshedge = logical(imadd(fshbin,nucedge)); %Adds nuclear periphery
            fshsel = fshedge - bwareaopen(fshedge,2000);
                 %removes periphery with attaching clusters
            fsh3f = bwareaopen(fshsel,minvol,18);
        else
            fsh3f = bwareaopen(fshbin,minvol,18);
        end
    else
        fsh3f = zmat;
    end
    clear fshsel fshedge fshbin im3 im2 im1 fshadj high low m d fshm fshs chb
    
    Lfsh = bwlabeln(fsh3f,18);
    Rfsh = regionprops3(Lfsh,fsh3,'Volume','MeanIntensity','MaxIntensity','Centroid','VoxelValues');
       data.fsh3v = single([Rfsh.Volume]);
       data.fsh3mean = single([Rfsh.MeanIntensity]);
       data.fsh3max = single([Rfsh.MaxIntensity]);
       data.fsh3c = single([Rfsh.Centroid]);
    Rhminfsh = regionprops3(Lfsh,hm1,'MeanIntensity','VoxelValues');
       data.hm1infsh3mean = single([Rhminfsh.MeanIntensity]);
      % IF channel 1 fluorescence in FISH 3 channel clusters
    Rhminfsh = regionprops3(Lfsh,hm1f,'MeanIntensity');
       data.hm1infsh2frov = single([Rhminfsh.MeanIntensity]);
       % IF channel 1 cluster fraction of overlap with FISH 3 channel clusters
    clear Rhminfsh hfsh Lhfsh Rhfsh

    if h2ch>0
        Rhminfsh = regionprops3(Lfsh,hm2,'MeanIntensity');
            data.hm2infsh3a = single([Rhminfsh.MeanIntensity]);
            % IF channel 2 fluorescence in FISH 3 channel clusters
        Rhminfsh = regionprops3(Lfsh,hm2f,'MeanIntensity');
            data.hm2infsh3frov = single([Rhminfsh.MeanIntensity]);
            % IF channel 2 cluster fraction of overlap with FISH 3 channel clusters
        clear Rhminfsh hfsh Lhfsh Rhfsh
    end
    clear Lfsh Rfsh fsh3 fsh3f
end
Lhm = bwlabeln(hm1f,18);
Rhm = regionprops3(Lhm,hm1,'Volume','MeanIntensity','Centroid');
    data.hm1v = single([Rhm.Volume]);
    data.hm1mean = single([Rhm.MeanIntensity]);
    data.hm1c = single(Rhm.Centroid);
    if h2ch>0
        Rhm2 = regionprops3(Lhm,hm2,'MeanIntensity');
        data.hm2inhm1mean = single([Rhm2.MeanIntensity]);
        % IF channel 2 fluorescence in IF channel 1 clusters
        Rhm2 = regionprops3(Lhm,hm2f,'MeanIntensity');
        data.hm2inhm1frov = single([Rhm2.MeanIntensity]);
        % IF channel 2 cluster fraction of overlap with IF channel 1 clusters
    end
clear Lhm Rhm2 rdhm1a rdhm2a Rhm

if h2ch>0
    Lhm = bwlabeln(hm2f,18);
    Rhm = regionprops3(Lhm,hm2,'Volume','MeanIntensity','Centroid');
        data.hm2v = single([Rhm.Volume]);
        data.hm2mean = single([Rhm.MeanIntensity]);
        data.hm2c = single(Rhm.Centroid);
    Rhm2 = regionprops3(Lhm,hm1,'MeanIntensity');
        data.hm1inhm2mean = single([Rhm2.MeanIntensity]);
        % IF channel 2 fluorescence in IF channel 1 clusters
    Rhm2 = regionprops3(Lhm,hm1f,'MeanIntensity');
    data.hm1inhm2frov = single([Rhm2.MeanIntensity]);
    % IF channel 2 cluster fraction of overlap with IF channel 1 clusters
end
data.zROI = zROId;
data.yROI = yROId;
data.xROI = xROId;

