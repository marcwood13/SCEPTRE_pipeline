function data = sceptre_rbgenerator(fname,maxlevel,nuclevel,choosefilter,xROId,yROId,zROId,ch,nuch,h1ch,h2ch,imorder,fstart,smooth,boxv,bn,bnd_nuc,imdimx,imdimy)

%% Step1: Obtain image info and XYZ boundaries

info = imfinfo(fname);
L = length(info);
imwidth = info(1).Width;
imheight = info(1).Height;

xrange = imdimx;
xROI = ((imwidth - imdimx)/2)+1:((imwidth + imdimx)/2);

yrange = imdimy;
yROI =  ((imheight - imdimy)/2)+1:((imheight + imdimy)/2);

zrange = (L/ch)+1-fstart;

zmat = zeros(xrange,yrange,zrange); %zero matrix for current ROI
laplace = fspecial('laplacian',0.2); %laplace filter

%% Selecting image Z indexes

if imorder>0
    
    zhm1 = h1ch+ch*(fstart-1):ch:L-(ch-h1ch);  %z index for h1ch
    
    if h2ch>0
        zhm2 = h2ch+ch*(fstart-1):ch:L-(ch-h2ch);  %z index for h2ch
    end
    
    if nuch>0
        znuc = nuch+ch*(fstart-1):ch:L-(ch-nuch); %z index for  nuch
    end
else

    zhm1 = (L*(h1ch-1)/ch)+fstart:L*h1ch/ch; %z index for  h1ch
    
    if h2ch>0
        zhm2 = (L*(h2ch-1)/ch)+fstart:L*h2ch/ch; %z index for  h2ch
    end
    
    if nuch>0
        znuc = (L*(nuch-1)/ch)+fstart:L*nuch/ch; %z index for  nuch
    end
end

%% Determining nuclear area and mask
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

%% Limiting the image regions to where the nucleus is

nucf = nucf0(xROId-min(xROI)+1,yROId-min(yROI)+1,zROId);
clear nucs nucadj nucbw nucbw nucp mnuc volnuc nucsel m d low high nucl thrsh nucf0 mip mx mx2 my my2 nucbin
zmat = zeros(length(xROId),length(yROId),length(zROId));
hm1 = zmat;
zhm1 = zhm1(zROId);

%% loding Immunofluorescence channel images
for k = 1:length(zhm1)
     im0 = imread(fname,'tif','Index',zhm1(k)); %load  IF channel image
     im1 = single(im0(xROId,yROId)); %crop image
     hm1(:,:,k) = im1;
end

if h2ch>0
    hm2 = zmat;
    zhm = zhm2(zROId);

    for k = 1:length(zhm)
         im0 = imread(fname,'tif','Index',zhm(k)); %load IF channel image
         im1 = single(im0(xROId,yROId)); %crop image
         hm2(:,:,k) = im1;
    end
end

%% Generating random boxes through the nucleus

srd = double(round((boxv)^(1/3)));
xr = length(xROId);
yr = length(yROId);
zr = length(zROId);

x1 = round((xr-srd)*rand(1,bn));
y1 = round((yr-srd)*rand(1,bn));
z1 = round((zr-srd)*rand(1,bn));

s1 = zeros(xr,yr,zr);

for h = 1:bn
    for i = 1:srd
        for j = 1:srd
            for k = 1:srd
                s1((x1(h)+i),(y1(h)+j),(z1(h)+k)) = 1;
            end
        end
    end
end

    rd0 = immultiply(s1,nucf); % Applying nuclear mask to random regions
    rd = bwareaopen(rd0,round((srd^3)-1),18); 
    % Selecting clusters that were not cropped by mask
    rd = rd - bwareaopen(rd0,round((srd^3)+1),18);
    % selecting clusters that are not overlapping

    Lrd = bwlabeln(rd,18);
Rrdhm1 = regionprops3(Lrd,hm1,'MeanIntensity','Volume','Centroid');
    data.rdhm1mean = single([Rrdhm1.MeanIntensity]);
    data.rdhm1v = single([Rrdhm1.Volume]);
    data.rdhm1c = single([Rrdhm1.Centroid]);
    
if h2ch>0
    Rrdhm2 = regionprops3(Lrd,hm2,'MeanIntensity');
        data.rdhm2mean = single([Rrdhm2.MeanIntensity]);

end

