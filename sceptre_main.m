
%main script to process images in a folder for SCEPTRE

projectname = 'Immunostained H3K4me3(AF568) and H3K27me3 (AF488) at FISH-labeled GAPDH (AT647N)(RPE)';

%% Input vaues for sceptre_processing
flist = struct2cell(dir('647*.*')); %acquires info list of images beginning
                                    %with assigned characters.
namelist = flist(1,:);
ch = 3; % number of channels in image stack
imorder = 0;    %if image stack order is z than channel, imorder = 0. 
                %If channel than z, imorder > 0.
fstart = 2; %choose which image in z to start with for processing.
smooth = 1; % # of standard deviations for 3x3 gaussian filter of channels
bnd_nuc = 3; %radius of sphere for morphological opening and closing.
edgedel = 1;    %If edgedel > 0, segmented FISH clusters contacting  
                %periphery of nuclear mask are removed.
maxlevel = 0.1*(2^16);  %Eliminates bias by hot pixels in image contrast 
                        %step e.g. contamination in image stack.
choosefilter = 'laplace'; %Selecte if 'Otsu' or 'laplace' binarizing method
nuch = 3;   %Channel used to generate nuclear mask. 
            %if 0, h1ch is used as nuclear mask channel
h1ch = 2;   %immunofluorescence channel 1, e.g., H3K4me3 
h2ch = 3;   %immunofluorescence channel 2, e.g., H3K27me3. 
            %If there is no immunofluorescence channel 3, h2ch = 0
f1ch = 1;   %DNA-FISH channel 1, e.g. GAPDH. 
            %If there is no FISH channel 1, f1ch = 0
f2ch = 0;   %DNA-FISH channel 2, e.g. HOXC. 
            %If there is no FISH channel 2, f2ch = 0
f3ch = 0;   %DNA-FISH channel 3, e.g. LINC-PINT. 
            %If there is no FISH channel 3, f3ch = 0
nuclevel = 3;   %threshold level for nuclear mask clip. 
hm1level = 3;   %threshold level for immunofluorescence channel 1 clip.  
hm2level = 2;   %threshold level for immunofluorescence channel 2 clip.
fsh1level = 15; %threshold level for DNA FISH channel 1 clip.
fsh2level = 0;  %threshold level for DNA FISH channel 2 clip. 
fsh3level = 0;  %threshold level for DNA FISH channel 3 clip. 
minvol = 20; %minimum volume for segmented clusters in channel
imdimx = 750; %x dimension of selected centered image stack region
imdimy = 750; %y dimension of selected centered image stack region 

%% sceptre_processing for batch of images

lnm = length(namelist);
h = waitbar(0,strcat('ENACTING SCEPTRE PROCESSING: ',' 1',' of ',{' '},num2str(lnm)));
setappdata(h,'canceling',0);
for k = 1:lnm
    if getappdata(h,'canceling')
        break
    end
    fname = namelist{k};
    r(k) = sceptre_processing(fname,ch,imorder,fstart,smooth,bnd_nuc,edgedel,maxlevel,choosefilter,nuch,h1ch,h2ch,f1ch,f2ch,f3ch,nuclevel,hm1level,hm2level,fsh1level,fsh2level,fsh3level,minvol,imdimx,imdimy);
    step = strcat('ENACTING SCEPTRE PROTOCOL: ',{' '},num2str(k+1),' of',{' '},num2str(lnm));
    waitbar(k/lnm,h,step)
    pause(0.1)

end
close(h)

