% Main development script for HN25 registration with pTVreg

%% 1. Load both of the autofluorescence channel images
clear fileNames t afImgs pemtAF immnAF;
fileNames = {'../HN25_slice1.tif', '../HN25_slice2.tif'};
for n = 1:2
	t = Tiff(fileNames{n});
	% MATLAB indexing starts at 1 but TIFF specification starts at 1
	% The AF channel is at IFD 8
	t.setDirectory(32);
	% Read the frame
	rawImg= read(t);
	% Frame type class is incorrect in the metadata, it should be single
	afImgs{n} = reshape(typecast(rawImg(:),'single'), size(rawImg));
end

% Images must be the same size: pad with 0s
maxRC = max([size(afImgs{1});size(afImgs{2})],[],1);

pemtAF = padarray(afImgs{1},floor((maxRC-size(afImgs{1}))/2),'replicate','pre');
pemtAF = padarray(pemtAF,ceil((maxRC-size(afImgs{1}))/2),'replicate','post');
immnAF = padarray(afImgs{2},floor((maxRC-size(afImgs{2}))/2),'replicate','pre');
immnAF = padarray(immnAF,ceil((maxRC-size(afImgs{2}))/2),'replicate','post');

% Normalize and compress intensities?

% Clear large objects
clear rawImg afImgs;

%% 2. Perform a rough rigid registration with phase correlation
% HN25 was manually rotated so not necessary this time

%% 3. Perform registration of autofluorescence channels with pTVreg
addpath(genpath('.'));

%  -- Registration options -- 

clear opts;
% Loss function metric: Local correlation coefficient. This seems to be
% the best. metric_param is the Gaussian kernal size. A size 7 or smaller
% is recommended. Larger kernels might be slower? 
opts.metric = 'loc_cc_fftn';
opts.metric_param = 7;
% Use an approximation to speed up solver (perhaps less accurate?)
opts.loc_cc_approximate = true;

% Downsampling factor for multi-level registration (0,1). A larger value
% downsamples by a smaller factor between levels... perhaps more accurate
% but slower. Typically use 0.5.
opts.k_down = 0.5;

% Spacing of displacement knots in (x,y) directions in pixels. Should be
% the same and a factor of 2.
opts.grid_spacing = [8 8];

% Maximum number of iterations at each level. More iterations are required
% for bigger deformations.
opts.max_iters = 50;

%  -- End of registration options -- 

% Run registration routine.
% The first image (pEMT AF channel) is deformed to match the second (Immun.
% AF channel).
% Tmin_out describes the deformation.
[pemtDef, Tmin_out] = ptv_register(pemtAF,immnAF, opts);

% -- Display registration -- 
% Show registered frames together
figure; imshowpair(pemtDef, immnAF,'falsecolor','Scaling','independent');

% Show displacment matrix
figure; plot_wiremesh(squeeze(Tmin_out(:,:,1,:)),4,'-r');

%% 4. Warp, pad and save full-size images
% Padding dimensions for full-size images
% scale factor between the scale used in the deformation and the full-size
% images
scale_factor = 4;
levelOffset = 8; 
targetRC = ceil(maxRC.*scale_factor);

% 4.1 Create full-size deformation matrix
fullTMin = zeros([targetRC(1),targetRC(2),1,2],'single');
fullTMin(:,:,1,1) = imresize(single(Tmin_out(:,:,1,1)),scale_factor).*scale_factor;
fullTMin(:,:,1,2) = imresize(single(Tmin_out(:,:,1,2)),scale_factor).*scale_factor;

% 4.2 Create new TIFF
tOut = Tiff(['HN25_composite_halfscale.tif'],'w8');

% TODO: figure out how to label the channels in a way HALO will recognize.
channelNames = {'DAPI','Opal520','Opal540','Opal570','Opal620','Opal650','Opal690','AF'};
	
tIn = Tiff(fileNames{1},'r');
clear tags;
tags.ImageLength = size(fullTMin,1);
tags.ImageWidth = size(fullTMin,2);
tags.SampleFormat = Tiff.SampleFormat.IEEEFP;
tags.Photometric = Tiff.Photometric.MinIsBlack;
tags.BitsPerSample = 32;
tags.SamplesPerPixel = 1;
tags.RowsPerStrip = 16;
tags.Compression = getTag(tIn,'Compression');
tags.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tags.Orientation = Tiff.Orientation.TopLeft;
tags.ResolutionUnit = getTag(tIn,'ResolutionUnit');
tags.XResolution = getTag(tIn,'XResolution')/4;
tags.YResolution = getTag(tIn,'YResolution')/4;

progressbar(0);
% 4.3 write warped Immne channels
for channel = 1:8
	% Load each channel
	tIn.setDirectory(channel+levelOffset);
	img = read(tIn);
	img = reshape(typecast(img(:),'single'), size(img));
	
	% Pad to correct size
	paddedImg = padarray(img, floor((targetRC-size(img))/2), 'replicate','pre');
	paddedImg = padarray(paddedImg, ceil((targetRC-size(img))/2), 'replicate','post');
	
	% Do warp
	deformedImg = ptv_deform(paddedImg, fullTMin);

	% Write this channel
	if channel ~= 1
		tOut.writeDirectory;
	end
	tOut.setTag(tags);
	write(tOut,deformedImg);
	
	progressbar(channel/16);
end

% 4.4 write unwarped pEMT channels
tIn = Tiff(fileNames{2},'r');
for channel = 1:8
	% Load each channel
	tIn.setDirectory(channel+levelOffset);
	img = read(tIn);
	img = reshape(typecast(img(:),'single'), size(img));
	
	% Pad to correct size
	paddedImg = padarray(img, floor((targetRC-size(img))/2), 'replicate','pre');
	paddedImg = padarray(paddedImg, ceil((targetRC-size(img))/2), 'replicate','post');

	% Write this channel
	tOut.writeDirectory;
	tOut.setTag(tags);
	write(tOut,paddedImg);
	
	progressbar(channel/16+.5);
end


close(tOut);
disp('Complete!');


%% We can modify Metadata aftewards if necessary...
% tOut = Tiff('HN25_immn_warped_test.tif','r+');

