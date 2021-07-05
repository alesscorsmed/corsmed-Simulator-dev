function [contrastSpace,contrastText,isContrast] = processMOLLI(...
    multiContrastSlice, imageData, MOLLIspecs, expControl)
%
% IMAGE.DICOM.PROCESSMOLLI
%
% Process the generated images so as to generate the T1 map:
%   single-shot bSSFP images
%   T1 map
% 
%
% INPUT
%
% OUTPUT
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'image.dicom.processMOLLI';

% info for debugging
if expControl.debug.debugMode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(expControl.debug.debugFlie,'a');
    catch
        fid = 1;
    end
    tTotal = tic();
    fprintf(fid, '\n%s : start', functionName);
end

for contrastNum=size(multiContrastSlice,3):-1:1
    contrastText{contrastNum} = imageData.contrast{contrastNum}.imageText.LTop{3};
end
contrastText{11} = ['T1map ',MOLLIspecs.textScheme];

%% is contrast or map
isContrast       = ones(11,1);
isContrast(11)    = 0;

contrastSpace = abs(multiContrastSlice);

%% Data preparation
numCO           = size(contrastSpace,3);
sizeImage       = size(contrastSpace(:,:,1));
totalPixels     = sizeImage(1)*sizeImage(2);
signalPerPixel  = zeros(totalPixels,numCO);

for iCO = 1:numCO
    imageCO                 = contrastSpace(:,:,iCO);
    signalPerPixel(:,iCO)   = imageCO(:);
end

IRtime = repmat(MOLLIspecs.TIs*1e3,size(signalPerPixel,1),1);
k0 = [signalPerPixel(:,size(signalPerPixel,2)), 1200*ones(size(signalPerPixel,1),1), 2*signalPerPixel(:,size(signalPerPixel,2))];

M = 1;
xdata = IRtime(:,M:end);
ydata = signalPerPixel(:,M:end);
%% Normalize data
max_ydata = max(ydata,[],1)';
for i = 1:size(ydata,2)
    ydata(:,i) =  ydata(:,i) ./ max_ydata(i);
end

%% Masking
% choose highest signal for mask
thres = 1e-1;
yy = ydata(:,end);
mask = find(yy > thres);

%% Initial values
% keeping in mind the fitting function:
% F(x) = abs(x(1) - x(3) .* exp(-xdata./x(2)))
ys = size(ydata,1);
x0 = [ydata(:,size(ydata,2)), 900*ones(size(ydata,1),1), 2*ydata(:,size(ydata,2))];
%% Optimization
options = optimoptions(@lsqcurvefit, 'Algorithm', 'trust-region-reflective',...
        'MaxIter',30,'Display','off','SpecifyObjectiveGradient',true);

% use box constraints
lb = [0,0,0];
ub = [2,2000, 4];

% allocation of solution vector
T1 = zeros(totalPixels,1);

% Use only masked pixels in parfor
x0Par = x0(mask,:);
xdataPar = xdata(mask,:);
ydataPar = ydata(mask,:);
T1par = zeros(numel(mask),1);

poolobj = gcp;
addAttachedFiles(poolobj,{'+image/+dicom/+processMOLLIfunctions/fitFun.m'})

tic
parfor ii = 1:numel(mask)
    x = lsqcurvefit(@image.dicom.processMOLLIfunctions.fitFun, x0Par(ii,:),...
        xdataPar(ii,:), ydataPar(ii,:),lb, ub, options);

    T1par(ii, 1) = x(2);
end
toc

T1(mask,1) = T1par(:);

% revert to original dimensions
contrastSpace(:,:,11) = reshape(T1,sizeImage);

%% final message
if  expControl.debug.debugMode
    figure();imagesc(contrastSpace(:,:,11),[300 2000])
    colormap jet
    
    fprintf(fid, ...
        '\n%s : elapsed time %.3fs',...
        functionName, toc(tTotal));    
end