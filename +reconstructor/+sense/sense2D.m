function [senseImage] = sense2D(kSpace, coilSens, encodingData, expControl)
%
% RECONSTRUCTOR.SENSE.SENSE2D
%
%   Applies SENSE reconstruction on 2D slice
%
%     Calculating coil sensitivities are the initial and most important step 
%     in the SENSE process. Low-resolution images are acquired separately from 
%     each surface coil at full field-of-view. These surface coil images are 
%     normalized by dividing them by a low-resolution body coil image. 
%     Filtering, thresholding, and point estimation are then applied to the 
%     data to generate coil sensitivity maps. These maps quantify 
%     the relative weighting of signals from different points of origin within 
%     the reception area of each coil.
%     (http://mriquestions.com/senseasset.html)
%     Implementation based on the paper Blaimer et al., SMASH, SENSE, PILS,
%     GRAPPA, Top Magn Reson Imaging 2004 Aug;15(4):223-36 and the
%     corresponding test code: 
%     https://github.com/veritas9872/Medical-Imaging-Tutorial
%
%
% INPUT
%     model         struct with slice anatomical model data 
%     rotMat        rotation matrix
%     refPoint      reference point for rotation
%     coilSystem    struct with coils data
%     expControl    control info
%
% OUTPUT
%     model    updated model struct with interpolated data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'reconstructor.sense.sense2D';
if (nargin < 4)
    ME = MException('Sense:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

%% open debugging
if expControl.debug.debugMode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(expControl.debug.debugFile,'a');
    catch
        fid = 1;
    end
    tTotal = tic();
    fprintf(fid, '\n%s : start', functionName);
end

%% extract encoding plan and FT operators
encodingPlan = encodingData.plan;
try
    iFTZ    = encodingData.operator.iFTZ;
    iFT     = encodingData.operator.iFT;    
catch
    iFTZ    = @(kSpace) fftshift(ifftn(ifftshift(kSpace),[],3));
    iFT     = @(kSpace) fftshift(ifftn(ifftshift(kSpace)));
end

%% Acceleration factor
Rfactor = encodingPlan.rFactorPE;

%% get compressed K-space from sparse one
% note, original K-Space is numFE, numPE, numSE, numCH, numCO
compressedKSpace = kSpace(:,1:Rfactor:end,:,:,:);
% get compressed dimensions:
[numFE,numPE,numSE,numCH,numCO] = size(compressedKSpace);

%% cropping image to reduce complexity
[coilX,coilY,coilZ,~] = size(coilSens);
kspacePadding = floor(([numFE,numPE,numSE]-[coilX,coilY,coilZ])/2);
% find the indexes in the extended coil maps
idxX = kspacePadding(1)+1:kspacePadding(1)+coilX;
idxZ = kspacePadding(3)+1:kspacePadding(3)+coilZ;

%% the sense image is going to be full-size
% take into account extension of foldover suppresion
%  scale size of coilY, and coil incidence indexes
if encodingPlan.factorY > 1
    idxStart = floor(coilY*(encodingPlan.factorY-1)/2);
    idxCoilY = idxStart+1:idxStart+coilY;
    coilY = encodingPlan.factorY*coilY;
else
    idxCoilY = 1:coilY;
end
%% use full (foldover extended) size: crop later
senseImage = zeros(coilX,coilY,coilZ,1,numCO);

%% loop on contrasts and slices
for ic = 1:numCO
    % note for these we swap the SE and CH for easier use in SENSE
    cMapsPerCH  = zeros(coilX,coilY,numCH,numSE);
    imagePerCH  = zeros(coilX,numPE,numCH,numSE);
    % cropp and assign
    for cc = 1:numCH
        % apply IFT to each slice independently
        for sl = 1:numSE
            if encodingPlan.is3D && numSE > 1
                % for 3D sequences
                % assuming cartesian encoding, undo Z dimension first FFT
                kSpace2D = iFTZ(compressedKSpace(:,:,:,cc,ic));
                % apply 2D iFT on each remaining k-Spaces
                imageCh = iFT(kSpace2D(:,:,idxZ(sl)));
            else
                % iFT to get image of channel from compressed K-space
                imageCh = iFT(compressedKSpace(:,:,idxZ(sl),cc,ic));
            end
            % crop image
            imagePerCH(:,:,cc,sl) = imageCh(idxX,:);
        end
        % corresponding coil map (3D)
        % NOTE: we reverse the Y-direction to be consistent with Image
        % cMapsPerCH(:,:,cc,:) = Cx(:,end:-1:1,:,cc) - 1j*Cy(:,end:-1:1,:,cc);
        cMapsPerCH(:,idxCoilY,cc,:) = coilSens(:,:,:,cc);

    end    

    %% recon each slice using SENSE
    rPeriod = numPE;
    peStart = round(Rfactor*numPE/2 - numPE/2);
    for sl = 1:numSE

        % % plot to illustrate
        % figure();
        % for cc = 1:numCH
        %     subplot(4,ceil(numCH/2),cc);
        %     imagesc(abs(cMapsPerCH(:,:,cc,sl)));
        %     colormap hot; colorbar; title(sprintf('Abs(C_%d)', cc));
        %     subplot(4,ceil(numCH/2),cc + numCH);
        %     imagesc(abs(imagePerCH(:,:,cc,sl)));
        %     colormap hot; colorbar; title(sprintf('Image(C_%d)', cc));
        % end

        % slice
        senseSlice = zeros(coilX,coilY);
        mapsPerCoilElement = cMapsPerCH(:,:,:,sl);
        imagesPerCoilElement = imagePerCH(:,:,:,sl);
        % get pixel image
        for x=1:coilX
            for y=1:rPeriod
                % find the indexes in the full domain 
                rIdx = (y-1)+peStart+(1:rPeriod:rPeriod*Rfactor);
                % adjust cyclic
                rIdx = mod(rIdx-1,coilY)+1; % make it cyclic
                cHat = reshape(mapsPerCoilElement(x,rIdx,:),[],numCH).';
                % solve inverse
                if nnz(cHat)
                    vectorI = reshape(imagesPerCoilElement(x,y,:),numCH,1);
                    senseSlice(x,rIdx) = pinv(cHat)*vectorI;
                end
            end
        end
        % assign to the overall SENSE image
        senseImage(:,:,idxZ(sl),1,ic) = senseSlice;
        
        % % plot slice for illustration
        % figure(); 
        % imagesc(abs(senseSlice));
        % colormap gray; colorbar; 
        % title(sprintf('SENSE IMAGE Slice %d)', sl));

    end

end

%% crop the image if fold over suppresion is applied
if encodingPlan.factorY > 1
    senseImage = senseImage(:,idxCoilY,:,:,:);
end

%% report
if expControl.debug.debugMode
    % correct sizes
    [numX,numY,numZ,~,numCO] = size(senseImage);
    fprintf(fid, ...
        '\n%s : 2D SENSE recon done, elapsed time %.3fs',...
        functionName, toc(tTotal));
    fprintf(fid, '\n  R-Factor       %d', Rfactor);
    fprintf(fid, '\n  I-Space size   %3d x %d x %d x %d',...
        numX, numY, numZ, numCO);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end
