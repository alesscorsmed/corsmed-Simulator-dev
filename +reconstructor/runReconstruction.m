function [kSpace, iSpace] = runReconstruction(...
    rawData, sensData, encodingData, expControl )
%
% RECONSTRUCTOR.RUNRECONSTRUCTION
%
%	Reconstruction interface.
%
% INPUT
%   timeSolution    struct with solution struct with initial data
%   acquisition     structure with acq info
%   coilSystem      struct with coil info (for SENSE)
%   expControl      experiment control struct
%
% OUTPUT
%   reconData       struct with K-space and image
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'reconstructor.runReconstruction';
if (nargin < 4)
    ME = MException('reconstructor:wrongArgCount',...
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

%% assemble k-space
% timeSolution.noise = zeros(size(timeSolution.Sx));

% try with mapping, if it fails, use old assemble
try
    [kSpace] = reconstructor.signal.mapKspace(rawData,...
        encodingData, expControl );
    kSpace = permute(kSpace, [1,2,3,5,4]);
catch
    disp('ERROR MAPPER');
    [kSpace] = reconstructor.signal.assembleKspace(rawData,...
        encodingData, expControl );
end

%% apply propper recon
switch lower(encodingData.reconstructor)
    
    case 'coremri_fft'
        % Sum-Of-Squares cartesian FFT based
        [iSpace] = reconstructor.cartesianFFT.sosFFT( kSpace,...
            encodingData, expControl );
        
    case 'fft'
        % basic FFT, generating an image per channel
        [iSpace] = reconstructor.cartesianFFT.basicFFT( kSpace,...
            encodingData, expControl );
        
    case 'sos'
        % basic FFT, generating an image per channel
        [iSpace] = reconstructor.cartesianFFT.sosFFT( kSpace,...
            encodingData, expControl );
        
    case 'sense'
        % apply SENSE 2D
        [iSpace] = reconstructor.sense.sense2D( kSpace, sensData, ...
            encodingData, expControl);
      
    case 'opt'
        % Apply sense recon w/ no acceleration: remove Coil Sens weighting
        [iSpace] = reconstructor.sense.sense2D( kSpace, sensData, ...
            encodingData, expControl);
        
    otherwise
        ME = MException('reconstructor:wrongMethod',...
            '%s : unknown reconstructor %s', ...
            functionName, encodingData.reconstructor );
        throw(ME);
        
end

%% report
if  expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : %s reconstruction done ',...
        functionName, encodingData.reconstructor);
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end