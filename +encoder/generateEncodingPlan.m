function [encoding] = generateEncodingPlan(acquisition, expControl)
%
% ENCODER.GENERATEENCODINGPLAN
%
%	Generates K-space encoding information from the acquisition.
%   This info is used to assemble the correct K-space.
%
% INPUT
%   acquisition         
%   expControl      
%
% OUTPUT
%   encoding   struct with K-space info to assemble the K-space
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'encoder.generateEncodingPlan';

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

%% generate the K-space info for the reconstruction
if acquisition.data.isCartesian
    
    %% generate encoding plan
    [encodingPlan] = encoder.cartesian.encodingPlanner(...
        acquisition.data, expControl);
    
    %% generate incidence matrices to map from time to k-space
    [encodingMap] = encoder.cartesian.encodingMapperKindex(...
        encodingPlan, expControl);
    
    %% create the Fourier Transform operators
    operators.iFT   = @(kSpace) fftshift(ifftn(ifftshift(kSpace)));
    operators.fFT   = @(iSpace) fftshift( fftn( fftshift(iSpace)));
   	operators.iFTZ  = @(kSpace) fftshift(ifftn(ifftshift(kSpace),[],3));
    operators.fFTZ  = @(iSpace) fftshift( fftn( fftshift(iSpace),[],3));
    
else
    % unknown sequence, throw
    ME = MException('sequence:wrongFamily',...
        '%s : non-cartesian or advanced encodings not supported',...
        functionName);
    throw(ME);
end

%% store in struct and return
encoding.plan           = encodingPlan;
encoding.map            = encodingMap;
encoding.operators      = operators;
encoding.reconstructor  = acquisition.data.reconstructor;

%% final message
if expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done, elapsed time %.3fs',...
        functionName, toc(tTotal));
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end
