function [acquisition] = generateEncodingInfo(acquisition, expControl)
%
% SEQUENCE.TOOLS.GENERATEENCODINGINFO
%
%	Generates K-space encoding information from the acquisition.
%   This info is used to assemble the correct K-space.
%
% INPUT
%   acquisition         
%   expControl      
%
% OUTPUT
%   kSpaceInfo   struct with K-space info to assemble the K-space
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'sequence.tools.generateEncodingInfo';

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


%% generate the K-Space info for the reconstruction
if acquisition.data.isCartesian
    [kSpaceInfo] = sequence.tools.cartesianEncodingInfo(...
        acquisition.data, expControl);
    acquisition.kSpaceInfo = kSpaceInfo;
else
    % unknown sequence, throw
    ME = MException('sequence:wrongFamily',...
        '%s : non-cartesian sequences not supported',...
        functionName);
    throw(ME);
end

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
