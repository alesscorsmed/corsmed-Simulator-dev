function [kSpace,iSpace] = sbrRecon(timeSolution,encoding,expControl)

% simJobs, spinModel, coilSystem, ...
%     acquisition, expControl, sessionData, pulseSequence, mrSystem)
%
% sbr.run.sbrRecon
%
%	Runs the reconstruction.
%
% INPUT
%
% OUTPUT
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'sbr.run.sbrRecon';
% if (nargin < 5)
%     ME = MException('eduTool:wrongArgCount',...
%         '%s : wrong argument count',functionName);
%     throw(ME);
% end
% 
% %% info for debugging
% if expControl.debug.debugMode
%     try % open file if possible, otherwise dump to stdout
%         fid = fopen(expControl.debug.debugFile,'a');
%     catch
%         fid = 1;
%     end
%     % time it
%     tTotal = tic();
%     fprintf(fid, '\n%s : start', functionName);
% end

%% if there is no DB connection
if ~isfield(expControl,'connLocalDB')
    expControl.connLocalDB = [];
end

%% create the Fourier Transform operators
reconData.operators.iFT = @(kSpace) fftshift(ifftn(ifftshift(kSpace)));
reconData.operators.FT  = @(iSpace) fftshift( fftn( fftshift(iSpace)));

%% assemble kSpace
[kSpace] = reconstructor.signal.assembleKspace(timeSolution,...
    encoding, expControl);

%% final message
if expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : done for experiment %d',...
        functionName, expControl.experimentID);
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n  Number of Slices  %d', spinModel.numSlices);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end