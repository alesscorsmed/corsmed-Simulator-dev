function [repetition] = applyPartialEcho(...
    repetition, startFE, numEchoes, expControl)
%
% SEQUENCE.TOOLS.APPLYPARTIALECHO
%
%	Applies partial echo by zeroing rx of the repetition.
%
% INPUT
%   repetition  original repetition  
%   startFE     index at which the actual FE starts (before this, RXs are zeroed)
%   numEchoes   number of repetitions of FE (for multi-echo or EPI)
%   expControl      
%
% OUTPUT
%   repetition  modified repetition
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'sequence.tools.applyPartialEcho';

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


% find the indexes of the positions of the Sampling Points
idxRx = find( repetition.rxSignal > 0);
% reshape the indexes: numFE x numEchoes
idxRx = reshape(idxRx,[],numEchoes);
% get the new indexes by cropping the first startFE-1
newIdxRx = idxRx(startFE:end,:);
% zero all previous rx Signals
repetition.rxSignal(:)  = 0;
% assign new values to the partial acquisition entries
repetition.rxSignal(newIdxRx(:)) = 1:length(newIdxRx(:));


%% final message
if  expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done for %s repetition, elapsed time %.3fs',...
        functionName, repetition.type, toc(tTotal));
    fprintf(fid, '\n  Original # samples  %d', length(idxRx(:)));
    fprintf(fid, '\n  Partial  # samples  %d', length(newIdxRx(:)));
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end
