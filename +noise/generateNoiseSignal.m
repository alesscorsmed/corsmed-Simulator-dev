function [raw] = generateNoiseSignal( raw, noiseData, expControl )
%
% NOISE.GENERATENOISESIGNAL
%
%	generates noise signal
%
% INPUT
%
%
% OUTPUT
%
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'noise.generateNoiseSignal';
if (nargin < 3)
    ME = MException('eduTool:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

%% info for debugging
if expControl.debug.debugMode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(expControl.debug.debugFile,'a');
    catch
        fid = 1;
    end
    % time it
    tTotal = tic();
    fprintf(fid, '\n%s : start', functionName);
end

%% extract data
numReads = raw.numReads;
numCoils = raw.numCoils;

%% get correct noise amplitude
% based on noise Std Dev and the simulation voxel volume
noiseSignal = noiseData.noiseAmp*( randn(numReads,numCoils) + 1j*randn(numReads,numCoils));

%% assign to raw struct
raw.noise = noiseSignal;

%% final message
if expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : done', functionName);
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n  Noise Std. Dev.   %.3g', noiseData.noiseSD);
    fprintf(fid, '\n  Noise Amplitude   %.3g', noiseData.noiseAmp);
    %% peak signal to noise
    if  noiseData.noiseAmp > 0
        SNR = max(abs(raw.Sx+1j*raw.Sy))/max(abs(noiseSignal));
        fprintf(fid, '\n  Signal    SNR     %.3g', SNR);
        fprintf(fid, '\n  Estimated SNR     %.3g', noiseData.SNR);
        fprintf(fid, '\n  Relative  SNR     %d%%', noiseData.relativeSNR);
    else
        fprintf(fid, '\n  Signal    SNR     INF');
        fprintf(fid, '\n  Estimated SNR     INF');
        fprintf(fid, '\n  Relative  SNR     INF');
    end

    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end

