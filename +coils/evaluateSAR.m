function [sarReport] = evaluateSAR(pulseSequence,coilSystem,acqData,expControl,...
    anatomicalModel)
%
% COILS.EVALUATESAR
%
%     Evaluates the SAR by weighting with pulse sequence RF data
%
%
% INPUT
%   coilSystem          struct with coils data
%   anatomicalModel     struct with model data
%
% OUTPUT
%   coilSystem          updated struct with correct active coils
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'coils.evaluateSAR';
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

%% compute time In-Real-Life and average RF power
% assume all repetitions, including slices and NEX
if acqData.is3D
    % number of averages
    numRep = acqData.NEX;
    numIL = 1; % force to no interleaving
else
    % number of averages and slices
    numRep = acqData.numSlices*acqData.NEX;
    % set the interleaving
    if acqData.interleaveSlice
        numIL = min(acqData.numSlices,pulseSequence.maxIL);
    else
        numIL = 1;
    end
end
% total sequence (1 rep, 1 slice) time
timeSeq = pulseSequence.time(end);

% peak signal
rfPeak  = max(abs(pulseSequence.rfmSignal(:)));
% cummulative RF deposition
rfSum   = sum( pulseSequence.timeDiff(:).*abs(pulseSequence.rfmSignal(:)) );
rfSum   = rfSum*numRep;

if strcmp(acqData.pulseSeqFamilyName,'cine-bssfp')
    timeIRL = sequence.tools.CINEcalculateTIRL(acqData,anatomicalModel,...
        pulseSequence);
elseif expControl.model.perfusion.apply % perfusion
    timeIRL = sequence.tools.perfusionCalculateTIRL(...
        expControl.model.perfusion,anatomicalModel);
else    
    % compute the total IRL time: scale by numRep, divided by IL factor
    timeIRL = timeSeq*ceil(numRep/numIL);
end

% average time (minimum is 1s) and 1s average RF
timeAvg = max(1,timeIRL);
rfAvg   = rfSum/timeAvg;

% store data
sarReport.B1rmsEst = rfAvg;
sarReport.B1rmsAvg = rfAvg;
sarReport.B1rmsSum = rfSum;
sarReport.B1rmsMax = rfPeak;
sarReport.timeAvg  = timeAvg;
sarReport.timeIRL  = timeIRL;

%% Scale SAR
sarData = coilSystem.coilModel{coilSystem.indexTx}.sar;
if isempty(sarData) || strcmpi(sarData.estSAR ,'n/a')
    %% no SAR available
    sarReport.sarAVG    = 'N/A';
    sarEst              = 0;
    sarReport.sarEST    = sarEst;   
    sarReport.sarCM3    = 0;
    %sarReport.sar10G = 0;
    sarReport.sarVOX    = 0;
else
    %% Compute SAR
    % B1+ and E fields scale linearly with I
    % But SAR = sigma/(2*dens) * Erms^2 scales with I^2
    rfScale             = rfAvg^2;
    sarEst              = rfScale*sarData.estSAR;
    sarReport.sarAVG    = num2str(round(100*sarEst)/100); % average 1 sec
    sarReport.sarEST    = sarEst;   
    sarReport.sarCM3    = rfScale*sarData.cm3SAR;
    %sarReport.sar10G    = rfScale*sarData.g10SAR;
    sarReport.sarVOX    = rfScale*sarData.voxSAR;
end
   

%% Report SAR and TIRL to front end
UIindicators = eduTool.frontend.reportSARtIRL(sarEst,timeIRL,expControl);

sarReport.UIindicators = UIindicators;

%% final message
if expControl.debug.debugMode
    fprintf(fid, ...
        '\n%s : done, elapsed time %.3fs',...
        functionName, toc(tTotal));
    fprintf(fid, '\n  Reported Est.   SAR   %s W/Kg/s', sarReport.sarAVG);
    fprintf(fid, '\n  Sequence Time         %.3f s', timeSeq);
    fprintf(fid, '\n  # Repetitions         %d (with %d Averages)',...
        numRep, acqData.NEX);
    fprintf(fid, '\n  # Interleaving    %d ',numIL);
    fprintf(fid, '\n  Time IRL              %.3f s', timeIRL);
    fprintf(fid, '\n  Peak         B1+rms   %.3f uT', rfPeak*1e6);
    fprintf(fid, '\n  Time Avg     B1+rms   %.3f uT/s', rfAvg*1e6);
    fprintf(fid, '\n  Time Avg Global SAR   %.3f W/Kg/s', sarReport.sarEST);
    %fprintf(fid, '\n  Time Avg    10g SAR   %.3f W/Kg/s', sarReport.sar10G);
    fprintf(fid, '\n  Time Avg   1cm3 SAR   %.3f W/Kg/s', sarReport.sarCM3);
    fprintf(fid, '\n  Time Avg  voxel SAR   %.3f W/Kg/s', sarReport.sarVOX);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end
