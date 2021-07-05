function [coilSystem] = precomputeMRSafety( coilSystem, anatomicalModel )
%
% COILS.PRECOMPUTEMRSAFETY
%
%     computes unitary SAR values for transmit coils in the system.
%
% INPUT
%   coilSystem          struct with coils data
%   anatomicalModel     struct with model data
%
% OUTPUT
%   coilSystem          updated struct with correct active coils
%
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'coils.precomputeMRSafety';
if (nargin < 2)
    ME = MException('Domain:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

%% report start
fid = 1;
tTotal = tic();
fprintf(fid, '\n%s : start', functionName);

%% compute SAR for transmit capable coils
numTx = 0;
for cc = 1:coilSystem.numModels
    
    %% compute SAR for transmit capable coils
    if ( coilSystem.coilModel{cc}.data.isTx )
  
        %vars
        tSAR  = tic();
        numTx = numTx + 1;
        
        %% make the coil transmit active
        coilSystem.indexTx = cc;
        
        %% compute SAR for the transmit coil with the defined isocenter
        [coilSystem] = coils.computeSAR(coilSystem, anatomicalModel);
        
        %% report
        fprintf(fid, '\n  Coil #%d: %20s -- SAR computed, Elapsed Time %.3fs',...
            cc, coilSystem.coilModel{cc}.data.name, toc(tSAR) );

    end
    
end

%% reset the default active transmit
coilSystem.indexTx = coilSystem.defaultIndex;

%% report
tTotal = toc(tTotal);
fprintf(fid, '\n%s : unitary SAR pre-computed for %d coils',...
    functionName, numTx );
fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
fprintf(fid, '\n');
