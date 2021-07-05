function [coilSystem] = computeSAR(coilSystem, anatomicalModel)
%
% COILS.COMPUTESAR
%
%     Calculates the SAR map and the Average SAR
%     given a Tx coil and an Anatomical Model
% 
%     Calculate average SAR - Check:
%       eq. 13 from Cao et al. (MRM 2014)
%       equation 6 from 
%       http://www.scielo.br/scielo.php?script=sci_arttext&pid=S2179-10742013000200010
%       and Christos' notes.
%
% NOTE: SAR values are unitary, unweighted by pulse Tx waveform
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
functionName = 'coils.computeSAR';
if (nargin < 2)
    ME = MException('Domain:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

need10g = 0; % avoid computing 10g SAR -- not used and expensive

%% default SAR values
sar.estSAR = 'N/A';
sar.avgSAR = 0;
sar.cm3SAR = 0;
sar.g10SAR = 0;
sar.voxSAR = 0;
sar.mapSAR = [];

% check if the model has electric properties
if  anatomicalModel.numTissues > 4
    hasElectricProperties = nnz(anatomicalModel.tissueValues(:,5:6));
else
    hasElectricProperties = 0;
end



%% load the required field data
try % load the E rms fields
    indexTx = coilSystem.indexTx;
    if isempty(coilSystem.coilModel{indexTx}.maps)
        coilSystem.coilModel{indexTx}.maps = load( ...
            coilSystem.coilModel{indexTx}.data.mapsFile,...
            'bx', 'by', 'sosB1m', 'sosB1p', 'sosErms', 'dim', 'spatial' );
    end
    coilStruct = coilSystem.coilModel{indexTx};
catch
    try
        %% if not possible to load the data, try the default Tx
        indexTx = coilSystem.defaultIndex;
        if isempty(coilSystem.coilModel{indexTx}.maps)
            coilSystem.coilModel{indexTx}.maps = load( ...
                coilSystem.coilModel{indexTx}.data.mapsFile,...
                'bx', 'by', 'sosB1m', 'sosB1p', 'sosErms', 'dim', 'spatial' );
        end
        coilStruct = coilSystem.coilModel{indexTx};
        message = sprintf('WARNING at %s: no EM fields for coil %s, default (%s) will be used for SAR',...
            functionName, coilSyste.activeTx, coilSystem.defaultTx );
        fprintf(1, '\n\n%s\n\n', message);
    catch
        %% no data available, not available SAR to report
        coilStruct = [];
        message = sprintf('WARNING at %s: no EM fields data available -- SAR N/A',...
            functionName );
        fprintf(1, '\n\n%s\n\n', message);
    end
end

%% compute the SAR given the txCoilData
if ~isempty(coilStruct) && hasElectricProperties
    
    if strcmpi(coilSystem.txMode, 'pTX')
        %% parallel transmit SAR calculation
        message = sprintf('WARNING at %s: parallel transmit not available',...
            functionName );
        fprintf(1, '\n\n%s\n\n', message);
    else
        %% use SOS map
        try
            %% coil data
            coilShift  = coilSystem.isocenter - coilStruct.data.coilIsocenter;
            rCoil      = coilStruct.maps.spatial + coilShift;
            dim        = coilStruct.maps.dim;
            
            %% interpolate map
            [Erms] = coils.interpolateMaps3D( ...
                anatomicalModel.spatial, coilStruct.maps.sosErms, rCoil, dim );
            
            %% Compute the SAR for the model
            
            %    SAR = sigma/(2*dens) * Erms^2
            densityMap  = anatomicalModel.tissueValues(anatomicalModel.tissueType,5);
            elecCondMap = anatomicalModel.tissueValues(anatomicalModel.tissueType,6);
            % scale the Erms by the B1+ value at the isocenter and cummulative excitation
            sar.mapSAR = (elecCondMap./(2*densityMap)).*...
                (Erms./(coilStruct.data.b1mIsocenter)).^2;
            
            % Compute average SAR over whole volume
            sar.avgSAR = mean(sar.mapSAR);
            sar.estSAR = sar.avgSAR; % use it as estimate
            
            % peak voxel
            [sar.voxSAR,idxPeak] = max(sar.mapSAR);
            
            % peak 1cm3-avg (find voxels within a sphere of 1cm3 around peak, and average)
            rsquare = (3/(4*pi)*0.01^3)^(2/3); % find r^2 from V = 4/3*PI*r^3 with V=1cm3
            idx1cm3 = sum((anatomicalModel.spatial - anatomicalModel.spatial(idxPeak,:)).^2,2) <= rsquare;
            sar.cm3SAR = mean(sar.mapSAR(idx1cm3));
            
            % 10g SAR
            if need10g
                stride = 10;
                rdelta = mean(anatomicalModel.resolution);
                voxelVolume = prod(anatomicalModel.resolution);
                for radius=rdelta:stride*rdelta:0.10
                    idx10g = sum((anatomicalModel.spatial - anatomicalModel.spatial(idxPeak,:)).^2,2) <= radius^2;
                    mass = sum(voxelVolume*densityMap(idx10g));
                    if mass >= 0.01
                        break;
                    end
                end
                for radius=radius-stride*rdelta:rdelta:radius
                    idx10g = sum((anatomicalModel.spatial - anatomicalModel.spatial(idxPeak,:)).^2,2) <= radius^2;
                    mass = sum(voxelVolume*densityMap(idx10g));
                    if mass >= 0.01
                        break;
                    end
                end
                sar.g10SAR = mean(sar.mapSAR(idx10g));
            else
                sar.g10SAR = 'N/A';
            end
            
        catch
            % something amiss
            message = sprintf('WARNING at %s: SAR computation failed',...
                functionName );
            fprintf(1, '\n\n%s\n\n', message);
        end
        
    end
    
else
    message = sprintf('WARNING at %s: no EM fields / Tissue data available -- SAR N/A',...
        functionName );
    fprintf(1, '\n\n%s\n\n', message);
end

%% assign sar to coil struct
coilSystem.coilModel{indexTx}.sar = sar;
