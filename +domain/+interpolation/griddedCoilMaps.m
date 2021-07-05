function [model] = griddedCoilMaps( model, rotMat, refPoint,...
    coilSystem, expControl)
%
% DOMAIN.INTERPOLATION.GRIDDEDCOILMAPS
%
%   Interpolates for the query points ( positions of the model )
%   using the original Coil grid and data.
%
%   Query (slice) points are rotated back to the coil coordinates.
%   The interpolation is applied on those original coordinates,
%   assuming a structured coordinate system (gridded interpolation).
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
functionName = 'domain.interpolation.griddedCoilMaps';
if (nargin < 5)
    ME = MException('Domain:wrongArgCount',...
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
    fprintf(fid, '\n\n%s : start', functionName);
end

%% basic info
% undo the transformations in the query points
% to bring them back to original coordinates
rSlice = model.r3D(model.nonZeroIndex,:); % only useful points
rSlice = rSlice * rotMat.'; % transpose of rotation to bring back
if ~isempty(refPoint)
    rSlice = rSlice + refPoint; % translate back to ref point
end

%% transmit coil interpolation
% get the active transmit coil structure
disp(['Coil index Tx',num2str(coilSystem.indexTx)])
coilStruct = coilSystem.coilModel{coilSystem.indexTx};
%% check transmit mode
switch lower(coilSystem.txMode)
    
    case lower('b1pSoS')
        %% Sum of Squares transmit: Use b1pSoS map
        txSensType = 'b1pSOS'; % transmit
        [txCoilMapsX,txCoilMapsY,txSensType] = coils.interpolateCoilSens( ...
            rSlice, coilStruct, txSensType, ...
            coilSystem.isocenter, coilSystem.b1pScaling );
        
    case lower('pTX')
        %% parallel transmit: Use the Bx and By maps to create the B1+ maps
        txSensType = 'pTX'; % transmit
        [txCoilMapsX,txCoilMapsY,txSensType] = coils.interpolateCoilSens( ...
            rSlice, coilStruct, txSensType, ...
            coilSystem.isocenter, coilSystem.b1pScaling );
        
    otherwise
        %% ideal mode, keep ones
        txSensType = 'ideal'; % transmit
        [txCoilMapsX,txCoilMapsY,txSensType] = coils.interpolateCoilSens( ...
            rSlice, coilStruct, txSensType, coilSystem.isocenter, 1.0 );
end

%% receive coil interpolation
% get the active transmit coil structure
coilStruct = coilSystem.coilModel{coilSystem.indexRx};
%% check receive mode
switch lower(coilSystem.rxMode)
    
    case lower('b1mSoS')
        %% Sum of Squares receiver: Use b1mSoS map
        rxSensType = 'b1mSOS'; % transmit
        [rxCoilMapsX,rxCoilMapsY,rxSensType] = coils.interpolateCoilSens( ...
            rSlice, coilStruct, rxSensType, ...
            coilSystem.isocenter, coilSystem.b1mScaling );
        
    case lower('pRX')
        %% parallel receive: Use the Bx and By maps to create the B1- maps
        rxSensType = 'pRX';
        [rxCoilMapsX,rxCoilMapsY,rxSensType] = coils.interpolateCoilSens( ...
            rSlice, coilStruct, rxSensType, ...
            coilSystem.isocenter, coilSystem.b1mScaling );
        
    otherwise
        %% ideal mode, keep ones
        rxSensType = 'ideal';
        [rxCoilMapsX,rxCoilMapsY,rxSensType] = coils.interpolateCoilSens( ...
            rSlice, coilStruct, rxSensType, coilSystem.isocenter, 1.0 );
end

%% assign
txNumCoils          = size(txCoilMapsX,2);
rxNumCoils          = size(rxCoilMapsX,2);
model.txCoilMapsX   = txCoilMapsX;
model.txCoilMapsY   = txCoilMapsY;
model.rxCoilMapsX   = rxCoilMapsX;
model.rxCoilMapsY   = rxCoilMapsY;
model.numRxCoils    = rxNumCoils;
model.numTxCoils    = txNumCoils;

%% report
if expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : coil maps generated for %s',...
        functionName, coilStruct.data.name );
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n  # Query Points    %d', size(rSlice,1));
    fprintf(fid, '\n  TX Coil (#%2d)     %s', ...
        coilSystem.indexTx, ...
        coilSystem.coilModel{coilSystem.indexTx}.data.name);
    fprintf(fid, '\n  TX op. mode       %s', txSensType);
    fprintf(fid, '\n  TX # Coil Maps    %d', txNumCoils);
    fprintf(fid, '\n  RX Coil (#%2d)     %s', ...
        coilSystem.indexRx, ...
        coilSystem.coilModel{coilSystem.indexRx}.data.name);
    fprintf(fid, '\n  RX op. mode       %s', rxSensType);
    fprintf(fid, '\n  RX # Coil Maps    %d', rxNumCoils);
    fprintf(fid, '\n\n');
    if fid ~=1
        fclose(fid);
    end
end

