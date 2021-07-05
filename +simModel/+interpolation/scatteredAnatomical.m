function [model] = scatteredAnatomical( model, rotMat, refPoint, ...
    anatomicalModel, dbgControl )
%
% SIMMODEL.INTERPOLATION.SCATTEREDANATOMICAL
%
%   Interpolates for the query points ( positions of the model )
%   using the original anatomicalModel grid and data.
%
%   Anatomical model points are rotated to the simulation plane.
%   The interpolation is applied on those rotated coordinates,
%   assuming a non-structured coordinate system (scattered interpolation).
%
% INPUT
%
%
%
% OUTPUT
%   model    updated model struct with interpolated data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'simModel:interpolation:scatteredAnatomical';
% check args
if (nargin < 4)
    ME = MException(['error:',functionName],...
        'Wrong arguments.');
    throw(ME);
end
% optional arguments
if (nargin < 5) || isempty(dbgControl)
    dbgControl.mode = 0;
    dbgControl.file = [];
end
% info for debugging
if dbgControl.mode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(dbgControl.file,'a');
    catch
        fid = 1;
    end
    tTotal = tic();
    fprintf(fid, '\n%s : start', functionName);
end

%% apply the transformations on the anatomical model points
% to bring them to the slice plane
rAnatomical = reshape( anatomicalModel.spatial, [], 3 );
if ~isempty(refPoint)
    rAnatomical = rAnatomical - refPoint; % translate to slice plane
end
rAnatomical = rAnatomical * rotMat; % rotate to the plane

%% query points from the model
rSlice = model.r3D;
nSlice = length(rSlice);

% cropping: find limits of query points to restrict locations in 3D
rSliceMax = max(rSlice) + anatomicalModel.resolution;
rSliceMin = min(rSlice) - anatomicalModel.resolution;

% valid indexes according to the cropping
idxSlice = find( ( rAnatomical(:,1) >= rSliceMin(1) ) ...
    & ( rAnatomical(:,1) <= rSliceMax(1) ) ...
    & ( rAnatomical(:,2) >= rSliceMin(2) ) ...
    & ( rAnatomical(:,2) <= rSliceMax(2) ) ...
    & ( rAnatomical(:,3) >= rSliceMin(3) ) ...
    & ( rAnatomical(:,3) <= rSliceMax(3) ) );

% use only required indexes
rAnatomical = rAnatomical(idxSlice,:);
nAnatomical = length(rAnatomical);

%% apply the interpolation
if isempty(idxSlice)
    
    %% slice is Out-of-bounds: Empty model
    model.nonZeroIndex   = []; % indexes with the entries of non zero
    model.numIsochromats = 0; % number of isochromats to simulate
    % external voxel properties (not from anatomical model)
    model.bi         = zeros(model.numIsochromats,1);
    model.pd         = zeros(model.numIsochromats,1);
    model.xDiffusion = zeros(model.numIsochromats,1);
    model.yDiffusion = zeros(model.numIsochromats,1);
    model.zDiffusion = zeros(model.numIsochromats,1);
    % tissue types
    model.tissueType = zeros(model.numIsochromats,1);
    
else
    %% slice is not Out-of-bounds
    try
        %% use Python functionality       
        % apply the interpolation function on indexes
        Fpy = py.scatteredNomat.scatteredPython( idxSlice, ...
            transpose(rAnatomical), transpose(rSlice), nAnatomical, nSlice);
        % get python result back into matlab format
        Fpy = double(py.array.array('d', py.numpy.nditer(Fpy)));
        % assign to the corresponding data
        vq = reshape(transpose(Fpy),[], 1);
        
        % vq has the indexes of the 'interpolated' entries (indexes of nearest)
        % get valid entries by removing NANs
        idxTissues = find(~isnan(vq));
        vq = vq(idxTissues);
        nIso = length(idxTissues);
        % assign to model
        model.nonZeroIndex      = idxTissues; % indexes with the entries of non zero
        model.numIsochromats    = nIso; % number of isochromats to simulate
        
        if anatomicalModel.numFrames > 1
            % find the entries in vq that are subject to frame motion
            [~,queryIdxFrame,origIdxFrame] = intersect(vq, anatomicalModel.idxFrame);
            % assign base tissue type
            vTissueType = reshape(anatomicalModel.tissueType(vq(:)),nIso,1);
            % create a matrix of tissueTypes for each frame
            model.tissueType = repmat( vTissueType, [1, anatomicalModel.numFrames] );
            % modify the entries subject to motion
            model.tissueType(queryIdxFrame,:) = anatomicalModel.tissueTypeFrame(origIdxFrame,:);
        else
            % non motion tissue
            model.tissueType = reshape(anatomicalModel.tissueType(vq(:)),nIso,1);
        end
        
        % interpolate the pdInhomogeneity if exists
        if ~isempty( anatomicalModel.pdInhomogeneity )
            % vAnatomical = anatomicalModel.pdInhomogeneity(idxSlice);
            % % apply the interpolation function
            % Fpy = py.scatteredNomat.scatteredPython( vAnatomical, ...
            %     rAnatomical(:).', rSlice(:).', nAnatomical, nSlice);
            % % get python result back into matlab format
            % Fpy = double(py.array.array('d', py.numpy.nditer(Fpy)));
            % % assign to the corresponding data
            % vq = reshape(transpose(Fpy),[], 1);
            % vq( isnan(vq) ) = 1.0;
            % model.pd = reshape(vq(idxTissues), [nIso,1]);
            
            % no need to interpolate, use indexes
            model.pd = reshape( anatomicalModel.pdInhomogeneity(vq(:)), [nIso,1]);
            model.pd(isnan(model.pd)) = 0.0;
            
        else
            model.pd = [];
        end
        
        % interpolate the b0Inhomogeneity if exists
        if ~isempty( anatomicalModel.b0Inhomogeneity )
            % vAnatomical = anatomicalModel.b0Inhomogeneity(idxSlice);
            % % apply the interpolation function
            % Fpy = py.scatteredNomat.scatteredPython( vAnatomical, ...
            %     rAnatomical(:).', rSlice(:).', nAnatomical, nSlice);
            % % get python result back into matlab format
            % Fpy = double(py.array.array('d', py.numpy.nditer(Fpy)));
            % % assign to the corresponding data
            % vq = reshape(transpose(Fpy),[], 1);
            % vq( isnan(vq) ) = 0.0;
            % model.bi = reshape(vq(idxTissues), [nIso,1]);
            
            % no need to interpolate for nearest, use indexes
            model.bi = reshape( anatomicalModel.b0Inhomogeneity(vq(:)), [nIso,1]);
            model.bi(isnan(model.bi)) = 0.0;
        else
            model.bi = zeros(nIso,1);
        end
        
    catch
        
        % MATLAB's native interpolator from the data
        % vAnatomical = anatomicalModel.tissueType(idxSlice);
        F = scatteredInterpolant(rAnatomical,idxSlice);
        F.Method                = 'nearest';
        F.ExtrapolationMethod   = 'none';
        vq = F(rSlice(:,1), rSlice(:,2), rSlice(:,3));
        
        % vq has the indexes of the 'interpolated' entries (indexes of nearest)
        % get valid entries by removing NANs
        idxTissues = find(~isnan(vq));
        vq = vq(idxTissues);
        nIso = length(idxTissues);
        % assign to model
        model.nonZeroIndex      = idxTissues; % indexes with the entries of non zero
        model.numIsochromats    = nIso; % number of isochromats to simulate
        
        if anatomicalModel.numFrames > 1
            % find the entries in vq that are subject to frame motion
            [~,queryIdxFrame,origIdxFrame] = intersect(vq, anatomicalModel.idxFrame);
            % assign base tissue type
            vTissueType = reshape(anatomicalModel.tissueType(vq(:)),nIso,1);
            % create a matrix of tissueTypes for each frame
            model.tissueType = repmat( vTissueType, [1, anatomicalModel.numFrames] );
            % modify the entries subject to motion
            model.tissueType(queryIdxFrame,:) = anatomicalModel.tissueTypeFrame(origIdxFrame,:);
        else
            % non motion tissue
            model.tissueType = reshape(anatomicalModel.tissueType(vq(:)),nIso,1);
        end

        % interpolate the pdInhomogeneity if exists
        if ~isempty( anatomicalModel.pdInhomogeneity )
            % vAnatomical = anatomicalModel.pdInhomogeneity(idxSlice);
            % F = scatteredInterpolant(rAnatomical,vAnatomical);
            % F.Method                = 'nearest';
            % F.ExtrapolationMethod   = 'none';
            % vq = F(rSlice(:,1), rSlice(:,2), rSlice(:,3));
            % vq( isnan(vq) ) = 1.0;
            % model.pd = reshape(vq(idxTissues), [nIso,1]);
            
            % no need to interpolate, use indexes
            model.pd = reshape( anatomicalModel.pdInhomogeneity(vq(:)), [nIso,1]);
            model.pd(isnan(model.pd)) = 0.0;
            
        else
            model.pd = [];
        end
        
        % interpolate the b0Inhomogeneity if exists
        if ~isempty( anatomicalModel.b0Inhomogeneity )
            % vAnatomical = anatomicalModel.b0Inhomogeneity(idxSlice);
            % F = scatteredInterpolant(rAnatomical,vAnatomical);
            % F.Method                = 'nearest';
            % F.ExtrapolationMethod   = 'none';
            % vq = F(rSlice(:,1), rSlice(:,2), rSlice(:,3));
            % vq( isnan(vq) ) = 0.0;
            % model.bi = reshape(vq(idxTissues), [nIso,1]);
            
            % no need to interpolate for nearest, use indexes
            model.bi = reshape( anatomicalModel.b0Inhomogeneity(vq(:)), [nIso,1]);
            model.bi(isnan(model.bi)) = 0.0;
        else
            model.bi = zeros(nIso,1);
        end
        
    end
end

%% report
if  dbgControl.mode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : scattered interpolation done', functionName);
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n  # Isochromats     %d', model.numIsochromats);
    fprintf(fid, '\n\n');
    if fid ~=1
        fclose(fid);
    end
end
