function [interpMaps] = interpolateMaps3D( rTarget,...
                   originalMaps, rOriginal, dim)
%
% COILS.INTERPOLATEMAPS3D
%
%     Interpolates a set 3D maps
%     into the query xq yq zq coordinates
%     returns intepolated maps
%
%	NOTE:   interp3 uses Matlab's meshgrid system, 
%           which is [y,x,z] coordinates
%           interpolation will be applied: interp3(Y,X,Z,V,Yq,Xq,Zq)
%           to be consistent with coil maps' coordinates
%
% INPUT
%   rTarget         query points as [xq, yq, zq]
%   originalMaps    original 3D maps
%   rOriginal       3D coordinates of the original maps
%   dim             dimensions of the 3D grid (nx, ny, nz)
%
% OUTPUT
%   interpMaps       interpolated maps
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'coils.interpolateMaps3D';
if (nargin < 4)
    ME = MException('Domain:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

%% reshape to 3D form
try
    rOriginal       = reshape(rOriginal, dim(1), dim(2), dim(3), 3);
    originalMaps    = reshape(originalMaps, dim(1), dim(2), dim(3), []);
catch
    ME = MException('Domain:wrongDimensions',...
        '%s : maps and coordinate dimensions do not match',functionName);
    throw(ME);
end

%% prepare query
numMaps = size(originalMaps,4);    
xq = reshape(rTarget(:,1),[],1);
yq = reshape(rTarget(:,2),[],1);
zq = reshape(rTarget(:,3),[],1);
interpMaps = zeros(size(xq,1),numMaps);

%% interpolator from the data
for mm = 1:numMaps
    % NOTE: interp3 uses Matlab's meshgrid system, which is [y,x,z] coordinates
    interpMaps(:,mm) = interp3( ...
        rOriginal(:,:,:,2), rOriginal(:,:,:,1), rOriginal(:,:,:,3),...
        originalMaps(:,:,:,mm), yq, xq, zq );
end
interpMaps(isnan(interpMaps)) = 0.0;

