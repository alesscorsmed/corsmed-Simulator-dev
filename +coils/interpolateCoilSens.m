function [coilmapX,coilmapY,sensType] = interpolateCoilSens( rSlice,...
                   coilStruct, sensType, isocenter, bmScale)
%
% COILS.INTERPOLATECOILSENS
%
%     Interpolates coil maps into a set of query points,
%     and returns the required sensitivity type
%
% INPUT
%   rSlice              query points [xq, yq, zq]
%   coilStruct          coilModel struct with data and maps
%   sensType            sensitivity to generate ideal/b1mSOS/b1pSOS/pRX/pTX
%   isocenter           system isocenter
%   bmScale             scale for the b map (usually b at Isocenter)
%
% OUTPUT
%   coilmapX,coilmapY   interpolated sensitivities (numIso,numCh)
%   sensType            returned sensType in case of fail
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'coils.interpolateCoilSens';
if (nargin < 5)
    ME = MException('Domain:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end
%% prepare the data
try
    if strcmpi(coilStruct.data.coilConfig,'optimal')
        sensType = 'ideal';
    else
        % coil valid FOV
        maxDomain   = max(abs(coilStruct.maps.spatial));
        maxRadius   = coilStruct.data.coilRadius - 0.005; % 0.5cm guard
        % move the coil to the system isocenter
        % take into account coil own isocenter ( position )
        coilShift   = isocenter - coilStruct.data.coilIsocenter;
        rCoil       = coilStruct.maps.spatial + coilShift;
        dim         = coilStruct.maps.dim;
    end
catch
    sensType = 'ideal';
    message = sprintf('WARNING at %s: coil %s with no field maps - ideal sensitivity returned',...
        functionName, coilStruct.data.name );
    fprintf(1, '\n\n%s\n\n', message);
end

% prepare as ideal (all 1) map
% assume zero phase (zero Cy)
%  received signal is (Cx -j Cy)*(Mx  +j My)
%  transmit signal is (Cx +j Cy)*(RFx +j RFy)
numIso      = size(rSlice,1);
coilmapX    =  ones(numIso,1);
coilmapY    = zeros(numIso,1);
%% check the kind of sensitivity to generate
switch lower(sensType)
    case lower('b1mSoS')
        numChannels = 1;
        %% Sum-of-Squares as receiver (B1m)
        if ~isempty(coilStruct.maps.b1mSOS)
            [coilmapX] = coils.interpolateMaps3D( ...
                rSlice, coilStruct.maps.b1mSOS, rCoil, dim );
            % for SoS assume zero phase (received signal is (Bx-jBy)*(Mx+jMy))
            coilmapY = zeros(size(coilmapX));
        else
            sensType = 'ideal';
            message = sprintf('WARNING at %s: coil %s with no b1mSOS fields - ideal sensitivity returned',...
                functionName, coilStruct.data.name );
            fprintf(1, '\n\n%s\n\n', message);
        end
        
    case lower('b1pSoS')
        numChannels = 1;
        %% Sum-of-Squares as transmitter (B1p)
        if ~isempty(coilStruct.maps.b1pSOS)
            [coilmapX] = coils.interpolateMaps3D( ...
                rSlice, coilStruct.maps.b1pSOS, rCoil, dim );
            % for SoS assume zero phase (transmitted signal is (Bx+jBy)*(Mx+jMy))
            coilmapY = zeros(size(coilmapX));
        else
            sensType = 'ideal';
            message = sprintf('WARNING at %s: coil %s with no b1pSOS fields - ideal sensitivity returned',...
                functionName, coilStruct.data.name );
            fprintf(1, '\n\n%s\n\n', message);
        end
        
    case lower('pRX')
        %% parallel receiver: Use the Bx and By maps to create the B1- map
        if ~isempty(coilStruct.maps.bx) && ~isempty(coilStruct.maps.by)
            numChannels = coilStruct.data.numCoils;
            % interpolate bx maps
            [interpMapX] = coils.interpolateMaps3D( ...
                rSlice, coilStruct.maps.bx, rCoil, dim );
            % interpolate by maps
            [interpMapY] = coils.interpolateMaps3D( ...
                rSlice, coilStruct.maps.by, rCoil, dim );
            % NOTE: things are a bit messed up here
            %    the received signal is         (Bx-jBy)*(Mx+jMy)
            %                                   where Bx and By can be complex
            %    the kernel is coded to return  (Cx-jCy)*(Mx+jMy)
            %                                   where Cx and Cy are no complex
            %    the coils come in Bx and By complex components
            %    we need to process them to maintain consistency
            %           generate receiving complex map Cmap = Bx-jBy
            %           Cx =  real(Cmap)
            %           Cy = -imag(Cmap)  <-- negative so that the kernel applies correct sign
            %    therefore
            %           Cx = real(interp_Bx) + imag(interp_By)
            %           Cy = real(interp_By) - imag(interp_Bx)
            %
            coilmapX = real(interpMapX) + imag(interpMapY);
            coilmapY = real(interpMapY) - imag(interpMapX);
        else
            sensType = 'ideal';
            numChannels = 1;
            message = sprintf('WARNING at %s: coil %s with no bx/by fields - ideal sensitivity returned',...
                functionName, coilStruct.data.name );
            fprintf(1, '\n\n%s\n\n', message);
        end
        
    case lower('pTX')
        %% parallel transmit: Use the Bx and By maps to create the B1+ map
        if ~isempty(coilStruct.maps.bx) && ~isempty(coilStruct.maps.by)
            numChannels = coilStruct.data.numCoils;
            % interpolate bx maps
            [interpMapX] = coils.interpolateMaps3D( ...
                rSlice, coilStruct.maps.bx, rCoil, dim );
            % interpolate by maps
            [interpMapY] = coils.interpolateMaps3D( ...
                rSlice, coilStruct.maps.by, rCoil, dim );
            % NOTE: things are a bit messed up here
            %    the transmit signal is          (Bx+jBy)*(Mx+jMy)
            %                                   where Bx and By can be complex
            %    the kernel is coded to generate (Cx+jCy)*(RFx+jRFy)
            %                                   where Cx and Cy are no complex
            %    the coils come in Bx and By complex components
            %    we need to process them to maintain consistency
            %           generate receiving complex map Cmap = Bx+jBy
            %           Cx = real(Cmap)
            %           Cy = imag(Cmap)
            %    therefore
            %           Cx = real(interp_Bx) - imag(interp_By)
            %           Cy = real(interp_By) + imag(interp_Bx)
            %
            coilmapX = real(interpMapX) - imag(interpMapY);
            coilmapY = real(interpMapY) + imag(interpMapX);
        else
            sensType = 'ideal';
            numChannels = 1;
            message = sprintf('WARNING at %s: coil %s with no bx/by fields - ideal sensitivity returned',...
                functionName, coilStruct.data.name );
            fprintf(1, '\n\n%s\n\n', message);
        end
        
    otherwise
        %% ideal case, sensitivities already computed
        sensType = 'ideal';
        numChannels = 1;
end

%% reshape
coilmapX = reshape(coilmapX, numIso, numChannels);
coilmapY = reshape(coilmapY, numIso, numChannels);
        
%% find and zero queries outside coil valid FOV domain
if ~strcmpi(sensType,'ideal')
    rSlice = rSlice - coilShift; % isocenter shift to be consistent with Coil
    idx = (rSlice(:,1).*rSlice(:,1) + rSlice(:,2).*rSlice(:,2) > maxRadius^2) ...
        | (abs(rSlice(:,3)) >= maxDomain(3));
    if ~isempty(idx)
        coilmapX(idx,:) = 0.0;
        coilmapY(idx,:) = 0.0;
    end
    % if it is not ideal, scale coilmaps with bmScale
    coilmapX = coilmapX / bmScale;
    coilmapY = coilmapY / bmScale;
end