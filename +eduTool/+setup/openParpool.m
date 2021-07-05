function p = openParpool()


% SETUP.OPENPARPOOL
%
%	new parallel pool
%
% INPUT
%
%
% OUTPUT
%
%   p          parallel pool object
%
%========================  CORSMED AB © 2020 ==============================
%


%% create Parpool depending on GPU number

if(gpuDeviceCount>1)
    if(isempty(gcp('nocreate')))
        p = parpool(gpuDeviceCount);
        p.IdleTimeout = Inf;
    else
        p = gcp();
    end
else
    p = gcp();
    p.IdleTimeout = Inf;
end
