function [rfPulse] = initializeRF(angle)
%
% DATA.ACQUISITIONDATA.INITIALIZERF
%
%	Function that initializes an RF pulse data structure.
%   Returns a rfPulse with empty fields.
%
%   This function is useful to define the fields.
%
% INPUT
%   None
%
% OUTPUT
%   rfPulse   rf pulse structure
%
%========================  CORSMED AB © 2020 ==============================
%

if (nargin < 1 || isempty(angle))
    angle=90;
end

% basic info
rfPulse.type            = 'sinc';
rfPulse.flipAngle       = angle; % in degrees
rfPulse.phase           = 0; % in radians
rfPulse.duration        = 3e-3;
rfPulse.cycles          = 2;
rfPulse.sliceThickness  = 6e-3;
rfPulse.doSliceSelect   = 0;
rfPulse.preRWScale      = 1.0;
rfPulse.postRWScale     = 1.0;
rfPulse.keepTimingSS    = 1;

