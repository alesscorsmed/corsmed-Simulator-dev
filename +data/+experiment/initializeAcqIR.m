function [prepIR] = initializeAcqIR()
%
% DATA.ACQUISITIONDATA.INITIALIZEACQIR
%
%	Function that initializes an IR (inversion recovery) data structure.
%   Returns a prepIR with empty fields.
%
%
% INPUT
%   None
%
% OUTPUT
%   prepIR   inversion recovery structure
%
%========================  CORSMED AB © 2020 ==============================
%

% basic info
prepIR.Apply            = 0; % do not apply by default
prepIR.type             = 'sinc';
prepIR.TI               = 100e-3; % inversion time
prepIR.flipAngle        = 180; % in degrees
prepIR.phase            = 0; % in radians
prepIR.duration         = 3e-3;
prepIR.cycles           = 2;
prepIR.sliceThickness   = 6e-3;
prepIR.doSliceSelect    = 0;
prepIR.preRWScale       = 1.0;
prepIR.postRWScale      = 1.0;
prepIR.keepTimingSS     = 1;
prepIR.postRFswc        = 0;  % software crusher at the end of the IR pulse