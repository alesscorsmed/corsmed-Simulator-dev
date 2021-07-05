function [reconData] = initialize()
%
% DATA.RECONDATA.INITIALIZE
%
%	Function that initializes the reconData data structure.
%   Returns a simulationControl with empty fields.
%
%   This function is useful to define the fields.
%
% INPUT
%   None
%
% OUTPUT
%   reconData   reconData structure
%
%========================  CORSMED AB © 2020 ==============================
%

% k-space and image space
reconData.kSpace = [];
reconData.iSpace = [];

% ISMRMRD header
reconData.ismrmrd.header.encoding.encodedSpace.matrixSize.x     = 128;
reconData.ismrmrd.header.encoding.encodedSpace.matrixSize.y     = 128;
reconData.ismrmrd.header.encoding.encodedSpace.matrixSize.z     = 1;

reconData.ismrmrd.header.encoding.reconSpace.matrixSize.x       = 128;
reconData.ismrmrd.header.encoding.reconSpace.matrixSize.y       = 128;
reconData.ismrmrd.header.encoding.reconSpace.matrixSize.z       = 1;

reconData.ismrmrd.header.encoding.reconSpace.fieldOfView_mm.x   = 200;
reconData.ismrmrd.header.encoding.reconSpace.fieldOfView_mm.y   = 200;
reconData.ismrmrd.header.encoding.reconSpace.fieldOfView_mm.z   = 6;