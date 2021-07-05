function [solution] = initializeSolution(numIsochromats,numRxs,numCoils)
%
% DATA.SIMULATION.INITIALIZESOLUTION
%
%	Function that initializes the solution with correct sizes from.
%   Returns a simulationControl with empty fields.
%
%   This function is useful to define the fields.
%
% INPUT
%   numIsochromats
%   numRxs
%   numCoils
%
% OUTPUT
%   expControl   expControl structure
%
%========================  CORSMED AB © 2020 ==============================
%

% basic info
solution.numCoils   = numCoils;
solution.numReads   = numRxs;
solution.indexes    = [];
solution.Mm         = zeros(numIsochromats,1);
solution.Mp         = zeros(numIsochromats,1);
solution.Mx         = zeros(numIsochromats,1);
solution.My         = zeros(numIsochromats,1);
solution.Mz         = zeros(numIsochromats,1);
solution.dMpDx      = zeros(numIsochromats,1);
solution.dMpDy      = zeros(numIsochromats,1);
solution.dMpDz      = zeros(numIsochromats,1);
solution.Sx         = zeros(numCoils*numRxs,1);
solution.Sy         = zeros(numCoils*numRxs,1);
solution.Sz         = zeros(numCoils*numRxs,1);

