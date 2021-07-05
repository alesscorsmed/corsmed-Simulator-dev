function [encodingPG] = initializePG()
%
% DATA.ACQUISITIONDATA.INITIALIZEPG
%
%	Function that initializes a Pulsed Gradient data structure,
%   to be used for diffusion sequences.
%   Returns a encodingPG with empty fields.
%
%   This function is useful to define the fields.
%
% INPUT
%   None
%
% OUTPUT
%   encodingPG   structure with PG encoding info
%
%========================  CORSMED AB © 2020 ==============================
%

% basic info
encodingPG.Apply    = 0; % do not apply by default
encodingPG.AG       = 30e-3; % gradient amplitude
encodingPG.TG       = 5e-3;  % gradient duration
encodingPG.TAU      = 150e-3; % echo time
encodingPG.Area     = [] ; % gradient area, T.B.D.
encodingPG.Time     = [] ; % gradient total time, T.B.D.
encodingPG.TP       = [] ; % plateau time, T.B.D. if applicable
encodingPG.Beta     = [] ; % effective Beta, T.B.D.
encodingPG.Dir      = 'x'; % encoding direction 'x' / 'y' / 'z'
