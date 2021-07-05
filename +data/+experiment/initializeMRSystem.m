function [mrSystem] = initializeMRSystem(scanner)
%
% DATA.EXPERIMENT.INITIALIZEMRSYSTEM
%
%	Function that initializes the mrSystem data structure,
%   with information of the scanner.
%   Returns a mrSystem with empty fields.
%
%   This function is useful to define the fields.
%
% INPUT
%   None
%
% OUTPUT
%   mrSystem   mrSystem structure
%
%========================  CORSMED AB © 2020 ==============================
%
if (nargin < 1)
    scanner = 'corsmed-Ideal-1p5T';
end

switch lower(scanner)
    
    case lower('Siemens-Avanto-1.5T')
        % basic info
        mrSystem.name   = scanner;
        mrSystem.b0     = 1.5; % T
        % gradient system
        mrSystem.maxGStrenght   = 0.045; % T/m
        mrSystem.SlewRate       = 200.0; % T/m/s
        
    case lower('Philips-Ingenia-1.5T')
        % basic info
        mrSystem.name   = scanner;
        mrSystem.b0     = 1.5; % T
        % gradient system
        mrSystem.maxGStrenght   = 0.045; % T/m
        mrSystem.SlewRate       = 200.0; % T/m/s
        
    case lower('GE-Optima450-1.5T')
        % basic info
        mrSystem.name   = scanner;
        mrSystem.b0     = 1.5; % T
        % gradient system
        mrSystem.maxGStrenght   = 0.044; % T/m
        mrSystem.SlewRate       = 200.0; % T/m/s
        
    case lower('Corsmed-1.5T')    
        % basic info
        mrSystem.name   = scanner;
        mrSystem.b0     = 1.5; % T
        % gradient system
        mrSystem.maxGStrenght   = 0.040; % T/m
        mrSystem.SlewRate       = 150.0; % T/m/s
        
    otherwise
        %% default is ideal 1.5T
        scanner = 'Corsmed-Ideal-1p5T'; 
        % basic info
        mrSystem.name   = scanner;
        mrSystem.b0     = 1.5; % T
        % gradient system
        mrSystem.maxGStrenght   = 0.040; % T/m
        mrSystem.SlewRate       = 0; % no ramps       
        % bandwidth limits, RF limits, others?

end
