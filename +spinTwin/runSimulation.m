function [timeSolution, stats] = runSimulation( ...
    simModel, pulseSequence, motionModel, simControl, dbgControl)
%
% SPINTWIN.RUNSIMULATION
%
%	Simulation interface.
%
% INPUT
%   pulseSequence   sub sequence to simulate
%   simModel        structure with slice model data from spinModel 
%   motionModel     struct with motion model
%   expControl      experiment control struct
%   timeSolution    struct with solution struct with initial data
%
% OUTPUT
%   timeSolution    solution struct filled with final data
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'spinTwin:runSimulation';
% check args
if (nargin < 2) || isempty(simModel) || isempty(pulseSequence)
    ME = MException(['error:',functionName],...
        'Wrong arguments.');
    throw(ME);
end
% optional arguments
if (nargin < 3) || isempty(simControl)
    % empty motionModel: will ignore motion
    motionModel = [];
end
if (nargin < 4) || isempty(simControl)
    % initialize the simulation control w/ defaults
    [simControl] = spinTwin.setup.initializeSimControl();
end
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

%% initialize the solution
timeSolution = data.simulation.initializeSolution(...
    simModel.numIsochromats, pulseSequence.numRxs, simModel.numRxCoils);
%% initialize to equilibrium using current PD
if ~isfield(simModel, 'pd') || isempty(simModel.pd)
    simModel.pd = reshape(simModel.tissueValues(simModel.tissueType(:),3),[],1);
    simModel.pd(isnan(simModel.pd)) = 0.0;
end
timeSolution.Mz(:) = simModel.mu*simModel.b0*simModel.pd(:);

%% call the simulator
if simModel.numIsochromats > 0
    switch lower(simControl.simulationEngine)
        case lower('phasorplus')
            simControl.simulationEngine = 'PhasorPlus';
            [timeSolution, stats] = spinTwin.fwdBloch.runPhasorPlus( ...
                timeSolution, pulseSequence, simModel, motionModel, ...
                simControl, dbgControl);
        case lower('phasor')
            simControl.simulationEngine = 'Phasor';
            [timeSolution, stats] = spinTwin.fwdBloch.runPhasor( ...
                timeSolution, pulseSequence, simModel, motionModel, ...
                simControl, dbgControl);
        case lower('diffusion')
            simControl.simulationEngine = 'Diffusion';
            [timeSolution, stats] = spinTwin.fwdBloch.runDiffusion(...
                timeSolution, pulseSequence, simModel, motionModel, ...
                simControl, dbgControl);
        otherwise % basic Bloch
            simControl.simulationEngine = 'Bloch';
            [timeSolution, stats] = spinTwin.fwdBloch.runBloch(...
                timeSolution, pulseSequence, simModel, motionModel, ...
                simControl, dbgControl);
    end
end

%% process solution
timeSolution.Sx = reshape(timeSolution.Sx, ...
    simModel.numRxCoils, pulseSequence.numRxs).';
timeSolution.Sy = reshape(timeSolution.Sy, ...
    simModel.numRxCoils, pulseSequence.numRxs).';
timeSolution.Sz = reshape(timeSolution.Sz, ...
    simModel.numRxCoils, pulseSequence.numRxs).';

%% report
if  dbgControl.mode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : %s simulation done for sequence %s',...
        functionName, simControl.simulationEngine, pulseSequence.name);
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end

