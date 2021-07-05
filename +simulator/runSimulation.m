function [timeSolution] = runSimulation(...
    simModel, pulseSequence, expControl, timeSolution , GPUindex)
%
% SIMULATOR.RUNSIMULATION
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
functionName = 'simulator.runSimulation';
if (nargin < 3)
    ME = MException('simulator:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end
if (nargin < 4)
    timeSolution = [];
end

%% open debugging
if expControl.debug.debugMode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(expControl.debug.debugFile,'a');
    catch
        fid = 1;
    end
    tTotal = tic();
    fprintf(fid, '\n%s : start', functionName);
end

%% initialize the solution
% generate the pd array if not available
if ~isfield(simModel, 'pd') || isempty(simModel.pd)
    simModel.pd = reshape(simModel.tissueValues(simModel.tissueType(:),3),[],1);
    simModel.pd(isnan(simModel.pd)) = 0.0;
end
% initialize solution with equilibrium
if isempty(timeSolution)
    %% allocate
    timeSolution = data.simulation.initializeSolution(...
        simModel.numIsochromats, pulseSequence.numRxs, simModel.numRxCoils);
    %% initialize to equilibrium
    timeSolution.Mz(:) = simModel.mu*simModel.b0*simModel.pd(:);
end

%% check if we have motion model, and the number of NEX is > 1
motionPattern = {'rotational','translational'};
if contains(lower(expControl.motionSpecs.pattern),motionPattern)
    % simulate all the NEX and add result to solution
    simulationNEX = pulseSequence.NEX;
else
    % single NEX is enough
    simulationNEX = 1;
end

%% copy current time Solution
currentTimeSolution = timeSolution;

%% loop on the NEX, generate the motion model, and simulate
for currentNEX = 1:simulationNEX
    
    %% generate the motion model for the current NEX
    [motionModel] = motion.generateMotionSequence( ...
        pulseSequence, currentNEX, expControl);
    
    %% call the simulator
    if simModel.numIsochromats > 0
        switch lower(expControl.simulation.simulationEngine)
            case lower('phasor')
                [currentTimeSolution] = simulator.bloch.kernel.runPhasor( ...
                    currentTimeSolution, pulseSequence, simModel, motionModel, expControl, GPUindex);
            case lower('diffusion')
                [currentTimeSolution] = simulator.diffusion.kernel.runSDW(...
                    currentTimeSolution, pulseSequence, simModel, motionModel, expControl, GPUindex);
            case lower('analytical')
                [currentTimeSolution] = simulator.bloch.kernel.runAnalytical(...
                    currentTimeSolution, pulseSequence, simModel, motionModel, expControl, GPUindex);
            otherwise % 'numerical' ODE for RF, analytical for RO/GR and DF
                expControl.simulation.simulationEngine = 'numerical';
                [currentTimeSolution] = simulator.bloch.kernel.runNumerical(...
                    currentTimeSolution, pulseSequence, simModel, motionModel, expControl, GPUindex);
        end
    end
    
    %% process solution
    timeSolution.Sx = timeSolution.Sx + currentTimeSolution.Sx;
    timeSolution.Sy = timeSolution.Sy + currentTimeSolution.Sy;
    timeSolution.Sz = timeSolution.Sz + currentTimeSolution.Sz;
    
end

%% average NEX and reshape
timeSolution.Sx = reshape( timeSolution.Sx./simulationNEX, ...
    simModel.numRxCoils, pulseSequence.numRxs).';
timeSolution.Sy = reshape( timeSolution.Sy./simulationNEX, ...
    simModel.numRxCoils, pulseSequence.numRxs).';
timeSolution.Sz = reshape( timeSolution.Sz./simulationNEX, ...
    simModel.numRxCoils, pulseSequence.numRxs).';

%% report
if  expControl.debug.debugMode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : %s simulation done for sequence %s',...
        functionName,expControl.simulation.simulationEngine,pulseSequence.name);
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end

