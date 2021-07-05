function [fx,JG] = runFun( x, returnType, pdActive, r1Active, r2Active, ...
    solution, sequence, simModel, simControl, dbgControl )
%
% SPINTWIN.SBRBLOCH.RUNFUN
%
%	Runs Forward function with derivatives
%      Given a sequence and model, uses x to populate pd, r1, r2 entries
%      Simulates the model and depending on the returnType,
%       returns the residual in fx, 
%               and either empty, the Gradient, or Jacobian in JG 
%
% INPUT
%   x               variables in flattened array
%   returnType      'residual', 'gradient', 'jacobian'
%   pdActive        flag to determine if this parameter is active
%   r1Active
%   r2Active   
%   solution        solution struct with initial data
%   sequence        sequence to simulate
%   simModel        structure with slice data from spinModel 
%   motionModel     struct with motion model
%   simControl      control struct for the simulation
%
% OUTPUT
%   fx              objective function: Sref - S
%   JG              either Gradient or Jacobian
%
% NOTE
%   solution        solution struct with filled data
%                    (Input / Output) Magnetizations at the center of the voxel
%                     solution.Mx  array with voxel's x component magnetization
%                     solution.My  array with voxel's y component magnetization
%                     solution.Mz  array with voxel's z component magnetization
%                    (Input / Output) Spatial derivatives of the Phase
%                     solution.dPx  d(Phase)/dx
%                     solution.dPy  d(Phase)/dy
%                     solution.dPz  d(Phase)/dz
%                    (Input / Output) Parameter derivatives: derivative w.r.t. R1 = 1/T1
%                     solution.dR1x    x component: d(Mx)/dR1
%                     solution.dR1y    y component: d(My)/dR1
%                     solution.dR1z    z component: d(Mz)/dR1
%                    (Input / Output) Parameter derivatives: derivative w.r.t. R2 = 1/T2
%                     solution.dR2x    x component: d(Mx)/dR2
%                     solution.dR2y    y component: d(My)/dR2
%                     solution.dR2z    z component: d(Mz)/dR2
%                    (Input / constant) Reference (acquired) signal: flatten array with order numRxCoils x numRxs
%                     solution.Srefx   x component of signal
%                     solution.Srefy   y component of signal
%                    (Output) Integrated signal: numRxCoils x numRxs
%                     solution.Sx   x component of signal
%                     solution.Sy   y component of signal
%                    (Output) Residual signal (Sref - S): numRxCoils x numRxs
%                     solution.Rx   x component
%                     solution.Ry   y component
%                    EMPTY: (Output) Gradients: numIsochromats  x numRxCoils
%                     solution.GPr  gradient for real part of Proton Density
%                     solution.GPi  gradient for imag part of Proton Density
%                     solution.GR1  gradient for R1 (1/T1)
%                     solution.GR2  gradient for R2 (1/T2)
%                    (Output) Jacobians: numIsochromats x numRxCoils x numRxs
%                     solution.JPDx  x jacobian for Proton Density, from d(Mx)/d(PD)
%                     solution.JPDy  y jacobian for Proton Density, from d(My)/d(PD)
%                     solution.JR1x  x jacobian for R1 (1/T1), from d(Mx)/d(R1)
%                     solution.JR1y  y jacobian for R2 (1/T2), from d(My)/d(R1)
%                     solution.JR2x  x jacobian for R2 (1/T1), from d(Mx)/d(R2)
%                     solution.JR2y  y jacobian for R2 (1/T2), from d(My)/d(R2)
%                     
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'spinTwin:sbrBloch:runFun';
% check args
if (nargin < 9) || isempty(x) ...
        || isempty(solution) || isempty(simModel) ...
        || isempty(sequence) || isempty(simControl)
    ME = MException(['error:',functionName],...
        'Wrong arguments.');
    throw(ME);
end
% optional arguments
if (nargin < 10) || isempty(dbgControl)
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

%% assign variables to model
numIsochromats = simModel.numIsochromats;
numRxCoils     = simModel.numRxCoils;
numRxs         = sequence.numRxs;
numVars = 0;
if pdActive
    simModel.pr(:) = x(numVars+1:numVars+numIsochromats);
    numVars = numVars + numIsochromats;
    simModel.pi(:) = x(numVars+1:numVars+numIsochromats);
    numVars = numVars + numIsochromats;
end
if r1Active
    simModel.r1(:) = x(numVars+1:numVars+numIsochromats);
    numVars = numVars + numIsochromats;
end
if r2Active
    simModel.r2(:) = x(numVars+1:numVars+numIsochromats);
    numVars = numVars + numIsochromats;
end

switch lower(returnType)
    
    case 'gradient'
        
        %% call the function to return the gradients
        [solution, ~] = spinTwin.sbrBloch.runGradients( solution, ...
            sequence, simModel, simControl, dbgControl );
        
        %% assemble the overall gradient
        % note that we assign transposed to make things easier to assemble
        JG = zeros(numRxCoils,numVars);
        % for each active parameter retrieve the gradient
        numVars = 0;
        if pdActive
            JG(:,numVars+1:numVars+numIsochromats) = transpose(solution.GPr);
            numVars = numVars + numIsochromats;
            JG(:,numVars+1:numVars+numIsochromats) = transpose(solution.GPi);
            numVars = numVars + numIsochromats;
        end
        if r1Active
            JG(:,numVars+1:numVars+numIsochromats) = transpose(solution.GR1);
            numVars = numVars + numIsochromats;
        end
        if r2Active
            JG(:,numVars+1:numVars+numIsochromats) = transpose(solution.GR2);
            numVars = numVars + numIsochromats;
        end
        % transpose and reshape to flattened format
        JG = reshape(transpose(JG),[],1);       
        
    case 'jacobian'
        
        %% call the function to return the jacobians
        [solution, ~] = spinTwin.sbrBloch.runJacobians( solution, ...
            sequence, simModel, simControl, dbgControl );
        
        
        %% assemble the overall jacobian
        %   in terms of Equations:
        %     | JPDx  -JPDy  JR1x   JR2x   | * | Pr | = | Sx |
        %     | JPDy   JPDx  JR1y   JR2y   |   | Pi |   | Sy |
        %                                      | R1 |
        %                                      | R2 |
        %
        %   individual jacobians come in numIso,numRxCoils,numRxs
        JG = zeros(numRxs,2,numRxCoils,numVars);
        % for each active parameter retrieve the Jacobian
        numVars = 0;
        if pdActive
            JG(:,1,:,numVars+1:numVars+numIsochromats) = permute(solution.JPDx,[3,2,1]);
            JG(:,2,:,numVars+1:numVars+numIsochromats) = permute(solution.JPDy,[3,2,1]);
            numVars = numVars + numIsochromats;
            JG(:,1,:,numVars+1:numVars+numIsochromats) = -permute(solution.JPDy,[3,2,1]);
            JG(:,2,:,numVars+1:numVars+numIsochromats) =  permute(solution.JPDx,[3,2,1]);
            numVars = numVars + numIsochromats;
        end
        if r1Active
            JG(:,1,:,numVars+1:numVars+numIsochromats) = permute(solution.JR1x,[3,2,1]);
            JG(:,2,:,numVars+1:numVars+numIsochromats) = permute(solution.JR1y,[3,2,1]);
            numVars = numVars + numIsochromats;
        end
        if r2Active
            JG(:,1,:,numVars+1:numVars+numIsochromats) = permute(solution.JR2x,[3,2,1]);
            JG(:,2,:,numVars+1:numVars+numIsochromats) = permute(solution.JR2y,[3,2,1]);
            numVars = numVars + numIsochromats;
        end
        % reshape to correct format
        JG = reshape(JG,[],numVars); 
        
    otherwise
        
        %% call the function to return the gradients for now
        [solution, ~] = spinTwin.sbrBloch.runGradients( solution, ...
            sequence, simModel, simControl, dbgControl );
        JG = [];
        
end

%% assemble the objective function
%  flattened array with [real_coil1; imag_coil1; real_coil2; imag_coil2; ...]
%  note that Rx and Ry are in numRxCoils x numRxs
fx = reshape( transpose([solution.Rx, solution.Ry]),[],1);

%% report
if  dbgControl.mode
    tTotal = toc(tTotal);
    fprintf(fid, '\n');
    fprintf(fid, '\n%s : done ', functionName);
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end
