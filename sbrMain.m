function [x,infoSbr] = sbrMain(jsonPath)

    %% Load the .json file and decode it
    disp('Loading experiment specs from file system')
    fid             = fopen(jsonPath);
    rawJson         = fread(fid,inf);
    strJson         = char(rawJson');
    fclose(fid);

    experimentData  = jsondecode(strJson);

    options         = experimentData.options;
    x0              = transpose(experimentData.x0);

    % Transpose the options.lb and options.ub
    options.lb      = transpose(options.lb);
    options.ub      = transpose(options.ub);
    
    %% Suppress warnings
    warning('off','parallel:gpu:device:DeviceDeprecated')
    
    %% Prepare the expControl
    % initialize expControl with sessionData
    sessionData = options.auxdata.sessionData;
    
    %% JV: CHANGED FOR COMPATIBILITY AFTER V3
    % experiment  = [];
    % expControl  = data.expControl.initialize(experiment,sessionData,'sbr');
    [expControl,~] = data.experiment.loadExperiment(sessionData,'sbr');
    
    options.auxdata.simulation.expControl = expControl;
    
    %% Prepare the anatomical model for the optimization process
    if isempty(options.auxdata.anatomicalModel.path)
        anatModelSpecs = options.auxdata.anatomicalModel.specs;
        if ~isfield(anatModelSpecs,'tissueValues')
            anatModelSpecs.tissueValues = ...
                zeros(anatModelSpecs.tissuesX*anatModelSpecs.tissuesY,6);
        end
        [spinModel,motionModel] = sbr.run.sbrSlicer(...
            anatModelSpecs.xMin,anatModelSpecs.xMax,...
            anatModelSpecs.yMin,anatModelSpecs.yMax,...
            anatModelSpecs.zMin,anatModelSpecs.zMax,...
            [anatModelSpecs.gridStepX,...
            anatModelSpecs.gridStepY,...
            anatModelSpecs.gridStepZ],...
            anatModelSpecs.tissuesX,anatModelSpecs.tissuesY);
    else
        load(options.auxdata.anatomicalModel.path)
        if ~exist('motionModel','var')
            motionModel.type = 'none';
        end            
    end

    options.auxdata.simulation.spinModel        = spinModel;
    options.auxdata.simulation.motionModel      = motionModel;
    
    %% Prepare the pulse sequence
    if isempty(options.auxdata.pulseSequence.path)
        
        error('The path for the pulse sequence is empty.')
        
%         % Generate the pulse sequence
%         [pulseSequence] = sequence.generatePulseSequence( ...
%             acquisition, mrSystem, expControl );
        
    else
        fprintf('%s','Loading anatomical model...')
        load(options.auxdata.pulseSequence.path)
        
        % If the pulse sequence is structured in the old format, convert it
        % so as to be used by the current framework
        if ~exist('pulseSequence','var')
            if exist('pulse_sequence','var')
                %% JV: CHANGED FOR COMPATIBILITY AFTER V3
                % pulseSequence = data.pulseSequence.initialize();
                [pulseSequence] = data.simulation.initializeSequence();
                % convert sequence
                pulseSequence = sequence.converter.oldToNewPulseSequence(...
                    pulseSequence,'N/A',pulse_sequence,isInKspace,...
                    soft_crushers,dt);
            else
                error(['The ',options.auxdata.pulseSequence.path,...
                ' does not include either the pulseSequence struct ',...
                '(new version) or the pulse_sequence variable (old version).'])
            end
        end
        fprintf('%s\n','DONE')
        
        % Prepare the encoding if it is not available at the 
        % options.auxdata.pulseSequence.path
        fprintf('%s','Prepare encoding...')
        if ~exist('kSpaceInfo','var')
            if exist('info','var')
                [encoding.map, encoding.plan] = ...
                    sequence.converter.oldKspaceToEncodingData(...
                    info.pulseSequence,'');
            else
                error(['If the pulse sequence is in an old format, the ',...
                    'info.pulseSequence.kspace is missing. If the pulse ',...
                    'sequence is in a new format, the kSpaceInfo struct ',...
                    'is missing.'])                
            end            
        end        
        fprintf('%s\n','DONE')
    end
    
    options.auxdata.simulation.pulseSequence    = pulseSequence;
    options.auxdata.simulation.encoding         = encoding;    
    
    
    %% Load the ground-truth kspace or run the ground truth simulation
    
    if (options.auxdata.runSimulationGT == 1)  % run the GT simulation
        
        fprintf('%s','Running the ground truth simulation...')
        % Load the anatomical model or create a new one
        if (options.auxdata.simulationGT.loadAnatModel == 1)  % load anatomical model
            load(options.auxdata.simulationGT.anatomicalPath)
            if ~exist('spinModelGT','var')
                error(['The spinModelGT is not available in the path: '...
                    options.auxdata.simulationGT.anatomicalPath])
            end
            if ~exist('motionModelGT','var')
                motionModelGT = 'none';
            end
        else  % create the anatomical model for the GT simulation
            anatModelGTSpecs = options.auxdata.simulationGT.anatomicalModel.specs;
            if ~isfield(anatModelGTSpecs,'tissueValues')
                anatModelGTSpecs.tissueValues = ...
                    zeros(anatModelGTSpecs.tissuesX*anatModelGTSpecs.tissuesY,6);
            end
            [spinModelGT,motionModelGT] = sbr.run.sbrSlicer(...
                anatModelGTSpecs.xMin,anatModelGTSpecs.xMax,...
                anatModelGTSpecs.yMin,anatModelGTSpecs.yMax,...
                anatModelGTSpecs.zMin,anatModelGTSpecs.zMax,...
                [anatModelGTSpecs.gridStepX,...
                anatModelGTSpecs.gridStepY,...
                anatModelGTSpecs.gridStepZ],...
                anatModelGTSpecs.tissuesX,anatModelGTSpecs.tissuesY,...
                anatModelGTSpecs.tissueValues);
        end
        
        % Run the GT simulation
        tic
        timeSolution    = sbr.run.sbrEngine(spinModelGT,pulseSequence,...
            motionModelGT,expControl,1);
        fprintf('%s\n','DONE')
        toc
        timeSolution.noise = zeros(size(timeSolution.Sy));    
        % Reshape the outcome of the simulator
        [kSpace] = reconstructor.signal.mapKspace(timeSolution,...
            encoding,expControl);
    else        
        fprintf('%s','Loading the ground-truth kspace...')
        load(options.auxdata.simulationGT.kSpace.path)

        if ~exist('kSpace','var')
            error(['The ',options.auxdata.simulationGT.kSpace.path,...
                ' does not include the kSpace var.'])
        end        
        fprintf('%s\n','DONE')
    end
    
    options.auxdata.simulation.kSpaceGT = kSpace;
    
    %% Define the size of the parpool
    p = parpool('local',gpuDeviceCount);
    
    %% Define the funcs for ipopt
    funcs.objective = @(x,auxdata) simulator(x,options.auxdata);
    funcs.gradient  = @(x,auxdata) gradient(x,options.auxdata);
    funcs.iterfunc  = @callback;
    

	%% Run IPOPT.
    tic
    [x, infoSbr]       = ipopt(x0,funcs,options);
    toc
    save(['outputs',filesep,datestr(now,'YYYYMMDD_hhmmss'),'_outputs.mat'],'x','infoSbr')
end

% ----------------------------------------------------------------------
function f = simulator(x,auxdata)

    fprintf('Iteration T1 = %d,T2 = %d,PD = %d - T1 = %d,T2 = %d,PD = %d - T1=%d,T2=%d,PD=%d - T1=%d,T2=%d,PD=%d\n',x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12))
    
    f = sbr.optFunctions.simulator3pNtissue(x,auxdata);
    
end

% ----------------------------------------------------------------------
function g = gradient(x,auxdata)

    fprintf('%s\n','CALLING GRADIENT OPERATOR:')
%     fprintf('Iteration T1 = %d,T2 = %d,PD = %d - T1 = %d,T2 = %d,PD = %d - T1=%d,T2=%d,PD=%d - T1=%d,T2=%d,PD=%d\n',x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12))
    gradientTimeWatcher = tic;
    g = sbr.optFunctions.gradientJacobian3pNtissue(x,auxdata);
    gradientTime = toc(gradientTimeWatcher);
    disp(['Gradient duration: ',num2str(gradientTime),'sec'])
    
end

% ----------------------------------------------------------------------
function H = hessianstructure() %#ok<DEFNU>
  H = sparse([ 1  0  0  0 
               1  1  0  0
               0  0  1  0
               0  1  1  1 ]);
end

% ----------------------------------------------------------------------
function H = hessian(x, sigma, lambda)
  H = [ 1200*x(1)^2-400*x(2)+2  0       0                          0
        -400*x(1)               220.2   0                          0
         0                      0       1080*x(3)^2- 360*x(4) + 2  0
         0                      19.8   -360*x(3)                   200.2 ];
  H = sparse(sigma*H);
end

% ----------------------------------------------------------------------
function b = callback(t, f, x)
  fprintf('%3d  %0.3g \n',t,f);
%   fmt = ['The vector P is: [', repmat('%g, ', 1, numel(x)-1), '%g]\n'];
%   fprintf(fmt, x)
  b = true;
end