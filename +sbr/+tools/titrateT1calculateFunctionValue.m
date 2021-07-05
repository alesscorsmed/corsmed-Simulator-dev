%% Prepare inputs

% Load the experiment specs
fid             = fopen(['+sbr',filesep,'testInputs',filesep,...
    'testInput_fromNewSimulator_20210110.json']);
rawJson         = fread(fid,inf);
strJson         = char(rawJson');
fclose(fid);

experimentData  = jsondecode(strJson);
options         = experimentData.options;

if 1
    % Load pulse sequence
    load(['+sbr',filesep,'testInputs',filesep,...
        'bSSFP_128x128_7_6_sineFA_k_10_singleLobe_withoutaby2_200kHz_CXedit.mat']);

    pulseSequence = data.pulseSequence.initialize();
    % convert sequence
    pulseSequence = sequence.converter.oldToNewPulseSequence(...
        pulseSequence,'N/A',pulse_sequence,isInKspace,...
        soft_crushers,dt);

    [encoding.map, encoding.plan] = ...
        sequence.converter.oldKspaceToEncodingData(...
        info.pulseSequence,'');
else
    load('+sbr/testInputs/bSSFP_EduToolv2.mat')
end

spinModelGT     = sbr.tools.createCircularSpinModelGT(1,0.25,1);
sessionData     = options.auxdata.sessionData;
experiment      = [];
expControl      = data.expControl.initialize(experiment,sessionData,'sbr');

expControl.simulation.kernelPtx = ...
   '/efs-mount-point-MATLAB/EduTool-Jorge/EduTool-CUDA/v26_sm70.ptx';

motionModel.type    = 'none';

%% Run the ground-truth experiment
fprintf('%s','Running GT simulation ... ')
timeSolutionGT  = sbr.run.sbrEngine(spinModelGT,pulseSequence,...
    motionModel,expControl,1);

timeSolutionGT.noise = zeros(size(timeSolutionGT.Sy));    
% Reshape the outcome of the simulator
[kSpaceGT] = reconstructor.signal.mapKspace(timeSolutionGT,...
    encoding,expControl);
fprintf('%s\n','DONE')

%% Calculate function value for T1 titration
index   = 0;
T1range = 0.7:0.01:1.7;
val1    = zeros(size(T1range));
tic
for T1new = T1range
    index = index + 1;
    fprintf('Running simulation %d ... ',index)    
    spinModelGPU = spinModelGT;
    
    % Modify the spinModel based on the T1new
    spinModelGPU.tissueValues(1,1)  = T1new;
    
    timeSolution        = sbr.run.sbrEngine(spinModelGPU,pulseSequence,...
        motionModel,expControl,1);
    
    timeSolution.noise  = zeros(size(timeSolution.Sy));
    
    % Reshape the outcome of the simulator
    [kSpace] = reconstructor.signal.mapKspace( timeSolution,...
        encoding, expControl );
    
    val1(index) = sum(abs((kSpace(:)-kSpaceGT(:))./max(kSpaceGT(:))).^2);
    fprintf('%s\n','DONE')
    
end
toc
plot(T1range,val1)
xlabel('T1 values')
ylabel('Function value')

%% Calculate function value for T2 titration
% index   = 0;
% T2range = 0.1:0.001:0.4;
% val2    = zeros(size(T2range));
% for T2new = T2range
%     index = index + 1;
%     fprintf('Running simulation %d ... ',index)    
%     spinModelGPU = spinModelGT;
%     
%     % Modify the spinModel based on the T1new
%     spinModelGPU.tissueValues(1,2)  = T2new;
%     
%     timeSolution        = sbr.run.sbrEngine(spinModelGPU,pulseSequence,...
%         motionModel,expControl,1);
%     
%     timeSolution.noise  = zeros(size(timeSolution.Sy));
%     
%     % Reshape the outcome of the simulator
%     [kSpace] = reconstructor.signal.mapKspace( timeSolution,...
%         encoding, expControl );
%     
%     val2(index) = sum(abs((kSpace(:)-kSpaceGT(:))./max(kSpaceGT(:))).^2);
%     fprintf('%s\n','DONE')
%     
% end
% 
% plot(T2range,val2)
% xlabel('T2 values')
% ylabel('Function value')

%% Calculate function value for PD titration
% index   = 0;
% PDrange = 0.9:0.001:1.1;
% val3    = zeros(size(PDrange));
% for PDnew = PDrange
%     index = index + 1;
%     fprintf('Running simulation %d ... ',index)    
%     spinModelGPU = spinModelGT;
%     
%     % Modify the spinModel based on the T1new
%     spinModelGPU.tissueValues(1,3)  = PDnew;
%     PDval                           = spinModelGPU.tissueValues(:,3);
%     spinModelGPU.pd                 = transpose(PDval(spinModelGPU.tissueType));
%     
%     timeSolution        = sbr.run.sbrEngine(spinModelGPU,pulseSequence,...
%         motionModel,expControl,1);
%     
%     timeSolution.noise  = zeros(size(timeSolution.Sy));
%     
%     % Reshape the outcome of the simulator
%     [kSpace] = reconstructor.signal.mapKspace( timeSolution,...
%         encoding, expControl );
%     
%     val3(index) = sum(abs((kSpace(:)-kSpaceGT(:))./max(kSpaceGT(:))).^2);
%     fprintf('%s\n','DONE')
%     
% end
% 
% plot(PDrange,val3)
% xlabel('PD values')
% ylabel('Function value')