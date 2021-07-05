function [pulseSequence,encoding] = oldPulseSequenceInterface( ...
         acquisition, encoding, mrSystem, expControl )
%
% SEQUENCE.OLDPULSESEQUENCEINTERFACE
%
%	Interface to generate pulse sequences based on.
%
% INPUT
%   acquisition         
%   mrSystem        
%   expControl      
%
% OUTPUT
%   pulseSequence   pulseSequence struct with data
%
%========================  CORSMED AB © 2020 ==============================
%

%% collect data to pass as it was
conn_localdb        = expControl.connLocalDB;
experiment_id       = expControl.experimentID;
pulseq_id           = expControl.pulseqID;
pulseSeqFamilyName  = acquisition.data.pulseSeqFamilyName;
noOfSlices          = acquisition.data.numSlices;

%% get data from DB into struct_conf
[struct_conf] = fetchGlobalParams( conn_localdb, 0 ,experiment_id);

%% extract struct_exper
sqlquery_user = ['SELECT pulseq_parameters.unique_name, ',...
    'pulseq_parameter_value.selected_value FROM `pulseq_parameter_value` ',...
    'RIGHT JOIN `pulseq_parameters` ON pulseq_parameter_value.parameter_id=',...
    'pulseq_parameters.id WHERE pulseq_parameter_value.pulseq_id=',...
    num2str(pulseq_id)];
sqlquery_results = exec(conn_localdb, sqlquery_user);
b_user = fetch(sqlquery_results);
cell_user = b_user.Data; % this is a cell array. Convert it to a struct
for i=1:size(cell_user,1)
    disp(['struct_exper.',cell_user{i,1},'=',cell_user{i,2}])
    struct_exper.(cell_user{i,1}) = cell_user{i,2};
end

if struct_conf.fastAlgorithm==1 && strcmp(struct_conf.motionPattern,'none')
    struct_exper.fast= '1';
else
    struct_exper.fast= '0';   
end

%% call the sequence generator
[pulse_sequence,isInKspace,soft_crushers,N_pulse,info,dt,gamma,...
    tirl,printForImage,desired_gridStep,slabThickness,dwPulseSequence] = ...
    pulseSequenceGenerator.createPulseSequence( ...
    pulseSeqFamilyName, struct_exper, conn_localdb, ...
    struct_conf.performance_gridZ_sliceThickness,...
    struct_conf.dwell_time,...
    experiment_id, pulseq_id, noOfSlices, ...
    struct_conf.advancedNotifications);

if (strcmp(pulseSeqFamilyName,'SE') || ...
        strcmp(pulseSeqFamilyName,'IR-SE')  || ...
        strcmp(pulseSeqFamilyName,'TSE') || ...
        strcmp(pulseSeqFamilyName,'IR-TSE') || ...
        strcmp(pulseSeqFamilyName,'SE-EPI') || ...
        strcmp(pulseSeqFamilyName,'SS-FSE'))
    pulse_sequence(6,:)=zeros(size(pulse_sequence(6,:)));
end

%% convert to new format
[pulseSequence] = data.simulation.initializeSequence();

%% convert sequence
pulseSequence = sequence.converter.oldToNewPulseSequence(pulseSequence,...
    acquisition.data.pulseSeqFamilyName,pulse_sequence,isInKspace,...
    soft_crushers,dt);

%% update encoding info
[encodingMap,encodingPlan] = sequence.converter.oldKspaceToEncodingData(...
    info,acquisition.data.pulseSeqFamilyName,...
    encoding);

%% add all extra info
pulseSequence.seqNum     = acquisition.data.pulseSeqNum;
pulseSequence.familyName = pulseSeqFamilyName;
pulseSequence.name      = pulseSeqFamilyName;
pulseSequence.type      = pulseSeqFamilyName;
pulseSequence.endEvent  = 'none'; % indicates what happens at the end
pulseSequence.totalTime = pulseSequence.time(end); % total time in seconds
pulseSequence.numSteps  = length(pulseSequence.time); % number of time steps
pulseSequence.numRxs    = nnz(pulseSequence.rxSignal); % number of readout points
pulseSequence.numReps   = encoding.plan.numPE;
pulseSequence.numEnc    = length(encoding.plan.peIncidence);
if isfield(info.pulseSequence,'EchoTrainLength')
    pulseSequence.numTL = info.pulseSequence.EchoTrainLength;
else
    pulseSequence.numTL = acquisition.data.ETL;
end
pulseSequence.numShots  = ceil(pulseSequence.numReps/pulseSequence.numTL);
if isfield(info.pulseSequence,'EchoSpacing')
    pulseSequence.ESP   = info.pulseSequence.EchoSpacing;
else
    pulseSequence.ESP   =  0;
end
pulseSequence.TE        = acquisition.data.TE; % repetition TE;
pulseSequence.TR        = acquisition.data.TR; % repetition TR;
pulseSequence.TI        = acquisition.prepIR.Apply*acquisition.prepIR.TI;
pulseSequence.effTE     = acquisition.data.TE; % effective echo time
pulseSequence.effTR     = acquisition.data.TR; % effective TR
pulseSequence.effTI     = pulseSequence.TI + pulseSequence.TE; % effective TI (to center of K space)

end

%% aux function
function [struct_conf] = fetchGlobalParams(conn_localdb,dummyRun,experiment_id)

sqlquery        = ['SELECT selected_value FROM ',...
    'edt_tool_local.global_configuration WHERE',...
    ' name IN (''performance_gridZ_sliceThickness'',',...
    '''PDinhomogeneity'',''isotropic_grid'',''analytical_simulation'',',...
    '''fast_algorithm'',''coil'',''motion'',''motion_rot_angle'',',...
    '''motion_rot_freq'',''motion_trans_magn'',''motion_trans_freq'',',...
    '''motion_trans_axis'',''kernel_threads'',''kernel_blocks'',',...
    '''gridSizeOption'',''activate_cs'',''dwell_time'',',...
    '''coil_basic_or_advanced'',''deactivateGx'',''deactivateGy'',',...
    '''deactivateGz'',''simulationMode'',''simulationKernel'',',...
    '''advanced_notifications'') ORDER BY FIELD ',...
    '(name,''performance_gridZ_sliceThickness'',''PDinhomogeneity'',',...
    '''isotropic_grid'',''analytical_simulation'',''fast_algorithm'',',...
    '''coil'',''motion'',''motion_rot_angle'',',...
    '''motion_rot_freq'',''motion_trans_magn'',''motion_trans_freq'',',...
    '''motion_trans_axis'',''kernel_threads'',''kernel_blocks'',',...
    '''gridSizeOption'',''activate_cs'',''dwell_time'',',...
    '''coil_basic_or_advanced'',''deactivateGx'',''deactivateGy'',',...
    '''deactivateGz'',''simulationMode'',''simulationKernel'',',...
    '''advanced_notifications'')'];
sqlquery_results    = exec(conn_localdb, sqlquery);
perfTypeInfo        = fetch(sqlquery_results);

struct_conf.performance_gridZ_sliceThickness = str2double(perfTypeInfo.Data{1,1});
struct_conf.PDinhomogeneity = str2double(perfTypeInfo.Data{2,1});
struct_conf.isotropicGrid   = str2double(perfTypeInfo.Data{3,1});
struct_conf.analyticalSim   = str2double(perfTypeInfo.Data{4,1});
struct_conf.fastAlgorithm   = str2double(perfTypeInfo.Data{5,1});
struct_conf.coilType        = perfTypeInfo.Data{6,1};
struct_conf.motionPattern   = perfTypeInfo.Data{7,1};

struct_conf.motionSpecs.rot_angle   = str2double(perfTypeInfo.Data{8,1});
struct_conf.motionSpecs.rot_freq    = str2double(perfTypeInfo.Data{9,1});
struct_conf.motionSpecs.trans_magn  = str2double(perfTypeInfo.Data{10,1});
struct_conf.motionSpecs.trans_freq  = str2double(perfTypeInfo.Data{11,1});
struct_conf.motionSpecs.trans_axis  = str2double(perfTypeInfo.Data{12,1});

struct_conf.threads     = str2double(perfTypeInfo.Data{13,1});
struct_conf.blocks      = str2double(perfTypeInfo.Data{14,1});

struct_conf.gridSizeOption  = perfTypeInfo.Data{15,1};
struct_conf.activateCS      = str2double(perfTypeInfo.Data{16,1});
struct_conf.dwell_time      = perfTypeInfo.Data{17,1};
struct_conf.coilMode        = perfTypeInfo.Data{18,1};  % It refers to the option coil_basic_or_advanced

struct_conf.deactivateGx    = str2double(perfTypeInfo.Data{19,1});
struct_conf.deactivateGy    = str2double(perfTypeInfo.Data{20,1});
struct_conf.deactivateGz    = str2double(perfTypeInfo.Data{21,1});

struct_conf.simulationMode      = perfTypeInfo.Data{22,1};
struct_conf.simulationKernel    = perfTypeInfo.Data{23,1};

struct_conf.advancedNotifications   = str2double(perfTypeInfo.Data{24,1});

%% Overwite global parameters if experiment is set through batch process
sqlquery        = ['SELECT batch_id FROM pulse_sequence WHERE exper_id=',num2str(experiment_id)];
sqlquery_results    = exec(conn_localdb, sqlquery);
qConf               = fetch(sqlquery_results);

if(~isnan(qConf.Data{1,1}))    
        sqlquery        = ['SELECT configurations FROM experiments WHERE id=',...
        num2str(experiment_id)];
        sqlquery_results    = exec(conn_localdb, sqlquery);
        params              = fetch(sqlquery_results);


        cell_qConf = split(params.Data,';');
        for i=1:size(cell_qConf,1)       
            attr{i} = split(cell_qConf{i,1},'=');
        end   
    
        struct_conf.performance_gridZ_sliceThickness = str2double(attr{1}{2});
        struct_conf.PDinhomogeneity = str2double(attr{2}{2});
        struct_conf.isotropicGrid   = str2double(attr{3}{2});
        struct_conf.analyticalSim   = str2double(attr{4}{2});
        struct_conf.fastAlgorithm   = str2double(attr{5}{2});
        struct_conf.coilType        = attr{6}{2};
        struct_conf.motionPattern   = attr{7}{2};

        struct_conf.motionSpecs.rot_angle   = str2double(attr{8}{2});
        struct_conf.motionSpecs.rot_freq    = str2double(attr{9}{2});
        struct_conf.motionSpecs.trans_magn  = str2double(attr{10}{2});
        struct_conf.motionSpecs.trans_freq  = str2double(attr{11}{2});
        struct_conf.motionSpecs.trans_axis  = str2double(attr{12}{2});

        struct_conf.threads     = str2double(attr{13}{2});
        struct_conf.blocks      = str2double(attr{14}{2});

        struct_conf.gridSizeOption  = attr{15}{2};
        struct_conf.activateCS      = attr{16}{2};
        struct_conf.dwell_time      = attr{17}{2};
        struct_conf.coilMode        = attr{18}{2};  % It refers to the option coil_basic_or_advanced

        struct_conf.deactivateGx    = str2double(attr{19}{2});
        struct_conf.deactivateGy    = str2double(attr{20}{2});
        struct_conf.deactivateGz    = str2double(attr{21}{2});

        struct_conf.simulationMode      = attr{22}{2};
        struct_conf.simulationKernel    = attr{23}{2};

        struct_conf.advancedNotifications   = str2double(attr{24}{2});
end

end



