%% This is the main framework for the educational tool

struct_tags         = aws.ec2.readInstanceTags();
struct_tags         = eduTool.modifyDefaultStructTags(struct_tags);

struct_application  = eduTool.initializeEduTool(struct_tags);

struct_application  = eduTool.findAnatomicalModelPath(struct_application,...
    str2num(struct_tags.Courseid));


%% Load the human anatomical model 
eduTool.frontend.runButtonEnabled(struct_application.conn_localdb,0);
eduTool.frontend.updateScannerStatus(struct_application.conn_localdb,...
    'The virtual MR scanner is booting (introducing anatomical model)');

struct_application.spinModel.anamDefaultType    = 1;  % It is the default model type of eduTool
struct_application.spinModel.anamSelectedType   = 1;  % It holds the previously-selected anatomical model type

timerStart_loadAnatModel = tic;
fprintf('Introducing anatomical model:\n            ');
struct_spinModel = load(struct_application.spinModel.defaultAnatModelPath);

timerEnd_loadAnatModel = toc(timerStart_loadAnatModel);
fprintf('-- %5.3f sec\n',timerEnd_loadAnatModel);


%% Load coil profiles