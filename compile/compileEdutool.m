function compileEdutool(edutoolVersion,env,versionMajor,versionMinor,versionBug)
% This function goes through the following steps:
%   1. Move older versions in the obsolete subfolder and keep up to 10
%       previous versions
%   2. Structure the mcc command. The mcc command inludes all the folders
%       within the edutool main folder (mainFolder) that start with the + 
%       sign and extra files  that are needed for the mcc
%   3. Compile edutool into the target folder (targetFolder)
%   4. Transfer files/folders from the main folder (mainFolder) into the 
%       target folder (targetFolder)

%
% INPUTS:
% edutoolVersion: Specifies the version of the edutool ('v1' for standalone
%   and 'v2' for jsonstandalone)
% env           : Specifies the environment ('production' or 'integration')
% versionMajor  : Specifies the major version number
% versionMinor  : Specifies the minor version number
% versionBug    : Specifies the bug fix maintenance release number
%
% EXAMPLE:
% compileEdutool('v1','integration','2','4','12')

%% Preparation

% Define the version of the edutool
timeStamp       = datestr(now,'yyyy-mm-dd HH:MM:SS');
timeStamp       = strrep(timeStamp,'-','');
timeStamp       = strrep(timeStamp,':','');
timeStamp       = strrep(timeStamp,' ','');

version.major   = versionMajor;
version.minor   = versionMinor;
version.bug     = versionBug;
version.build   = timeStamp;

% Folders
mainFolder      = '/efs-mount-point-MATLAB/cxanthis/corsmed-Simulator/';
targetFolder    = ['/efs-mount-point/compiled_projects/edutool2/',...
    edutoolVersion,'_',env];

% Find the latest active .sh file
latestEdutool   = dir([targetFolder,filesep,'run_*.sh']);

% Find the number of the latest active version
latestVerTxt    = dir([targetFolder,filesep,'version_*.txt']);
if size(latestVerTxt,1)
    latestVersion = latestVerTxt(1).name(9:end-4);
else
    latestVersion = '0_0_0_0';
end

disp(['Latest version: ',latestVersion])
%% Obsolete folder
% Transfer previous version into the obsolete folder. Keep up to 10 
% previous versions plus the previous one
versionsToKeep = 10;
obsoleteFolder = [targetFolder,filesep,'obsolete'];

fprintf(1,'%s   ','Transfer the latest version into the obsolete folder')

if ~isfolder(obsoleteFolder)
    mkdir(obsoleteFolder)
end

previousEdutool = dir([obsoleteFolder,filesep,'run_*.sh']);
[~,idx]         = sort([previousEdutool.datenum]);
previousEdutool = previousEdutool(idx);

% Delete previous versions
if size(previousEdutool,1)>versionsToKeep
    for i=1:(size(previousEdutool,1)-versionsToKeep)
        filename = previousEdutool(i).name;
        delete([obsoleteFolder,filesep,filename])
        delete([obsoleteFolder,filesep,filename(5:end-3)])
    end
end

% Transfer the latest one into the obsolete folder
for j=1:size(latestEdutool,1)
    copyfile([targetFolder,filesep,latestEdutool(j).name],...
        [obsoleteFolder,filesep,latestEdutool(j).name(1:end-3),'_',latestVersion,'.sh'])
    copyfile([targetFolder,filesep,latestEdutool(j).name(5:end-3)],...
        [obsoleteFolder,filesep,latestEdutool(j).name(5:end-3),'_',latestVersion])
end
fprintf(1,'%s\n','DONE')
%% Add all the + folders into the mcc command
k               = dir(mainFolder);

% find all the entries with the + symbol in their name
idx             = ~cellfun('isempty',strfind({k.name},'+'));  %#ok
plusFolders     = k(idx);

stringPlusFolders = '';
for i=1:size(plusFolders,1)
    stringPlusFolders = strjoin({stringPlusFolders,strjoin({'-a',...
        [mainFolder,plusFolders(i).name]})});
end

%% Add extra files needed for the mcc
extraString = '-a /efs-mount-point-MATLAB/cxanthis/corsmed-Simulator/test.dcm';

%% Structure the mcc command
currentVersion = [version.major,'_',version.minor,'_',version.bug,'_',version.build];

compileMatlabFunction = ['mcc ',...
    '-o mainEdutool_',edutoolVersion,'_',env,' ',...
    '-W main:mainEdutool2 ',...
    '-T link:exe ', ...
    '-d ',targetFolder,' ',...
    '-v /efs-mount-point-MATLAB/cxanthis/corsmed-Simulator/mainEdutool2.m',...
    strjoin({stringPlusFolders,extraString})];


%% Execute mcc
fprintf(1,'mcc command:\n%s\n',compileMatlabFunction)
fprintf(1,'Compilation just started:\n\n')
[~,cmdout] = system(compileMatlabFunction);

fprintf(1,'%s\n\n',cmdout)

%% Store the version number as the name of a .txt file
fclose(fopen([targetFolder,filesep,'version_',currentVersion,'.txt'], 'w'));

%% Delete the previous .txt file
if size(latestVerTxt,1)
    delete([targetFolder,filesep,latestVerTxt(1).name]);
end

%% Copy files/folders into the target folder 
copyfile([mainFolder,'inputs'],[targetFolder,filesep,'inputs'])
fprintf(1,'%s\n','Files transferred successfully from the main folder to the target folder')