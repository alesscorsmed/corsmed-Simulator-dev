function [sessionData] = connectDatabase(sessionData,tagsStruct)
%
% EDUTOOL.BACKEND.CONNECTDATABASE
%
%	builds the connection with the Database 
%
% INPUT
%
%   sessionData    struct with session and instance info
%   tagsStruct     struct with necessary data from AWS tags
%
% OUTPUT
%
%   sessionData     struct with Data necessary for start of Edutool
%
%========================  CORSMED AB © 2020 ==============================
%
%
functionName = 'eduTool.backend.connectDatabase';
% time it
tTotal = tic();
fprintf(1, '\n%s : start', functionName);
if nargin<2
    tagsStruct.Key = '';
end

%% start collecting
if ispc
    fid = fopen('C:\My_Documents\SHARED_VM_UBUNTU_EDUCATIONAL_TOOL\IMPORTANT\edt_details.txt','r');
    urlstring = textscan(fid,'%s');
    fclose(fid);
else
    credentialsFileFilePath = '/efs-mount-point/IMPORTANT/';
    struct_index_txtFile   = find(strcmp({tagsStruct.Key}, 'CredentialsFile')==1);
    if ~isempty(struct_index_txtFile)
        edtDetailsFile   = [credentialsFileFilePath,...
            tagsStruct(struct_index_txtFile).Value];
    else
        if sessionData.developmentUse == 1
            % specific file for developers
            edtDetailsFile = sprintf('%sedt_details_%d.txt',...
                credentialsFileFilePath,sessionData.userID);
        else
            edtDetailsFile = [credentialsFileFilePath,'edt_details.txt'];
        end
    end

    fprintf(1, '\n%s : reading %s', functionName, edtDetailsFile);
    fid = fopen(edtDetailsFile,'r');
    urlstring = textscan(fid,'%s');
    fclose(fid);
end

IndexLocalDB = strfind(urlstring{1}, '#local_db');
isOne = cellfun(@(x)isequal(x,1),IndexLocalDB);
row_LocalDB = find(isOne);

IndexRemoteDB = strfind(urlstring{1}, '#remote_db');
isOne = cellfun(@(x)isequal(x,1),IndexRemoteDB);
row_RemoteDB = find(isOne);

IndexTestMode = strfind(urlstring{1}, '#test_mode');
isOne = cellfun(@(x)isequal(x,1),IndexTestMode);
row_TestMode = find(isOne);

IndexCommonFolder = strfind(urlstring{1}, '#common_folder');
isOne = cellfun(@(x)isequal(x,1),IndexCommonFolder);
row_CommonFolder = find(isOne);

IndexKernelsFolder = strfind(urlstring{1}, '#kernels_ptx_folder');
isOne = cellfun(@(x)isequal(x,1),IndexKernelsFolder);
row_KernelFolder = find(isOne);

IndexResultsFolder = strfind(urlstring{1}, '#results_folder');
isOne = cellfun(@(x)isequal(x,1),IndexResultsFolder);
row_ResultsFolder = find(isOne);

IndexGadgetronFolder = strfind(urlstring{1}, '#gadgetron_folder');
isOne = cellfun(@(x)isequal(x,1),IndexGadgetronFolder);
row_GadgetronFolder = find(isOne);

IndexGadgetronResultsFolder = strfind(urlstring{1}, '#gadgetron_results_folder');
isOne = cellfun(@(x)isequal(x,1),IndexGadgetronResultsFolder);
row_GadgetronResultsFolder = find(isOne);

IndexGadgetronISMRMRDFolder = strfind(urlstring{1}, '#gadgetron_ISMRMRD_folder');
isOne = cellfun(@(x)isequal(x,1),IndexGadgetronISMRMRDFolder);
row_GadgetronISMRMRDFolder = find(isOne);

IndexPulseSequenceFolder = strfind(urlstring{1}, '#pulseSequenceFolder');
isOne = cellfun(@(x)isequal(x,1),IndexPulseSequenceFolder);
row_PulseSequenceFolder = find(isOne);

IndexAnalyticalReconstructionFolder = strfind(urlstring{1}, '#reconstruction_analytical_results_folder');
isOne = cellfun(@(x)isequal(x,1),IndexAnalyticalReconstructionFolder);
row_AnalyticalReconstructionFolder = find(isOne);    

%% local DB
localdb_url     = char(urlstring{1}(row_LocalDB+1,1));
localdb_db      = char(urlstring{1}(row_LocalDB+2,1));
localdb_user    = char(urlstring{1}(row_LocalDB+3,1));
localdb_pass    = char(urlstring{1}(row_LocalDB+4,1));

if strcmp(localdb_db,'-')
    localdb_db = '';
end

if strcmp(localdb_pass,'-')
    localdb_pass = '';
end

%% remote DB
remotedb_url     = char(urlstring{1}(row_RemoteDB+1,1));
remotedb_db      = char(urlstring{1}(row_RemoteDB+2,1));
remotedb_user    = char(urlstring{1}(row_RemoteDB+3,1));
remotedb_pass    = char(urlstring{1}(row_RemoteDB+4,1));

if strcmp(remotedb_db,'-')
    remotedb_db = '';
end

if strcmp(remotedb_pass,'-')
    remotedb_pass = '';
end

%% store remote DB info
sessionData.remoteDB.db          = remotedb_db;
sessionData.remoteDB.user        = remotedb_user;
sessionData.remoteDB.pass        = remotedb_pass;
sessionData.remoteDB.driver      = [];
sessionData.remoteDB.dblocalurl	 = remotedb_url;

%% get folders
sessionData.folderSystem.testMode                = char(urlstring{1}(row_TestMode+1,1));
sessionData.folderSystem.commonFolder            = char(urlstring{1}(row_CommonFolder+1,1));
sessionData.folderSystem.kernelFolder            = char(urlstring{1}(row_KernelFolder+1,1));
sessionData.folderSystem.resultsFolder           = char(urlstring{1}(row_ResultsFolder+1,1));
sessionData.folderSystem.gadgetronFolder         = char(urlstring{1}(row_GadgetronFolder+1,1));
sessionData.folderSystem.gadgetronResultsFolder  = char(urlstring{1}(row_GadgetronResultsFolder+1,1));
sessionData.folderSystem.gadgetronISMRMRDFolder  = char(urlstring{1}(row_GadgetronISMRMRDFolder+1,1));
sessionData.folderSystem.pulseSequenceFolder     = char(urlstring{1}(row_PulseSequenceFolder+1,1));
sessionData.folderSystem.analyticReconFolder     = char(urlstring{1}(row_AnalyticalReconstructionFolder+1,1));    

%% CONNECT TO DATABASES
driver = 'com.mysql.cj.jdbc.Driver';

% javaclasspath([pwd,filesep,'mysql-connector-java-8.0.12',filesep,'mysql-connector-java-8.0.12.jar']);

if ispc
    javaclasspath('mysql-connector-java-8.0.12/mysql-connector-java-8.0.12.jar');
else
    javaclasspath('/efs-mount-point/IMPORTANT/mysql-connector-java-8.0.12.jar');
end

%% Local DB info
dblocalurl          = ['jdbc:mysql://',localdb_url];
sessionData.connLocalDB = ...
    database(localdb_db,localdb_user,localdb_pass,driver,dblocalurl);

sessionData.localDB.db          = localdb_db;
sessionData.localDB.user        = localdb_user;
sessionData.localDB.pass        = localdb_pass;
sessionData.localDB.driver      = driver;
sessionData.localDB.dblocalurl	= dblocalurl;

%% report
fprintf(1, ...
    '\n%s : done, elapsed time %.3fs',...
    functionName, toc(tTotal));
fprintf(1, '\n  User       : %s', sessionData.localDB.user);
fprintf(1, '\n  Driver     : %s', sessionData.connLocalDB.driver);
fprintf(1, '\n  URL        : %s', sessionData.connLocalDB.URL);
fprintf(1, '\n  DataSource : %s', sessionData.connLocalDB.DataSource);
fprintf(1, '\n  Message    : %s', sessionData.connLocalDB.Message);
fprintf(1, '\n  Type       : %s', sessionData.connLocalDB.Type);
fprintf(1, '\n');
