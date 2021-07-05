function expControl = fetchExperimentInfo(experiment,expControl,sessionData)
%
% EDUTOOL.FRONTEND.EXPERIMENTINFO
%
%	creates the experimentStruct with info for the upcoming experiment
%
% INPUT
%
%   sessionData         solution struct with initial data
%   experiment          info of experiment taken from the frontEnd
%   expControl          expControl structure
%
% OUTPUT
%   expControl          expControl structure
%
%========================  CORSMED AB © 2020 ==============================
%


expControl.experimentID    = experiment.Data{1,1}; 
expControl.remotedbID      = experiment.Data{1,2};
expControl.status          = experiment.Data{1,3};
expControl.experInfo       = experiment.Data{1,4};
expControl.reconstructor   = experiment.Data{1,5};
expControl.reconInfo       = experiment.Data{1,6};
expControl.pulseqMatFile   = experiment.Data{1,7};
expControl.pulseqID        = experiment.Data{1,8};

if ~isempty(expControl.pulseqMatFile)
    expControl.pathPulseSeq = [sessionData.folderSystem.pulseSequenceFolder,filesep,expControl.pulseqMatFile];
else
    expControl.pathPulseSeq = '';
end

sqlquery = ['SELECT user_id FROM',...
    ' edt_tool_local.experiments WHERE id=',num2str(expControl.experimentID)];
sqlqueryResults = exec(expControl.connLocalDB, sqlquery);
b = fetch(sqlqueryResults);
expControl.userID  = b.Data{1,1};
        
% Creation of new empty Struct
for i=1:4 
    for j=1:4
        expControl.stack_struct.slice.image.corners(i).line(j).text = "";
    end
end
            
sqlquery2           = ['SELECT selected_value FROM',...
    ' edt_tool_local.global_configuration WHERE name=''outer_fov_ratio'''];
sqlqueryResults2   = exec(expControl.connLocalDB, sqlquery2);
sqlqueryResults2   = fetch(sqlqueryResults2);

expControl.outerFOVratio       = str2num(sqlqueryResults2.Data{1,1});
        
        