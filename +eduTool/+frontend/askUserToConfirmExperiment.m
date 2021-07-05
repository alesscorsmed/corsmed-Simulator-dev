function askUserToConfirmExperiment(expControl,sessionData,msg)
% This function will pause the experiment and ask the user to confirm that
% he agrees with the msg. If the user has already accepted this message
% once, this will not appear again for the current experiment.
%
% Future fixes:
% If the workflow has more than one messages to present to the user, only
% the first message will appear and the rest will be omitted due to the 
% current design of the process. A potential fix could involve the
% introduction of different ids for every confirm message that we send to
% the user.

if expControl.experimentID ~= sessionData.latestConfirmedExperimentID
    eduTool.frontend.updateExperimentProgress(expControl,'','confirm',msg);
end