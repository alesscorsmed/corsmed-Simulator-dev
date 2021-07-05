function [] = updateExperimentProgress(expControl,progress,status,message,...
    messageAdmins)
%
% edutool.frontend.updateExperimentProgress
%
% updates the progress of the experiment
%
% INPUT
%   expControl          structure that holds info for the experiment
%   progress            progress of the experiment
%   status              status we want to change into
%                           status may be one of the following:
%                             - error
%                             - cancelled-error
%                             - confirm
%                             - info
%                             - started
%                             - finished
%   message             message sent to the frontend
%
% OUTPUT
%   none
%
% ========================= CORSMED AB @ 2020 ============================
%

functionName = 'edutool.frontend.updateExperimentProgress';

if nargin==1 
    if ~isfield(expControl,'progress')
        ME = MException('eduTool:wrongArgCount',...
            '%s : wrong argument count',functionName);
        throw(ME);
    else
        progress        = num2str(expControl.progress);
        status          = '';
        message         = '';
        messageAdmins   = '';
    end
elseif nargin==2
    status          = 'update';
    message         = '';
    messageAdmins   = '';
elseif nargin==3
    message         = '';
    messageAdmins   = '';
elseif nargin==4
    messageAdmins   = '';
end

if isfield(expControl,'approach') && strcmp(expControl.approach,'jsonstandalone')
    if ~isempty(progress)
        status = 'update';
    end
    messages.message        = message;
    messages.messageAdmins  = messageAdmins;
    tools.updateJsonExperimentProgress(expControl,status,messages,progress);
else
    if isfield(expControl,'connLocalDB') && ~isempty(expControl.connLocalDB)
        experimentID    = expControl.experimentID;
        conn            = expControl.connLocalDB;
        
        if ~isempty(status)
            if any(strcmp(status,{'cancelled-error','confirm','info'}))
                eduTool.frontend.errorAndDBentry(conn,message,status,...
                    expControl.experimentID,expControl.pulseqID)
            else
                exec(conn, ...
                    ['UPDATE experiments SET status=''',status,''' WHERE id=',...
                    num2str(experimentID)]);
            end
        end

        if ~isempty(progress)
            eduTool.frontend.progressUpdate(expControl);
        end
    end
end