function instanceData = readInstanceTags(instance_id)

%% get the info
% Read the assigned tags
command2        = ['sudo aws ec2 describe-instances --instance-ids ',...
    instance_id,' --query "Reservations[].Instances[].Tags[?Key==''simulator-instance''].Value[]"'];
fprintf('\n%s\n',command2)
[~,tags]        = system(command2);
% decode
tagsStruct      = jsondecode(tags);

%% assign info to structure
% for i=1:size(tagsStruct,2)
%     strName = regexprep(tagsStruct(i).Key,'-','');
%     instanceData.(strName) = tagsStruct(i).Value;
% end

instanceData.simulatorinstance = tagsStruct{1,1};