function [targetStruct] = deepCopyStruct(targetStruct, originStruct)
%
% TOOLS.MISC.deepCopyStruct
%
%	recursive function to assign all fields of a struct (deepCopy)
%
% INPUT
%
%
% OUTPUT
%
%
%========================  CORSMED AB © 2020 ==============================
%

%% get all fields of originStruct as cell array
originFields = fieldnames(originStruct);

%% loop on the field names and call recursive for structs, otherwise assign
for ff = 1:length(originFields)
    if isstruct( originStruct.(originFields{ff}) )
        % is a structure in origin
        if isfield(targetStruct,originFields{ff})
            % the same exists in the target: recursive call
            targetStruct.(originFields{ff}) = tools.misc.deepCopyStruct(...
                targetStruct.(originFields{ff}), originStruct.(originFields{ff}) );
        else
            % we do not have this field in target: copy
            targetStruct.(originFields{ff}) = originStruct.(originFields{ff});
        end
    else
        % is just a field: copy in target
        targetStruct.(originFields{ff}) = originStruct.(originFields{ff});
    end
end


