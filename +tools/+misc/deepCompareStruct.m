function [errorStruct] = deepCompareStruct(targetStruct, originStruct, ...
    errorStruct, throwError )
%
% TOOLS.MISC.deepCompareStruct
%
%	recursive function to compare all fields of a struct
%
% INPUT
%
%
% OUTPUT
%
%
%========================  CORSMED AB © 2020 ==============================
%
if (nargin < 4) || isempty(throwError)
    throwError = 0;
end

%% get all fields of originStruct as cell array
originFields = fieldnames(originStruct);
%% loop on the field names and call recursive for structs, otherwise assign
for ff = 1:length(originFields)
    if isstruct( originStruct.(originFields{ff}) )
        % is a structure in origin
        if isfield(targetStruct,originFields{ff})
            % the same exists in the target: recursive call
            if ~isfield(errorStruct,originFields{ff})
                errorStruct.(originFields{ff}) = [];
            end
            errorStruct.(originFields{ff}) = tools.misc.deepCompareStruct(...
                    targetStruct.(originFields{ff}), originStruct.(originFields{ff}),errorStruct.(originFields{ff}));
        else
            % we do not have this field in target: report
             msg = sprintf('Field %s of Origin structure does not exists in Target', originFields{ff});
             if throwError
                 ME = MException('error:verification', '%s', msg);
                 throw(ME);
             else
                fprintf(1, '\n %s', msg);
                %targetStruct.(originFields{ff}) = originStruct.(originFields{ff});
            end
        end
    else
        % is just a field: compare
        if isfield(targetStruct,originFields{ff})
            if ~isequal(targetStruct.(originFields{ff}),originStruct.(originFields{ff}))
                % error in field
                if throwError
                    msg = sprintf('Field %s differs', originFields{ff});
                    ME = MException('error:verification', '%s', msg);
                    throw(ME);
                else
                    fprintf(1, '\n Field %s differs between structures', originFields{ff});
                    errorStruct.(originFields{ff}).V1 = originStruct.(originFields{ff});
                    errorStruct.(originFields{ff}).V2 = targetStruct.(originFields{ff});
                end
            end
        else
            % we do not have this field in target: report
             msg = sprintf('Field %s of Origin structure does not exists in Target', originFields{ff});
             if throwError
                 ME = MException('error:verification', '%s', msg);
                 throw(ME);
             else
                fprintf(1, '\n %s', msg);
            end
        end
    end
end

