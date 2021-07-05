function [targetStruct] = deepCellRestore(originStruct)
%
% TOOLS.MISC.deepCellRestore
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
    currentData = originStruct.(originFields{ff});
    % apply recursion if the entry is a struct
    if isstruct( currentData )
        currentData = tools.misc.deepCellRestore( currentData );
    end
    % check if is a cell array to restore
    if contains(originFields{ff},'_cell_')
        % get name and number
        nameParts = strsplit(originFields{ff},'_cell_');
        % generate the entry as cell
        targetStruct.(nameParts{1}){str2num(nameParts{2})} = ...
            currentData;
    else
        % not cell, assign directly
        targetStruct.(originFields{ff}) = currentData;
    end
end

