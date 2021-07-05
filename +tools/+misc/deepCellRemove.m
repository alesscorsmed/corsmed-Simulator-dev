function [targetStruct] = deepCellRemove(originStruct)
%
% TOOLS.MISC.deepCellRemove
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
    if iscell( originStruct.(originFields{ff}) )
        currentCell = originStruct.(originFields{ff});
        for ii = 1:length( originStruct.(originFields{ff}) )
            % apply recursion if the entry is a struct
            if isstruct( currentCell{ii} )
                currentCell{ii} = tools.misc.deepCellRemove( currentCell{ii} );
            end
            % else
            newName = sprintf('%s_cell_%d',originFields{ff},ii);
            targetStruct.(newName) = currentCell{ii};
        end
    else
        targetStruct.(originFields{ff}) = originStruct.(originFields{ff});
    end
end

