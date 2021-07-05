function [ok] = saveStruct2Json( myStruct, jsonFileName )
%
% TOOLS.IO.SAVESTRUCT2JSON
%
%	saves a structure into a json file, removing cell arrays
%
% INPUT
%
%
% OUTPUT
%
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'tools.IO.saveStruct2Json';
try
    % remove cell arrays
    [myStruct] = tools.misc.deepCellRemove(myStruct);
    % json Encode
    jsonData = jsonencode(myStruct);
    % json file write
    fid = fopen(jsonFileName,'w');
    fwrite(fid, jsonData, 'char');
    fclose(fid);
    % all fine
    ok = 1;
catch
    % fail
    ok = 0;
    % throw error
    ME = MException('json:Error',...
        '%s : Structure could not be saved to json file %s', ...
        functionName, jsonFileName);
    throw(ME);
end