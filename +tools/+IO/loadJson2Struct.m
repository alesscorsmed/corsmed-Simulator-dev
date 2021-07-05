function [myStruct] = loadJson2Struct( jsonFileName )
%
% TOOLS.IO.LOADJSON2STRUCT
%
%	loads a structure from a json file, restoring cell arrays
%
% INPUT
%
%
% OUTPUT
%
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'tools.IO.loadJson2Struct';
try
    % load json and decode
    fid = fopen(jsonFileName,'r');
    myStruct = jsondecode(fread(fid,inf,'*char').');
    fclose(fid);
    % restore
    [myStruct] = tools.misc.deepCellRestore(myStruct);
catch

    % throw error
    ME = MException('json:Error',...
        '%s : Structure could not be loaded from json file %s', ...
        functionName, jsonFileName);
    throw(ME);
end