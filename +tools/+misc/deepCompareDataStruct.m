function [fail,msg] = deepCompareDataStruct(targetStruct, originStruct, ...
    relTol, throwError )
%
% TOOLS.MISC.deepCompareDataStruct
%
%	recursive function to compare all fields of a struct, 
%   and for array fields verifies that the norm of the error is small
%
% INPUT
%
%
% OUTPUT
%
%
%========================  CORSMED AB © 2020 ==============================
%
if (nargin < 3) || isempty(relTol)
    relTol = 1e-6;
end
if (nargin < 4) || isempty(throwError)
    throwError = 0;
end

% initial
fail    = 0;
msg     = 'PASS';

%% check if it is a structure:
if isstruct( originStruct )
       %% it is a structure: loop on the fields and recursive call
       
       originFields = fieldnames(originStruct);
       for ff = 1:length(originFields)
           %% check if exist on the target
           if isfield(targetStruct,originFields{ff})
               % recursive call
               [recfail,msg] = tools.misc.deepCompareDataStruct( ...
                   targetStruct.(originFields{ff}), ...
                   originStruct.(originFields{ff}), ...
                   relTol, throwError );
               % if fail comparison, throw
               fail = fail + recfail;
               if fail
                   msg = sprintf('Field %s : %s', originFields{ff},msg);
                   if throwError
                       ME = MException('error:verification', '%s', msg);
                       throw(ME);
                   else
                       fprintf(1, '\n %s', msg);
                   end
               end
               
           else
               %% we do not have this field in target: report
               msg = sprintf('Origin structure does not exists in Target');
               fail = fail + 1;
               return;
           end
               
       end
    
else
    if iscell(originStruct)
        %% cell array: loop on entries and recursive call
        
        % check number of elements
        numCellEntries = numel(originStruct);
        if numel(targetStruct) ~= numCellEntries
            fail = fail +1;
            msg = sprintf('Cell varible with different size');
            return;
        end
        
        % loop and call
        for kk = 1:numCellEntries
            [recfail,msg] = tools.misc.deepCompareDataStruct(...
                targetStruct{kk} ,...
                originStruct{kk} ,...
                relTol, throwError );
            fail = fail + recfail;
            if fail
                return;
            end
        end
        
    elseif isnumeric(originStruct) && numel(originStruct) > 1
        %% numeric vector: compare with tolerance
        refNorm = originStruct;
        refNorm = norm(refNorm(:));
        refMax  = max(abs(originStruct(:)));
        errData = targetStruct - originStruct;
        errNorm = norm(errData(:));
        errMax  = max(abs(errData(:)));
        if (errNorm > relTol*refNorm) && (errMax > 1e-2*relTol*refMax)
            msg = sprintf('Numeric array differs (error Max %g, Norm %g , reference Max %g, Norm %g)',...
                errMax, errNorm, refMax, refNorm);
            fail = fail + 1;
            return;
        end
        
    else
        
        %% not numeric vector: compare with isequal (unless function or NaN)
        if ~strcmpi(class(targetStruct),'function_handle')
            if isnan(targetStruct)
                if ~isnan(originStruct)
                    msg = sprintf('Variable differs (NaN)');
                    fail = fail +1;
                    return;
                end
            else
                if ~isequal(targetStruct,originStruct)
                    msg = sprintf('Variable differs');
                    fail = fail +1;
                    return;
                end
            end
        end
    end
end

if fail
    msg = sprintf(' Structures are different ');
    if throwError
        ME = MException('error:verification', '%s', msg);
        throw(ME);
    else
        fprintf(1, '\n %s', msg);
    end
    fail = 0;
else
    msg = sprintf(' Structures are equivalent ');
end
