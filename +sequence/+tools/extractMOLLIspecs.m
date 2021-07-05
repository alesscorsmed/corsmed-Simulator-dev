function [schemeMOLLI,pauseMOLLI,TIs] = ...
    extractMOLLIspecs(acquisitionMOLLI,cc_duration)

schemeCellArray = strsplit(acquisitionMOLLI.scheme,'(');

schemeMOLLI     = zeros(1,size(schemeCellArray,2));
pauseMOLLI    	= zeros(1,size(schemeCellArray,2)-1);

LLexper = 0;
for iLL = 1:size(schemeCellArray,2)
    LLexper = LLexper + 1;
    schemeMOLLI(1,LLexper) = str2double(schemeCellArray{1,iLL}(end));
    if iLL>1
        pauseMOLLI(1,LLexper-1) = str2double(schemeCellArray{1,iLL}(1));
    end
end

TIsCellArray = strsplit(acquisitionMOLLI.TIs,',');

% Remove ( and ) from the first and last cell of TIsCellArray
TIsCellArray{1,1} = TIsCellArray{1,1}(2:end);
TIsCellArray{1,end} = TIsCellArray{1,end}(1:end-1);

TIs = str2double(TIsCellArray)/1000;  % convert it to sec