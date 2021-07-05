function [MOLLI_scheme,pause_cc,TI_initial] = ...
    extractMOLLIspecs(struct_pulseq,cc_duration)

schemeCellArray = strsplit(struct_pulseq.molli_scheme,'(');

MOLLI_scheme    = zeros(1,size(schemeCellArray,2));
pause_cc        = zeros(1,size(schemeCellArray,2)-1);

LLexper = 0;
for iLL = 1:size(schemeCellArray,2)
    LLexper = LLexper + 1;
    MOLLI_scheme(1,LLexper) = str2double(schemeCellArray{1,iLL}(end));
    if iLL>1
        pause_cc(1,LLexper-1) = str2double(schemeCellArray{1,iLL}(1));
    end
end

TIsCellArray = strsplit(struct_pulseq.tis,',');

% Remove ( and ) from the first and last cell of TIsCellArray
TIsCellArray{1,1} = TIsCellArray{1,1}(2:end);
TIsCellArray{1,end} = TIsCellArray{1,end}(1:end-1);

TIs = str2double(TIsCellArray);

% Create the array of the TIs
TI_initial = [];
for iTI = 1:size(TIs,2)
    TI_initial = [TI_initial,TIs(1,iTI),((1:MOLLI_scheme(1,iTI)-1)*cc_duration)+TIs(1,iTI)];
end