function [contrastSpace,contrastText,isContrast] = processSWI(multiContrastSlice)

functionName = 'image.dicom.processSWI';

%% allocate
[nX,nY,~]       = size(multiContrastSlice);
contrastSpace   = zeros(nX,nY,3);

% SWI post-processing with threshold pi (assumes zero phase-background):
iMag    = conj(multiContrastSlice) .* multiContrastSlice;
iPhase  = angle(multiContrastSlice);
% phaseMaskplus = ones(size(iPhase));
phaseMaskMinus = ones(size(iPhase));
phaseMaskMinus(iPhase<0) = (iPhase(iPhase<0) + pi)/pi;
% phaseMaskplus(iPhase>0) = (iPhase(iPhase>0) - pi)/pi;

contrastSpace(:,:,1) = iMag;
contrastSpace(:,:,2) = iPhase;
contrastSpace(:,:,3) = (iMag.*phaseMaskMinus);

%% text
contrastText{1} = 'Magnitude';
contrastText{2} = 'Phase';
contrastText{3} = 'SWI';

%% is contrast or map
isContrast       = ones(3,1);
isContrast(2:3)  = 0;
