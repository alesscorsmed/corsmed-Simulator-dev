function [dicomInfo] = findDicomOrientationPosition(slicePlane, dicomInfo)
%
% IMAGE.DICOM.findDicomOrientationPosition
%
%     updates dicomInfo for the slice
%
% INPUT
%
% OUTPUT
%
%========================  CORSMED AB © 2020 ==============================

% pixel size
nX      = slicePlane.sizeX;
nY      = slicePlane.sizeY;
fovXmm  = slicePlane.fovX*1e3;
fovYmm  = slicePlane.fovY*1e3;
PixelSpacing = [fovXmm/nX; fovYmm/nY];

% orientation
point_tl_mm = slicePlane.LTop*1e3;
point_tr_mm = slicePlane.RTop*1e3;
point_bl_mm = slicePlane.LBot*1e3;
Dh = point_tr_mm - point_tl_mm;  % direction cosine of the first row
Dv = point_bl_mm - point_tl_mm;  % direction cosine of the first column

% Convert them to unit vectors
Dh = Dh./norm(Dh);
Dv = Dv./norm(Dv);

ImagePositionPatient = point_tl_mm + PixelSpacing(1)/2*Dh + PixelSpacing(2)/2*Dv;
ImagePositionPatient = ImagePositionPatient';
ImageOrientationPatient = [Dh';Dv'];

% store data
dicomInfo.foldoverDir               = slicePlane.foldoverDir;
dicomInfo.PixelSpacing              = PixelSpacing;
dicomInfo.ImagePositionPatient      = ImagePositionPatient;
dicomInfo.ImageOrientationPatient   = ImageOrientationPatient;
