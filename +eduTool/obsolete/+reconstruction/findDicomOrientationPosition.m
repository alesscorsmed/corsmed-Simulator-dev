function [ImagePositionPatient,ImageOrientationPatient,PixelSpacing] = ...
    findDicomOrientationPosition(outputImage)

matrix_x    = size(outputImage.image,2);  % defines the size of k-space along x direction
matrix_y    = size(outputImage.image,1);  % defines the size of k-space along y direction
FOV_x_mm    = outputImage.fovFE*1000;
FOV_y_mm    = outputImage.fovPE*1000;

PixelSpacing = [FOV_x_mm/matrix_x;FOV_y_mm/matrix_y];

point_tl_mm = outputImage.LTop*1000;
point_tr_mm = outputImage.RTop*1000;
point_bl_mm = outputImage.LBot*1000;

Dh = point_tr_mm - point_tl_mm;  % direction cosine of the first row
Dv = point_bl_mm - point_tl_mm;  % direction cosine of the first column

% Convert them to unit vectors
Dh = Dh./norm(Dh);
Dv = Dv./norm(Dv);

ImagePositionPatient = point_tl_mm + PixelSpacing(1,1)/2*Dh + PixelSpacing(2,1)/2*Dv;
ImagePositionPatient = ImagePositionPatient';

ImageOrientationPatient = [Dh';Dv'];