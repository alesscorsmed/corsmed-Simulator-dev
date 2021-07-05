function iMapCropped = cropImageAlongFE(iMap,reconstruction,parallelImaging)
% Crop image along the FE so as to avoid the appearance of foldover
% artifacts in this direction. This is a simulator-based workaround for the
% application of low-pass filtering in the acquired signal in true MR
% systems.

reconNx = reconstruction.ismrmrd.header.encoding.reconSpace.matrixSize.x;
reconNy = reconstruction.ismrmrd.header.encoding.reconSpace.matrixSize.y;

if strcmp(parallelImaging,'sense')
    iMapCropped = iMap(:,round(size(iMap,2)/2 - reconNx/2)+1:...
        round(size(iMap,2)/2 + reconNx/2));
else
    iMapCropped = iMap(round(size(iMap,1)/2 - reconNy/2)+1:...
        round(size(iMap,1)/2 + reconNy/2),...
        round(size(iMap,2)/2 - reconNx/2)+1:...
        round(size(iMap,2)/2 + reconNx/2));
end