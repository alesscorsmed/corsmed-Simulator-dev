clear all; close all; clc;

%% create structure with nested cell arrays
imageData.name          = 'test';
imageData.bodyPartName  = 'head';
imageData.is3D          = 0;
imageData.numSlices     = 3;
imageData.numContrasts  = 2;

% info of image coordinates for each slice
for ss = 1:imageData.numSlices
    
    % dimensions in plane (image FOV)
    plane.fovX    = 0.2;
    plane.fovY    = 0.2;
    plane.sizeX   = 128;  % image pixels
    plane.sizeY   = 128;  % image pixels
    % coordinates
    plane.LTop = 1*[0.1; 0.2; 0.3]; % Top-Left
    plane.RTop = 2*[0.1; 0.2; 0.3]; % Top-Right
    plane.LBot = 3*[0.1; 0.2; 0.3]; % Bottom-Left
    plane.RBot = 4*[0.1; 0.2; 0.3]; % Bottom-ROrientFOV
    % extended FOV points
    plane.LTopExtFOV = 1*[0.1; 0.2; 0.3]; % Top-Left
    plane.RTopExtFOV = 2*[0.1; 0.2; 0.3]; % Top-Right
    plane.LBotExtFOV = 3*[0.1; 0.2; 0.3]; % Bottom-Left
    plane.RBotExtFOV = 4*[0.1; 0.2; 0.3]; % Bottom-ROrientFOV
    % orientation
    plane.BOrient = 'A';
    plane.TOrient = 'P';
    plane.LOrient = 'L';
    plane.ROrient = 'R';
    % flip and rotation of the image
    plane.flip = 0; % 0: no flip / 1: flip first dimension / 2 flip second
    plane.rot  = 90; % angle in degrees to flip the image (90, -90, 180...)
    % foldover direction
    plane.foldoverDir = 'X';
    
    imageData.slice{ss}.plane = plane;
    
    for cc = 1:imageData.numContrasts
        % this is filled in SLICER
        imageData.slice{ss}.contrast{cc}.imageText  = ['stuff here'];
        
        % the next are allocated, but not filled: will be filled after RECON, in DICOM
        imageData.slice{ss}.contrast{cc}.uniqueID = 1;
        imageData.slice{ss}.contrast{cc}.image = rand(128,128);
        imageData.slice{ss}.contrast{cc}.mask  = rand(128,128);
        imageData.slice{ss}.contrast{cc}.resizeFactor = 1;
        imageData.slice{ss}.contrast{cc}.bmpName = 'bmp';
        imageData.slice{ss}.contrast{cc}.txtName = 'txt';
        imageData.slice{ss}.contrast{cc}.dcmName = 'dcm';
        imageData.slice{ss}.contrast{cc}.tissueMapBmpName = 'tissueMap';
        imageData.slice{ss}.contrast{cc}.kspaceOutputComplex = 'kSpace.cplx';
        imageData.slice{ss}.contrast{cc}.kspaceTxtName = 'kSpace.txt';
        imageData.slice{ss}.contrast{cc}.kspaceBmpName = 'kSpace.bmp';
        
    end
    
end

 % structure with info for the dicom
imageData.dicomInfo = [];

jsonFileName = './testIO.json';

%% save to json
[ok] = tools.IO.saveStruct2Json(imageData, jsonFileName);

%% restore from json
[imageDataRestored] = tools.IO.loadJson2Struct(jsonFileName);

%% check correctness
isequal(imageDataRestored,imageData)




