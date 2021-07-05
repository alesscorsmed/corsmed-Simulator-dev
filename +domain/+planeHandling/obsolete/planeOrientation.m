function [b_orient,t_orient,l_orient,r_orient, planeType] = ...
    planeOrientation(point_bl_recon,point_br_recon,point_tl_recon,...
    point_tr_recon)

RefPoint = point_bl_recon;

point_bl_recon2 = point_bl_recon-RefPoint;
point_br_recon2 = point_br_recon-RefPoint;
point_tl_recon2 = point_tl_recon-RefPoint;
point_tr_recon2 = point_tr_recon-RefPoint;

% The nornal vector of the plane
planeNormalVector = cross(point_bl_recon2-point_br_recon2,point_bl_recon2-point_tl_recon2);
%disp(planeNormalVector);

% Check whether the plane is 1. not oblique, 2. single oblique or 3. double
% oblique. Check the dot product of the plane's normal vector against the
% normal vectors of the 3 main axes.
dotX = dot(planeNormalVector,[1,0,0]);
dotY = dot(planeNormalVector,[0,1,0]);
dotZ = dot(planeNormalVector,[0,0,1]);

dotArray = [dotX,dotY,dotZ];
nonZeroDotProd = find(dotArray);

% Identify the orientation of the plane and make the necessary rotations
if size(nonZeroDotProd,2) == 3
    planeType = 'oblique';
    % Center point of plane
    planeCenter     = (point_tl_recon + point_br_recon)/2;
    point_bl_recon1 = point_bl_recon-planeCenter;
    %disp(point_bl_recon1);
    point_br_recon1 = point_br_recon-planeCenter;
    %disp(point_br_recon1);
    point_tl_recon1 = point_tl_recon-planeCenter;
    %disp(point_tl_recon1);
    point_tr_recon1 = point_tr_recon-planeCenter;
    %disp(point_tr_recon1);

    % The nornal vector of the plane
    planeNormalVector = cross(point_bl_recon1-point_br_recon1,point_bl_recon1-point_tl_recon1);
    %disp(planeNormalVector);
    
    % Find its closest axis vector
    absNormalV  = abs(planeNormalVector);
    index1      = find(absNormalV == max(absNormalV));
    
    closestAxisUnitVector           = zeros(1,3);
    closestAxisUnitVector(1,index1) = sign(planeNormalVector(1,index1));
    %disp(closestAxisUnitVector);
    
    % Find the rotation matrix
    r = vrrotvec(planeNormalVector,closestAxisUnitVector);
    rotationMatrix = vrrotvec2mat(r);
    %disp(rotationMatrix);

    point_bl_recon2 = rotationMatrix*point_bl_recon1';
    %disp(point_bl_recon2);
    point_br_recon2 = rotationMatrix*point_br_recon1';
    %disp(point_br_recon2);
    point_tl_recon2 = rotationMatrix*point_tl_recon1';
    %disp(point_tl_recon2);
    point_tr_recon2 = rotationMatrix*point_tr_recon1';
    %disp(point_tr_recon2);
    
    % The NEW nornal vector of the plane
    planeNormalVector2 = cross(point_bl_recon2-point_br_recon2,point_bl_recon2-point_tl_recon2);
    %disp(planeNormalVector2);
    
    dotX2 = dot(planeNormalVector2,[1,0,0]);
    dotY2 = dot(planeNormalVector2,[0,1,0]);
    dotZ2 = dot(planeNormalVector2,[0,0,1]);

    dotArray2 = [dotX2,dotY2,dotZ2];
    index2 = find(abs(dotArray2)<0.000001);
    dotArray2(1,index2) = 0;
    
    rotationMatrix = domain.planeHandling.findRotationMatrix(dotArray2,...
        planeNormalVector,point_bl_recon2,point_br_recon2,point_tl_recon2,point_tr_recon2);
    %disp(rotationMatrix);
    
    point_bl_recon3 = rotationMatrix*point_bl_recon2;
    %disp(point_bl_recon3);
    point_br_recon3 = rotationMatrix*point_br_recon2;
    point_tl_recon3 = rotationMatrix*point_tl_recon2;
    point_tr_recon3 = rotationMatrix*point_tr_recon2;
    %disp(point_tr_recon3);
    
    RefPoint = point_bl_recon3;

    point_bl_recon3 = point_bl_recon3-RefPoint;
    point_br_recon3 = point_br_recon3-RefPoint;
    point_tl_recon3 = point_tl_recon3-RefPoint;
    point_tr_recon3 = point_tr_recon3-RefPoint;
    %disp(point_tr_recon3);
    
    % Bring the center of the plane at the center of the coordinate system
    % Find its closest axis vector
    % Find the rotation matrix that brings the normal vector to the closest
    % axis vector
    % https://stackoverflow.com/questions/25825464/get-closest-cartesian-axis-aligned-vector-in-javascript
    % https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
elseif size(nonZeroDotProd,2) == 2
    planeType = 'single oblique';
    rotationMatrix = domain.planeHandling.findRotationMatrix(dotArray,...
        planeNormalVector,point_bl_recon2,point_br_recon2,point_tl_recon2,point_tr_recon2);
    point_bl_recon3 = rotationMatrix*point_bl_recon2';
    point_br_recon3 = rotationMatrix*point_br_recon2';
    point_tl_recon3 = rotationMatrix*point_tl_recon2';
    point_tr_recon3 = rotationMatrix*point_tr_recon2';
elseif size(nonZeroDotProd,2) == 1
    planeType = 'not oblique';
    point_bl_recon3 = point_bl_recon2';
    point_br_recon3 = point_br_recon2';
    point_tl_recon3 = point_tl_recon2';
    point_tr_recon3 = point_tr_recon2';
end

% After rotation, check the dot product of the rotated-plane's normal 
% vector against the normal vectors of the 3 main axes and identify its
% position
% The nornal vector of the plane
planeNormalVector3 = cross(point_bl_recon3-point_br_recon3,point_bl_recon3-point_tl_recon3);
%disp(planeNormalVector3);

% Check whether the plane is 1. not oblique, 2. single oblique or 3. double
% oblique. Check the dot product of the plane's normal vector against the
% normal vectors of the 3 main axes.
dotX3 = dot(planeNormalVector3,[1,0,0]);
dotY3 = dot(planeNormalVector3,[0,1,0]);
dotZ3 = dot(planeNormalVector3,[0,0,1]);

dotArray3   = [dotX3,dotY3,dotZ3];
index3      = find(abs(dotArray3)<0.000001);
dotArray3(1,index3) = 0;
nonZeroDotProd3     = find(dotArray3);
if nonZeroDotProd3 == 1    
    planeType = 'Sagittal';
elseif nonZeroDotProd3 == 2    
    planeType = 'Coronal';
elseif nonZeroDotProd3 == 3   
    planeType = 'Transversal';
end

% Center point of plane
planeCenter     = (point_tl_recon3 + point_br_recon3)/2;

pBottomCenter   = (point_bl_recon3 + point_br_recon3)/2;
pLeftCenter     = (point_bl_recon3 + point_tl_recon3)/2;
pTopCenter      = (point_tl_recon3 + point_tr_recon3)/2;
pRightCenter    = (point_br_recon3 + point_tr_recon3)/2;

% Relative position of the points to the center of the plane
relpBottomCenter   = pBottomCenter - planeCenter;
%disp(relpBottomCenter);
relpLeftCenter     = pLeftCenter - planeCenter;
relpTopCenter      = pTopCenter - planeCenter;
relpRightCenter    = pRightCenter - planeCenter;
%disp(relpRightCenter);

% TEMP - Zero the small components of the vector
absrelpBottomCenter         = abs(relpBottomCenter);
index3                      = find(absrelpBottomCenter ~= max(absrelpBottomCenter));
relpBottomCenter(index3,1)  = 0;

absrelpLeftCenter           = abs(relpLeftCenter);
index4                      = find(absrelpLeftCenter ~= max(absrelpLeftCenter));
relpLeftCenter(index4,1)    = 0;

absrelpTopCenter            = abs(relpTopCenter);
index5                      = find(absrelpTopCenter ~= max(absrelpTopCenter));
relpTopCenter(index5,1)     = 0;

absrelpRightCenter          = abs(relpRightCenter);
index6                      = find(absrelpRightCenter ~= max(absrelpRightCenter));
relpRightCenter(index6,1)   = 0;

b_orient = domain.planeHandling.pointOrientation(relpBottomCenter,'bottom of the image');
t_orient = domain.planeHandling.pointOrientation(relpTopCenter,'top of the image');
l_orient = domain.planeHandling.pointOrientation(relpLeftCenter,'left of the image');
r_orient = domain.planeHandling.pointOrientation(relpRightCenter,'right of the image');