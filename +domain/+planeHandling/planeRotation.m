function [t_matrixX,t_matrixY,t_matrixZ,p2_rot_final,p3_rot_final] = ...
    planeRotation(p1,p2,p3)
% This function takes as inputs the coordinates of the three points 
% (p1, p2, p3) and returns the rotation matrices about the X, Y and Z axes.
% We assume that p1 is the reference point, p2 must lie across the X axis 
% after rotation and p3 must lie across the Y axis after rotation.
dotProduct = dot(p1-p2,p1-p3);
if dotProduct~=0
    disp(['WARNING: the points of the plane are not exactly perpendicular',...
        '(dotProduct = ',num2str(dotProduct),')'])
end
RefPoint = p1;

p1 = p1-RefPoint;
p2 = p2-RefPoint;
p3 = p3-RefPoint;

%% for justification, I compute the angle between the two vectors
% atan2d(norm(cross(v1,v2)),dot(v1,v2));

%% Rotate across Z axis
% find the angle between the projection of p2 and the X axis
p2_proj=p2;
p2_proj(3)=0;

% if p2_proj(1)<=0 && p2_proj(2)<=0
%     ZangleInDegrees = 90+abs(atand(p2_proj(1)/p2_proj(2)));
% elseif p2_proj(1)<=0 && p2_proj(2)>=0
%     ZangleInDegrees = 180+abs(atand(p2_proj(2)/p2_proj(1)));
% elseif p2_proj(1)>=0 && p2_proj(2)<=0
%     ZangleInDegrees = abs(atand(p2_proj(2)/p2_proj(1)));
% elseif p2_proj(1)>=0 && p2_proj(2)>=0
%     ZangleInDegrees = -abs(atand(p2_proj(2)/p2_proj(1)));
% end
% 
% if (isnan(ZangleInDegrees))
%     ZangleInDegrees=0;
% end
% 
% t_matrix=[cosd(ZangleInDegrees),-sind(ZangleInDegrees),0;...
%     sind(ZangleInDegrees),cosd(ZangleInDegrees),0;...
%     0,0,1]; %Z axis

zAngle   = atan2(p2_proj(2),p2_proj(1));
t_matrix = [ ...
    cos(zAngle),  sin(zAngle), 0;...
    -sin(zAngle),  cos(zAngle), 0;...
    0, 0, 1 ];

t_matrixZ=t_matrix;
p2_rot=t_matrix*p2';
p3_rot=t_matrix*p3';

%% Rotate across Y axis
% find the angle between p2_rot and the X axis

% if p2_rot(1)<=0 && p2_rot(3)<=0
%     YangleInDegrees = 180+abs(atand(p2_rot(3)/p2_rot(1)));
% elseif p2_rot(1)<=0 && p2_rot(3)>=0
%     YangleInDegrees = 90+abs(atand(p2_rot(1)/p2_rot(3)));
% elseif p2_rot(1)>=0 && p2_rot(3)<=0
%     YangleInDegrees = -abs(atand(p2_rot(3)/p2_rot(1)));
% elseif p2_rot(1)>=0 && p2_rot(3)>=0
%     YangleInDegrees = abs(atand(p2_rot(3)/p2_rot(1)));
% end
% if (isnan(YangleInDegrees))
%     YangleInDegrees=0;
% end
% 
% t_matrix=[cosd(YangleInDegrees),0, sind(YangleInDegrees);...
%     0,1,0;...
%     -sind(YangleInDegrees),0,cosd(YangleInDegrees)]; %Y axis

yAngle   = atan2(p2_rot(3),p2_rot(1));
t_matrix = [ ...
    cos(yAngle), 0, sin(yAngle);...
    0, 1, 0;
    -sin(yAngle), 0, cos(yAngle)];


t_matrixY=t_matrix;
p2_rot2=t_matrix*p2_rot;
p3_rot2=t_matrix*p3_rot;

%% Rotate across X axis
% find the angle between p3 and the Y axis

% if p3_rot2(2)<=0 && p3_rot2(3)<=0
%     XangleInDegrees = 180-abs(atand(p3_rot2(3)/p3_rot2(2)));
% elseif p3_rot2(2)<=0 && p3_rot2(3)>=0
%     XangleInDegrees = -(90+abs(atand(p3_rot2(2)/p3_rot2(3))));
% elseif p3_rot2(2)>=0 && p3_rot2(3)<=0
%     XangleInDegrees = -(-90+abs(atand(p3_rot2(2)/p3_rot2(3))));
% elseif p3_rot2(2)>=0 && p3_rot2(3)>=0
%     XangleInDegrees = -(90-abs(atand(p3_rot2(2)/p3_rot2(3))));
% end
% if (isnan(XangleInDegrees))
%     XangleInDegrees=0;
% end
% 
% t_matrix=[1,0,0;...
%     0,cosd(XangleInDegrees), -sind(XangleInDegrees);...
%     0,sind(XangleInDegrees),cosd(XangleInDegrees)]; %X axis

xAngle   = atan2(p3_rot2(3),p3_rot2(2));
t_matrix = [ 1, 0, 0;...
    0,  cos(xAngle), sin(xAngle);...
    0, -sin(xAngle), cos(xAngle)];

t_matrixX=t_matrix;
p2_rot3=t_matrix*p2_rot2;
p3_rot3=t_matrix*p3_rot2;

%% Resolve rounding issues with decimal digits 
p2_rot3=p2_rot3*1e5;
p3_rot3=p3_rot3*1e5;
p2_rot3=round(p2_rot3);
p3_rot3=round(p3_rot3);
p2_rot_final=p2_rot3*1e-5;
p3_rot_final=p3_rot3*1e-5;

% disp(['p2_rot_final = [',num2str(p2_rot_final(1)),',',num2str(p2_rot_final(2)),',',num2str(p2_rot_final(3)),']'])
% disp(['p3_rot_final = [',num2str(p3_rot_final(1)),',',num2str(p3_rot_final(2)),',',num2str(p3_rot_final(3)),']'])