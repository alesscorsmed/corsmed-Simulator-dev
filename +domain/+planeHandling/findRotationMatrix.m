function rotationMatrix = findRotationMatrix(dotArray,planeNormalVector,...
    point_bl_recon,point_br_recon,point_tl_recon,point_tr_recon)

nonZeroDotProd = find(dotArray);
ZeroDotProd = find(dotArray==0);

if size(nonZeroDotProd,2) == 3  % The plane is double oblique
    pBottomCenter   = (point_bl_recon + point_br_recon)/2;
    pLeftCenter     = (point_bl_recon + point_tl_recon)/2;
    pTopCenter      = (point_tl_recon + point_tr_recon)/2;
    pRightCenter    = (point_br_recon + point_tr_recon)/2;
    
elseif size(nonZeroDotProd,2) == 2 % The plane is single oblique
    if ZeroDotProd==3  % Normal vector perpendicular to Z-axis
        v1 = [1,0,0];
        Theta1 = atan2d(norm(cross(planeNormalVector,v1)),...
            dot(planeNormalVector,v1));
        if Theta1 > 90
            Theta1 = 180 - Theta1;
        end
        if Theta1 > 45
            Theta1 = Theta1 - 90;
        end
        rotationMatrix = [cosd(Theta1),-sind(Theta1),0;...
            sind(Theta1),cosd(Theta1),0;...
            0,0,1];
    elseif ZeroDotProd==2  % Normal vector perpendicular to Y-axis
        v1 = [0,0,1];
        Theta1 = atan2d(norm(cross(planeNormalVector,v1)),...
            dot(planeNormalVector,v1));
        if Theta1 > 90
            Theta1 = 180 - Theta1;
        end
        if Theta1 > 45
            Theta1 = Theta1 - 90;
        end
        rotationMatrix = [cosd(Theta1),0, sind(Theta1);...
            0,1,0;...
            -sind(Theta1),0,cosd(Theta1)]; %Y axis
    elseif ZeroDotProd==1  % Normal vector perpendicular to X-axis
        v1 = [0,1,0];
        Theta1 = atan2d(norm(cross(planeNormalVector,v1)),...
            dot(planeNormalVector,v1));
        if Theta1 > 90
            Theta1 = 180 - Theta1;
        end
        if Theta1 > 45
            Theta1 = Theta1 - 90;
        end
        rotationMatrix = [1,0,0;...
            0,cosd(Theta1), -sind(Theta1);...
            0,sind(Theta1),cosd(Theta1)]; %X axis
    end
else
    rotationMatrix = [1,0,0;0,1,0;0,0,1];
end

