function [p1_zero,p2_FE,p3_PE,point_bl_recon,point_br_recon,...
    point_tl_recon,point_tr_recon] = ...
    planePosition(point_tl_frontend,point_tr_frontend,...
    point_bl_frontend,point_br_frontend,point_tl,point_tr,...
            point_bl,point_br)
% When images were acquired from the same plane but with opposite FE and PE
% directions, the anatomy within the reconstructed images were slightly on 
% different posisions. In order to fix this issue (Jira issue ED-393), this
% function has been developed. Regardless of where the FE and PE directions
% are, this function will always return the correct orientation of the slice.
% Check notes on Onenote (Educational Tool -> DOCUMENTATION
% -> fixPositionOfPlane).

pointsOfPlane = [point_tl_frontend;point_tr_frontend;point_bl_frontend;point_br_frontend];
pointsOfPlane_2norm = vecnorm(pointsOfPlane,2,2);

closestToZeroPoint = find(pointsOfPlane_2norm==min(pointsOfPlane_2norm));

% Find the two vectors that define the plane and calculate its normal
% vector
vector1 = point_tl_frontend - point_tr_frontend;
vector2 = point_br_frontend - point_tr_frontend;
crossVector = cross(vector1,vector2);

% Calculate the testVector so as to find out if the plane "looks" toward
% one direction or the opposite
testVector = point_tr_frontend + crossVector;

% If the normal vector points towards the center of the system execute the
% commands in the if statement. If not, the plane is flipped and you should
% execute the commands in the else statement.
if norm(testVector)<norm(point_tr_frontend)
    if closestToZeroPoint==2 || closestToZeroPoint==3

        p1_zero = point_br;
        p2_FE   = point_bl;
        p3_PE   = point_tr;

        point_bl_recon = point_br;
        point_br_recon = point_bl;
        point_tl_recon = point_tr;
        point_tr_recon = point_tl;

    else

        p1_zero = point_tl;
        p2_FE   = point_tr;
        p3_PE   = point_bl;

        point_bl_recon = point_tl;
        point_br_recon = point_tr;
        point_tl_recon = point_bl;
        point_tr_recon = point_br;

    end
else
    if closestToZeroPoint==1 || closestToZeroPoint==4

        p1_zero = point_br;
        p2_FE   = point_bl;
        p3_PE   = point_tr;

        point_bl_recon = point_br;
        point_br_recon = point_bl;
        point_tl_recon = point_tr;
        point_tr_recon = point_tl;

    else

        p1_zero = point_tl;
        p2_FE   = point_tr;
        p3_PE   = point_bl;

        point_bl_recon = point_tl;
        point_br_recon = point_tr;
        point_tl_recon = point_bl;
        point_tr_recon = point_br;

    end
end

