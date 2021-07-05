function [point_tl_recon,point_tr_recon,point_bl_recon,point_br_recon] = ...
    findReconPoints(point_tl,point_tr,point_bl,point_br,outerFOVratio)

Dh = point_tr - point_tl;
Dv = point_bl - point_tl;

point_tl_recon = point_tl - outerFOVratio*Dh - outerFOVratio*Dv;
point_tr_recon = point_tr + outerFOVratio*Dh - outerFOVratio*Dv;
point_bl_recon = point_bl - outerFOVratio*Dh + outerFOVratio*Dv;
point_br_recon = point_br + outerFOVratio*Dh + outerFOVratio*Dv;