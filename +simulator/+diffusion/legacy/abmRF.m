function [mm, mp, mz] = abmRF(mm, mp, mz, ...
    pulseTdiff, pulseRx, pulseSWC, ...
    pulseGx, pulseGy, pulseGz,...
    pulseRfMag, pulseRfPhase,pulseRfFreq,...
    x,y,z,bi,pd,...
    tissueProperties, tissueType,...
    nIsochromats, nTimeStart, nTimeEnd,...
    gamma, b0, mu)
%
%
%

NTISSUEPROPERTIES = 6;
nTimeStart = nTimeStart+1;
wgamma = 2*pi*gamma;
mm = mm + 1e-4*mu*(b0 + bi).*pd;

% tissue properties
r1 =  1.0./( tissueProperties( NTISSUEPROPERTIES*( tissueType(:)-1 ) +1 ) );
r2 =  1.0./( tissueProperties( NTISSUEPROPERTIES*( tissueType(:)-1 ) +2 ) );
cs =  tissueProperties( NTISSUEPROPERTIES*( tissueType(:)-1 ) +4 )*1.5*0.000001;
%relaxation term
lz =  r1.*(mu*(b0 + bi)).*pd;

% ------------------------------------------------------
% step 1: ABM order 1 (FE)

% evaluate time delta for next time step
dt = pulseTdiff(nTimeStart)*1e-4;

% evaluate fields at current timestep
sphase = sin(pulseRfPhase(nTimeStart));
cphase = cos(pulseRfPhase(nTimeStart));

wbx =  wgamma*pulseRfMag(nTimeStart)*cphase;
wby =  wgamma*pulseRfMag(nTimeStart)*sphase;
wbz =  wgamma*(cs + bi + x*pulseGx(nTimeStart) + y*pulseGy(nTimeStart) + z*pulseGz(nTimeStart));

% evaluate function
% feval( prevF, mm, mp, mz, wbx, wby, wbz, r1, r2, lz, dt );
% convert to mx, my
sphase = sin(mp);
cphase = cos(mp);
mx = mm.*cphase;
my = mm.*sphase;
% evaluate the matrix and function
mxtemp = - r2.*mx + wbz.*my - wby.*mz;
mytemp = -wbz.*mx -  r2.*my + wbx.*mz;
mztemp =  wby.*mx - wbx.*my -  r1.*mz + lz;
% correct back to polar
prevFm  =  mxtemp.*cphase + mytemp.*sphase;
prevFp  =  mm;
prevFp( prevFp == 0 ) = dt*prevFm( prevFp == 0 );
prevFp( prevFp == 0 ) = dt;
prevFp  = (mytemp.*cphase - mxtemp.*sphase)./ prevFp;
prevFz  =  mztemp;

% FE update of solution for nextStep
mm = mm + dt*prevFm;
mp = mp + dt*prevFp;
mz = mz + dt*prevFz;

% ------------------------------------------------------
% rest of steps: ABM order 2
for timestep = nTimeStart:nTimeEnd
    
    % evaluate time delta for next time step
    dt = pulseTdiff(timestep);
    h = 0.5*dt;
    
    % evaluate fields at current timestep
    sphase = sin(pulseRfPhase(timestep));
    cphase = cos(pulseRfPhase(timestep));
    
    wbx =  wgamma*pulseRfMag(timestep)*cphase;
    wby =  wgamma*pulseRfMag(timestep)*sphase;
    wbz =  wgamma*(cs + bi + x*pulseGx(timestep) + y*pulseGy(timestep) + z*pulseGz(timestep));
    
    % evaluate function
    % feval( prevF, mm, mp, mz, wbx, wby, wbz, r1, r2, lz, dt );
    % convert to mx, my
    sphase = sin(mp);
    cphase = cos(mp);
    mx = mm.*cphase;
    my = mm.*sphase;
    % evaluate the matrix and function
    mxtemp = - r2.*mx + wbz.*my - wby.*mz;
    mytemp = -wbz.*mx -  r2.*my + wbx.*mz;
    mztemp =  wby.*mx - wbx.*my -  r1.*mz + lz;
    % correct back to polar
    currFm  =  mxtemp.*cphase + mytemp.*sphase;
    currFp  =  mm;
    currFp( currFp == 0 ) = dt*currFm( currFp == 0 );
    currFp( currFp == 0 ) = dt;
    currFp = (mytemp.*cphase - mxtemp.*sphase)./ currFp;
    currFz  =  mztemp;
    
    deltaP = 3*h*currFp - h*prevFp;
    repeats = 1;
    while any( abs(deltaP) > 1e-1 )
        repeats = repeats+1;
        h = 0.5*dt/repeats;
        deltaP = 3*h*currFp - h*prevFp;
    end
    
    for rr = 1:repeats
        % ABM2 update of solution for nextStep
        mm = mm + 3*h*currFm - h*prevFm;
        mp = mp + 3*h*currFp - h*prevFp;
        mz = mz + 3*h*currFz - h*prevFz;
    end
    
%     figure(10);clf;
%     plot(3*h*currFp - h*prevFp,'LineWidth', 2 );
%     
%     figure(20);clf; hold on;
%     plot(dt*single(1:nTimeEnd),pulseRfMag, 'LineWidth', 2 );
%     plot(dt*single(timestep),pulseRfMag(timestep), 'v', 'LineWidth', 2 );
%     xlabel(sprintf(' time (repeats = %d, h=%.3f ns, current: t = %.3f ms, s = %.3f nT)', repeats, h*1e9, dt*single(timestep)*1e3, pulseRfMag(timestep)*1e6));
%     
%     pause();
%     
    % update data for next step
    prevFm = currFm;
    prevFp = currFp;
    prevFz = currFz;
    
end
