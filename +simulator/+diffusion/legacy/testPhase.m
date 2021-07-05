clc
Bpeak = 2.4303e-05;
alpha = Bpeak*dt*2*pi*gamma;
theta = 0.010*40e-3*dt*2*pi*gamma;
phi = 0;

cphi = cos(phi);
sphi = sin(phi);
calpha = cos(alpha);
salpha = sin(alpha);
ctheta = cos(theta);
stheta = sin(theta);

mx=0; my=0; mz=1; mm = 0; mp = 0;

aux1 = cphi*cphi + sphi*sphi*calpha;
aux2 = sphi*sphi + cphi*cphi*calpha;
aux3 = cphi*sphi*(1- calpha);
aux4 = sphi*salpha;
aux5 = cphi*salpha;


tempx = aux1*mx + aux3*my -   aux4*mz;
tempy = aux3*mx + aux2*my +   aux5*mz;
mz    = aux4*mx - aux5*my + calpha*mz;
mx = tempx;
my = tempy;

mp = atan2(my,mx);
mm = abs(mx+1j*my);

fprintf(1, '\n mx = %.4f, my = %.4f, mm = %.5f, mp = %.2f \n', mx,my,mm,mp);


for ii = 1:600

mp = mp - theta;
mx = mm*cos(mp);
my = mm*sin(mp);

fprintf(1, '\n %d Rot: mx = %.4f, my = %.4f, mm = %.5f, mp = %.2f', ii, mx,my,mm,mp);

tempx = aux1*mx + aux3*my -   aux4*mz;
tempy = aux3*mx + aux2*my +   aux5*mz;
mz    = aux4*mx - aux5*my + calpha*mz;
mx = tempx;
my = tempy;

localmp = atan2(my,mx);
mm = abs(mx+1j*my);
deltaP = localmp - rem(mp,2*pi);
% if deltaP > pi
%     deltaP = deltaP - 2*pi;
% end
% if deltaP < -pi
%     deltaP = deltaP + 2*pi;
% end
mp = mp + deltaP;

fprintf(1, '\n %d RF:  mx = %.4f, my = %.4f, mm = %.5f, mp = %.2f', ii, mx,my,mm,localmp);
fprintf(1, '\n %d IT:  mx = %.4f, my = %.4f, mm = %.5f, mp = %.2f (delta = %.3f) \n', ii, mx,my,mm,mp,deltaP);

if deltaP > pi
     break
end

end

