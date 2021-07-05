function [elementsize_new,alarm] = testAccumPhaseGRE(N_pulse,N_TR,dt,...
    G_pulse,elementsize)
% This function tests for each gradient during the pulse if the phase
% difference between contiguous elements along each axis is greater than
% 180 degrees (0.5 cycle). If this happens it changes the dimensions of the
% element of the object.

accum_G = zeros(3,N_TR);
alarm=zeros(1,3);
elementsize_new = zeros(1,3);

% CHECK FOR Gx------------------------------------------------------------%
accum_G(1,1) = 42.56*10^6*G_pulse(1,1)*dt*elementsize(1,1);

for i=2:N_TR
    accum_G(1,i) = accum_G(1,i-1)+42.56*10^6*G_pulse(1,i)*dt*elementsize(1,1);
end

[r1,c1] = find(accum_G(1,:)>=0.5);

if isempty(r1)
    alarm(1,1)=0;
    elementsize_new(1,1) = elementsize(1,1);
else
    alarm(1,1)=1;
    elementsize_new(1,1) = elementsize(1,1)*0.5;
end

%%% END OF Gx CHECK-------------------------------------------------------%

% CHECK FOR Gy------------------------------------------------------------%
accum_G(2,1) = 42.56*10^6*G_pulse(2,1)*dt*elementsize(1,2);

for i=2:N_TR
    accum_G(2,i) = accum_G(2,i-1)+42.56*10^6*G_pulse(2,i)*dt*elementsize(1,2);
end

[r2,c2] = find(accum_G(2,:)>=0.5);

if isempty(r2)
    alarm(1,2)=0;
    elementsize_new(1,2) = elementsize(1,2);
else
    alarm(1,2)=1;
    elementsize_new(1,2) = elementsize(1,2)*0.5;
end

%%% END OF Gy CHECK-------------------------------------------------------%

% CHECK FOR Gz------------------------------------------------------------%
accum_G(3,1) = 42.56*10^6*G_pulse(3,1)*dt*elementsize(1,3);

for i=2:N_TR
    accum_G(3,i) = accum_G(3,i-1)+42.56*10^6*G_pulse(3,i)*dt*elementsize(1,3);
end

[r3,c3] = find(accum_G(3,:)>=0.5);

if isempty(r3)
    alarm(1,3)=0;
    elementsize_new(1,3) = elementsize(1,3);
else
    alarm(1,3)=1;
    elementsize_new(1,3) = elementsize(1,3)*0.5;
end

%%% END OF Gz CHECK-------------------------------------------------------%