function [signal_final,BW,total_timesteps]=addRFsinc(total_dur,cycles,...
    angle,dt,gamma)

sidelobes=(cycles-1)*2;
total_lobes = sidelobes + 2;  % because the main lobe has width of 2 times 
                              % the width of the sidelobe
total_timesteps = round(total_dur/dt);

sinc_step = total_lobes/total_timesteps;  % because the sinc function is 
                              % created in cycles. So, since I have 4
                              % cycles in total, I divide this by the
                              % timesteps that I want my sinc function to
                              % have
                              
t1 = 0:sinc_step:total_lobes-sinc_step;

signal = zeros(1,1,total_timesteps);
signal_new = zeros(1,total_timesteps);

for i=1:total_timesteps
    signal(1,1,i) = sinc(t1(1,i)-(2+sidelobes)/2);
end

for i=1:total_timesteps
    signal_new(1,i) = signal(1,1,i);
end

% t2=(0:total_timesteps-1)*dt;
% figure();  % plot1
% plot(t2,signal_new)

% ------------------------ RF magnitude --------------------------------- %
sum_s = sum(signal_new);
B1_magn = angle/(360*gamma*sum_s*dt);

signal_final = B1_magn.*signal_new;


% ---------------- Accumulated angle of the RF pulse -------------------- %
angle_per_timestep = 360*signal_final*gamma*dt;

accum_angle = zeros(1,total_timesteps);
for i=2:total_timesteps
    accum_angle(1,i)=accum_angle(1,i-1)+angle_per_timestep(1,i);
end

%----------------------- Calculate the RF BW ---------------------------- %
%More info on how to calculate the bandwidth (BW) of the rf pulse (sinc)
%see pg.37 on the book "Handbook of MRI Pulse Sequences"

t0 = total_dur/(cycles*2);
BW = 1/t0;

disp(['RF pulse bandwidth: ', int2str(BW),'Hz'])