function [B1,Dph,Dfr] = addRF_AdiabaticTanTanh_HalfPassage(B1max,fmax,z,k,rf_duration,dt)
% It returns the rf amplitude B1(t) and phase Dph modulation functions
% of the adiabatic pulse as it is described in the paper by Garwood M, 
% JMR, 1991 "Symmetric Pulses to induce arbitrary flip angles with 
% compensation for RF inhomogeneity and Resonance Offsets"
% It returns also the freq. modulation matrix which is zero-filled since
% this pulse is phase modulated only.

timesteps = round(rf_duration/dt);
angle = 90;
angle_step = angle/timesteps;

t=0:dt:rf_duration-dt;

B1 = zeros(1,timesteps);
Dph = zeros(1,timesteps);

phmax = (-2*pi*fmax*rf_duration/(k*tan(k)))*log(cos(k));

B1 = B1max*tanh(z*t/rf_duration);
Dph = phmax-(2*pi*fmax*rf_duration/(k*tan(k)))*log(cos(k*(1-t/rf_duration))/cos(k));

figure();
subplot(1,2,1), plot(t,B1), title('Magnitude modulation');
axis tight
subplot(1,2,2), plot(t,Dph,'r'), title('Phase modulation');
axis tight

figure();
[haxes,hline1,hline2] = plotyy(t,B1/B1max,t,Dph/fmax,'plot','plot');
% grid on

axes(haxes(1))
% ylabel('B1(t)/B1max','rot',0)
set(hline1,'LineWidth',3)
set(haxes(1),'FontSize',12)
ylim([-0.1 1.1])

axes(haxes(2))
% ylabel('Äù(t)/Äùmax','rot',0)
set(hline2,'LineWidth',3)
% 'LineStyle','--',
set(haxes(2),'FontSize',12)
ylim([-1.2 1.2])

% xlabel('Time') 
% title('Adiabatic Amplitude/Frequency modulation')
% legend('RF frequency modulation','RF amplitude modulation') 

Dfr = zeros(size(Dph));
end

