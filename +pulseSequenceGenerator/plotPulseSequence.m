timeaxis = [1:N_pulse]*dt;

a1 = subplot(7,1,1);
plot(timeaxis,pulse_sequence(1,:),'r')
ylabel('RF magnitude') % in Tesla

a2 = subplot(7,1,2);
plot(timeaxis,pulse_sequence(2,:),'r')
ylabel('RF phase') % in radians

a3 = subplot(7,1,3);
plot(timeaxis,pulse_sequence(3,:),'r')
ylabel('RF Freq.') % in Hz

a4 = subplot(7,1,4);
plot(timeaxis,pulse_sequence(4,:))
if exist('times_fitting_single_point')
    hold on
    plot(times_fitting_single_point*dt,pulse_sequence(4,times_fitting_single_point),'go')
end
ylabel('Gx') % in T/m

a5 = subplot(7,1,5);
plot(timeaxis,pulse_sequence(5,:))
if exist('times_fitting_single_point')
    hold on
    plot(times_fitting_single_point*dt,pulse_sequence(5,times_fitting_single_point),'go')
end
ylabel('Gy') % in T/m

a6 = subplot(7,1,6);
plot(timeaxis,pulse_sequence(6,:))
ylabel('Gz') % in T/m

a7 = subplot(7,1,7);
plot(timeaxis,isInKspace,'k')
ylabel('isInKspace')

set(gcf,'Name','Pulse Sequence')
linkaxes([a1 a2 a3 a4 a5 a6 a7],'x');
set(gca,'XTick',[])