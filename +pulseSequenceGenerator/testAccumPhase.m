function [elementsize_new,alarm,alarm_count] = testAccumPhase(TRmin,dt,...
    G_pulse,pulse_sequence_timings,elementsize,gamma)
% This function tests for each gradient during the pulse if the phase
% difference between contiguous elements along each axis is greater than
% 180 degrees (0.5 cycle). If this happens it changes the dimensions of the
% element of the object.

accum_G         = zeros(3,size(pulse_sequence_timings,2));
diffTime        = (pulse_sequence_timings(2,:)-pulse_sequence_timings(1,:)+1)*dt;
accum_G(1,:)    = cumsum(diffTime.*G_pulse(1,:).*(gamma*elementsize(1,1)));
accum_G(2,:)    = cumsum(diffTime.*G_pulse(2,:).*(gamma*elementsize(1,2)));
accum_G(3,:)    = cumsum(diffTime.*G_pulse(3,:).*(gamma*elementsize(1,3)));

maxPhase        = max(accum_G');
reductFactor    = ceil(maxPhase/0.5);
alarm_count     = sum(ceil(maxPhase/0.5));

elementsize_new = elementsize./reductFactor;

alarm           = reductFactor;