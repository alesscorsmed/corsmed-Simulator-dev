function plotPulseSequence(pulseSequence,MRparts)
% MRparts is a cell array with the following inputs
% MRparts = {'RFm','Gx','Gy','Gz','RO','RX','PX','SWC'};
% RFm is for RF magnitude
% Gx, Gy, Gz are for gradients along the x,y and z direction
% RO is for readout points (when the receiver is on)
% RX is for the limits of the RF pulses
% PX is for the limits of the other components
% SWC is for software crushers
%
% If MRparts = {'all'}, all the abovementioned components will be printed

if nargin<2
    MRparts = {'all'};
end

legendCell = {};

figure();
if any(strcmp(MRparts,'RFm')) || any(strcmp(MRparts,'all'))
    plot(pulseSequence.time,pulseSequence.rfmSignal*1e3);
    hold on
    legendCell(end+1) = {'RF amp'};
end

if any(strcmp(MRparts,'Gx')) || any(strcmp(MRparts,'all'))
    plot(pulseSequence.time,pulseSequence.feSignal);
    hold on
    legendCell(end+1) = {'Gx'};
end

if any(strcmp(MRparts,'Gy')) || any(strcmp(MRparts,'all'))
    plot(pulseSequence.time,pulseSequence.peSignal);
    hold on
    legendCell(end+1) = {'Gy'};
end

if any(strcmp(MRparts,'Gz')) || any(strcmp(MRparts,'all'))
    plot(pulseSequence.time,pulseSequence.ssSignal);
    hold on
    legendCell(end+1) = {'Gz'};
end

if any(strcmp(MRparts,'RO')) || any(strcmp(MRparts,'all'))
    plot(pulseSequence.time(pulseSequence.rxSignal>0), ...
        pulseSequence.gxSignal(pulseSequence.rxSignal>0), 'o');
    hold on
    legendCell(end+1) = {'Readout'};
end

if any(strcmp(MRparts,'RX')) || any(strcmp(MRparts,'all'))
    plot(pulseSequence.time(pulseSequence.rxLimits(:,1)), ...
        pulseSequence.gxSignal(pulseSequence.rxLimits(:,1)), '>');
    hold on
    plot(pulseSequence.time(pulseSequence.rxLimits(:,2)), ...
        pulseSequence.gxSignal(pulseSequence.rxLimits(:,2)), '<');
    hold on
    legendCell(end+1) = {'RF start'};
    legendCell(end+1) = {'RF end'};
end

if any(strcmp(MRparts,'PX')) || any(strcmp(MRparts,'all'))
    plot(pulseSequence.time(pulseSequence.partLimits(:,1)),...
        zeros(pulseSequence.numParts,1), '^');
    hold on
    plot(pulseSequence.time(pulseSequence.partLimits(:,2)),...
        zeros(pulseSequence.numParts,1), 'v');
    hold on
    legendCell(end+1) = {'Part start'};
    legendCell(end+1) = {'Part end'};
end

if any(strcmp(MRparts,'SWC')) || any(strcmp(MRparts,'all'))
    if nnz(pulseSequence.swcSignal) > 0
        plot(pulseSequence.time(pulseSequence.swcSignal>0), ...
            zeros(nnz(pulseSequence.swcSignal),1), 's',...
            'LineWidth', 2, 'MarkerSize', 10);
    end
    legendCell(end+1) = {'SWC'};
end

if isfield(pulseSequence,'name')
    title(pulseSequence.name)
else
    title(pulseSequence.type)
end
xlabel('time (s)');
legend(legendCell)