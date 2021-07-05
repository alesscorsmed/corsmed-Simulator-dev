function [f] = splot(seq)
% Plot the sequence trace for qc
%   plot(...,'TimeRange',[start stop]) within start and stop times

tres = median(seq.time(2:end)-seq.time(1:end-1)); % grid time for resampling on a full view

fig=figure;
if nargout>0
    f=fig;
end
ax=zeros(1,6);
for i=1:6
    ax(i)=subplot(3,2,i);
end
ax=ax([1 3 5 2 4 6]);   % Re-order axes
arrayfun(@(x)hold(x,'on'),ax);
arrayfun(@(x)grid(x,'on'),ax);
labels={'ADC/labels','RF mag (T)','RF/ADC ph (rad)','Gx (T/m)','Gy (T/m)','Gz (T/m)'};
arrayfun(@(x)ylabel(ax(x),labels{x}),1:6);

plot(seq.time(seq.rxSignal>0), seq.rxSignal(seq.rxSignal>0), "*r", 'Parent',ax(1));
plot(seq.time, abs(seq.rfmSignal), 'Parent',ax(2));
plot(seq.time, angle(exp(1i*angle(seq.rfmSignal)).*exp(1i*seq.rfpSignal).*exp(1i*2*pi*seq.time.*seq.rffSignal)), 'Parent',ax(3));
xlabel(ax(3),'t (s)');


vt = [seq.time(1)];
Gx = [seq.gxSignal(1),seq.gySignal(1),seq.gzSignal(1)];
for ii = 2:length(seq.time)
    if abs(seq.time(ii) - seq.time(ii-1)) > 2*tres
        vt = [vt; seq.time(ii-1)+tres];
        Gx = [Gx; 0.0,0.0,0.0];
    end
    vt = [vt; seq.time(ii)];
    Gx = [Gx; seq.gxSignal(ii),seq.gySignal(ii),seq.gzSignal(ii)];
end

plot(vt, Gx(:,1),'Parent',ax(4));
plot(vt, Gx(:,2),'Parent',ax(5));
plot(vt, Gx(:,3),'Parent',ax(6));
xlabel(ax(6),'t (s)');

end








                 