function [pulse_sequence_new,kspace_times_new,isInKspace,soft_crushers_new] = ...
    pulseOptimization(pulse_sequence,N,kspace,soft_crushers)
% It is a copy of the NEW_pulse_optimization_20150118_fetal.m but it is
% changed to be used in the coreMRI framework
fastAlgTime = tic;
kspace_times = [kspace(:,2),kspace(:,2)+kspace(:,3)-1];

l = 1;
l_temp = 1;
start = pulse_sequence(7,1);
test1 = pulse_sequence(1:6,1);
kspace_line = 1;

% PREALLOCATE MEMORY
% for t=1:N
%     test2 = pulse_sequence(1:6,t);
%     k = isequal(test1, test2);
%     
%     if ~(k==1 && (t<kspace_times(kspace_line,1) || ...
%             t>kspace_times(kspace_line,2)) && (soft_crushers(1,t)~=1))
%         l_temp=l_temp+1;
%     end
%     
%     if t==kspace_times(kspace_line,2) && kspace_line~=size(kspace_times,1)
%         kspace_line = kspace_line + 1;
%     end
%     
%     test1 = test2;
%     
% end
% 
% pulse_sequence_new = zeros(8,l_temp);
% soft_crushers_new = zeros(1,l_temp);
% END PREALLOCATE MEMORY

% develop matrix opt
test1 = pulse_sequence(1:6,1);
kspace_line = 1;

kspace_times_new = zeros(size(kspace_times));

for t=1:N    
    test2 = pulse_sequence(1:6,t);
    k = isequal(test1, test2);
    
    % IF THERE IS NO CHANGE
    % only if the test1 is equal to test2 AND t doesn't belong to 
    % the kspace times AND it is not the end of TR
    if k==1 && (t<kspace_times(kspace_line,1) || ...
            t>kspace_times(kspace_line,2)) && (soft_crushers(1,t)~=1) 
        
       pulse_sequence_new(1:8,l) = [test1;start;pulse_sequence(8,t)];
    
    else  % IF THERE IS A CHANGE
        l=l+1;
        start = pulse_sequence(7,t);
        test1 = test2;
        pulse_sequence_new(1:8,l) = [test1;start;pulse_sequence(8,t)];
        
        soft_crushers_new(1,l) = soft_crushers(1,t);  % 20140417
    end
    
    if t==kspace_times(kspace_line,1)
        kspace_times_new(kspace_line,1)=l;
    elseif t==kspace_times(kspace_line,2)
        kspace_times_new(kspace_line,2)=l;
        if kspace_line~=size(kspace_times,1)
            kspace_line = kspace_line + 1;
        end
    end
end

% develop isInKspace matrix
kspace_line = 1;
isInKspace = zeros(1,size(pulse_sequence_new,2));
index = 0;
for i=1:size(pulse_sequence_new,2)
    if i>kspace_times_new(kspace_line,2) && kspace_line<size(kspace_times,1)
        kspace_line = kspace_line + 1;
    end
    if i>=kspace_times_new(kspace_line,1) && i<=kspace_times_new(kspace_line,2) 
        index = index + 1;
        isInKspace(1,i) = index;
    end   
end
tEnd = toc(fastAlgTime);
disp(['Fast algorith: ',num2str(tEnd),'sec'])