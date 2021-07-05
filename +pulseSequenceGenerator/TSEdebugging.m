load('C:\Users\chris\Downloads\RawData_user_539_exper_4113.mat')
% load('C:\Users\chris\OneDrive\Desktop\tse_ETL5_CORRECT_try2.mat')
load('C:\Users\chris\OneDrive\Desktop\tse_ETL5_CORRECT_try3_TE_0p05.mat')

ETL = 5;
rBWstar         = 1/info.pulseSequence.dt;
rBW             = info.pulseSequence.rBW;
factor_BW       = round(rBWstar/rBW);

PEsteps = info.pulseSequence.Ny;

Mx = reshape(Mx_per_timestep,[round(size(Mx_per_timestep,2)/PEsteps),PEsteps]);
My = reshape(My_per_timestep,[round(size(My_per_timestep,2)/PEsteps),PEsteps]);
Mz = reshape(Mz_per_timestep,[round(size(Mz_per_timestep,2)/PEsteps),PEsteps]);

Mxy = Mx_per_timestep+My_per_timestep*exp(sqrt(-1)*pi/2);
kspaceMatrix = reshape(Mxy,[round(size(Mx_per_timestep,2)/PEsteps),PEsteps]);
kspaceMatrix = kspaceMatrix(1:factor_BW:end,:);

orderOfKspacelines      = info.pulseSequence.kspace(:,6)';
kspaceMatrix_reordered  = zeros(size(kspaceMatrix));

for i = 1:size(orderOfKspacelines,2)
    
    kspace_line = orderOfKspacelines(1,i);
    kspaceMatrix_reordered(:,kspace_line) = kspaceMatrix(:,i);    
    
end

I = fftshift(fft2(kspaceMatrix_reordered));
Inew = abs(I);
figure(10)
subplot(1,2,1)
imagesc([0 info.pulseSequence.FOVx], [0 info.pulseSequence.FOVy],Inew);
colormap gray
axis image    
title('Reconstructed image')
subplot(1,2,2)
waterfall(abs(kspaceMatrix_reordered))
title('Reordered')
xlabel('FE')
ylabel('PE')