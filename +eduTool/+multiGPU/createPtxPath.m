function ptxPath = createPtxPath(...
    pdInhomogeneity,motionPattern,blocksPerIteration,threadsPerBlock,simulationKernel)
%
%DYNAMICRESOURCES.createPtxPath
%
%   create ptx path based on kernel
%   and the model of GPU which leads to 
%   different sm cuda architecture
%
%
% INPUT
%   pdInhomogeneity     boolean value for use of PDinhomogeneity
%   motionPattern       choice of Motion
%   threadsPerBlock     the default threads extracted from Global Configs
%   blockerPerIteration the default blocks extracted from Global Configs
%   gpuName             the name of GPU
%   simulationKernel    the kernel used
%
% OUTPUT
%   ptxPath             the final ptx path for the compiled cuda code
%                       that will be used
%
% ========================= CORSMED AB @ 2020 ============================
%

useKernel25 =1;

%% For now use only kernel_25
if ~useKernel25 

    gpuName = gpuDevice(1).Name;

    %ptx based on the kernel we use
    if strcmp(simulationKernel,'stable')

        %cuda architecture based on GPU
        if(contains(gpuName,'V100'))
            ptxPath = ['/efs-mount-point/kernels/sm70/cuda_kernel_16_Rx_',...
                num2str(blocksPerIteration),'_',num2str(threadsPerBlock),'.ptx'];

        elseif(contains(gpuName,'K80'))
            ptxPath = ['/efs-mount-point/kernels/sm37/cuda_kernel_16_Rx_',...
               num2str(blocksPerIteration),'_',num2str(threadsPerBlock),'.ptx'];

        else
            ptxPath = ['/efs-mount-point/kernels/cuda_kernel_16_Rx_',...
               num2str(blocksPerIteration),'_',num2str(threadsPerBlock),'.ptx'];
        end

        if pdInhomogeneity
            disp('Kernel file has changed. PD inhom is taken into account.')
            ptxPath = ptxPath(1:end-4);
            ptxPath = [ptxPath,'_PDinhom.ptx'];
        end

        if strcmp(motionPattern,'translational')
            disp('Kernel file has changed. Translational motion is taken into account.')
            ptxPath = ptxPath(1:end-4);
            ptxPath = [ptxPath,'_motionTranslational.ptx'];

        elseif strcmp(motionPattern,'rotational')
            disp('Kernel file has changed. Rotational motion is taken into account.')
            ptxPath = ptxPath(1:end-4);
            ptxPath = [ptxPath,'_motionRotational.ptx'];
        end

    else

        if(contains(gpuName,'V100'))
            ptxPath = '/efs-mount-point/kernels/sm70/cuda_kernel_20_universal.ptx';

        elseif(contains(gpuName,'K80'))
            ptxPath = '/efs-mount-point/kernels/sm37/cuda_kernel_20_universal.ptx';

        else
            ptxPath = '/efs-mount-point/kernels/cuda_kernel_20_universal.ptx';
        end    
    end
    
else
    
   ptxPath  = '/efs-mount-point-MATLAB/EduTool-Jorge/EduTool-CUDA/cuda_kernel_25.ptx';

end
    


