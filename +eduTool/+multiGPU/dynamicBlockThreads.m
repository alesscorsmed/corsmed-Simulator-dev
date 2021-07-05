function [ptxPath,threadsPerBlock,blocksPerIteration] = ...
    dynamicBlockThreads(threadsPerBlock,blocksPerIteration,...
    isochromats,currentCPU,pdInhomogeneity,motionPattern,...
    simulationKernel,runParallel,workers)
%
%DYNAMICRESOURCES.dynamicThreadsBlocks
%
%   Select Threads / Blocks based on
%   GPU and Kernel
%
%
% INPUT
%   threadsPerBlock     the default threads extracted from Global Configs
%   blockerPerIteration the default blocks extracted from Global Configs
%   isochromats         the isochromats of corresponding slice
%   currentGPU          the  designated GPU used
%   pdInhomogeneity     boolean value for the use of PDinhomogeneity
%   motionPattern       choice of Motion
%   simulationKernel    the kernel used
%   runParallel         boolean value for use of parallel run
%   workers             number of workers corresponnding GPUs available
%
% OUTPUT
%   finalThreads        optimal number of threads
%   finalBlocks         optimal number of blocks
% 
% ========================= CORSMED AB @ 2020 ============================
%

%case we use the stable kernel
if strcmp(simulationKernel,'stable')
    
    %case we use the K80 GPU (p2) or V100 (p3)
    if(currentCPU.Name == "Tesla K80")
        threadsPerBlock = 128;
        blocksPerIteration = 128;
    else    
        num_of_iterations = ceil(isochromats / (threadsPerBlock * ...
            blocksPerIteration));    
        if(num_of_iterations > 6)
            threadsPerBlock = 512;
            blocksPerIteration = 256;
        else
            threadsPerBlock = 256;
            blocksPerIteration = 256; 
        end
    end
    
else
    
    if(currentCPU.Name == "Tesla K80")
        threadsPerBlock = 128;
    else
        threadsPerBlock = 512;
    end
    
    %case we have parallel run
    if runParallel
        blocksPerIteration = ceil(isochromats/(threadsPerBlock*workers));
    else
        blocksPerIteration = ceil(isochromats/threadsPerBlock);
    end
end

ptxPath = kernel.find_ptxFile(pdInhomogeneity,motionPattern,...
        blocksPerIteration,threadsPerBlock,currentCPU.Name,...
        simulationKernel);

