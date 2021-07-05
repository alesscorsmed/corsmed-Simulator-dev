function [currentGPU] = chooseGPU(parfevalTag,spmdTag,sliceNum)
%
% EDUTOOL.AWS.EC2.READINSTANCETAGS
%
%	chooses the GPU that will be used in this Slice/Experiment
%   based on parallel run and the multiGPU Design
%
% INPUT
%
%   parfevalTag       The boolean value for parallel run (if 0 then serial)
%   spmdTag           The boolean value for SPMD way ( if 0 then GPU/slice)
%   sliceNum          The slice number of the experiment
%
% OUTPUT
%   currentGPU        The chosen GPU of the aws instance
%
%========================  CORSMED AB © 2020 ==============================
%
if(nargin == 1)    
    sliceNum = parfevalTag;
    if(gpuDeviceCount > 1)

            %in case slices > GPUs           
            if(sliceNum > gpuDeviceCount)
                if(rem(sliceNum,gpuDeviceCount) ~=0)
                    currentGPU = rem(sliceNum,gpuDeviceCount);
                else
                    currentGPU = gpuDeviceCount;
                end
            else
                currentGPU = sliceNum;
            end
    else
            currentGPU = 1;
    end
                   
else
           
    %Choose GPU 
    if(parfevalTag)
        if(spmdTag)

            %if we have a spmd design (mostly 3D)
            currentGPU = gpuDevice(labindex);
        else
            if(gpuDeviceCount > 1)

                %in case slices > GPUs           
                if(sliceNum > gpuDeviceCount)
                    if(rem(sliceNum,gpuDeviceCount) ~=0)
                        currentGPU = gpuDevice(rem(sliceNum,gpuDeviceCount));
                    else
                        currentGPU = gpuDevice(gpuDeviceCount);
                    end
                else
                    currentGPU = gpuDevice(sliceNum);
                end
            else
            currentGPU = gpuDevice(labindex);
            end
        end
    else
        currentGPU = gpuDevice(1);
    end

end
    