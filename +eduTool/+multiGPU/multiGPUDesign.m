function [spmdTag,parfevalTag] = multiGPUDesign(acquisition,parfevalTag,expControl,tags_struct)
%
% BACKEND.MULTIGPUDESIGN
%
%   select multiGPU design based on experiment	
%
% INPUT
%
%   acquisition          struct of acquired experiment Data
%   parfevalTag          boolean value for parallel design
%   expControl           struct of various experiment Data
%   tagsStruct           struct with necessary data from AWS tags
%
% OUTPUT
%
%   spmdTag        boolean value for spmd design
%   parfevalTag    boolean value for parallel design
%
%========================  CORSMED AB © 2020 ==============================
%


%% serial approach if we have light experiment
if(parfevalTag ==1)
    
    if(acquisition.data.pulseSeqFamilyName ~= "GRE-3D" && expControl.simulation.simulationKernel =="latest")
       
       parfevalTag = 0;
    
       %% cases the experiment is considered heavy or light  
       if(acquisition.data.numSlices >= 4)
           
           %% small matrix and default isotropic grid 
           if(acquisition.data.matrixX >= 128 && acquisition.data.matrixY >= 128 ...
           && expControl.model.gridSizeOption ~= "1mm")            
                parfevalTag =1;
           end

           %% if even one matrix is maxed out
           if(acquisition.data.matrixX >= 384 || acquisition.data.matrixY >= 384)
               parfevalTag = 1;
           end
           
           %% Big Matrix and high slice Thickness
           if(acquisition.data.matrixX >= 192 && acquisition.data.matrixY >= 192 && acquisition.data.sliceThickness >= 0.015)
               parfevalTag = 1;
           end

           %% Big Matrix and heavy coil types
           if acquisition.data.matrixX >= 192 && acquisition.data.matrixY >= 192
               switch expControl.model.coilType

                   case '4xConformalRx'
                       parfevalTag=1;
                   case '8xConformalRx'
                       parfevalTag=1;
                   case '8xPlanarRx'
                       parfevalTag=1;
               end
           end
           
           %% in case of Rotational Motion
           if(strcmp(expControl.motionSpecs.pattern,'rotational'))
               parfevalTag=1;
           end
           
       end 
    end
end

%Choose SPMD way for 3D
spmdTag = 0;
if(gpuDeviceCount>1 && parfevalTag==1 && acquisition.data.pulseSeqFamilyName == "GRE-3D")
        spmdTag=1;
end