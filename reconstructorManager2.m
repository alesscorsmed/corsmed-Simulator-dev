function reconstructorManager2(redisAddress,redisPort,...
    redisRequestKey,redisResponseKey,redisUpdatesKey,uniqueID,...
    dummyData)
%
% SERVICES.RECONSTRUCTORMANAGER
%
%	Mock up of the SIMULATOR manager service
%
% INPUT
%
% OUTPUT
%
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'reconstructorManager';
if (nargin < 1)
    ME = MException('eduTool:wrongArgCount',...
        '%s : wrong argument count',functionName);
    throw(ME);
end

if (nargin < 7)
    dummyData = 0;
end

useMatlabBuiltinRedisFunction = 0;

try

    %% Establish connection with Redis
    % Connect to a redis server
    redis.address   = redisAddress;
    redis.port      = redisPort;

    R = tools.redis.redisEstablishConnection(redis.address,redis.port);
    
    %%
    % Check if the key exists and it is not empty
    [redisRequest,~,redisReport] = tools.redis.redisGet(R,redisRequestKey);

    % If there is no key, or the key is empty, continue
    if isempty(redisRequest) || strcmp(redisReport,'ERROR - NONEXISTANT KEY')
        ME = MException('eduTool:emptySimReqRedisKey',...
            '%s : this key is either empty or does not exist',redisRequestKey);
        throw(ME);
    elseif iscell(redisRequest)
        if isempty(redisRequest{1,1})
            ME = MException('eduTool:emptySimReqRedisKey',...
            '%s : this key is either empty',redisRequestKey);
        throw(ME);
        end
    end
    
    tTotal = tic();
    %% GET DATA FROM REDIS
    fprintf(1, '\n Getting data from redis for %s experiment',...
        uniqueID);
    fprintf(1, '\n');

    totalJobsKey        = strcat(uniqueID,'_totalJobs');
    reconDataKey        = strcat(uniqueID,'_reconData');
    expControlKey       = strcat(uniqueID,'_expControl');
    imageDataKey        = strcat(uniqueID,'_imageData');
    
    if useMatlabBuiltinRedisFunction
    
        % Establish redis connection using the MATLAB built-in function
        c = tools.redisMatlab.redisEstablishConnectionMatlab(...
            redis.address,redis.port);
        
%         v = get(c,{pulseSequenceKey,motionModelKey,...
%             expControlKey,simJobKey});
% 
%         pulseSequence   = v{1,1};
%         motionModel     = v{1,2};
%         expControl      = v{1,3};
%         simJob          = v{1,4};
    
    else
        % Flush data from input and output buffer of R
        flushinput(R)
        flushoutput(R)
    
        [expControlValue,~,~] = tools.redis.redisGet(R,expControlKey);
        expControl = jsondecode(expControlValue{1,1});
        
        % @@@ TEMP
        expControl = rmfield(expControl,'connLocalDB');        

        if dummyData
            load('20210125_reconInputs.mat')
        else
            [reconDataValue,~,~] = tools.redis.redisGet(R,reconDataKey);
            reconData = jsondecode(reconDataValue{1,1});
        
            [imageDataValue,~,~] = tools.redis.redisGet(R,imageDataKey);
            imageData = jsondecode(imageDataValue{1,1});
        end
        
        % Add the operators in the reconData structure 
        % create the Fourier Transform operators
        operators.iFT   = @(kSpace) fftshift(ifftn(ifftshift(kSpace)));
        operators.fFT   = @(iSpace) fftshift( fftn( fftshift(iSpace)));
        operators.iFTZ  = @(kSpace) fftshift(ifftn(ifftshift(kSpace),[],3));
        operators.fFTZ  = @(iSpace) fftshift( fftn( fftshift(iSpace),[],3));
        reconData.encoding.operators = operators;

        [totalJobsValue,~,~] = tools.redis.redisGet(R,totalJobsKey);
        totalJobs = str2num(totalJobsValue{1,1});
        
        % Flush data from input and output buffer of R
        flushinput(R)
        flushoutput(R)

        for i = 1:totalJobs    
            simSignalKey = strcat(uniqueID,'_simSignal_',num2str(i));

            [simSignalValue,~,~] = tools.redis.redisGet(R,simSignalKey);
            simSignalPerJob{1,i} = jsondecode(simSignalValue{1,1});
        end    
    end
    
    simSignal.numJobs = totalJobs;
    simSignal.numSlices = simSignalPerJob{1,1}.numSlices;
    simSignal.numCoils = simSignalPerJob{1,1}.numCoils;
    simSignal.numReads = simSignalPerJob{1,1}.numReads;
    
    for i = 1:totalJobs
        simSignal.timeSolution{i}=simSignalPerJob{1,i}.timeSolution;
    end
    
    %% RECONSTRUCTION
    [reconData] = eduTool.run.recon(...
        simSignal, reconData, expControl );
    
    %% DICOM GENERATION
    [imageData,jsonStr] = eduTool.run.dicom(...
        reconData, imageData, expControl );
    
    %% EXPORT OUTPUT TO REDIS
    [~, statusSetResponse] = tools.redis.redisSetJson(R,redisResponseKey,jsonStr);

    %% Delete redisRequestKey from redis
    [~, statusDel] = tools.redis.redisDEL(R,redisRequestKey);
    if ~strcmp(statusDel,'OK')
        ME = MException('eduTool:delKeyRedisFailure',...
            '%s : could not be deleted from redis',redisRequestKey);
        throw(ME);
    end

    %% REPORT OK
    fprintf(1, '\n');
    fprintf(1, '\n RECONSTRUCTOR MANAGER on %s DONE', redisRequestKey);
    fprintf(1, '\n  Elapsed time  : %.3fs', toc(tTotal));
    fprintf(1, '\n');
    
    % Flush data from input and output buffer of R
    flushinput(R)
    flushoutput(R)
    
catch ME
    %% send error back to UI
    errorMessage = sprintf(['Error in function %s() at line %d.',...
        '\n Error Message: %s'], ....
        ME.stack(1).name,ME.stack(1).line,...
        ME.message);
    % errorMessage = tools.printErrorMessage(expControl,ME); 

    tools.updateJsonProgress(R,redisUpdatesKey,...
        'error',errorMessage);
    
    % Flush data from input and output buffer of R
    flushinput(R)
    flushoutput(R)

    tools.redis.redisDisconnect(R);

    disp(errorMessage)
    
    
end
