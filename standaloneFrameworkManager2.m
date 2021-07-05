redis.address   = 'cache.internal.integration.corsmed.com';
redis.port      = '6379';

% Prepare the redis input that holds the json string with all experiment'
% specs
uniqueID = 'CX14';
if 1
    R = tools.redis.redisEstablishConnection(redis.address,redis.port);
    inputs.test1 = 1;
    jsonStr=jsonencode(inputs);
%     [~, statusSet1] = tools.redis.redisSet(R,'testRequest1',jsonStr);
    [~, statusSet2] = tools.redis.redisSet(R,'simreq_1',jsonStr);    
    [~, statusSet3] = tools.redis.redisSet(R,'recreq_1',jsonStr);
    [~, statusSet4] = tools.redis.redisSet(R,'testCalc1',uniqueID);
    tools.redis.redisDisconnect(R);
end

%% CALCULATOR
tTotalSlicer = tic;
% breakLoop = 1, if we want to terminate the while loop after the 1st 
% successful iteration 
breakLoop = 1;
loadDummyData = 1;
slicerManager3(redis.address,redis.port,...
    'calcreq_13','testResponse1','testUpdate1','testCalc1',...
    '[4,6]','/efs-mount-point/S20/INPUTS/anatomical/edutool',...
    '/efs-mount-point/S20/INPUTS/coils/edutool',breakLoop,loadDummyData);
totalSlicer = toc(tTotalSlicer);

%% SIMULATOR
tTotalSimulator = tic;
% dummyData = 1, if we want to skip loading all data from Redis (due to
% instability issues) and load dummy data from the local file system
dummyData = 0;
simulatorManager3(redis.address,redis.port,...
    'simreq_1','simresp_1','testUpdate1',uniqueID,'1',dummyData);
totalSimulator = toc(tTotalSimulator);

%% CHECK REDIS AVAILABLE KEYS
R = tools.redis.redisEstablishConnection(redis.address,redis.port);
[Response, R, S] = tools.redis.redisCommand(R,...
    tools.redis.redisCommandString(sprintf('keys *')));
tools.redis.redisDisconnect(R);

ResponseCell    = regexp(Response, '\s', 'split');
redisKeys       = sort(ResponseCell(1,5:4:end));

fprintf(1, '\n  REDIS KEYS : \n');
for i=1:size(redisKeys,2)
    fprintf(1,'%s\n',redisKeys{1,i});
end

%% RECONSTRUCTOR
tTotalReconstructor = tic;
% simDummyData = 1, if we want to skip loading all data from Redis (due to
% instability issues) and load dummy data from the local file system
simDummyData = 1;
reconstructorManager2(redis.address,redis.port,...
    'recreq_1','recresp_1','testUpdate1',uniqueID,...
    simDummyData)
totalReconstructor = toc(tTotalReconstructor);

%%
fprintf(1, '\n  SLICE MANAGER TOTAL TIME : %.3fs\n', totalSlicer);
fprintf(1, '\n  SIMULATOR MANAGER TOTAL TIME : %.3fs\n', totalSimulator);
fprintf(1, '\n  RECONSTRUCTOR MANAGER TOTAL TIME : %.3fs\n', totalReconstructor);