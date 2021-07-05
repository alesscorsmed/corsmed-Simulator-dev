%% Load anatomical model
anatomicalModel = ...
    load('/efs-mount-point/S20/INPUTS/anatomical/edutool/MIDA_20210110.mat');

disp('Anatomical model loaded.')

%% Connect to redis
try
    ctrl = mps.cache.control('myRedisConnection','Redis',...
    'Host','redis-nodeport.internal.integration.corsmed.com',...
    'Port',30007);
    start(ctrl)
catch ME
    
end
c = mps.cache.connect('myCache','Connection','myRedisConnection');

disp('Redis connection established.')

%% Prepare matrices and upload to redis
spatialModel = anatomicalModel.spatial;

batchSize = 1000000;
noOfFiles = ceil(size(spatialModel,1)/batchSize);

for i=1:noOfFiles
    startIndex = (i-1)*batchSize+1;
    if i == noOfFiles
        endIndex = size(spatialModel,1);
    else
        endIndex = i*batchSize;
    end
    
    spatial = spatialModel(startIndex:endIndex,:);
    put(c,['anatomicalMIDAspatial_',num2str(i)],spatial)
    clear spatial
    disp(['Part ',num2str(i),' uploaded successfully'])
end

%% Print the redis keys
keys(c)