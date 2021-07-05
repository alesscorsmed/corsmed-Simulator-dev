function [jsonStructure,R,redisReportGet] = redisGetJsonWrapper(R,keyInput)

redisAddress        = R.RemoteHost;
redisPort           = num2str(R.RemotePort);
inputBufferSize     = R.inputBufferSize;
outputBufferSize    = R.outputBufferSize;
    
redisReportGet = '';
while ~strcmp(redisReportGet,'OK')

    if strcmp(R.Status,'closed')
        R = tools.redis.redisEstablishConnection(redisAddress,redisPort,...
            inputBufferSize,outputBufferSize);
    end

    [jsonRedis,~,redisReportGet] = tools.redis.redisGet(R,keyInput);
    
    % If there is no key, or the key is empty, continue
    if strcmp(redisReportGet,'ERROR - NONEXISTANT KEY')
        redisReportGet = 'ERROR - NONEXISTANT KEY';
        jsonStructure = '';
        return;
    elseif iscell(jsonRedis)
        if isempty(jsonRedis{1,1})
            redisReportGet = 'EMPTY KEY';
            jsonStructure = '';
            return;
        end        
    end

    if contains(redisReportGet,'WrongBufferSize')
        redisReportWords = regexp(redisReportGet,'-','split');
        if strcmp(redisReportWords{1,2},'GET')
            inputBufferSize = str2num(redisReportWords{1,3});
            disp('Update input buffer size')
        elseif strcmp(redisReportWords{1,2},'SET')
            outputBufferSize = str2num(redisReportWords{1,3});
            disp('Update output buffer size')
        end
        
        flushinput(R)
        flushoutput(R)
        tools.redis.redisDisconnect(R);
        disp('Redis connection stopped due to improper buffersize.')
        disp('Restarting redis connection...')

    end
    
end

% Find the last occurence of the character } in the response from redis
occChar = strfind(jsonRedis{1,1},'}');

jsonStructure = jsondecode(jsonRedis{1,1}(1:occChar(1,end)));