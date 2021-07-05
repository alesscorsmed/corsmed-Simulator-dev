function R = redisSetJsonWrapper(R,keyOutput,jsonStructure,applyBracketsFix)
% applyBracketsFix: if 1, then brackets are added at the start and end of 
% the json string (if they are missing). This is mainly for storing the
% redis keys in the correct format for the UI to be able to run properly.

if nargin<4
    applyBracketsFix = 0;
end

redisAddress        = R.RemoteHost;
redisPort           = num2str(R.RemotePort);
inputBufferSize     = R.inputBufferSize;
outputBufferSize    = R.outputBufferSize;

redisReportSet = '';
while ~strcmp(redisReportSet,'OK')

    if strcmp(R.Status,'closed')
        R = tools.redis.redisEstablishConnection(redisAddress,redisPort,...
            inputBufferSize,outputBufferSize);
    end
    
    keyValue = jsonencode(jsonStructure);
    
    % Workaround for passing json strings to redis with brackets at the start
    % and end of the json string (if they are missing)
    if applyBracketsFix
        if ~strcmp(keyValue(1,1),'[') && ~strcmp(keyValue(1,1),']')
            keyValue = strcat('[',keyValue,']');
        end
    end

    [~, redisReportSet] = tools.redis.redisSetJson(R,keyOutput,...
                keyValue);

    if contains(redisReportSet,'WrongBufferSize')
        redisReportWords = regexp(redisReportSet,'-','split');
        if strcmp(redisReportWords{1,2},'GET')
            inputBufferSize = str2num(redisReportWords{1,3});
            disp('Update input buffer size')
        elseif strcmp(redisReportWords{1,2},'SET')
            outputBufferSize = str2num(redisReportWords{1,3});
            disp('Update output buffer size')
        end

        tools.redis.redisDisconnect(R);
        disp('Redis connection stopped due to improper buffersize.')
        disp('Restarting redis connection...')

    end
    
end