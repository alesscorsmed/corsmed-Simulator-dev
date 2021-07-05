function [Output, R, S] = redisCommand(R, Cmd)

%% Check if the method is a GET or SET
wordsOfCommand = regexp(Cmd, '\s', 'split');
if strcmp(wordsOfCommand{1,3},'GET')
    cmdType = 'GET';
elseif strcmp(wordsOfCommand{1,3},'SET')
    cmdType = 'SET';
    if (size(Cmd,2)+100)>R.OutputBufferSize
        Output = strcat('WrongBufferSize-',cmdType,'-',num2str(size(Cmd,2)+100));
        S = 'ERROR';
        return
    end
else
    cmdType = '';
end    

%% Execute command
fprintf(R, Cmd);

%% wait for bytes to show up
timeout = 1.0;
tic
while R.BytesAvailable == 0
  pause(0.005)
  if toc >= timeout
    Output = '';
    S = 'ERROR - REDIS TIMEOUT';
    return;
  end
end

%% Read response
response = '';
tic
while R.BytesAvailable > 0
  [chunk, ~, msg] = fread(R, R.BytesAvailable);
  response = [response char(chunk')];
  
  % if the response is not +OK, make the following tests
  responceFirstCell = response;
  if ~strcmp(responceFirstCell(1:3),'+OK') && ~strcmp(responceFirstCell(1:5),'+PONG')
      wordsOfResponse   = regexp(response, '\s', 'split');  
      sizeOfResponse    = str2num(wordsOfResponse{1,1}(2:end));
      if isnumeric(sizeOfResponse)
          if (strcmp(cmdType,'GET') && sizeOfResponse>R.InputBufferSize) || ...
                  (strcmp(cmdType,'SET') && sizeOfResponse>R.OutputBufferSize)

              response = strcat('WrongBufferSize-',cmdType,'-',num2str(sizeOfResponse+100));
              break

          end
      end
%       if (strcmp(cmdType,'GET'))
%           response = strcat(wordsOfResponse{3:end});
%       end
  end
  if toc >= timeout
    break
  end
end

Output = response;
S = 'OK';