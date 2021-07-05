function [r, status] = redisConnection(host,port,...
    inputBufferSize,outputBufferSize)

if nargin == 0
  host = '127.0.0.1';
end

if nargin < 2
  port = 6379;
end

if nargin < 3
    inputBufferSize     = 10000000; % 10MB
    outputBufferSize    = 10000000; % 10MB 
end

r = tcpip(host, port);
r.InputBufferSize   = inputBufferSize;
r.OutputBufferSize  = outputBufferSize;
r.terminator = 'CR/LF';
fopen(r);

if ~strcmp(r.Status, 'open')
  error('Unable to connect to Redis server at host %s on port %d', host, port)
  r = 0;
  status = 'ERROR';
  return
end

% there shouldn't be any info on the port, but if there is, let's get rid of it
if r.BytesAvailable > 0
  fread(r, r.BytesAvailable);
end

status = 'OK';
