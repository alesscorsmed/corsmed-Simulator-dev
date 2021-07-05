function parpoolConn = parpoolConnection(conn_localdb,localdb)

%
% BACKEND.PARPOOLCONNECTION
%
%	builds the  parpool connection 
%
% INPUT
%
%   tagsStruct           struct with necessary data from AWS tags
%   instanceAtts         attributes from the AWs instance
%
% OUTPUT
%
%   parpoolConn          parallel pool connection object
%========================  CORSMED AB © 2020 ==============================
%



parpoolDBtimer = tic;
fprintf('Opening JDBC data source connection for parpool ...\n');
%Create ParPool DB Connection
urlsplit    = split(conn_localdb.URL,'/');
opts_url    = split(urlsplit{3},':');
opts_port   = opts_url{2};
opts_server = opts_url{1};

opts                    = configureJDBCDataSource('Vendor','MySQL');
opts.DataSourceName     = 'ParpoolConnection';
opts.DatabaseName       = conn_localdb.DataSource;
opts.PortNumber         = str2num(opts_port);
opts.JDBCDriverLocation = '/efs-mount-point/IMPORTANT/mysql-connector-java-8.0.12.jar';
opts.Server             = opts_server;
opts                    = addConnectionOptions(opts,'useSSL','YES');
opts                    = addConnectionOptions(opts,'verifyServerCertificate','false');

% Test connection
[status,message]        = testConnection(opts,localdb.user,localdb.pass);  %#ok

p = gcp();
saveAsJDBCDataSource(opts)
parpoolConn = createConnectionForPool(p,opts.DataSourceName,localdb.user,localdb.pass);
fprintf(' done -- %.3f sec\n',toc(parpoolDBtimer));