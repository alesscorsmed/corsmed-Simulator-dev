function [instanceAtts] = initialize()
%
% DATA.INSTANCEATTS.INITIALIZE
%
%	Function that initializes elements of the instance
%
%   This function is useful to define the fields.
%
% INPUT
%   None
%
% OUTPUT
%   Instance Attributes Structure
%
%========================  CORSMED AB © 2020 ==============================
%



%All info goes to instanceAtts struct
instanceAtts.versionNum = 'v20200930a';
instanceAtts.instanceId = '';
instanceAtts.courseId           = 0;
instanceAtts.AWStagUserID       = 0;

instanceAtts.cudaVersion     = 8;
instanceAtts.Parfeval    = 0;
instanceAtts.pythonVersion   = '3.7';
instanceAtts.developmentUse = 0;
instanceAtts.NSIM    = 0;

if instanceAtts.NSIM
    setup_python_libs;
end



