function appSpecs = initializeApplication(APP,APPROACH,MODE,VERSION,...
    jsontestfile,testPauseTime)

appSpecs.application    = APP;
appSpecs.approach       = APPROACH;
appSpecs.mode           = MODE;
appSpecs.version        = VERSION;

%%
% use inputs 2 and 3 of main function only if edutool-jsonstandalone-test
if strcmp(APP,'edutool') && strcmp(APPROACH,'jsonstandalone')
    if strcmp(MODE,'test')
        if nargin < 2
            jsontestfile    = '';
            testPauseTime   = '3';
        elseif nargin == 2        
            testPauseTime   = '3';
        end
    else
        jsontestfile    = 'none';
        testPauseTime   = '0';
    end
    
%     % if the code is not deployed, execute the script that copies files
%     % from local folder (/home/ubuntu/edutoolTransferToS3) to s3
%     if ~isdeployed
%         tcmdS3   = tic();
%         % If in dev mode and the folder edutoolTransferToS3 not available,
%         % create thus foldr
%         if strcmp(MODE,'dev') && exist('/home/ubuntu/edutoolTransferToS3','dir')~=7
%             [~,cmdout1] = system('sudo mkdir /home/ubuntu/edutoolTransferToS3');  %#ok
%         end
%         disp('Calling the .sh for establishing connection with s3')
%         [~,cmdout2] = system('sudo ./+tools/s3bucket/s3TransferFilesService.sh');
%         disp(cmdout2)
%         [~,cmdout3] = system('sudo ps -ef | grep s3-upload');
%         disp(cmdout3)
%         tcmdS3Time = toc(tcmdS3);
%         disp(['Elapsed time: ',num2str(tcmdS3Time),'sec'])
%     end
    
end

appSpecs.jsontestfile   = jsontestfile;
appSpecs.testPauseTime  = testPauseTime;

%%
if strcmp(APP,'edutool') && strcmp(APPROACH,'centralized')
    if nargin<4
        error(['The centralized solution requires specs for the ',...
            'simulator to be defined.'])
    end
    simulatorInputs	= strsplit(simulatorSpecs,'-');
    userID          = str2num(simulatorInputs{1,1});
    if size(simulatorInputs,2)==1        
        courseID        = 6;
        anatomicalID    = 6;
    elseif size(simulatorInputs,2)==2
        courseID        = str2num(simulatorInputs{1,2});
        anatomicalID    = 6;
    else        
        courseID        = str2num(simulatorInputs{1,2});
        anatomicalID    = str2num(simulatorInputs{1,3});
    end
    
    appSpecs.userID         = userID;
    appSpecs.courseID       = courseID;
    appSpecs.anatomicalID   = anatomicalID;
    
end