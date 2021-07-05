function [tissueValues] = perfusionUpdateTissueProp( tissueValues,...
    perfusionData, contrNum, dbgControl)
%
% MODELS.PERFUSION.PERFUSIONUPDATETISSUEPROP
%
%	Function that updates the tissue properties based on the perfusion
%	model selected for the experiment
%
% INPUTS
%   tissueValues    original tissue values
%   perfusionData   perfusion data, see below details
%
% OUTPUT
%   tissueValues    modified tissue values
%
%
% INPUT Details
%   perfusionData   struct with parameters:
%     perfusionData.GdRelaxivity = 5.6 / 1000;   % Gd relaxivity         [l/(mmol*ms)]
%     perfusionData.MBFrest      = 1.0 / 60;     % Rest MBF              [ml/g/s]
%     perfusionData.MBFstress    = 3.5 / 60;     % Stress MBF            [ml/g/s]
%     perfusionData.Fermi_alpha  = 0.25;         % Fermi model parameter alpha
%     perfusionData.Fermi_beta   = 0.25;         % Fermi model parameter beta
%     perfusionData.Tshift       = 3;            % temporal LV-myo shift [s]
%     perfusionData.contrastDose = 0.075;        % reference dose [mmol/kg]
%     perfusionData.RestStress   = 2;            % 1=rest; 2=stress
%     % ids of tissues subject to perfusion
%     perfusionData.idMyoLV      = [];
%     perfusionData.idMyoLA      = [];
%     perfusionData.idMyoRV      = [];
%     perfusionData.idMyoRA      = [];
%     perfusionData.idBloodLV    = [];
%     perfusionData.idBloodLA    = [];
%     perfusionData.idBloodRV    = [];
%     perfusionData.idBloodRA    = [];
%     %  ids of tissues to keep intact (non zero)
%     perfusionData.idSteady     = [];
%     % timing
%     perfusionData.contrasts    = 32;  % Number of Dynamics
%     perfusionData.Tcc          = 1.0; % Cardiac cycle duration [s]
%
% INFO
% - Volume of Gd = D x W / C, where D = dose by weight (mmol/kg), W =
% weight (kg) and C = concentration (mmol/ml)
%
%========================  CORSMED AB © 2020 ==============================
%
functionName = 'simModel:perfusion:perfusionUpdateTissueProp';
% check args
if (nargin < 3)
    ME = MException(['error:',functionName],...
        'Wrong arguments.');
    throw(ME);
end
% optional arguments
if (nargin < 4) || isempty(dbgControl)
    dbgControl.mode = 0;
    dbgControl.file = [];
end
% info for debugging
if dbgControl.mode
    try % open file if possible, otherwise dump to stdout
        fid = fopen(dbgControl.file,'a');
    catch
        fid = 1;
    end
    tTotal = tic();
    fprintf(fid, '\n%s : start', functionName);
end

%% Convert relaxation times to ms
tissueValues(:,1:2) = tissueValues(:,1:2)*1000;

%% apply perfusion
try
    % T1Muscle        = modelTissues(muscleID,1); % T1 muscle [ms]
    % T1Fat           = modelTissues(fatID,1);    % T1 fat    [ms]
    % T1Blood         = modelTissues(bloodID,1);  % T1 blood  [ms]
    % T1Liver         = modelTissues(liverID,1);  % T1 liver  [ms]
    % T1Bone          = modelTissues(boneID,1);   % T1 bone   [ms]
    
    % extra T1 properties of tissues subject to perfusion
    T1myoLV         = tissueValues(perfusionData.idMyoLV,1);
    T1myoLA         = tissueValues(perfusionData.idMyoLA,1);
    T1myoRV         = tissueValues(perfusionData.idMyoRV,1);
    T1myoRA         = tissueValues(perfusionData.idMyoRA,1);
    T1bloodLV       = tissueValues(perfusionData.idBloodLV,1);
    T1bloodLA       = tissueValues(perfusionData.idBloodLA,1);
    T1bloodRV       = tissueValues(perfusionData.idBloodRV,1);
    T1bloodRA       = tissueValues(perfusionData.idBloodRA,1);
    
    GdRelaxivity    = perfusionData.GdRelaxivity;	% Gd relaxivity         [l/(mmol*ms)]
    MBFrest         = perfusionData.MBFrest;      % Rest MBF              [ml/g/s]
    MBFstress       = perfusionData.MBFstress;    % Stress MBF            [ml/g/s]
    Fermi_alpha     = perfusionData.Fermi_alpha;  % Fermi model parameter alpha
    Fermi_beta      = perfusionData.Fermi_beta ;  % Fermi model parameter beta
    Tshift          = perfusionData.Tshift;       % temporal LV-myo shift [s]
    RestStress      = perfusionData.RestStress;   % 1=rest; 2=stress
    contrastDose    = perfusionData.contrastDose; % reference dose [mmol/kg]
    
    contrasts       = perfusionData.contrasts;             % Number of Dynamics
    Tcc             = perfusionData.Tcc;                   % Cardiac cycle duration [s]
    
    % aif template
    aif = [0,0,0,0,0.00088825,0.017793,0.13429,0.5466,1.453,2.8276,4.3365,...
        5.5112,6.0165,5.7934,5.0205,3.9772,2.9161,1.9988,1.2913,0.79165,...
        0.46315,0.25983,0.14035,0.073254,0.037056,0.018216,0.0087227,...
        0.0040768,0.0018632,0.00083405,0.00036621,0.00015793,6.6972e-05,...
        2.7956e-05,1.1499e-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    
    % ca is the arterial input at dose specified in contrastDose [mmol/l]
    aif01   = aif*12/max(aif);
    ca      = aif01*contrastDose/0.1;   % scale to desired dose
    t       = linspace(0,(contrasts-1)*Tcc,length(ca));
    
    % Upsample AIF (by a factor of 100)
    tinf    = 0:1/100:t(end);
    cainf   = interp1(t,ca,tinf,'pchip','extrap');
    
    % Scale flow
    qfl     = MBFrest;          % rest flow
    if (RestStress==2) % stress flow
        qfl = MBFstress;
    end
    qfl     = qfl*length(ca)/round(contrasts*Tcc);
    
    % Calculate impulse residue function
    irf	= (1+Fermi_beta)./(1+Fermi_beta*exp(Fermi_alpha*tinf));
    
    % Calculate cm (myocardial tissue concentration [mmol/l]), apply Tshift and
    % downsample
    cm = convolve(tinf,cainf,irf,qfl);                  % myocardial conc
    cm = interp1(tinf+Tshift,cm,tinf,'pchip','extrap'); % apply Tshift
    ca = interp1(tinf,cainf,t);                         % downsample AIF (discrete measurement)
    cm = interp1(tinf,cm,t);                            % downsample MYO
    
    %% Contrast concentration for different tissues
    % Get contrast concentration in right atrium (ra), right ventricle (rv),
    % left atrium (la) and left ventricle (lv) at specific time t (times of
    % contrasts)
    idynamic    = 1:contrasts;
    ra          = min(round((idynamic+4)*length(ca)/contrasts)+1,length(ca));
    rv          = min(round((idynamic+3)*length(ca)/contrasts)+1,length(ca));
    la          = min(round((idynamic+0)*length(ca)/contrasts)+1,length(ca));
    lv          = min(round((idynamic-1)*length(ca)/contrasts)+1,length(ca));
    
    contrastBP_ra   = ca(ra(contrNum));   % RA blood pool
    contrastBP_rv   = ca(rv(contrNum)); 	% RV blood pool
    contrastBP_la   = ca(la(contrNum));  	% LA blood pool
    contrastBP_lv   = ca(lv(contrNum));  	% LV blood pool
    
    contrastMYO_ra  = cm(ra(contrNum)); 	% RA myocardium
    contrastMYO_rv  = cm(rv(contrNum)); 	% RV myocardium
    contrastMYO_la  = cm(la(contrNum)); 	% LA myocardium
    contrastMYO_lv  = cm(lv(contrNum));  	% LV myocardium
    
    
    %% Update T1 values
    % Myocardium
    r1_myoLV  = 1/T1myoLV+contrastMYO_lv*GdRelaxivity;
    r1_myoRV  = 1/T1myoRV+contrastMYO_rv*GdRelaxivity;
    r1_myoLA  = 1/T1myoLA+contrastMYO_la*GdRelaxivity;
    r1_myoRA  = 1/T1myoRA+contrastMYO_ra*GdRelaxivity;
    
    % Blood
    r1_bldplLV  = 1/T1bloodLV+contrastBP_lv*GdRelaxivity;
    r1_bldplRV  = 1/T1bloodRV+contrastBP_rv*GdRelaxivity;
    r1_bldplLA  = 1/T1bloodLA+contrastBP_la*GdRelaxivity;
    r1_bldplRA  = 1/T1bloodRA+contrastBP_ra*GdRelaxivity;
    
    % re-assign
    tissueValues(perfusionData.idMyoLV,1)    = 1/r1_myoLV;
    tissueValues(perfusionData.idMyoLA,1)    = 1/r1_myoLA;
    tissueValues(perfusionData.idMyoRV,1)    = 1/r1_myoRV;
    tissueValues(perfusionData.idMyoRA,1)    = 1/r1_myoRA;
    tissueValues(perfusionData.idBloodLV,1)  = 1/r1_bldplLV;
    tissueValues(perfusionData.idBloodLA,1)  = 1/r1_bldplLA;
    tissueValues(perfusionData.idBloodRV,1)  = 1/r1_bldplRV;
    tissueValues(perfusionData.idBloodRA,1)  = 1/r1_bldplRA;
    
    %     % Body
    %     r1_fat_body	= 1/T1Fat;    % that covers body and pericardium
    %     r1_muscle   = 1/T1Muscle;   % that covers muscle
    %     r1_liver    = 1/T1Liver;    % that covers liver
    %     r1_bone     = 1/T1Bone;     % that covers bones (ribs, cortical bones, spine, bone marrow, etc)
    
    % Keep fat, Muscle, Liver and Bone intact and set all other tissues PD equal to zero
    idToZero = setdiff(1:size(tissueValues,1),...
        union(perfusionData.idSteady(:), ...
        [perfusionData.idMyoLV; ...
        perfusionData.idMyoLA; ...
        perfusionData.idMyoRV; ...
        perfusionData.idMyoRA; ...
        perfusionData.idBloodLV; ...
        perfusionData.idBloodLA; ...
        perfusionData.idBloodRV; ...
        perfusionData.idBloodRA ] ) );
    
    tissueValues(idToZero,3) = 0.0;
    tissueValues(perfusionData.idSteady(:),3) = 0.2;
    
    if  dbgControl.mode
        fprintf(fid, '\n%s : perfusion applied', functionName);
    end
    
catch
    if  dbgControl.mode
        fprintf(fid, '\nWARNING  %s : perfusion could not be applied', functionName);
    end
end

%% report
if  dbgControl.mode
    tTotal = toc(tTotal);
    fprintf(fid, '\n%s : perfusion applied', functionName);
    fprintf(fid, '\n  Elapsed Time      %.3fs', tTotal);
    fprintf(fid, '\n');
    if fid ~=1
        fclose(fid);
    end
end

%% Plot
% if 0
%     plot(1./r1_myoLV,'b')
%     hold on
%     plot(1./r1_myoRV,'r')
%     hold on
%     plot(1./r1_myoLA,'g')
%     hold on
%     plot(1./r1_myoRA,'k')
%     hold on
%     plot(1./r1_bldplLV,'b--')
%     hold on
%     plot(1./r1_bldplRV,'r--')
%     hold on
%     plot(1./r1_bldplLA,'g--')
%     hold on
%     plot(1./r1_bldplRA,'k--')
%     legend('Myo LV','Myo RV','Myo LA','Myo RA','Blood LV','Blood RV',...
%         'Blood LA','Blood RA')
% end

%% Convert relaxation times to sec
tissueValues(:,1:2) = tissueValues(:,1:2)/1000;

end


%% Other functions
function c = convolve(t, c, h, q)
dt = (t(end)-t(1)) / (length(t)-1);

c_ = zeros(length(c));
for i=1:length(c)
    for j =1:length(c)
        if(i+1-j>0)
            c_(i,j) = c(i+1-j);
        end
    end
end
if size(h,1) == 1
    h=h';
end
c = q*c_*h*dt;
end