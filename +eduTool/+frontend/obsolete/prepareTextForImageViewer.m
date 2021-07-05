function textForImage = prepareTextForImageViewer(acquisition,pulseSequence,...
    reconData,spinModel,expControl,sliceNum,conNum)

pulseSeqFamilyName 	= acquisition.data.pulseSeqFamilyName;
pulseSeqNum         = acquisition.data.pulseSeqNum;

%% Bottom Right
textForImage.RBot{3} = datestr(now,'yyyy-mm-dd HH:MM:SS');

%% Top Left
textForImage.LTop{1} = [num2str(pulseSeqNum),' - ',pulseSeqFamilyName];
textForImage.LTop{2} = ['Slice ',num2str(sliceNum),'/',num2str(size(reconData.slice,2))];

if reconData.numC > 1
    textForImage.LTop{3} = ['Contrast ',num2str(conNum),'/',...
        num2str(reconData.numC)];
end

%% Return if the old pulse sequence generator is utilized
if expControl.useOldSequence
    textForImage.RTop{1} = '';
    textForImage.LBot{1} = '';
    return;
end

%% Top Right

gridStep                = expControl.model.gridStep;

textForImage.RTop{1}    = spinModel.bodyPartName;
textForImage.RTop{2}    = ['Grid: ',num2str(gridStep(1,1)*1000),...
    'x',num2str(gridStep(1,2)*1000),...
    'x',num2str(gridStep(1,3)*1000),' mm'];
textForImage.RTop{3}    = expControl.model.coilType;


%% Bottom Left
FOV(1,1)            = acquisition.reconstruction.ismrmrd.header.encoding.reconSpace.fieldOfView_mm.x;
FOV(1,2)            = acquisition.reconstruction.ismrmrd.header.encoding.reconSpace.fieldOfView_mm.y;
kspace(1,1)         = acquisition.reconstruction.ismrmrd.header.encoding.reconSpace.matrixSize.x;
kspace(1,2)         = acquisition.reconstruction.ismrmrd.header.encoding.reconSpace.matrixSize.y;

textForImage.LBot{1} = ['FOV: ',num2str(FOV(1,1)),'x',num2str(FOV(1,2)),'mm'];
textForImage.LBot{2} = ['MATRIX: ',num2str(kspace(1,1)),'x',num2str(kspace(1,2))];

if ~strcmp(acquisition.data.partialFourier,'no')
    if strcmp(acquisition.data.partialFourier,'readConjugate')
        partialFourierText = [' (',num2str(round(acquisition.data.fFactor*kspace(1,1))),...
            'x',num2str(kspace(1,2)),')'];
    elseif strcmp(acquisition.data.partialFourier,'phaseConjugate')
        partialFourierText = [' (',num2str(kspace(1,1)),...
            'x',num2str(round(acquisition.data.fFactor*kspace(1,2))),')'];
    else
        partialFourierText = '';
    end
else
    partialFourierText = '';
end


if strcmp(pulseSeqFamilyName,'TSE') || strcmp(pulseSeqFamilyName,'IR-TSE') || ...
        strcmp(pulseSeqFamilyName,'SS-FSE')
    
    textForImage.LBot{3} = sprintf('TR/TEeff: %d/%.1fms - ETL: %d',...
        pulseSequence.TR*1000,pulseSequence.effTE*1000,pulseSequence.numTL);
    
    if strcmp(acquisition.data.parallelImaging,'sense')
        if strcmp(pulseSeqFamilyName,'IR-TSE') || ...
            (strcmp(pulseSeqFamilyName,'TSE') && ...
            strcmp(acquisition.data.fatsat,'lipidir'))
            textForImage.LBot{4} = sprintf('TI= %.1fms - SENSE: %d',...
                pulseSequence.TI*1e3, acquisition.data.rFactor);
        else
            textForImage.LBot{4} = sprintf('SENSE: %d',...
                acquisition.data.rFactor);
        end
    else
        if strcmp(pulseSeqFamilyName,'IR-TSE') || ...
                (strcmp(pulseSeqFamilyName,'TSE') && ...
                strcmp(acquisition.data.fatsat,'lipidir'))
            textForImage.LBot{4} = sprintf('TI= %.1fms',...
                pulseSequence.TI*1e3);
        else
            textForImage.LBot{4} = '';
        end
    end
elseif strcmp(pulseSeqFamilyName,'SE') || strcmp(pulseSeqFamilyName,'IR-SE')

    textForImage.LBot{2}        = [textForImage.LBot{2},' ',partialFourierText];
    textForImage.LBot{3}        = ['TR/TE: ',num2str(acquisition.data.TR*1000),...
        '/',num2str(acquisition.data.TE*1000),'ms'];
    if strcmp(acquisition.data.parallelImaging,'sense')
        if strcmp(pulseSeqFamilyName,'IR-SE')
            textForImage.LBot{4} = ['TI: ',num2str(acquisition.data.TI*1000),...
                'ms',' - SENSE: ',struct_pulseq.rfactor];
        else
            textForImage.LBot{4} = sprintf('SENSE: %d',...
                acquisition.data.rFactor);
        end
    else
        if strcmp(pulseSeqFamilyName,'IR-SE')
            textForImage.LBot{4} = sprintf('TI= %.1fms',...
                acquisition.data.TI*1e3);
        else
            textForImage.LBot{4} = '';
        end
    end
    
    
elseif strcmp(pulseSeqFamilyName,'bSSFP') ||...
        strcmp(pulseSeqFamilyName,'IR-bSSFP') ||...
        strcmp(pulseSeqFamilyName,'MOLLI')
    
    textForImage.LBot{3} = ['TR/TE: ',num2str(acquisition.data.TR*1000),'/',...
        num2str(acquisition.data.TE*1000),'ms'];    
    if strcmp(acquisition.data.parallelImaging,'sense')
        if strcmp(pulseSeqFamilyName,'IR-bSSFP') || strcmp(acquisition.data.fatsat,'lipidir')
            textForImage.LBot{4} = sprintf('TI= %.1fms - SENSE: %d',...
                acquisition.data.TI*1e3, acquisition.data.rFactor);
        else
            textForImage.LBot{4} = sprintf('SENSE: %d',...
                acquisition.data.rFactor);
        end
    else
        if strcmp(pulseSeqFamilyName,'IR-bSSFP') || strcmp(acquisition.data.fatsat,'lipidir')
            textForImage.LBot{4} = sprintf('TI= %.1fms',...
                acquisition.data.TI*1e3);
        else
            textForImage.LBot{4} = '';
        end
    end
    
    
elseif strcmp(pulseSeqFamilyName,'PG-SE-EPI')
    
	% DW information
	if acquisition.data.encPG.Beta > 0
		betaScale = floor(log10(acquisition.data.encPG.Beta));
	else
		betaScale = 0;
	end
	textForImage.LBot{3}        = sprintf('B (%s): %#1.2fe%d s/m2',...
		upper(acquisition.data.encPG.Dir), BETA/(10^betaScale), betaScale);
    
    if strcmp(acquisition.data.parallelImaging,'sense')
        textForImage.LBot{4} = sprintf('SENSE: %d',...
            acquisition.data.rFactor);
    else
		textForImage.LBot{4}    = '';
    end
    
elseif strcmp(pulseSeqFamilyName,'SE-EPI') || ...
        strcmp(pulseSeqFamilyName,'EPI')
    
    textForImage.LBot{3}        = ['ESP: ',num2str(acquisition.data.effESP*1000),'ms'];
    if strcmp(acquisition.data.parallelImaging,'sense')
        textForImage.LBot{4} = sprintf('SENSE: %d',...
            acquisition.data.rFactor);
    else
        textForImage.LBot{4}    = '';
    end
    
elseif strcmp(pulseSeqFamilyName,'MP-RAGE') || ...
        strcmp(pulseSeqFamilyName,'GRE-3D')
        
    textForImage.LBot{3}        = ['TR/TE: ',num2str(acquisition.data.TR*1000),'/',...
        num2str(acquisition.data.TE*1000),'ms'];
    textForImage.LBot{4}        = '';
    
elseif strcmp(pulseSeqFamilyName,'GRE-OOP') || ...
        strcmp(pulseSeqFamilyName,'GRE-conc') || ...
        strcmp(pulseSeqFamilyName,'GRE')    
    
    textForImage.LBot{3}        = ['TR/TE: ',num2str(acquisition.data.TR*1000),'/',...
        num2str(acquisition.data.TE*1000),'ms'];
    if strcmp(acquisition.data.parallelImaging,'sense')
        textForImage.LBot{4} = sprintf('SENSE: %d',...
            acquisition.data.rFactor);
    else
        textForImage.LBot{4}    = '';
    end
    
end
