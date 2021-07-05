function [acquisition] = initializeAcquisition()
%
% DATA.EXPERIMENT.INITIALIZEACQUISITION
%
%	Function that initializes acquisition data structure.
%   Returns a acquisition with empty fields, or filed if connDB is passed.
%
%   This function is useful to define the fields.
%
% INPUT
%   expControl   struct with experiment info
%
% OUTPUT
%   acquisition   acquisition structure
%
%========================  CORSMED AB © 2020 ==============================
%

%% generate basic default data
acquisition.data 	= data.experiment.initializeAcqData();
%% main RF pulse
acquisition.mainRF  = data.experiment.initializeAcqRF(90);
%% refocusing pulse
acquisition.refRF 	= data.experiment.initializeAcqRF(180);
%% IR preparation info
acquisition.prepIR 	= data.experiment.initializeAcqIR();
%% PG encoding info
acquisition.encPG  	= data.experiment.initializeAcqPG();
%% noise info
acquisition.noise   = data.experiment.initializeAcqNoise();
