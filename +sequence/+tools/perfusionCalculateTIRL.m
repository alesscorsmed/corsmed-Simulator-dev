function TIRL = perfusionCalculateTIRL(expControlPerfusion,anatomicalModel)

cardiacCycleDuration = 60/anatomicalModel.HR; % [sec]

TIRL = expControlPerfusion.contrasts * cardiacCycleDuration;