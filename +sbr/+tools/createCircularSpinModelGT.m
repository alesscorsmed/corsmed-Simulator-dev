function spinModelGT = createCircularSpinModelGT(T1,T2,PDvalue)
% T1      = 1;
% T2      = 0.25;
% PDvalue = 1;

% FOV of the image
FOVx    = 0.3;
FOVy    = 0.3;

% Grid size
pointsx = 128; 
pointsy = 128;

resolution = [FOVx/pointsx,FOVy/pointsy,0.001];

circleRadius = 0.1;

[X,Y] = meshgrid(-FOVx/2:resolution(1,1):FOVx/2,-FOVy/2:resolution(1,2):FOVy/2);

x = X(:);
y = Y(:);

radius  = sqrt(x.^2+y.^2);

% Find spins that are within a circle with radius equal to circleRadius
inds    = find(radius<=circleRadius);

tissueType          = 2*ones(size(x));
tissueType(inds)    = 1;

pd          = zeros(size(x));
pd(inds)    = PDvalue;

spinModelGT.numIsochromats	= size(x,1);
spinModelGT.resolution      = resolution;
spinModelGT.mu           	= 1;
spinModelGT.b0            	= 1;
spinModelGT.x              	= x;
spinModelGT.y              	= y;
spinModelGT.z              	= zeros(size(x));
spinModelGT.bi             	= zeros(size(x));
spinModelGT.pd           	= pd;
spinModelGT.tissueDiff     	= zeros(size(x,1),3);
spinModelGT.tissueType    	= tissueType;
spinModelGT.rxCoilMapsX    	= ones(size(x));
spinModelGT.rxCoilMapsY   	= zeros(size(x));
spinModelGT.numRxCoils     	= 1;
spinModelGT.tissueValues    = [T1,T2,1,0,0,0;zeros(1,6)];

% temp = reshape(tissueType,129,129);
% imshow(temp)