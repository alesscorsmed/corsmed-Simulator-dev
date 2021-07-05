clear all; close all; clc;

numX = 256;
numY = 128;
numZ = 3;
numCH = 2;
numCI = 2;

x = linspace(-1,1,numX);
y = linspace(-1,1,numY);
if numZ == 1
z = 1;
else
    z = linspace(-1,1,numZ);
end
[yq,xq,zq] = meshgrid(y,x,z);

% circle
circleImage = zeros(size(xq));
circleImage( abs((xq+0.1) + 1j*(yq-0.2)) <= 0.5 ) = 1;

% images and kSpace
circleDSpace = zeros(numX,numY,numZ,numCH,numCI);
circleKSpace = zeros(numX,numY,numZ,numCH,numCI);
for ic=1:numCI
    for cc=1:numCH
        currentImage = sqrt(1/numCH)*ic*circleImage;
        circleDSpace(:,:,:,cc,ic) = currentImage;
        circleKSpace(:,:,:,cc,ic) = fftshift(fftn(fftshift(currentImage)));
    end
end
circleISpace = sqrt(sum(conj(circleDSpace).*circleDSpace,4));

% recon
kSpaceInfo.xPadd = 256;
kSpaceInfo.yPadd = 256;
kSpaceInfo.zPadd = 1;
expControl.debug.debugMode = 1;
expControl.debug.debugFile = [];

[iSpace] = reconstruction.coreFFT.sosFFT( circleKSpace, kSpaceInfo, expControl );

for ic=1:numCI
    for zz=1:numZ
    figure();
    % original image
    subplot(2,2,1);
    imagesc(abs(ic*circleImage(:,:,zz)));
    axis image; colorbar;
    title(sprintf('Original slice %d contrast %d', zz, ic));
    
    % original image from coils
    subplot(2,2,3);
    imagesc(abs(circleISpace(:,:,zz,1,ic)));
    axis image; colorbar;
    title(sprintf('Original SOS slice %d contrast %d', zz, ic));
    
    % original image
    subplot(2,2,2);
    imagesc(abs(iSpace(:,:,zz,1,ic)));
    axis image; colorbar;
    title(sprintf('Recon Abs slice %d contrast %d', zz, ic));
    
    % original image
    subplot(2,2,4);
    imagesc(angle(iSpace(:,:,zz,1,ic)));
    axis image; 
    title(sprintf('Recon Phase slice %d contrast %d', zz, ic));
    
    end
end
