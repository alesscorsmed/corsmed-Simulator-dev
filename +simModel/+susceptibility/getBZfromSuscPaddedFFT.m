function  Bz = getBZfromSuscPaddedFFT(Svol,B0, Xe, Msz)
%% Method:
%	Marques, J.P. and Bowtell, R.
%	Concepts in Magnetic Resonance Part B: MR engineering 25, pages 65-78
%   2005
%

N = size(Svol);
Nhalf=round(N/2);
%Msz = 2^nextpow2(Msz);
M = [Msz, Msz, Msz];
 
indice1 = (M(1)/2)-Nhalf(1)+1:(M(1)/2)-Nhalf(1)+N(1);
indice2 = (M(2)/2)-Nhalf(2)+1:(M(2)/2)-Nhalf(2)+N(2);
indice3 = (M(3)/2)-Nhalf(3)+1:(M(3)/2)-Nhalf(3)+N(3);
stemp=Svol; 
Svol=Xe*ones(M); 
Svol(indice1, indice2, indice3)=stemp;


% Kspace coordinates 
[Nx,Ny,Nz]=size(Svol);
szS = [Nx,Ny,Nz];
fx=(0:1/Nx:.5); fx=[fx -fliplr(fx(2:end-1))];   %freq. axis suitable for FFT
fy=(0:1/Ny:.5); fy=[fy -fliplr(fy(2:end-1))];   %freq. axis suitable for FFT
fz=(0:1/Nz:.5); fz=[fz -fliplr(fz(2:end-1))];   %freq. axis suitable for FFT
[kx,ky,kz]=ndgrid(fx,fy,fz);
Dk=(kz.^2)./(kx.^2+ky.^2+kz.^2);    %oops, devide by zero
Dk(1,1,1)=0;                        %we don't care and set the NaN to zero

% FFT:
Sk=fftn(Svol,szS);
Bz = (1/3-Dk).*Sk;
Bz=real(ifftn(Bz,szS));

% rescale to corrected average:
Bz = (B0*Xe/3)*Bz/mean(Bz(:));

% Remove padding:
Bz=Bz(indice1, indice2, indice3);