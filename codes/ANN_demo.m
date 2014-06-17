%% This script demonstrates the use of ANN for simulation of real images
%% Before simulation, you should finish each training first
%   Yongdeng Zhang & Lusheng Gu
%   6th of May, 2012
%**************************************************************************

clear
clc

% parameters for initialization
self.a = 100;            % pixel size(nm)      
self.wl = 580;           % emission wavelength(nm)
self.n = 1.518;          % glass/oil
self.n0 = 1.33;          % sample/water
self.M = 160;            % magnification
self.NA = 1.45;          % numerical aperture 
self.tirf = 1;           % TIRF illumination
self.approx = 1;         % small pixel approximation
self.npix = 7;           % the size of image
self.thetai = 1.2;       % incidence angle

% parameters for data processing
conversion = 10.5;       % conversion factor (photons per AD count at unit gain)
gain = 300;              % gain of EMCCD

%% load real data
%% data1 is an image of biotin-conjugated mEos2, hence it represents free dipole
I = double(imread('free.tif'));
I_raw = I;
I = reshape(I,self.npix^2,1);
maxI = max(I(:));
I = I/maxI;

% load net
P = importdata('net_free.mat');

% simulation
info(1:2) = sim(P.Net.position,I);
info(3:4) = 0;
info(5) = 180;
info(6) = 0;

% get fit image
I_fit = ANN_getPSF(info,self);
outputs = sim(P.Net.photon,I);
outputs(1) = outputs(1).*maxI*self.npix^2;
outputs(2) = outputs(2).*maxI;
outputs(1) = outputs(1)-outputs(2)*self.npix^2;
info(7) = outputs(1);
I_fit = I_fit*info(7);

figure; 
hold on
I_res = I_raw-I_fit;
subplot(3,3,1);imshow(I_raw,[]);ylabel('Raw');title('free dipole');
subplot(3,3,4);imshow(I_fit,[]);ylabel('Simulation')
subplot(3,3,7);imshow(I_res,[]);ylabel('Residual')

%% data2 is an image of fixed mEos2 immobilized in PMMA layer, hence it represents fixed dipole
I = double(imread('fixed.tif'));
I_raw = I;
I = reshape(I,self.npix^2,1);
maxI = max(I(:));
I = I/maxI;

% load net
P = importdata('net_fixed.mat');

% simulation
info(1:2) = sim(P.Net.position,I);
out = sim(P.Net.phi,I);    
k = sqrt(out(1)^2+out(2)^2);
k=[k
    k];
out(1:2) = out(1:2)./k;
if out(1)>=0 && out(2)>=0       % 0-90
    outputs(3) = asind(out(1));
elseif out(1)>=0 && out(2)<0    % 90-180
    outputs(3) = 180-asind(out(1));
elseif out(1)<0 && out(2)>=0    % 270-360
    outputs(3) = 360+asind(out(1));
elseif out(1)<0 && out(2)<0     % 180-270
    outputs(3) = 180-asind(out(1));
end
outputs(4) = sim(P.Net.theta,I);   
outputs(4,outputs(4,:)>1) = 1;
outputs(4,outputs(4,:)<0) = 0;
outputs(4) = outputs(4)*90;
info(3:4) = outputs(3:4)*pi/180;
info(5:6) = 0;

% get fit image
I_fit = ANN_getPSF(info,self);
outputs = sim(P.Net.photon,I);
outputs(1) = outputs(1).*maxI*self.npix^2;
outputs(2) = outputs(2).*maxI;
outputs(1) = outputs(1)-outputs(2)*self.npix^2;
info(7) = outputs(1);
I_fit = I_fit*info(7);

I_res = I_raw-I_fit;
subplot(3,3,2);imshow(I_raw,[]);title('fixed dipole');
subplot(3,3,5);imshow(I_fit,[]);
subplot(3,3,8);imshow(I_res,[]);

%% data3 is an image of mEos2 from fixed cell PALM experiment, hence it represents restricted dipole
I = double(imread('restricted.tif'));
I_raw = I;
I = reshape(I,self.npix^2,1);
maxI = max(I(:));
I = I/maxI;

% load net
P = importdata('net_restricted.mat');

% simulation
info(1:2) = sim(P.Net.position,I);
out = sim(P.Net.phi,I);    
out(1) = out(1)./out(3);
out(2) = out(2)./out(3);
out(1) = min(max(out(1),-1),1);
out(2) = min(max(out(2),-1),1);
if out(1)>=0 && out(2)>=0       % 0-90
    outputs(3) = asind(out(1));
elseif out(1)>=0 && out(2)<0    % 90-180
    outputs(3) = 180-asind(out(1));
elseif out(1)<0 && out(2)>=0    % 270-360
    outputs(3) = 360+asind(out(1));
elseif out(1)<0 && out(2)<0     % 180-270
    outputs(3) = 180-asind(out(1));
end
out = sim(P.Net.theta,I); 
k = sqrt(out(1)^2+out(2)^2);
out(1) = out(1)/k;
out(2) = out(2)/k;
if out(1)>=0 && out(2)>=0       % 0-90
    outputs(4) = asind(out(1))/2;
elseif out(1)>=0 && out(2)<0    % 90-180
    outputs(4) = (180-asind(out(1)))/2;
end
k = min(max(k,0),1);
outputs(5) = acosd(sqrt(2*k+1/4)-1/2);
info(3:5) = outputs(3:5);
info(6) = 0;

% get fit image
I_fit = ANN_getPSF(info,self);
outputs = sim(P.Net.photon,I);
outputs(1) = outputs(1).*maxI*self.npix^2;
outputs(2) = outputs(2).*maxI;
outputs(1) = outputs(1)-outputs(2)*self.npix^2;
info(7) = outputs(1);
I_fit = I_fit*info(7);

I_res = I_raw-I_fit;
subplot(3,3,3);imshow(I_raw,[]);title('restricted dipole');
subplot(3,3,6);imshow(I_fit,[]);
subplot(3,3,9);imshow(I_res,[]);
hold off