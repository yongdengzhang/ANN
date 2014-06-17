%% This script demonstrates the generation of training data for ANN
%   Yongdeng Zhang & Lusheng Gu
%   5th of May, 2012
%%*************************************************************************

clear
clc

% parameters for initialization
prompt={'Enter pixel size (nm):',...
        'Enter emission wavelength (nm):'...
        'Enter glass/oil refractive index:'...
        'Enter sample/water refractive index:'...
        'Enter magnification:'...
        'Enter numerical aperture:'...
        'Enter flag for tirf illumination:'...
        'Enter the size of image (pixel):'...
        'Enter incidence angle (rad):'...
        'Enter image type (free,fixed,restricted):'...
        'Enter training data number:'};
name='Input for initialization';
numlines=1;
defaultanswer={'100','580','1.518','1.33','160','1.45','1','7','1.2','1','1000'};
answer=inputdlg(prompt,name,numlines,defaultanswer);
if ~isempty(answer)
    self.a = str2double(answer(1));             % pixel size(nm)      
    self.wl = str2double(answer(2));            % emission wavelength(nm)
    self.n = str2double(answer(3));             % glass/oil
    self.n0 = str2double(answer(4));            % sample/water
    self.M = str2double(answer(5));             % magnification
    self.NA = str2double(answer(6));            % numerical aperture 
    self.tirf = str2double(answer(7));          % TIRF illumination
    self.npix = str2double(answer(8));          % the size of image
    self.thetai = str2double(answer(9));        % incidence angle
    type = str2double(answer(10));              % 1:free dipole;2: fixed dipole;3:restricted dipole
    N = str2double(answer(11));                 % number of training data
else
    return
end

% generate training data                    
cx = (rand(N,1)-0.5);                   % x position of molecule(pixel)
cy = (rand(N,1)-0.5);                   % y position of molecule(pixel)
phi = zeros(N,1);                       % azimuthal angle(0-360бу)
theta = zeros(N,1);                     % polar angle(0-90бу)
delta = zeros(N,1);                     % half cone angle(0-180бу)
defocus = zeros(N,1);                   % defocus (um)
switch type
    case 1  % free dipole
        delta = ones(N,1)*180;  
        filename = 'data_free.mat';
    case 2  % fixed dipole
        phi = 360*rand(N,1);      
        theta = 90*rand(N,1);   
        filename = 'data_fixed.mat';
    case 3  % restricted dipole
        phi = 360*rand(N,1);     
        theta = 90*rand(N,1);     
        delta = 90*rand(N,1);
        filename = 'data_restricted.mat';
end

% generate fixed data
info = [cx,cy,phi,theta,delta,defocus];
data = ANN_getPSF(info,self,1);

% save data
info = info';
save(filename,'data','info','type'); 