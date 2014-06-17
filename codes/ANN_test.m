%% This script demostrates basic the basic steps of ANN for single molecule localization
%   Yongdeng Zhang & Lusheng Gu
%   6th of May, 2012
%**************************************************************************

clear
clc

% load images
[filename, pathname] = uigetfile( ...
    {'*.tif;*.tiff', 'All TIF-Files (*.tif,*.tiff)'; ...
        '*.*','All Files (*.*)'}, ...
    'Open Image File');
if isequal([filename,pathname],[0,0])
    return
else
    FileStr = fullfile(pathname,filename);
end
A = tiffread(FileStr);
N = length(A(1,1,:));
[row column] = size(A(:,:,1));

% find single molecules
V = [];
ROI = [];
n = 1;
for i=1:N
    I = single(A(:,:,i));
    [out,M] = Gauss2D(I,1);
    ImgMax = locmax2d(out,[5 5]);
    [y,x] = find(ImgMax>500);
    for j=1:length(y)
        if x(j)>3 && x(j)<=column-3 && y(j)>3 && y(j)<=row-3
           ROI(:,:,n)=I(y(j)-3:y(j)+3,x(j)-3:x(j)+3);
           V(n,1) = y(j);
           V(n,2) = x(j);
           n = n+1;
        end
    end
end

% simulation using trained neural network
load net_free.mat;
outputs = ANN_simulation(ROI,Net,type,[1,1,0,0],10.5,300);

% reconstruction of image
amp = 2;
Inew = zeros(row*amp,column*amp);
for i=1:length(V(:,1))
    cy = (V(i,1)+outputs(1,i))*amp;
    cx = (V(i,2)+outputs(2,i))*amp;
    Inew(round(cy),round(cx)) = 1;
end
figure;imshow(Inew,[]);