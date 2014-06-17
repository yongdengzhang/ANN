%% This script demonstrates the training of ANN
%   Note: you should generate traing data first!
%   Using dipimage 2.3 (http://www.diplib.org/)
%   Yongdeng Zhang & Lusheng Gu
%   5th of May, 2012
%%*************************************************************************

clear
clc

% load training data
[filename, pathname] = uigetfile( ...
    {'*.mat', 'All MAT-Files (*.mat)'; ...
        '*.*','All Files (*.*)'}, ...
    'Open Data File',cd);
if isequal([filename,pathname],[0,0])
    return
else
    FileStr = fullfile(pathname,filename);
end
load(FileStr);      

% parameters for data processing
bg = 150;                       % constant background(intensity)
var = 42;                       % variance of Gaussian noise(intensity)
conversion = 10.5;              % conversion factor (photons per AD count at unit gain)
gain = 300;                     % gain of EMCCD
npix = length(data(:,:,1));     % the size of image
N = length(data(1,1,:));        % the number of images
photonN = rand(N,1)*10000+10;   % photon number for each image

% photon and background information
info(7,:) = photonN;
info(8,:) = bg;

% convert photons into intensity
for i=1:N
    data(:,:,i) = data(:,:,i)*photonN(i);
end
data = data*gain/conversion;

% corrupted with Poisson and Gaussian noise
dipstart;   % start dipimage
if conversion>0
    data_noise = noise(data,'poisson',conversion/gain);
    data_noise = dip_array(data_noise);
end
data_noise = data_noise+bg;
if var>0
    data_noise = noise(data_noise,'gaussian',var);
    data_noise = dip_array(data_noise);
end

% preprocessing of inout/target pairs
data_noise = double(reshape(data_noise,npix^2,N));
maxI=zeros(1,N);
for i=1:N
    maxI(1,i)=max(data_noise(:,i));
    data_noise(:,i)=data_noise(:,i)/maxI(1,i); 
end

% training flag
flag1 = 1;
flag2 = 1;
flag3 = 1;
flag4 = 1;
Net.position = [];
Net.photon = [];
Net.phi = [];
Net.theta = [];

if flag1
    % create network for position
    numHiddenNeurons = 10;                      % number of hidden neurons
    targets = info(1:2,:);                      % training target      
    net = newfit(data_noise,targets,numHiddenNeurons);
    net.trainFcn = 'trainlm';                   % training function
    net.divideFcn = 'divideint';                % division function              
    net.performFcn = 'msne';                    % performance function
    net.divideParam.trainRatio = 70/100;        % adjust as desired
    net.divideParam.valRatio = 15/100;          % adjust as desired
    net.divideParam.testRatio = 15/100;         % adjust as desired
    net.trainParam.max_fail = 6;                % maximum validation failures
    net.trainParam.epochs = 100;                % maximum number of epochs to train
    % training
    [Net.position,tr] = train(net,data_noise,targets); 
end

if flag2
    % create network for photons
    numHiddenNeurons = 10;                      
    targets(1,:) = info(7,:)*gain/conversion./maxI/(npix^2)+info(8,:)./maxI;
    targets(2,:) = info(8,:)./maxI;
    net = newfit(data_noise,targets,numHiddenNeurons);
    net.trainFcn = 'trainlm';                  
    net.divideFcn = 'divideint';               
    net.performFcn = 'msne';                    
    net.divideParam.trainRatio = 70/100;        
    net.divideParam.valRatio = 15/100;          
    net.divideParam.testRatio = 15/100;         
    net.trainParam.max_fail = 6;                
    net.trainParam.epochs = 100;               
    [Net.photon,tr] = train(net,data_noise,targets); 
end

if type==1  % free dipole
    flag3 = 0;
    flag4 = 0;
end

if flag3
    % create network for phi
    numHiddenNeurons = 10;                     
    if type==3      % restricted
        targets(1,:)=sind(info(3,:)).*sind(info(4,:)*2).*abs(cosd(info(5,:))).*(1+cosd(info(5,:)))/2;
        targets(2,:)=cosd(info(3,:)).*sind(info(4,:)*2).*abs(cosd(info(5,:))).*(1+cosd(info(5,:)))/2;
        targets(3,:)=sind(info(4,:)*2).*abs(cosd(info(5,:))).*(1+cosd(info(5,:)))/2;
    elseif type==2  % fixed
        targets(1,:)=sind(info(3,:)).*sind(info(4,:)*2);
        targets(2,:)=cosd(info(3,:)).*sind(info(4,:)*2);  
    end
    net = newfit(data_noise,targets,numHiddenNeurons);
    net.divideFcn = 'divideint';               
    net.performFcn = 'msne';                    
    net.divideParam.trainRatio = 70/100;        
    net.divideParam.valRatio = 15/100;          
    net.divideParam.testRatio = 15/100;        
    net.trainFcn='trainlm';
    net.trainParam.max_fail = 6;
    net.trainParam.epochs = 100;
    [Net.phi,tr] = train(net,data_noise,targets); 
end

if flag4
    % create network for theta and delta
    targets = [];
    numHiddenNeurons = 10;                      
    if type==3      % restricted
        targets(1,:)=sind(info(4,:)*2).*abs(cosd(info(5,:))).*(1+cosd(info(5,:)))/2;
        targets(2,:)=cosd(info(4,:)*2).*abs(cosd(info(5,:))).*(1+cosd(info(5,:)))/2;
    elseif type==2  % fixed
        targets(1,:)=info(4,:)/90;  
    end
    net = newfit(data_noise,targets,numHiddenNeurons);
    net.divideFcn = 'divideint';                
    net.performFcn = 'msne';                   
    net.divideParam.trainRatio = 70/100;        
    net.divideParam.valRatio = 15/100;          
    net.divideParam.testRatio = 15/100;        
    net.trainFcn='trainlm';
    net.trainParam.max_fail = 6;
    net.trainParam.epochs = 100;
    [Net.theta,tr] = train(net,data_noise,targets); 
end

% simulation with training data
outputs = zeros(8,N);
if flag1
    outputs(1:2,:) = sim(Net.position,data_noise);
end
if flag2
    outputs(7:8,:) = sim(Net.photon,data_noise);
    outputs(7,:) = outputs(7,:).*maxI*npix^2;
    outputs(8,:) = outputs(8,:).*maxI;
    outputs(7,:) = (outputs(7,:)-outputs(8,:)*npix^2)*conversion/gain;
end
if flag3
    if type==3 % restricted
        out = sim(Net.phi,data_noise);    
        n = length(out(1,:));
        out(1,:) = out(1,:)./out(3,:);
        out(2,:) = out(2,:)./out(3,:);
        for i=1:N
            out(1,i) = min(max(out(1,i),-1),1);
            out(2,i) = min(max(out(2,i),-1),1);
        end
        for i = 1:n
            if out(1,i)>=0 && out(2,i)>=0       % 0-90
                outputs(3,i) = asind(out(1,i));
            elseif out(1,i)>=0 && out(2,i)<0    % 90-180
                outputs(3,i) = 180-asind(out(1,i));
            elseif out(1,i)<0 && out(2,i)>=0    % 270-360
                outputs(3,i) = 360+asind(out(1,i));
            elseif out(1,i)<0 && out(2,i)<0     % 180-270
                outputs(3,i) = 180-asind(out(1,i));
            end
        end
    elseif type==2  % fixed  
        out = sim(Net.phi,data_noise);    
        n = length(out(1,:));
        k = sqrt(out(1,:).^2+out(2,:).^2);
        k=[k
            k];
        out(1:2,:) = out(1:2,:)./k;
        for i=1:n
            if out(1,i)>=0 && out(2,i)>=0       % 0-90
                outputs(3,i) = asind(out(1,i));
            elseif out(1,i)>=0 && out(2,i)<0    % 90-180
                outputs(3,i) = 180-asind(out(1,i));
            elseif out(1,i)<0 && out(2,i)>=0    % 270-360
                outputs(3,i) = 360+asind(out(1,i));
            elseif out(1,i)<0 && out(2,i)<0     % 180-270
                outputs(3,i) = 180-asind(out(1,i));
            end
        end
    end
end

if flag4
    if type==3 % restricted
        out = sim(Net.theta,data_noise); 
        n = length(out(1,:));
        k = sqrt(out(1,:).^2+out(2,:).^2);
        out(1,:) = out(1,:)./k;
        out(2,:) = out(2,:)./k;
        for i = 1:n
            if out(1,i)>=0 && out(2,i)>=0       % 0-90
                outputs(4,i) = asind(out(1,i))/2;
            elseif out(1,i)>=0 && out(2,i)<0    % 90-180
                outputs(4,i) = (180-asind(out(1,i)))/2;
            end
        end
        for i = 1:n
            k(1,i) = min(max(k(1,i),0),1);
        end
        outputs(5,:)=acosd(sqrt(2*k+1/4)-1/2);
    elseif type==2  % fixed  
        outputs(4,:)=sim(Net.theta,data_noise);   
        outputs(4,outputs(4,:)>1)=1;
        outputs(4,outputs(4,:)<0)=0;
        outputs(4,:)=outputs(4,:)*90;
    end
end

% calculate errors
error=outputs-info;
for i=1:N
    if error(3,i)>180
        error(3,i)=error(3,i)-360;
    end
    if error(3,i)<-180
        error(3,i)=360+error(3,i);
    end
end

save('results.mat','outputs','info','error');
if flag1
    error(8,:) = sqrt(error(1,:).^2+error(2,:).^2);
    e = sqrt(sum(error(8,:).^2)/N);
    disp(['ANN position error RMSE (pixel):',num2str(e)]);
end
if flag2
    e = sqrt(sum(error(7,:).^2)/N);
    disp(['ANN photon error RMSE (photon):',num2str(e)]);
end
if flag3
    e = sqrt(sum(error(3,:).^2)/N);
    disp(['ANN phi error RMSE (degree):',num2str(e)]);
end
if flag4
    e = sqrt(sum(error(4,:).^2)/N);
    disp(['ANN theta error RMSE (degree):',num2str(e)]);
    if type==3
        e = sqrt(sum(error(5,:).^2)/N);
        disp(['ANN delta error RMSE (degree):',num2str(e)]);
    end
end

% save results
switch type
    case 1
        filestr = 'net_free.mat';
    case 2
        filestr = 'net_fixed.mat';
    case 3
        filestr = 'net_restricted.mat';
end
save(filestr,'Net','type');