function outputs = ANN_simulation(data,Net,type,flag,conversion,gain)
% Input:    data,raw images for simulation (npix*npix*N)
%           Net, trained neural networks
%           type, 1:free dipoe;2,fixed dipole;3:restricted dipole
%           flag, [1,1,1,1],determine which net for simulation
%           conversion,conversion factor (photons per AD count at unit
%           gain,gain of EMCCD
% Output:   outputs,7*N array,[cx,cy,phi,theta,delta,photonN,bg]

% preprocessing of inout/target paris
npix = length(data(:,1,1));
N = length(data(1,1,:));
data_noise = double(reshape(data,npix^2,N));
maxI=zeros(1,N);
for i=1:N
    maxI(1,i)=max(data_noise(:,i));
    data_noise(:,i)=data_noise(:,i)/maxI(1,i); 
end

outputs = zeros(8,N);
if flag(1)
    outputs(1:2,:) = sim(Net.position,data_noise);
end
if flag(2)
    outputs(6:7,:) = sim(Net.photon,data_noise);
    outputs(6,:) = outputs(6,:).*maxI*npix^2;
    outputs(7,:) = outputs(7,:).*maxI;
    outputs(6,:) = (outputs(6,:)-outputs(7,:)*npix^2)*conversion/gain;
end
if flag(3)
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

if flag(4)
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