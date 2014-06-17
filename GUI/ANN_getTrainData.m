function [data,info]=ANN_getTrainData(handles)
   
    % parameters for initialization
    self.a = str2double(get(handles.edit3,'string'));               % pixel size(nm)      
    self.wl = str2double(get(handles.edit6,'string'));              % emission wavelength(nm)
    self.n = str2double(get(handles.edit8,'string'));               % glass/oil
    self.n0 = str2double(get(handles.edit7,'string'));              % sample/water
    self.M = str2double(get(handles.edit2,'string'));               % magnification
    self.NA = str2double(get(handles.edit4,'string'));              % numerical aperture 
    self.tirf = str2double(get(handles.edit16,'string'));           % tirf illumination         
    self.npix = str2double(get(handles.edit5,'string'));            % the size of image
    self.thetai = 1.2;                                              % incidence angle
    
    n=str2double(get(handles.edit1,'string'));
    cx=(2*rand(n,1)-1)*0.5;
    cy=(2*rand(n,1)-1)*0.5;

    fixed=get(handles.radiobutton3,'value');
    free=get(handles.radiobutton4,'value');
    restricted=get(handles.radiobutton5,'value');

    if fixed || restricted
        phi=rand(n,1)*360;
        theta=rand(n,1)*90;
    end
    if free
       phi=zeros(n,1);
       theta=zeros(n,1);
    end   
    if fixed
        delta=zeros(n,1);  
    elseif free
        delta=zeros(n,1)+180;  
    elseif restricted
        delta=rand(n,1)*90;
    end

    focus=(2*rand(n,1)-1)*str2double(get(handles.edit15,'string'));
    
    N=zeros(n,1);
    B=zeros(n,1);
    
    info=[cx,cy,phi,theta,delta,focus,N,B]';
    info1=[cx,cy,phi,theta,delta,focus];
    data=ANN_getPSF(info1,self,1);           
end