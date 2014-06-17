function varargout = ANN_gui(varargin)
% ANN_GUI M-file for ANN_gui.fig
%      ANN_GUI, by itself, creates a new ANN_GUI or raises the existing
%      singleton*.
%
%      H = ANN_GUI returns the handle to a new ANN_GUI or the handle to
%      the existing singleton*.
%
%      ANN_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANN_GUI.M with the given input arguments.
%
%      ANN_GUI('Property','Value',...) creates a new ANN_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ANN_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ANN_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ANN_gui

% Last Modified by GUIDE v2.5 09-Jul-2012 21:18:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ANN_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @ANN_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ANN_gui is made visible.
function ANN_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ANN_gui (see VARARGIN)

% Choose default command line output for ANN_gui
handles.output = hObject;
dipstart
clc
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ANN_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ANN_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {'*.mat', 'All MAT-Files (*.mat)'
        '*.*','All Files (*.*)'}, ...
    'Load Data File');
if isequal([filename,pathname],[0,0])
    return
else
    filestr=fullfile(pathname,filename);
end
P=importdata(filestr);
handles.data=P.data;
handles.info=P.info;
handles.type=P.type;
handles.npix=size(handles.data,1);
clear P
guidata(hObject,handles) 

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n=length(handles.info(1,:));
N=str2double(get(handles.edit12,'string')); 
PhotonN=round(rand(n,1)*N)+10;
for i=1:n
    handles.data(:,:,i)=handles.data(:,:,i)*PhotonN(i);
end
handles.info(7,:)=PhotonN;

handles.gain=str2double(get(handles.edit17,'string'));   
handles.factor=str2double(get(handles.edit18,'string')); 
handles.data_noise=handles.data*handles.gain/handles.factor;
conversion=handles.factor/handles.gain;
handles.data_noise=noise(handles.data_noise,'poisson',conversion);

handles.data_noise=dip_array(handles.data_noise);
b=str2double(get(handles.edit13,'string'));   
handles.info(8,:)=b;
handles.data_noise=handles.data_noise+b;
var=str2double(get(handles.edit18,'string')); 
handles.data_noise=noise(handles.data_noise,'gaussian',var);
handles.data_noise=dip_array(handles.data_noise);

npix=handles.npix;
handles.data_noise=double(reshape(handles.data_noise,npix^2,n));
maxI=zeros(1,n);
for i=1:n
    maxI(1,i)=max(handles.data_noise(:,i));
    handles.data_noise(:,i)=handles.data_noise(:,i)/maxI(1,i); 
end
handles.maxI=maxI;
guidata(hObject,handles) 

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.position=get(handles.checkbox9,'value');
handles.photon=get(handles.checkbox10,'value');
handles.phi=get(handles.checkbox11,'value');
handles.theta=get(handles.checkbox12,'value');
handles.focus=get(handles.checkbox13,'value');

if handles.type==2  % free
    handles.phi=0;
    handles.theta=0;
end

npix=handles.npix;
maxI=handles.maxI;
info=handles.info;

numHiddenNeurons=str2double(get(handles.edit20,'string')); 
mem_reduc=str2double(get(handles.edit20,'string')); 
epochs=str2double(get(handles.edit20,'string')); 
if handles.position==1
    targets=info(1:2,:);  
    net = newfit(handles.data_noise,targets,[numHiddenNeurons numHiddenNeurons]);
    net.divideFcn = 'divideint';
    net.performFcn = 'msne';
    net.divideParam.trainRatio = 70/100;  % Adjust as desired
    net.divideParam.valRatio = 15/100;  % Adjust as desired
    net.divideParam.testRatio = 15/100;  % Adjust as desired
    net.trainFcn = 'trainlm';
    net.trainParam.mem_reduc = mem_reduc;
    net.trainParam.max_fail = 6;
    net.trainParam.epochs = epochs;
    [net,tr] = train(net,handles.data_noise,targets); 
    handles.net_position=net;
end
if handles.photon==1
    info(7,:)=info(7,:)*handles.gain/handles.factor;
    info(7,:)=info(7,:)./maxI/(npix^2)+info(8,:)./maxI;
    info(8,:)=info(8,:)./maxI;
    targets=info(7:8,:);    
    net = newfit(handles.data_noise,targets,[numHiddenNeurons numHiddenNeurons]);
    net.divideFcn = 'divideint';
    net.performFcn = 'msne';
    net.divideParam.trainRatio = 70/100;  % Adjust as desired
    net.divideParam.valRatio = 15/100;  % Adjust as desired
    net.divideParam.testRatio = 15/100;  % Adjust as desired    
    net.trainFcn='trainlm';
    net.trainParam.mem_reduc = mem_reduc;
    net.trainParam.max_fail = 6;
    net.trainParam.epochs = epochs;
    [net,tr] = train(net,handles.data_noise,targets); 
    handles.net_photon=net;
end
if handles.phi==1
    if handles.type==3 % restricted
        targets(1,:)=sind(info(3,:)).*sind(info(4,:)*2).*abs(cosd(info(5,:))).*(1+cosd(info(5,:)))/2;
        targets(2,:)=cosd(info(3,:)).*sind(info(4,:)*2).*abs(cosd(info(5,:))).*(1+cosd(info(5,:)))/2;
        targets(3,:)=sind(info(4,:)*2).*abs(cosd(info(5,:))).*(1+cosd(info(5,:)))/2;
    elseif handles.type==1  % fixed
        targets(1,:)=sind(info(3,:)).*sind(info(4,:)*2);
        targets(2,:)=cosd(info(3,:)).*sind(info(4,:)*2);  
    end
    net = newfit(handles.data_noise,targets,[numHiddenNeurons numHiddenNeurons]);
    net.divideFcn = 'divideint';
    net.performFcn = 'msne';
    net.divideParam.trainRatio = 70/100;  % Adjust as desired
    net.divideParam.valRatio = 15/100;  % Adjust as desired
    net.divideParam.testRatio = 15/100;  % Adjust as desired
    net.trainFcn='trainlm';
    net.trainParam.mem_reduc = mem_reduc;
    net.trainParam.max_fail = 6;
    net.trainParam.epochs = epochs;
    [net,tr] = train(net,handles.data_noise,targets); 
    handles.net_phi=net;
end
if handles.theta==1
    targets=[];
    if handles.type==3 % restricted
        targets(1,:)=sind(info(4,:)*2).*abs(cosd(info(5,:))).*(1+cosd(info(5,:)))/2;
        targets(2,:)=cosd(info(4,:)*2).*abs(cosd(info(5,:))).*(1+cosd(info(5,:)))/2;
    elseif handles.type==1  % fixed
        targets(1,:)=info(4,:)/90;  
    end
    net = newfit(handles.data_noise,targets,[numHiddenNeurons numHiddenNeurons]);
    net.divideFcn = 'divideint';
    net.performFcn = 'msne';
    net.divideParam.trainRatio = 70/100;  % Adjust as desired
    net.divideParam.valRatio = 15/100;  % Adjust as desired
    net.divideParam.testRatio = 15/100;  % Adjust as desired
    net.trainFcn='trainlm';
    net.trainParam.mem_reduc = mem_reduc;
    net.trainParam.max_fail = 6;
    net.trainParam.epochs = epochs;
    [net,tr] = train(net,handles.data_noise,targets); 
    handles.net_theta=net;
end
if handles.focus==1
    targets=[];
    targets(1,:)=info(6,:);  
    net = newfit(handles.data_noise,targets,[numHiddenNeurons numHiddenNeurons]);
    net.divideFcn = 'divideint';
    net.performFcn = 'msne';
    net.divideParam.trainRatio = 70/100;  % Adjust as desired
    net.divideParam.valRatio = 15/100;  % Adjust as desired
    net.divideParam.testRatio = 15/100;  % Adjust as desired
    net.trainFcn='trainlm';
    net.trainParam.mem_reduc = mem_reduc;
    net.trainParam.max_fail = 6;
    net.trainParam.epochs = epochs;
    [net,tr] = train(net,handles.data_noise,targets); 
    handles.net_focus=net;
end
guidata(hObject,handles) 

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.type==2
    handles.phi=0;
    handles.theta=0;
end
npix=handles.npix;
maxI=handles.maxI;
info=handles.info;
% info(8,:)=sind(info(6,:)*2).*(1-sind(info(7,:)/2));
n=length(handles.data_noise(1,:));
outputs=zeros(8,n);

if handles.position
    outputs(1:2,:) = sim(handles.net_position,handles.data_noise);
end
if handles.photon
    outputs(7:8,:) = sim(handles.net_photon,handles.data_noise);
    outputs(7,:)=outputs(7,:).*maxI*npix^2;
    outputs(8,:)=outputs(8,:).*maxI;
    outputs(7,:)=(outputs(7,:)-outputs(8,:)*npix^2)*handles.factor/handles.gain;
end
if handles.phi
    if handles.type==3 % restricted
        out=sim(handles.net_phi,handles.data_noise);    
        n=length(out(1,:));
        out(1,:)=out(1,:)./out(3,:);
        out(2,:)=out(2,:)./out(3,:);
        for i=1:n
            out(1,i)=min(max(out(1,i),-1),1);
            out(2,i)=min(max(out(2,i),-1),1);
        end
    %     k=sqrt(out(1,:).^2+out(2,:).^2);
    %     k=[k
    %         k];
    %     out(1:2,:)=out(1:2,:)./k;
        for i=1:n
            if out(1,i)>=0 && out(2,i)>=0       % 0-90
                outputs(3,i)=asind(out(1,i));
            elseif out(1,i)>=0 && out(2,i)<0    % 90-180
                outputs(3,i)=180-asind(out(1,i));
            elseif out(1,i)<0 && out(2,i)>=0    % 270-360
                outputs(3,i)=360+asind(out(1,i));
            elseif out(1,i)<0 && out(2,i)<0     % 180-270
                outputs(3,i)=180-asind(out(1,i));
            end
        end
    elseif handles.type==1  % fixed  
        out = sim(handles.net_phi,handles.data_noise);    
        n=length(out(1,:));
        k=sqrt(out(1,:).^2+out(2,:).^2);
        k=[k
            k];
        out(1:2,:)=out(1:2,:)./k;
        for i=1:n
            if out(1,i)>=0 && out(2,i)>=0       % 0-90
                outputs(3,i)=asind(out(1,i));
            elseif out(1,i)>=0 && out(2,i)<0    % 90-180
                outputs(3,i)=180-asind(out(1,i));
            elseif out(1,i)<0 && out(2,i)>=0    % 270-360
                outputs(3,i)=360+asind(out(1,i));
            elseif out(1,i)<0 && out(2,i)<0     % 180-270
                outputs(3,i)=180-asind(out(1,i));
            end
        end
    end
end
if handles.theta
    if handles.type==3 % restricted
        out=sim(handles.net_theta,handles.data_noise); 
%         outputs(8,:)=out(1,:);
        n=length(out(1,:));
        k=sqrt(out(1,:).^2+out(2,:).^2);
        out(1,:)=out(1,:)./k;
        out(2,:)=out(2,:)./k;
        for i=1:n
            if out(1,i)>=0 && out(2,i)>=0       % 0-90
                outputs(4,i)=asind(out(1,i))/2;
            elseif out(1,i)>=0 && out(2,i)<0    % 90-180
                outputs(4,i)=(180-asind(out(1,i)))/2;
            end
        end
        for i=1:n
            k(1,i)=min(max(k(1,i),0),1);
        end
        outputs(5,:)=acosd(sqrt(2*k+1/4)-1/2);
    %     k1=ones(1,n)-k;
    %     for i=1:n
    %         k1(1,i)=min(max(k1(1,i),0),1);
    %     end
    %     outputs(7,:)=2*asind(k1);
    elseif handles.type==1  % fixed  
        outputs(4,:)=sim(handles.net_theta,handles.data_noise);   
        outputs(4,outputs(4,:)>1)=1;
        outputs(4,outputs(4,:)<0)=0;
        outputs(4,:)=outputs(4,:)*90;
    end
end
if handles.focus
    outputs(6,:)=sim(handles.net_focus,handles.data_noise); 
end
error=outputs-info;
for i=1:n
    if error(3,i)>180
        error(3,i)=error(3,i)-360;
    end
    if error(3,i)<-180
        error(3,i)=360+error(3,i);
    end
end
save('results.mat','outputs','info','error');
error(8,:)=sqrt(error(1,:).^2+error(2,:).^2);
e=sqrt(sum(error(8,:).^2)/n);
disp(['ANN xy position error (RMSE):',num2str(e)]);
e=sqrt(sum(error(7,:).^2)/n);
disp(['ANN photon error (RMSE):',num2str(e)]);
if handles.phi
    e=sqrt(sum(error(3,:).^2)/n);
    disp(['ANN phi error (RMSE):',num2str(e)]);
end
if handles.theta
    e=sqrt(sum(error(4,:).^2)/n);
    disp(['ANN theta error (RMSE):',num2str(e)]);
    if handles.type==3
        e=sqrt(sum(error(5,:).^2)/n);
        disp(['ANN delta error (RMSE):',num2str(e)]);
%         e=sqrt(sum(error(8,:).^2)/n);
%         disp(['ANN RMSE:',num2str(e)]);
    end
end
if handles.focus
    e=sqrt(sum(error(6,:).^2)/n);
    disp(['ANN z postion error (RMSE):',num2str(e)]);
end

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {'*.mat', 'All MAT-Files (*.mat)'
        '*.*','All Files (*.*)'}, ...
    'Load Data File');
if isequal([filename,pathname],[0,0])
    return
else
    filestr=fullfile(pathname,filename);
end
P=importdata(filestr);
handles.position=0;
handles.photon=0;
handles.phi=0;
handles.theta=0;
handles.focus=0;
handles.type=P.type;
if isfield(P,'net_position');
    handles.position=1;
    handles.net_position=P.net_position;
end
if isfield(P,'net_photon');
    handles.photon=1;
    handles.net_photon=P.net_photon;
end
if isfield(P,'net_phi');
    handles.phi=1;
    handles.ne_phi=P.net_phi;
end
if isfield(P,'net_theta');
    handles.theta=1;
    handles.net_theta=P.net_theta;
end
if isfield(P,'net_focus');
    handles.focus=1;
    handles.net_focus=P.net_focus;
end
P=[];
guidata(hObject,handles) 

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
npix=handles.npix;
n=length(handles.data_noise(1,:));
data_noise=reshape(handles.data_noise,npix,npix,n);
for i=1:n
    data_noise(:,:,i)=data_noise(:,:,i)*handles.maxI(1,i);
end
tiffwrite(uint16(data_noise));

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uiputfile( ...
    {'*.mat', 'All MAT-Files (*.mat)'
        '*.*','All Files (*.*)'}, ...
    'Save Data File');
if isequal([filename,pathname],[0,0])
    return
else
    filestr=fullfile(pathname,filename);
end
data=handles.data;
info=handles.info;
type=handles.type;
save(filestr,'data','info','type');
clear data info type

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.radiobutton3,'value')
    handles.type=1;
elseif get(handles.radiobutton5,'value')
    handles.type=3;
else
    handles.type=2;
end
handles.npix = str2double(get(handles.edit5,'string'));
[handles.data,handles.info]=ANN_getTrainData(handles);
guidata(hObject,handles) 

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uiputfile( ...
    {'*.mat', 'All MAT-Files (*.mat)'
        '*.*','All Files (*.*)'}, ...
    'Save Data File');
if isequal([filename,pathname],[0,0])
    return
else
    filestr=fullfile(pathname,filename);
end
P.type=handles.type;
if isfield(handles,'net_position');
    P.net_position=handles.net_position;
end
if isfield(handles,'net_photon');
    P.net_photon=handles.net_photon;
end
if isfield(handles,'net_phi');
    P.net_phi=handles.net_phi;
end
if isfield(handles,'net_theta');
    P.net_theta=handles.net_theta;
end
if isfield(handles,'net_focus');
    P.net_focus=handles.net_focus;
end
save(filestr,'P');
clear P;


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7


% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9


% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox10


% --- Executes on button press in checkbox11.
function checkbox11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox11


% --- Executes on button press in checkbox12.
function checkbox12_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox12


% --- Executes on button press in checkbox13.
function checkbox13_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox13



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
