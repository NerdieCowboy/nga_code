function varargout = HSIDataViewer(varargin)
% HSIDATAVIEWER MATLAB code for HSIDataViewer.fig
%      HSIDATAVIEWER, by itself, creates a new HSIDATAVIEWER or raises the existing
%      singleton*.
%
%      H = HSIDATAVIEWER returns the handle to a new HSIDATAVIEWER or the handle to
%      the existing singleton*.
%
%      HSIDATAVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HSIDATAVIEWER.M with the given input arguments.
%
%      HSIDATAVIEWER('Property','Value',...) creates a new HSIDATAVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HSIDataViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HSIDataViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HSIDataViewer

% Last Modified by GUIDE v2.5 17-Nov-2013 19:27:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HSIDataViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @HSIDataViewer_OutputFcn, ...
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


% --- Executes just before HSIDataViewer is made visible.
function HSIDataViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HSIDataViewer (see VARARGIN)

% Choose default command line output for HSIDataViewer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes HSIDataViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = HSIDataViewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

p3 = str2num(get(handles.p3,'String')); %#ok<ST2NM>
p2 = str2num(get(handles.p2,'String')); %#ok<ST2NM>
a0 = str2num(get(handles.a0,'String')); %#ok<ST2NM>
a1 = str2num(get(handles.a1,'String')); %#ok<ST2NM>
a2 = str2num(get(handles.a2,'String'))*(10^(-3)); %#ok<ST2NM>
a3 = str2num(get(handles.a3,'String'))*(10^(-6)); %#ok<ST2NM>
lambda3 = a0 + a1*p3 + a2*p3^2 + a3*p3^3;
lambda2 = a0 + a1*p2 + a2*p2^2 + a3*p2^3;
set(handles.lambda3,'String',num2str(lambda3));
set(handles.lambda2,'String',num2str(lambda2));
set(handles.slider1,'Min',p3); %#%#ok<MSNU> ok<ST2NM>
set(handles.slider1,'Max',p2);
ss = 1/(p2-p3+1);
set(handles.slider1,'SliderStep',[ss 20*ss]);

p1 = round(get(handles.slider1,'Value'));
set(handles.p1,'string',num2str(p1));
lambda1 = a0 + a1*p1 + a2*p1^2 + a3*p1^3;
set(handles.lambda1,'string',num2str(lambda1));

n0 = str2num(get(handles.n0,'String')); %#ok<ST2NM>
m0 = str2num(get(handles.m0,'String')); %#ok<ST2NM>
p0 = str2num(get(handles.p0,'String')); %#ok<ST2NM>
rot = get(handles.rot,'Value');
ff = get(handles.ff,'Value');

img_path = get(handles.img_path,'String');
trial = get(handles.trial,'String');
trial_out = img_path;
if (rot == 1)
    trial_out = [trial_out 'rot_'];
end
if (ff == 1)
    trial_out = [trial_out 'ff_'];
end
if ((ff==1) || (rot==1))
    trial_out = [trial_out trial];
else
    trial_out = [trial_out trial '1'];
end  
%cube = single(multibandread(trial_out,[m0,n0,p0],'uint16',0,'bsq','ieee-le',{'Band','Direct',p1}));
cube = multibandread(trial_out,[m0,n0,p0],'uint16',0,'bsq','ieee-le',{'Band','Direct',p1});

axes(handles.axes1)
%imagesc(single(cube))
imagesc(cube)
colormap gray
%axis image
axis off
clear cube


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function x0_Callback(hObject, eventdata, handles)
% hObject    handle to x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function x0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y0_Callback(hObject, eventdata, handles)
% hObject    handle to y0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function y0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function n0_Callback(hObject, eventdata, handles)
% hObject    handle to n0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function n0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function header_sz_Callback(hObject, eventdata, handles)
% hObject    handle to header_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function header_sz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to header_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rot.
function rot_Callback(hObject, eventdata, handles)
% hObject    handle to rot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in ff.
function ff_Callback(hObject, eventdata, handles)
% hObject    handle to ff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in file_info.
function file_info_Callback(hObject, eventdata, handles)
% hObject    handle to file_info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fn,path] = uigetfile('*.img','Select the cube file.');
set(handles.img_path,'String',path);
ind = regexp(fn, '\.');
set(handles.img_ext,'String',fn((ind+1):end));
fn(ind:end) = [];
set(handles.trial,'String',fn);


% --- Executes on button press in file_info2.
function file_info2_Callback(hObject, eventdata, handles)
% hObject    handle to file_info2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tmp = get(handles.img_path,'String');
[fn,path] = uigetfile('*.img','Select the white cube file.',tmp);
set(handles.w_path,'String',path);
ind = regexp(fn, '\.');
set(handles.w_ext,'String',fn((ind+1):end));
fn(ind:end) = [];
set(handles.w_fn,'String',fn);

function img_path_Callback(hObject, eventdata, handles)
% hObject    handle to img_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function img_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to img_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trial_Callback(hObject, eventdata, handles)
% hObject    handle to trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function trial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function img_ext_Callback(hObject, eventdata, handles)
% hObject    handle to img_ext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function img_ext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to img_ext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function w_path_Callback(hObject, eventdata, handles)
% hObject    handle to w_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function w_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to w_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function w_fn_Callback(hObject, eventdata, handles)
% hObject    handle to w_fn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function w_fn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to w_fn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function w_ext_Callback(hObject, eventdata, handles)
% hObject    handle to w_ext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function w_ext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to w_ext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run_button.
function run_button_Callback(hObject, eventdata, handles)
% hObject    handle to run_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

p3 = str2num(get(handles.p3,'String')); %#ok<ST2NM>
p2 = str2num(get(handles.p2,'String')); %#ok<ST2NM>
a0 = str2num(get(handles.a0,'String')); %#ok<ST2NM>
a1 = str2num(get(handles.a1,'String')); %#ok<ST2NM>
a2 = str2num(get(handles.a2,'String'))*(10^(-3)); %#ok<ST2NM>
a3 = str2num(get(handles.a3,'String'))*(10^(-6)); %#ok<ST2NM>
lambda3 = a0 + a1*p3 + a2*p3^2 + a3*p3^3;
lambda2 = a0 + a1*p2 + a2*p2^2 + a3*p2^3;
set(handles.lambda3,'String',num2str(lambda3));
set(handles.lambda2,'String',num2str(lambda2));
set(handles.slider1,'Enable','off');
set(handles.slider1,'Visible','on');
set(handles.slider1,'Min',p3); %#%#ok<MSNU> ok<ST2NM>
set(handles.slider1,'Max',p2);
ss = 1/(p2-p3+1);
set(handles.slider1,'Value',p3);
set(handles.slider1,'SliderStep',[ss 20*ss]);
set(handles.lambda1,'Visible','on');
set(handles.lambda1,'string',lambda3);
set(handles.p1,'Visible','off');
set(handles.p1,'string',p3);
rot = get(handles.rot,'Value');
useeasel = get(handles.useeasel,'Value');
read_only = get(handles.read_only,'Value');
ff = get(handles.ff,'Value');

img_path = get(handles.img_path,'String');
trial = get(handles.trial,'String');
img_ext = get(handles.img_ext,'String');

if (ff == 1)
    w_path = get(handles.w_path,'String');
    w_fn = get(handles.w_fn,'String');
    w_ext = get(handles.w_ext,'String');
    w_full = [w_path w_fn '.' w_ext];

    b_path = get(handles.b_path,'String');
    b_fn = get(handles.b_fn,'String');
    b_ext = get(handles.b_ext,'String');
    b_full = [b_path b_fn '.' b_ext];
    if (exist(b_full,'file') == 2)
        b_on = 1;
    else
        b_on = 0;
    end
end

trial_out = img_path;
if (rot == 1)
    trial_out = [trial_out 'rot_'];
end
if (ff == 1)
    trial_out = [trial_out 'ff_'];
end
if ((ff==1) || (rot == 1))
    trial_out = [trial_out trial];
else
    trial_out = [trial_out trial '1'];
end
fnh = [trial_out '.hdr'];

if (read_only == 1)
    %{
    fid = fopen(fnh);
    while ~feof(fid)
        line = fgetl(fid);
        msk = isspace(line);
        line(msk==1) = '';
        [~,~,e] = regexp(line,'^samples=','match','start','end');
        if (e>0)
            n0 = str2double(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^lines=','match','start','end');
        if (e>0)
            m0 = str2double(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^bands=','match','start','end');
        if (e>0)
            p0 = str2double(line((e+1):end));
        end
    end
    fclose(fid);
    set(handles.n0,'string',num2str(n0));
    set(handles.m0,'string',num2str(m0));
    set(handles.p0,'string',num2str(p0));
    %}
elseif (read_only == 0)
    fn_full = [img_path trial '.' img_ext];
    
    fid = fopen(fn_full,'r');
    tmp = fread(fid,1,'uint16',0,'ieee-le');
    y0 = fread(fid,1,'uint16',0,'ieee-le');
    set(handles.y0,'String',num2str(y0));
    x0 = fread(fid,1,'uint16',0,'ieee-le');
    set(handles.x0,'String',num2str(x0));
    tmp = fread(fid,18,'uint16',0,'ieee-le');
    n0 = fread(fid,1,'uint16',0,'ieee-le');
    set(handles.n0,'String',num2str(n0));
    %header_sz = str2num(get(handles.header_sz,'String')); %#ok<ST2NM>
    header_sz = y0*x0*10 + 512;
    set(handles.header_sz,'String',num2str(header_sz));
    header = fread(fid,(header_sz-44));
    clear tmp header
    
    if (ff == 1)
       fidw = fopen(w_full,'r');
       header = fread(fidw,header_sz);  %header
       clear header
       
       if (b_on == 1)
           fidb = fopen(b_full,'r');
           header = fread(fidb,header_sz);  %header
           clear header
       end
    end
    
    if (rot == 1)
        %{
        if (ff == 1)
            tmp2 = fread(fidw,[x0,y0],'uint16',0,'ieee-le');
            tmp2 = tmp2';
            tmp2 = single(tmp2);
            fclose(fidw);
            
            if (b_on == 1)
                tmp3 = fread(fidb,[x0,y0],'uint16',0,'ieee-le');
                tmp3 = tmp3';
                tmp3 = single(tmp3);
                fclose(fidb);
            end
        end
        %}
        %{
        n1 = 0;
        cube = uint16(zeros(y0,x0,n0));
        while ~feof(fid)
            tmp = fread(fid,[x0,y0],'uint16',0,'ieee-le');
            if (ff == 1)
                tmp2 = fread(fidw,[x0,y0],'uint16',0,'ieee-le');
                if (b_on == 1)
                    tmp3 = fread(fidb,[x0,y0],'uint16',0,'ieee-le');
                end
            end
            if (~isempty(tmp))
                n1 = n1 + 1;
                tmp = tmp';
                tmp = single(tmp);
                if (ff == 1)
                    tmp2 = tmp2';
                    tmp2 = single(tmp2);
                    if (b_on == 1)
                        tmp3 = tmp3';
                        tmp3 = single(tmp3);
                    end
                end

                if (ff == 1)
                    if (b_on == 0)
                        cube(:,:,n1) = uint16((tmp./tmp2)*(2^14));
                    elseif (b_on == 1)
                        cube(:,:,n1) = uint16(((tmp-tmp3)./(tmp2-tmp3))*(2^14));
                    end
                else
                    cube(:,:,n1) = uint16(tmp);
                end
                clear tmp
            end
        end
        fclose(fid);
        if (ff == 1)
            fclose(fidw);
            clear tmp
            if (b_on == 1)
                fclose(fidb);
                clear tmp3
            end
        end
        %}
    elseif (rot == 0)
        %{
        if (ff == 1)
            tmp2 = fread(fidw,[y0,x0],'uint16',0,'ieee-le');
            tmp2 = single(tmp2);
            fclose(fidw);
        end
        %}
        
        if (useeasel == 1)
            if (ff == 1)
                tmp2 = fread(fidw,[y0,x0],'uint16',0,'ieee-le');
                if (b_on == 1)
                    tmp3 = fread(fidb,[y0,x0],'uint16',0,'ieee-le');
                end
            end
        end
        
        n1 = 0;
        cube = uint16(zeros(y0,x0,n0));
        while ~feof(fid)
            tmp = fread(fid,[y0,x0],'uint16',0,'ieee-le');
            if (useeasel == 0)
                if (ff == 1)
                    tmp2 = fread(fidw,[y0,x0],'uint16',0,'ieee-le');
                    if (b_on == 1)
                        tmp3 = fread(fidb,[y0,x0],'uint16',0,'ieee-le');
                    end
                end
            end
            if (~isempty(tmp))
                n1 = n1 + 1;
                tmp = single(tmp);
                if (useeasel == 0)
                    if (ff == 1)
                        tmp2 = single(tmp2);
                        if (b_on == 1)
                            tmp3 = single(tmp3);
                        end
                    end
                end

                if (ff == 1)
                    if (b_on == 0)
                        cube(:,:,n1) = uint16((tmp./tmp2)*(2^14));
                    elseif (b_on == 1)
                        cube(:,:,n1) = uint16(((tmp-tmp3)./(tmp2-tmp3))*(2^14));
                    end
                else
                    cube(:,:,n1) = uint16(tmp);
                end
                clear tmp
            end
        end
        fclose(fid);
        if (ff == 1)
            fclose(fidw);
            clear tmp2
            if (b_on == 1)
                fclose(fidb);
                clear tmp3
            end
        end
    end

    for n = 1:256:n0
        nend = n+256-1;
        if (nend > n0)
            nend = n0;
        end
        tmp = cube(:,:,n:nend);
        tmp = shiftdim(tmp,1);
        [m0,~,p0] = size(tmp);
        set(handles.m0,'string',num2str(m0));
        set(handles.p0,'string',num2str(p0));
        tmp = tmp(:,:,p0:-1:1);    
        multibandwrite(uint16(tmp),trial_out,'bsq',[1 n 1],[m0 n0 p0],'machfmt','ieee-le')
        clear tmp
    end
    
    fid1 = fopen(fnh,'w');
    fprintf(fid1,'ENVI\n');
    fprintf(fid1,'description = {}\n');
    fprintf(fid1,'samples = %u\n',n0);
    fprintf(fid1,'lines = %u\n',m0);
    fprintf(fid1,'bands = %u\n',p0);
    fprintf(fid1,'header offset = 0\n');
    fprintf(fid1,'file type = ENVI Standard\n');
    fprintf(fid1,'data type = 12\n');
    fprintf(fid1,'interleave = bsq\n');
    fprintf(fid1,'byte order = 0\n');
    fprintf(fid1,'Wavelength = {');
    for i = 1:(p0-1)
        lambda1 = a0 + a1*i + a2*i^2 + a3*i^3;
        fprintf(fid1,'%f, ',lambda1);
    end
    fprintf(fid1,'%f}\n',a0 + a1*p0 + a2*p0^2 + a3*p0^3);
    fclose(fid1);
end
clear cube

axes(handles.axes1)
p1 = str2num(get(handles.p1,'string')); %#ok<ST2NM>
lambda1 = a0 + a1*p1 + a2*p1^2 + a3*p1^3;
set(handles.lambda1,'string',num2str(lambda1));
%tmp= single(multibandread(trial_out,[m0,n0,p0],'uint16',0,'bsq','ieee-le',{'Band','Direct',p1}));
tmp= multibandread(trial_out,[m0,n0,p0],'uint16',0,'bsq','ieee-le',{'Band','Direct',p1});
%imagesc(single(tmp))
imagesc(tmp)
colormap gray
%axis image
axis off
clear tmp

lambda = a0 + a1*(p3:p2) + a2*((p3:p2).^2) + a3*((p3:p2).^3);
axes(handles.axes2)
plot(lambda,zeros(1,(p2-p3+1)))
xlabel('wavelength (nm)');
ylabel('reflectance');
xlim([lambda(1) lambda(end)]);
ylim([0 1]);

set(handles.axes1,'Visible','on');
set(handles.axes2,'Visible','on');

set(handles.axes1,'Units','pixels');
set(handles.axes2,'Units','points');

set(gcf,'WindowButtonDownFcn',@mouseMove);
set(gcf,'Units','pixels');
datacursormode off

set(handles.slider1,'Enable','on');

function p1_Callback(hObject, eventdata, handles)
% hObject    handle to p1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function p1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in useeasel.
function useeasel_Callback(hObject, eventdata, handles)
% hObject    handle to useeasel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




function a3_Callback(hObject, eventdata, handles)
% hObject    handle to a3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function a3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a2_Callback(hObject, eventdata, handles)
% hObject    handle to a2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function a2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a1_Callback(hObject, eventdata, handles)
% hObject    handle to a1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function a1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a0_Callback(hObject, eventdata, handles)
% hObject    handle to a0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function a0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lambda2_Callback(hObject, eventdata, handles)
% hObject    handle to lambda2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function lambda2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambda2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lambda3_Callback(hObject, eventdata, handles)
% hObject    handle to lambda3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function lambda3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambda3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p2_Callback(hObject, eventdata, handles)
% hObject    handle to p2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function p2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p3_Callback(hObject, eventdata, handles)
% hObject    handle to p3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function p3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in read_only.
function read_only_Callback(hObject, eventdata, handles)
% hObject    handle to read_only (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function mouseMove (hobject, eventdata)

handles = guidata(HSIDataViewer);
pos = get(handles.axes1,'Position');
pos2 = get(handles.axes2,'Position');
C = get(gcf,'CurrentPoint');
x = C(1,1) - pos(1) + 1;
y = pos(4) - C(1,2) + pos(2);
x2 = C(1,1) - pos2(1) + 1;
y2 = pos(4) - C(1,2) + pos2(2);
if ((x>=1) && (x<=pos(3)) && (y>=1) && (y<=pos(4)))
    datacursormode off
    img_path = get(handles.img_path,'String');
    trial = get(handles.trial,'String');

    n0 = str2num(get(handles.n0,'String')); %#ok<ST2NM>
    m0 = str2num(get(handles.m0,'String')); %#ok<ST2NM>
    p0 = str2num(get(handles.p0,'String')); %#ok<ST2NM>
    p2 = str2num(get(handles.p2,'String')); %#ok<ST2NM>
    p3 = str2num(get(handles.p3,'String')); %#ok<ST2NM>
    coordn = round(n0*x/pos(3));
    coordm = round(m0*y/pos(4));
    set(handles.x1,'String',num2str(coordn));
    set(handles.y1,'String',num2str(coordm));

    a0 = str2num(get(handles.a0,'String')); %#ok<ST2NM>
    a1 = str2num(get(handles.a1,'String')); %#ok<ST2NM>
    a2 = str2num(get(handles.a2,'String'))*(10^(-3)); %#ok<ST2NM>
    a3 = str2num(get(handles.a3,'String'))*(10^(-6)); %#ok<ST2NM>
    lambda3 = a0 + a1*p3 + a2*p3^2 + a3*p3^3;
    lambda2 = a0 + a1*p2 + a2*p2^2 + a3*p2^3;
    lambda = a0 + a1*(p3:p2) + a2*((p3:p2).^2) + a3*((p3:p2).^3);
    set(handles.lambda3,'string',num2str(lambda3));
    set(handles.lambda2,'string',num2str(lambda2));

    ff = get(handles.ff,'Value');
    rot = get(handles.rot,'Value');
    trial_out = img_path;
    if (rot == 1)
        trial_out = [trial_out 'rot_'];
    end
    if (ff == 1)
        trial_out = [trial_out 'ff_'];
    end
    if ((ff==1) || (rot == 1))
        trial_out = [trial_out trial];
    else
        trial_out = [trial_out trial '1'];
    end

    cube = single(multibandread(trial_out,[m0,n0,p0],'uint16',0,'bsq','ieee-le',{'Row','Range',[y,y]},{'Column','Range',[x,x]},{'Band','Range',[p3,p2]}))/(2^14);
    axes(handles.axes2)
    plot(lambda,cube(:));
    ylim([0 max(cube(:))]);
    xlabel('wavelength (nm)');
    ylabel('reflectance');
    xlim([lambda(1) lambda(end)]);
elseif ((x2>=1) && (x2<=pos2(3)) && (y2>=1) && (y2<=pos2(4)))
    set(handles.axes2,'Units','points');
    fig = figure(HSIDataViewer);
    %datacursormode on
    %dcm_obj = datacursormode(fig);
    %set(dcm_obj,'UpdateFcn',@myupdatefcn)
%else
    %datacursormode(HSIDataViewer)
end

function txt = myupdatefcn(empt,event_obj)
% Customizes text of data tips

ptpos = get(event_obj,'Position');
txt = {['Wavelength: ',num2str(ptpos(1)),' nm'],...
       ['Reflectance: ',num2str(ptpos(2))]};

function eid_Callback(hObject, eventdata, handles)
% hObject    handle to eid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of eid

eid_state = get(handles.eid,'UserData');
if (eid_state == 0)
    set(handles.eid,'UserData',1);
    set(handles.eid,'String','Wavelength');
    
    fig = figure(HSIDataViewer);
    datacursormode on
    dcm_obj = datacursormode(fig);
    set(dcm_obj,'UpdateFcn',@myupdatefcn)
else
    set(handles.eid,'UserData',0);
    set(handles.eid,'String','Spatial pts');
    
    fig = figure(HSIDataViewer);
    datacursormode off
end



function p0_Callback(hObject, eventdata, handles)
% hObject    handle to p0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function p0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function m0_Callback(hObject, eventdata, handles)
% hObject    handle to m0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function m0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to m0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in file_info3.
function file_info3_Callback(hObject, eventdata, handles)
% hObject    handle to file_info3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tmp = get(handles.b_path,'String');
[fn,path] = uigetfile('*.img','Select the dark cube file.',tmp);
set(handles.b_path,'String',path);
ind = regexp(fn, '\.');
set(handles.b_ext,'String',fn((ind+1):end));
fn(ind:end) = [];
set(handles.b_fn,'String',fn);

function b_path_Callback(hObject, eventdata, handles)
% hObject    handle to b_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function b_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to b_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function b_fn_Callback(hObject, eventdata, handles)
% hObject    handle to b_fn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function b_fn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to b_fn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function b_ext_Callback(hObject, eventdata, handles)
% hObject    handle to b_ext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function b_ext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to b_ext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
