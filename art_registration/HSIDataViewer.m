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

% Last Modified by GUIDE v2.5 15-Jul-2015 11:53:52

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

pc = ~ismac;
set(handles.pc,'String',num2str(pc));
ncores = feature('NumCores');
set(handles.ncores,'String',num2str(ncores));
set(handles.axes1,'Units','pixels');
set(handles.axes2,'Units','pixels');

if (pc == 0)
    fnt = '/Applications/RISDataViewer/application/RISDataViewer.txt';
elseif (pc == 1)
    fnt = 'RISDataViewer.txt';
end
vnir_on = get(handles.VNIR,'Value');
xnir_on = get(handles.xNIR,'Value');
nir_on = get(handles.NIR,'Value');
if (xnir_on == 1)
    cam0 = 'xnir';
elseif (vnir_on == 1)
    cam0 = 'vnir';
elseif (nir_on == 1)
    cam0 = 'nir';
end
if (exist(fnt,'file') == 2)
    fid = fopen(fnt);
    while ~feof(fid)
        line = fgetl(fid);
        msk = isspace(line);
        line(msk==1) = '';
        [~,~,e] = regexp(line,['^' cam0 '_a0='],'match','start','end');
        if (e>0)
            a0 = str2num(line((e+1):end));
        end
        [~,~,e] = regexp(line,['^' cam0 '_a1='],'match','start','end');
        if (e>0)
            a1 = str2num(line((e+1):end));
        end
        [~,~,e] = regexp(line,['^' cam0 '_a2='],'match','start','end');
        if (e>0)
            a2 = str2num(line((e+1):end));
        end
        [~,~,e] = regexp(line,['^' cam0 '_a3='],'match','start','end');
        if (e>0)
            a3 = str2num(line((e+1):end));
        end
        [~,~,e] = regexp(line,['^deltan='],'match','start','end');
        if (e>0)
            deltan = str2num(line((e+1):end));
        end
    end
    fclose(fid);
    set(handles.a0,'String',num2str(a0));
    set(handles.a1,'String',num2str(a1));
    set(handles.a2,'String',num2str(a2));
    set(handles.a3,'String',num2str(a3));
    set(handles.deltan,'String',num2str(deltan));
end

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

read_only = get(handles.read_only,'Value');
img_path = get(handles.img_path,'String');
trial = get(handles.trial,'String');
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
if (exist(w_full,'file') == 2)
    w_on = 1;
else
    w_on = 0;
end
if ((b_on == 1) && (w_on == 1))
    ff = 1;
    div_only = 0;
elseif ((w_on == 1) && (b_on == 0))
    ff = 1;
    div_only = 1;
else
    ff = 0;
    div_only = 0;
end
trial_out = img_path;
if (read_only == 1)
    trial_out = [trial_out trial];
elseif (ff == 1) 
    trial_out = [trial_out 'ff_' trial];
elseif (ff == 0)
    trial_out = [trial_out trial '1'];
end   
n0 = str2num(get(handles.n0,'string'));
m0 = str2num(get(handles.m0,'string'));
p0 = str2num(get(handles.p0,'string'));
datatype = get(handles.datatype,'string');
interleave = get(handles.interleave,'string');
byteorder = get(handles.byteorder,'string');
headeroffset = str2num(get(handles.headeroffset,'string'));

p3 = str2num(get(handles.p3,'String')); %#ok<ST2NM>
p2 = str2num(get(handles.p2,'String')); %#ok<ST2NM>
a0 = str2num(get(handles.a0,'String')); %#ok<ST2NM>
a1 = str2num(get(handles.a1,'String')); %#ok<ST2NM>
a2 = str2num(get(handles.a2,'String'))*(10^(-3)); %#ok<ST2NM>
a3 = str2num(get(handles.a3,'String'))*(10^(-6)); %#ok<ST2NM>
%lambda3 = a0 + a1*p3 + a2*p3^2 + a3*p3^3;
%lambda2 = a0 + a1*p2 + a2*p2^2 + a3*p2^3;
%set(handles.lambda3,'String',num2str(lambda3));
%set(handles.lambda2,'String',num2str(lambda2));
set(handles.slider1,'Min',p3); %#%#ok<MSNU> ok<ST2NM>
set(handles.slider1,'Max',p2);
ss = 1/(p2-p3+1);
set(handles.slider1,'SliderStep',[ss 20*ss]);

p1 = round(get(handles.slider1,'Value'));
set(handles.p1,'string',num2str(p1));
lambda1 = a0 + a1*p1 + a2*p1^2 + a3*p1^3;
set(handles.lambda1,'string',num2str(lambda1));

axes(handles.axes1)
cube = multibandread(trial_out,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Band','Direct',p1});

b0 = get(handles.axes1,'Position');
m1 = b0(4);
n1 = b0(3);
cube = imresize(cube,[m1 n1]);
imagesc(cube)
colormap gray
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

%vnir_on = get(handles.VNIR,'Value');
%xnir_on = get(handles.xNIR,'Value');
%nir_on = get(handles.NIR,'Value');
%read_only = get(handles.read_only,'Value');
[fn,path] = uigetfile({'*.*'},'Select the first cube in the set.');
%{
if (read_only == 1)
    [fn,path] = uigetfile({'*.*'},'File Selector');
elseif (vnir_on == 1)
    [fn,path] = uigetfile({'*.cube';'*.img';'*.bsq';'*.bil',;'*.bip'},'File Selector');
elseif (xnir_on == 1)
    [fn,path] = uigetfile({'*.img'},'File Selector');
elseif (nir_on == 1)
    [fn,path] = uigetfile({'*.cube';'*.img';'*.bsq';'*.bil',;'*.bip'},'File Selector');
end
%}
set(handles.img_path,'String',path);
ind = regexp(fn, '\.');
if (~isempty(ind))
    set(handles.img_ext,'String',fn((ind+1):end));
    fn(ind:end) = [];
else
    set(handles.img_ext,'String','');
end
ind = regexp(fn, '_[0-9]$');
if (~isempty(ind))
    set(handles.trial,'String',fn(1:(ind-1)));
    set(handles.num1,'String',fn((ind+1):end));
else
    set(handles.trial,'String',fn);
end

% --- Executes on button press in file_info2.
function file_info2_Callback(hObject, eventdata, handles)
% hObject    handle to file_info2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tmp = get(handles.img_path,'String');
%vnir_on = get(handles.VNIR,'Value');
%xnir_on = get(handles.xNIR,'Value');
%nir_on = get(handles.NIR,'Value');
[fn,path] = uigetfile({'*.*'},'File Selector',tmp);
%{
if (vnir_on == 1)
    [fn,path] = uigetfile({'*.cube';'*.img';'*.bsq';'*.bil',;'*.bip'},'File Selector',tmp);
elseif (xnir_on == 1)
    [fn,path] = uigetfile({'*.img'},'File Selector',tmp);
elseif (nir_on == 1)
    [fn,path] = uigetfile({'*.cube';'*.img';'*.bsq';'*.bil',;'*.bip'},'File Selector',tmp);
end
%}
set(handles.w_path,'String',path);
ind = regexp(fn, '\.');
set(handles.w_ext,'String',fn((ind+1):end));
fn(ind:end) = [];
set(handles.w_fn,'String',fn);

% --- Executes on button press in file_info3.
function file_info3_Callback(hObject, eventdata, handles)
% hObject    handle to file_info3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tmp = get(handles.img_path,'String');
%vnir_on = get(handles.VNIR,'Value');
%xnir_on = get(handles.xNIR,'Value');
%nir_on = get(handles.NIR,'Value');
[fn,path] = uigetfile({'*.*'},'File Selector',tmp);
%{
if (vnir_on == 1)
    [fn,path] = uigetfile('*.dark',tmp);
elseif (xnir_on == 1)
    [fn,path] = uigetfile({'*.img'},'File Selector',tmp);
elseif (nir_on == 1)
    [fn,path] = uigetfile({'*'},'File Selector',tmp);
end
%}
set(handles.b_path,'String',path);
ind = regexp(fn, '\.');
set(handles.b_ext,'String',fn((ind+1):end));
fn(ind:end) = [];
set(handles.b_fn,'String',fn);

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

ncores = str2num(get(handles.ncores,'String'));
deltan = str2num(get(handles.deltan,'String'));
pc = str2num(get(handles.pc,'String'));
p3 = str2num(get(handles.p3,'String'));
p2 = str2num(get(handles.p2,'String'));
a0 = str2num(get(handles.a0,'String'));
a1 = str2num(get(handles.a1,'String'));
a2 = str2num(get(handles.a2,'String'))*(10^(-3));
a3 = str2num(get(handles.a3,'String'))*(10^(-6));
read_only = get(handles.read_only,'Value');
vnir_on = get(handles.VNIR,'Value');
xnir_on = get(handles.xNIR,'Value');
nir_on = get(handles.NIR,'Value');
if (get(handles.bip_button,'value') == 1)
    interleave1 = 'bip';
elseif (get(handles.bil_button,'value') == 1)
    interleave1 = 'bil';
elseif (get(handles.bsq_button,'value') == 1)
    interleave1 = 'bsq';
end

startr = zeros(ncores,1);
startr(1) = 1;
stopr = zeros(ncores,1);
stopr(1) = ceil(deltan/ncores);
for k = 2:ncores
    startr(k) = (k-1)*ceil(deltan/ncores) + 1;
    stopr(k) = k*ceil(deltan/ncores);
end
stopr(ncores) = deltan;
tmp = stopr-startr+1;
rmx = max(tmp);
clear tmp

img_path = get(handles.img_path,'String');
trial = get(handles.trial,'String');
num1 = get(handles.num1,'String');
num2 = get(handles.num2,'String');
if (~isempty(num2))
    num2 = str2num(num2);
elseif (~isempty(num1))
    num2 = str2num(num1);
end
if (~isempty(num1))
    num1 = str2num(num1);
    mult_cubes = 1;
else
    mult_cubes = 0;
    num1 = 1;
    num2 = 1;
end
img_ext = get(handles.img_ext,'String');

if (mult_cubes == 1)
    trial_out = [img_path trial '_' num2str(1)];
elseif (mult_cubes == 0)
    trial_out = [img_path trial];
end
fnh = [trial_out '.hdr'];
fid = fopen(fnh);
while ~feof(fid)
    line = fgetl(fid);
    msk = isspace(line);
    line(msk==1) = '';
    [~,~,e] = regexp(line,'^samples=','match','start','end');
    if (e>0)
        n0 = str2num(line((e+1):end));
    end
    [~,~,e] = regexp(line,'^lines=','match','start','end');
    if (e>0)
        m0 = str2num(line((e+1):end));
    end
    [~,~,e] = regexp(line,'^bands=','match','start','end');
    if (e>0)
        p0 = str2num(line((e+1):end));
    end
    [~,~,e] = regexp(line,'^datatype=','match','start','end');
    if (e>0)
        datatype0 = str2num(line((e+1):end));
        if (datatype0 == 1)
            datatype = 'uint8';
        elseif (datatype0 == 4)
            datatype = 'single';
        elseif (datatype0 == 5)
            datatype = 'double';
        elseif (datatype0 == 12)
            datatype = 'uint16';
        end
    end
    [~,~,e] = regexp(line,'^interleave=','match','start','end');
    if (e>0)
        interleave = line((e+1):end);
    end
    [~,~,e] = regexp(line,'^byteorder=','match','start','end');
    if (e>0)
        byteorder = str2num(line((e+1):end));
        if (byteorder == 0)
            byteorder = 'ieee-le';
        elseif (byteorder == 1)
            byteorder = 'ieee-be';
        end
    end
    [~,~,e] = regexp(line,'^headeroffset=','match','start','end');
    if (e>0)
        headeroffset = str2num(line((e+1):end));
    end
end
fclose(fid);
set(handles.p3,'String',num2str(1));
set(handles.p2,'String',num2str(p0));

if (read_only == 0)
    if (ncores > 1)
        isOpen = matlabpool('size') > 0;
        if (isOpen == 1)
            matlabpool close
        end
        matlabpool(ncores)
    end
    
    useeasel = get(handles.VNIR,'Value');
    
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
    if (exist(w_full,'file') == 2)
        w_on = 1;
    else
        w_on = 0;
    end
    if ((b_on == 1) && (w_on == 1))
        ff = 1;
        div_only = 0;
    elseif ((w_on == 1) && (b_on == 0))
        ff = 1;
        div_only = 1;
    else
        ff = 0;
        div_only = 0;
    end
    set(handles.ff,'Value',ff);
    set(handles.div_only,'Value',div_only);
    
    for i = num1:num2
        if (mult_cubes == 1)
            if (~isempty(img_ext))
                fn_full = [img_path trial '_' num2str(i) '.' img_ext];
            else
                fn_full = [img_path trial '_' num2str(i)];
            end
        elseif (mult_cubes == 0)
            if (~isempty(img_ext))
                fn_full = [img_path trial '.' img_ext];
            else
                fn_full = [img_path trial];
            end
        end
            
        if ((ff == 1) || (div_only == 1))
            if (mult_cubes == 1)
                trial_out = [img_path 'ff_' trial '_' num2str(i)];
            elseif (mult_cubes == 0)
                trial_out = [img_path 'ff_' trial];
            end
        else
            if (mult_cubes == 1)
                trial_out = [img_path trial '1_' num2str(i)];
            elseif (mult_cubes == 0)
                trial_out = [img_path trial '1'];
            end
        end
        fnh = [trial_out '.hdr'];
        
        % flat_field
        if (nir_on == 1)
            d_header = 32768;
            md = 25;
            m0 = 640;
            n0 = 640;
            p0 = 256;
            if (ff == 1)
                w = single(multibandread(w_full,[m0,n0,p0],'uint16',0,'bsq','ieee-le'));
                mu_w = mean(w,1);
                mu_w = ones(m0,n0,p0,'single');
            elseif (ff == 0)
                mu_w = ones(m0,n0,p0,'single');
            end
            if (div_only == 0)
                d = single(multibandread(d_full,[m0,n0,p0],'uint16',0,'bsq','ieee-le'));
                mu_d = mean(d,1);
                mu_d = ones(m0,n0,p0,'single');
            elseif (div_only == 1)
                mu_d = zeros(m0,n0,p0,'single');
            end
            
            for m = 1:deltan:m0
                parfor j = 1:ncores
                    m11 = m+startr(j)-1;
                    m12 = m+stopr(j)-1;
                    if (m12 > m0)
                        m12 = m0;
                    end
                    cube = single(multibandread(fn_full,[m0,n0,p0],'uint16',0,'bsq','ieee-le'),{'Rows','Range',[m11 m12]});
                    dtot2 = mu_d(m11:m12,:,:);
                    wtot2 = mu_w(m11:m12,:,:);
                    cube = (cube - dtot2)./(wtot2 - dtot2);
                    multibandwrite(single(cube),trial_out,'bsq',[m11 1 1],[m0 n0 p0],'machfmt','ieee-le')
                end
            end
            clear cube mu_d mu_w dtot2 wtot2
        elseif (vnir_on == 1)
            if (ff == 1)
                fnh = [wfull '.hdr'];
                fid = fopen(fnh);
                while ~feof(fid)
                    line = fgetl(fid);
                    msk = isspace(line);
                    line(msk==1) = '';
                    [~,~,e] = regexp(line,'^lines=','match','start','end');
                    if (e>0)
                        mw = str2num(line((e+1):end));
                    end
                end
                fclose(fid);
                w = single(multibandread(w_full,[mw,n0,p0],'uint16',0,'bsq','ieee-le'));
                mu_w = mean(w,1);
                mu_w = repmat(mu_w,[m0 1 1]);
            elseif (ff == 0)
                mu_w = ones(m0,n0,p0,'single');
            end
            if (div_only == 0)
                d_header = 32768;
                md = 25;
                d = single(multibandread(d_full,[md,n0,p0],'uint16',d_header,'bsq','ieee-le'));
                mu_d = mean(d,1);
                mu_d = repmat(mu_d,[m0 1 1]);
            elseif (div_only == 1)
                mu_d = zeros(m0,n0,p0,'single');
            end
            
            for m = 1:deltan:m0
                parfor j = 1:ncores
                    m11 = m+startr(j)-1;
                    m12 = m+stopr(j)-1;
                    if (m12 > m0)
                        m12 = m0;
                    end
                    fn_vnir = [img_path trial '_' num2str(j) '.' img_ext];
                    cube = single(multibandread(fn_vnir,[m0,n0,p0],'uint16',0,'bsq','ieee-le'),{'Rows','Range',[m11 m12]});
                    dtot3 = mu_d(m11:m12,:,:);
                    wtot3 = mu_w(m11:m12,:,:);
                    cube = (cube - dtot3)./(wtot3 - dtot3);
                    multibandwrite(single(cube),trial_out,'bsq',[m11 1 1],[m0 n0 p0],'machfmt','ieee-le')
                end
            end
            clear cube mu_d mu_w dtot3 wtot3
        elseif (xnir_on == 1)
            if (ff == 1)
                fidw = fopen(w_full,'r');
                header = fread(fidw,256,'uint16',0,'ieee-le');
                p0 = header(2);
                m0 = header(3);
                nw = header(22);
                header = fread(fidw,5*p0*m0,'uint16',0,'ieee-le');
                clear header
                tmpw = zeros(p0,m0,nw,'single');
                for n = 1:nw
                    tmpw(:,:,n) = single(fread(fidw,[p0,m0],'uint16',0,'ieee-le'));
                end
                fclose(fidw);
                w_tot = mean(tmpw,3);
                clear tmpw
            elseif (ff == 0)
                tmpw = ones(p0,m0,1,'single');
                w_tot = mean(tmpw,3);
                clear tmpw
            end

            if (div_only == 0)
                fidd = fopen(b_full,'r');
                header = fread(fidd,256 + 5*p0*m0,'uint16',0,'ieee-le');
                nd = header(22);
                clear header
                tmpd = zeros(p0,m0,nd,'single');
                for n = 1:nd
                    tmpd(:,:,n) = single(fread(fidd,[p0,m0],'uint16',0,'ieee-le'));
                end
                fclose(fidd);
                d_tot = mean(tmpd,3);
                clear tmpd
            elseif (div_only == 1)
                tmpd = zeros(p0,m0,1,'single');
                d_tot = mean(tmpd,3);
                clear tmpd
            end

            fid = fopen(fn_full,'r');
            header = fread(fid,256 + 5*p0*m0,'uint16',0,'ieee-le');
            n0 = header(22);
            clear header
            cube = zeros(p0,m0,deltan,'single');
            cube1 = zeros(p0,m0,rmx,ncores,'single');

            n1 = 0;
            for n = 1:n0
                n1 = n1 + 1;
                cube(:,:,n1) = single(fread(fid,[p0,m0],'uint16',0,'ieee-le'));
                if ((n1 == deltan) || (n == n0))
                    for j = 1:ncores
                        cube1(:,:,1:(stopr(j)-startr(j)+1),j) = cube(:,:,startr(j):stopr(j));
                    end
                    parfor j = 1:ncores
                        dtot1 = repmat(d_tot,[1 1 (stopr(j)-startr(j)+1)]);
                        wtot1 = repmat(w_tot,[1 1 (stopr(j)-startr(j)+1)]);
                        cube1(:,:,:,j) = (cube1(:,:,:,j) - dtot1)./(wtot1 - dtot1);
                    end
                    clear dtot1 wtot1
                    for j = 1:nores
                        cube(:,:,startr(j):stopr(j)) = cube1(:,:,1:(stopr(j)-startr(j)+1),j);
                    end
                    clear cube1
                    
                    tmp = shiftdim(cube(:,:,1:n1),1);
                    tmp = tmp(:,:,p0:-1:1);
                    multibandwrite(single(tmp),trial_out,interleave1,[1 (n-n1+1) 1],[m0 n0 p0],'machfmt','ieee-le')
                    clear tmp
                    n1 = 0;
                end
            end
            fclose(fid);
            clear cube w_tot d_tot
        end
        fid1 = fopen(fnh,'w');
        fprintf(fid1,'ENVI\n');
        fprintf(fid1,'description = {}\n');
        fprintf(fid1,'samples = %u\n',n0);
        fprintf(fid1,'lines = %u\n',m0);
        fprintf(fid1,'bands = %u\n',p0);
        fprintf(fid1,'header offset = 0\n');
        fprintf(fid1,'file type = ENVI Standard\n');
        fprintf(fid1,'data type = 4\n');
        fprintf(fid1,['interleave = ' interleave1 '\n']);
        fprintf(fid1,'byte order = 0\n');
        fprintf(fid1,'Wavelength = {');
        for p = 1:(p0-1)
            lambda1 = a0 + a1*p + a2*p^2 + a3*p^3;
            fprintf(fid1,'%f, ',lambda1);
        end
        fprintf(fid1,'%f}\n',a0 + a1*p0 + a2*p0^2 + a3*p0^3);
        fclose(fid1);
    end
    
    if (ncores > 1)
        isOpen = matlabpool('size') > 0;
        if (isOpen == 1)
            matlabpool close
        end
    end
end

set(handles.m0,'String',m0);
set(handles.n0,'String',n0);
set(handles.p0,'String',p0);
set(handles.datatype,'String',datatype);
set(handles.byteorder,'String',byteorder);
set(handles.interleave,'String',interleave);
set(handles.headeroffset,'String',headeroffset);
lambda3 = a0 + a1*p3 + a2*p3^2 + a3*p3^3;
lambda2 = a0 + a1*p2 + a2*p2^2 + a3*p2^3;
set(handles.lambda3,'String',num2str(lambda3));
set(handles.lambda2,'String',num2str(lambda2));
set(handles.slider1,'Enable','off');
set(handles.slider1,'Visible','on');
set(handles.slider1,'Min',p3);
set(handles.slider1,'Max',p2);
ss = 1/(p2-p3+1);
set(handles.slider1,'Value',p3);
set(handles.slider1,'SliderStep',[ss 20*ss]);
set(handles.lambda1,'Visible','on');
set(handles.lambda1,'string',lambda3);
set(handles.p1,'Visible','off');
set(handles.p1,'string',p3);

axes(handles.axes1)
tmp= multibandread(trial_out,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Band','Direct',p3});

b0 = get(handles.axes1,'Position');
sc1 = b0(3)/n0;
sc2 = b0(4)/m0;
if (sc2 < sc1)
    sc = sc2;
else
    sc = sc1;
end
n1 = floor(n0*sc);
m1 = floor(m0*sc);
tmp = imresize(tmp,[m1 n1]);
b0(3) = n1;
b0(4) = m1;
set(handles.axes1,'Position',b0)

imagesc(tmp)
colormap gray
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

pc = str2num(get(handles.pc,'String'));
if (pc == 0)
    fnt = '/Applications/RISDataViewer/application/RISDataViewer.txt';
elseif (pc == 1)
    fnt = 'RISDataViewer.txt';
end

if (exist(fnt,'file') == 2)
    fid = fopen(fnt);
    while ~feof(fid)
        line = fgetl(fid);
        msk = isspace(line);
        line(msk==1) = '';
        [~,~,e] = regexp(line,'^xnir_a0=','match','start','end');
        if (e>0)
            xnir_a0 = str2num(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^xnir_a1=','match','start','end');
        if (e>0)
            xnir_a1 = str2num(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^xnir_a2=','match','start','end');
        if (e>0)
            xnir_a2 = str2num(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^xnir_a3=','match','start','end');
        if (e>0)
            xnir_a3 = str2num(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^vnir_a0=','match','start','end');
        if (e>0)
            vnir_a0 = str2num(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^vnir_a1=','match','start','end');
        if (e>0)
            vnir_a1 = str2num(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^vnir_a2=','match','start','end');
        if (e>0)
            vnir_a2 = str2num(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^vnir_a3=','match','start','end');
        if (e>0)
            vnir_a3 = str2num(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^nir_a0=','match','start','end');
        if (e>0)
            nir_a0 = str2num(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^nir_a1=','match','start','end');
        if (e>0)
            nir_a1 = str2num(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^nir_a2=','match','start','end');
        if (e>0)
            nir_a2 = str2num(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^nir_a3=','match','start','end');
        if (e>0)
            nir_a3 = str2num(line((e+1):end));
        end
    end
    fclose(fid);
end
a0 = str2num(get(handles.a0,'String'));
a1 = str2num(get(handles.a1,'String'));
a2 = str2num(get(handles.a2,'String'));
a3 = str2num(get(handles.a3,'String'));

if (pc == 0)
    fn = '/Applications/XRFMapping/application/RISDataViewer.txt';
elseif (pc == 1)
    fn = 'RISDataViewer.txt';
end
fid = fopen(fn,'w+t','n');

if (xnir_on == 1)
    fprintf(fid,'xnir_a0 = %s\n',num2str(a0));
    fprintf(fid,'xnir_a1 = %s\n',num2str(a1));
    fprintf(fid,'xnir_a2 = %s\n',num2str(a2));
    fprintf(fid,'xnir_a3 = %s\n',num2str(a3));
else
    fprintf(fid,'xnir_a0 = %s\n',num2str(xnir_a0));
    fprintf(fid,'xnir_a1 = %s\n',num2str(xnir_a1));
    fprintf(fid,'xnir_a2 = %s\n',num2str(xnir_a2));
    fprintf(fid,'xnir_a3 = %s\n',num2str(xnir_a3));
end
if (vnir_on == 1)
    fprintf(fid,'vnir_a0 = %s\n',num2str(a0));
    fprintf(fid,'vnir_a1 = %s\n',num2str(a1));
    fprintf(fid,'vnir_a2 = %s\n',num2str(a2));
    fprintf(fid,'vnir_a3 = %s\n',num2str(a3));
else
    fprintf(fid,'vnir_a0 = %s\n',num2str(vnir_a0));
    fprintf(fid,'vnir_a1 = %s\n',num2str(vnir_a1));
    fprintf(fid,'vnir_a2 = %s\n',num2str(vnir_a2));
    fprintf(fid,'vnir_a3 = %s\n',num2str(vnir_a3));
end
if (nir_on == 1)
    fprintf(fid,'nir_a0 = %s\n',num2str(a0));
    fprintf(fid,'nir_a1 = %s\n',num2str(a1));
    fprintf(fid,'nir_a2 = %s\n',num2str(a2));
    fprintf(fid,'nir_a3 = %s\n',num2str(a3));
else
    fprintf(fid,'nir_a0 = %s\n',num2str(nir_a0));
    fprintf(fid,'nir_a1 = %s\n',num2str(nir_a1));
    fprintf(fid,'nir_a2 = %s\n',num2str(nir_a2));
    fprintf(fid,'nir_a3 = %s\n',num2str(nir_a3));
end
fprintf(fid,'deltan = %s\n',num2str(deltan));
fclose(fid);

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

handles = guidata(RISDataViewer);
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

    n0 = str2num(get(handles.n0,'String'));
    m0 = str2num(get(handles.m0,'String'));
    datatype = get(handles.datatype,'String');
    interleave = get(handles.interleave,'String');
    byteorder = get(handles.byteorder,'String');
    headeroffset = str2num(get(handles.headeroffset,'String'));
    p2 = str2num(get(handles.p2,'String'));
    p3 = str2num(get(handles.p3,'String'));
    coordn = round(n0*x/pos(3));
    coordm = round(m0*y/pos(4));
    set(handles.x1,'String',num2str(coordn));
    set(handles.y1,'String',num2str(coordm));

    a0 = str2num(get(handles.a0,'String'));
    a1 = str2num(get(handles.a1,'String'));
    a2 = str2num(get(handles.a2,'String'))*(10^(-3));
    a3 = str2num(get(handles.a3,'String'))*(10^(-6));
    lambda3 = a0 + a1*p3 + a2*p3^2 + a3*p3^3;
    lambda2 = a0 + a1*p2 + a2*p2^2 + a3*p2^3;
    lambda = a0 + a1*(p3:p2) + a2*((p3:p2).^2) + a3*((p3:p2).^3);
    set(handles.lambda3,'string',num2str(lambda3));
    set(handles.lambda2,'string',num2str(lambda2));

    read_only = get(handles.read_only,'Value');
    ff = get(handles.ff,'Value');
    trial_out = img_path;
    if (read_only == 1)
        trial_out = [trial_out trial];
    else
        if (ff == 1)
            trial_out = [trial_out 'ff_' trial];
        else
            trial_out = [trial_out trial '1'];
        end
    end

    cube = single(multibandread(trial_out,[m0,n0,p2-p3+1],datatype,headeroffset,interleave,byteorder,{'Row','Direct',coordm},{'Column','Direct',coordn},{'Band','Range',[p3,p2]}));
    
    cube(cube<0) = 0;
    cube(isnan(cube)) = 0;
    cube(isinf(cube)) = 0;
    axes(handles.axes2)
    plot(lambda,cube(:));
    ylim([0 max(cube(:))]);
    xlabel('wavelength (nm)');
    ylabel('reflectance');
    xlim([lambda(1) lambda(end)]);
    clear cube
elseif ((x2>=1) && (x2<=pos2(3)) && (y2>=1) && (y2<=pos2(4)))
    set(handles.axes2,'Units','points');
    fig = figure(RISDataViewer);
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
    
    fig = figure(RISDataViewer);
    datacursormode on
    dcm_obj = datacursormode(fig);
    set(dcm_obj,'UpdateFcn',@myupdatefcn)
else
    set(handles.eid,'UserData',0);
    set(handles.eid,'String','Spatial pts');
    
    fig = figure(RISDataViewer);
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


% --- Executes on button press in parcomp.
function parcomp_Callback(hObject, eventdata, handles)
% hObject    handle to parcomp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of parcomp


function n0_Callback(hObject, eventdata, handles)
% hObject    handle to n0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n0 as text
%        str2num(get(hObject,'String')) returns contents of n0 as a double


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



function datatype_Callback(hObject, eventdata, handles)
% hObject    handle to datatype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of datatype as text
%        str2num(get(hObject,'String')) returns contents of datatype as a double


% --- Executes during object creation, after setting all properties.
function datatype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datatype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function byteorder_Callback(hObject, eventdata, handles)
% hObject    handle to byteorder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of byteorder as text
%        str2num(get(hObject,'String')) returns contents of byteorder as a double


% --- Executes during object creation, after setting all properties.
function byteorder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to byteorder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function interleave_Callback(hObject, eventdata, handles)
% hObject    handle to interleave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of interleave as text
%        str2num(get(hObject,'String')) returns contents of interleave as a double


% --- Executes during object creation, after setting all properties.
function interleave_CreateFcn(hObject, eventdata, handles)
% hObject    handle to interleave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function headeroffset_Callback(hObject, eventdata, handles)
% hObject    handle to headeroffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of headeroffset as text
%        str2num(get(hObject,'String')) returns contents of headeroffset as a double


% --- Executes during object creation, after setting all properties.
function headeroffset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to headeroffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function ff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ncores_Callback(hObject, eventdata, handles)
% hObject    handle to ncores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ncores as text
%        str2num(get(hObject,'String')) returns contents of ncores as a double


% --- Executes during object creation, after setting all properties.
function ncores_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ncores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in div_only.
function div_only_Callback(hObject, eventdata, handles)
% hObject    handle to div_only (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of div_only



function num2_Callback(hObject, eventdata, handles)
% hObject    handle to num2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num2 as text
%        str2num(get(hObject,'String')) returns contents of num2 as a double


% --- Executes during object creation, after setting all properties.
function num2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num1_Callback(hObject, eventdata, handles)
% hObject    handle to num1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num1 as text
%        str2num(get(hObject,'String')) returns contents of num1 as a double


% --- Executes during object creation, after setting all properties.
function num1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function deltan_Callback(hObject, eventdata, handles)
% hObject    handle to deltan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of deltan as text
%        str2num(get(hObject,'String')) returns contents of deltan as a double


% --- Executes during object creation, after setting all properties.
function deltan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to deltan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in update_coeffs.
function update_coeffs_Callback(hObject, eventdata, handles)
% hObject    handle to update_coeffs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pc = str2num(get(handles.pc,'String'));

if (pc == 0)
    fnt = '/Applications/RISDataViewer/application/RISDataViewer.txt';
elseif (pc == 1)
    fnt = 'RISDataViewer.txt';
end
vnir_on = get(handles.VNIR,'Value');
xnir_on = get(handles.xNIR,'Value');
nir_on = get(handles.NIR,'Value');
if (xnir_on == 1)
    cam0 = 'xnir';
elseif (vnir_on == 1)
    cam0 = 'vnir';
elseif (nir_on == 1)
    cam0 = 'nir';
end
if (exist(fnt,'file') == 2)
    fid = fopen(fnt);
    while ~feof(fid)
        line = fgetl(fid);
        msk = isspace(line);
        line(msk==1) = '';
        [~,~,e] = regexp(line,['^' cam0 '_a0='],'match','start','end');
        if (e>0)
            a0 = line((e+1):end);
        end
        [~,~,e] = regexp(line,['^' cam0 '_a1='],'match','start','end');
        if (e>0)
            a1 = line((e+1):end);
        end
        [~,~,e] = regexp(line,['^' cam0 '_a2='],'match','start','end');
        if (e>0)
            a2 = line((e+1):end);
        end
        [~,~,e] = regexp(line,['^' cam0 '_a3='],'match','start','end');
        if (e>0)
            a3 = line((e+1):end);
        end
    end
    fclose(fid);
    set(handles.a0,'String',a0);
    set(handles.a1,'String',a1);
    set(handles.a2,'String',a2);
    set(handles.a3,'String',a3);
end
