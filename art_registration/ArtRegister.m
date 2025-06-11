function varargout = ArtRegister(varargin)
% ARTREGISTER M-file for ArtRegister.fig
%      ARTREGISTER, by itself, creates a new ARTREGISTER or raises the existing
%      singleton*.
%
%      H = ARTREGISTER returns the handle to a new ARTREGISTER or the handle to
%      the existing singleton*.
%
%      ARTREGISTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ARTREGISTER.M with the given input arguments.
%
%      ARTREGISTER('Property','Value',...) creates a new ARTREGISTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ArtRegister_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ArtRegister_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ArtRegister

% Last Modified by GUIDE v2.5 20-Jul-2015 15:10:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ArtRegister_OpeningFcn, ...
                   'gui_OutputFcn',  @ArtRegister_OutputFcn, ...
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

% --- Executes just before ArtRegister is made visible.
function ArtRegister_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ArtRegister (see VARARGIN)

% Choose default command line output for ArtRegister
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

pc = ~ismac;
set(handles.pc,'String',num2str(pc));
ncores = feature('NumCores');
set(handles.ncores,'String',num2str(ncores));

% UIWAIT makes XRFDataViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ArtRegister_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in turnonoutput.
function turnonoutput_Callback(hObject, eventdata, handles)
% hObject    handle to turnonoutput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function Npoly_Callback(hObject, eventdata, handles)
% hObject    handle to Npoly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function Npoly_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Npoly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Ntran_Callback(hObject, eventdata, handles)
% hObject    handle to Ntran (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function Ntran_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ntran (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.stop,'Value',0);
pc = get(handles.pc,'Value');
fn_path = get(handles.fn_path,'String');
trial = get(handles.trial,'String');
N1 = get(handles.N1,'String');
N2 = get(handles.N2,'String');
img_ext = get(handles.img_ext,'String');
if (pc == 1)
    if (~isempty(N1))
        ir_fn = [fn_path '\' trial '_' N1 '.' img_ext];
    else
        ir_fn = [fn_path '\' trial '.' img_ext];
    end
elseif (pc == 0)
    if (~isempty(N1))
        ir_fn = [fn_path '/' trial '_' N1 '.' img_ext];
    else
        ir_fn = [fn_path '/' trial '.' img_ext];
    end
end
rgb_fn = get(handles.rgb_fn,'String');
sc = str2double(get(handles.sc,'String'));
j01 = str2double(get(handles.j01,'String'));
j02 = str2double(get(handles.j02,'String'));
repcnt_max = str2double(get(handles.repcnt_max,'String'));
max_shift = str2double(get(handles.max_shift,'String'));
minpts = str2double(get(handles.minpts,'String'));
Npoly = str2double(get(handles.Npoly,'String'));
Ntran = str2double(get(handles.Ntran,'String'));
memory_limited = get(handles.memory_limited,'Value');
morepts = get(handles.morepts,'Value');
useresults = get(handles.useresults,'Value');
filteron = get(handles.filteron,'Value');
addsharp = get(handles.addsharp,'Value');
debugmode = get(handles.debugmode,'Value');
usegui = 1;
feather_width = str2double(get(handles.feather_width,'String'));
mxd = str2double(get(handles.mxd,'String'));
art_register(rgb_fn,sc,ir_fn,j01,j02,repcnt_max,max_shift,minpts,Npoly,Ntran,memory_limited,morepts,useresults,filteron,addsharp,debugmode,usegui,feather_width,mxd)

if (get(handles.shift_ref,'Value') == 1)
    sc = 1;
    if (pc == 1)
        if (~isempty(N1))
            rgb_fn = [fn_path '\' trial '_' N1 '_IR.tif'];
        else
            rgb_fn = [fn_path '\' trial '_IR.tif'];
        end
    elseif (pc == 0)
        if (~isempty(N1))
            rgb_fn = [fn_path '/' trial '_' N1 '_IR.tif'];
        else
            rgb_fn = [fn_path '/' trial '_IR.tif'];
        end
    end
end
%{
if (get(handles.buildon,'Value') == 1)
    tmp = get(handles.ir_fn,'String');
    ir_last = [tmp(1:(end-4)) '_IR.tif'];
    tmp = get(handles.ir_fn2,'String');
    ir_next = [tmp(1:(end-4)) '_IR.tif'];
    copyfile(ir_last,ir_next);
end
ir_fn2 = get(handles.ir_fn2,'String');
addsharp = get(handles.addsharp2,'Value');
if (exist(ir_fn2,'file') == 2)
    art_register(rgb_fn,sc,ir_fn2,j01,j02,repcnt_max,max_shift,minpts,Npoly,Ntran,memory_limited,morepts,useresults,filteron,addsharp,debugmode,usegui,feather_width,mxd)
end
%}
%{
if (get(handles.shift_ref,'Value') == 1)
    sc = 1;
    tmp = get(handles.ir_fn2,'String');
    rgb_fn = [tmp(1:(end-4)) '_IR.tif'];
end
if (get(handles.buildon,'Value') == 1)
    tmp = get(handles.ir_fn2,'String');
    ir_last = [tmp(1:(end-4)) '_IR.tif'];
    tmp = get(handles.ir_fn3,'String');
    ir_next = [tmp(1:(end-4)) '_IR.tif'];
    copyfile(ir_last,ir_next);
end
ir_fn3 = get(handles.ir_fn3,'String');
addsharp = get(handles.addsharp3,'Value');
if (exist(ir_fn3,'file') == 2)
    art_register(rgb_fn,sc,ir_fn3,j01,j02,repcnt_max,max_shift,minpts,Npoly,Ntran,memory_limited,morepts,useresults,filteron,addsharp,debugmode,usegui,feather_width,mxd)
end
%}
function j01_Callback(hObject, eventdata, handles)
% hObject    handle to j01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function j01_CreateFcn(hObject, eventdata, handles)
% hObject    handle to j01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function rgb_fn_Callback(hObject, eventdata, handles)
% hObject    handle to rgb_fn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rgb_fn as text
%        str2double(get(hObject,'String')) returns contents of rgb_fn as a double

% --- Executes during object creation, after setting all properties.
function rgb_fn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rgb_fn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fn,path] = uigetfile('*.tif','Select the RGB file.');
set(handles.rgb_fn,'String',[path fn]);

function j02_Callback(hObject, eventdata, handles)
% hObject    handle to j02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of j02 as text
%        str2double(get(hObject,'String')) returns contents of j02 as a double

% --- Executes during object creation, after setting all properties.
function j02_CreateFcn(hObject, eventdata, handles)
% hObject    handle to j02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sc_Callback(hObject, eventdata, handles)
% hObject    handle to sc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sc as text
%        str2double(get(hObject,'String')) returns contents of sc as a double

% --- Executes during object creation, after setting all properties.
function sc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function repcnt_max_Callback(hObject, eventdata, handles)
% hObject    handle to repcnt_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of repcnt_max as text
%        str2double(get(hObject,'String')) returns contents of repcnt_max as a double

% --- Executes during object creation, after setting all properties.
function repcnt_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to repcnt_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function max_shift_Callback(hObject, eventdata, handles)
% hObject    handle to max_shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_shift as text
%        str2double(get(hObject,'String')) returns contents of max_shift as a double

% --- Executes during object creation, after setting all properties.
function max_shift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function minpts_Callback(hObject, eventdata, handles)
% hObject    handle to minpts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minpts as text
%        str2double(get(hObject,'String')) returns contents of minpts as a double

% --- Executes during object creation, after setting all properties.
function minpts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minpts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in filteron.
function filteron_Callback(hObject, eventdata, handles)
% hObject    handle to filteron (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of filteron

% --- Executes on button press in useresults.
function useresults_Callback(hObject, eventdata, handles)
% hObject    handle to useresults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of useresults

% --- Executes on button press in morepts.
function morepts_Callback(hObject, eventdata, handles)
% hObject    handle to morepts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of morepts

% --- Executes on button press in memory_limited.
function memory_limited_Callback(hObject, eventdata, handles)
% hObject    handle to memory_limited (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of memory_limited

% --- Executes on button press in pushbutton_extract1.
function pushbutton_extract1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_extract1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pc = get(handles.pc,'Value');
fn_path = get(handles.fn_path,'String');
trial = get(handles.trial,'String');
N1 = get(handles.N1,'String');
N2 = get(handles.N2,'String');
sc1 = str2double(get(handles.sc1,'String'));
dir1 = str2double(get(handles.extract_dir1,'String'));

if (~isempty(N1))
    N1 = str2double(N1);
    N2 = str2double(N2);
    for i = N1:N2
        if (pc == 1)
            fn = [fn_path '\' trial '_' num2str(i) '.img'];
        elseif (pc == 0)
            fn = [fn_path '/' trial '_' num2str(i) '.img'];
        end
        if (exist(fn,'file') == 2)
            ir_extract(fn,dir1,1,sc1)
        end
    end
    if (pc == 1)
        fn1 = [fn_path '\' trial '_' num2str(N1) '\' trial '_' num2str(N1) '_' num2str(sprintf('%03.0f',1)) '.tif'];
    elseif (pc == 0)
        fn1 = [fn_path '/' trial '_' num2str(N1) '/' trial '_' num2str(N1) '_' num2str(sprintf('%03.0f',1)) '.tif'];
    end
else
    if (pc == 1)
        fn = [fn_path '\' trial '.img'];
    elseif (pc == 0)
        fn = [fn_path '/' trial '.img'];
    end
    if (exist(fn,'file') == 2)
        ir_extract(fn,dir1,1,sc1)
    end
    if (pc == 1)
        fn1 = [fn_path '\' trial '\' trial '_' num2str(sprintf('%03.0f',1)) '.tif'];
    elseif (pc == 0)
        fn1 = [fn_path '/' trial '/' trial '_' num2str(sprintf('%03.0f',1)) '.tif'];
    end
end

img = imread(fn1,'tif');
[~,~,p0] = size(img);
if (p0 == 1)
    set(handles.rgb_data,'Value',0);
elseif (p0 == 3)
    set(handles.rgb_data,'Value',1);
end
clear img

message = sprintf('Extraction has completed.');
uiwait(msgbox(message));


% --- Executes on button press in pushbutton_extract2.
function pushbutton_extract2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_extract2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pc = get(handles.pc,'Value');
fn_path = get(handles.fn_path,'String');
trial = get(handles.trial,'String');
img_ext = get(handles.img_ext,'String');
N1 = str2double(get(handles.N1,'String'));
N2 = str2double(get(handles.N2,'String'));
sc1 = str2double(get(handles.sc1,'String'));
sz = str2double(get(handles.sz,'String'));

if (~isempty(N1))
    for j = N1:N2
        if (pc == 1)
            fn = [fn_path '\'  trial '_' num2str(j) '.' img_ext];
        elseif (pc == 0)
            fn = [fn_path '/'  trial '_' num2str(j) '.' img_ext];
        end
        if (exist(fn,'file') == 2)
            xray_extract(fn_path,trial,j,img_ext,sc1,sz,pc);
        end
    end
    p0 = 3;
else
    if (pc == 1)
        fn = [fn_path '\' trial '.' img_ext];
    elseif (pc == 0)
        fn = [fn_path '/' trial '.' img_ext];
    end
    if (exist(fn,'file') == 2)
        xray_extract(fn_path,trial,[],img_ext,sc1,sz,pc);
    end
    if (pc == 1)
        fn1 = [fn_path '\' trial '\' trial '_' num2str(sprintf('%03.0f',1)) '.tif'];
    elseif (pc == 0)
        fn1 = [fn_path '/' trial '/' trial '_' num2str(sprintf('%03.0f',1)) '.tif'];
    end
    img = imread(fn1,'tif');
    [~,~,p0] = size(img);  
    clear img
end

if (p0 == 1)
    set(handles.rgb_data,'Value',0);
elseif (p0 == 3)
    set(handles.rgb_data,'Value',1);
end

message = sprintf('Segmentation has completed.');
uiwait(msgbox(message));


function sc1_Callback(hObject, eventdata, handles)
% hObject    handle to sc1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sc1 as text
%        str2double(get(hObject,'String')) returns contents of sc1 as a double

% --- Executes during object creation, after setting all properties.
function sc1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sc1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addsharp.
function addsharp_Callback(hObject, eventdata, handles)
% hObject    handle to addsharp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of addsharp

% --- Executes on button press in debugmode.
function debugmode_Callback(hObject, eventdata, handles)
% hObject    handle to debugmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of debugmode


% --- Executes on button press in shift_ref.
function shift_ref_Callback(hObject, eventdata, handles)
% hObject    handle to shift_ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of shift_ref


% --- Executes on button press in buildon.
function buildon_Callback(hObject, eventdata, handles)
% hObject    handle to buildon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of buildon

% --- Executes on button press in stop.
function stop_Callback(hObject, eventdata, handles)
% hObject    handle to stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of stop



function feather_width_Callback(hObject, eventdata, handles)
% hObject    handle to feather_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of feather_width as text
%        str2double(get(hObject,'String')) returns contents of feather_width as a double


function mxd_Callback(hObject, eventdata, handles)
% hObject    handle to mxd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mxd as text
%        str2double(get(hObject,'String')) returns contents of mxd as a double


% --- Executes during object creation, after setting all properties.
function mxd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mxd (see GCBO)
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
%        str2double(get(hObject,'String')) returns contents of ncores as a double


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



function fn_path_Callback(hObject, eventdata, handles)
% hObject    handle to fn_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fn_path as text
%        str2double(get(hObject,'String')) returns contents of fn_path as a double


% --- Executes during object creation, after setting all properties.
function fn_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fn_path (see GCBO)
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

% Hints: get(hObject,'String') returns contents of trial as text
%        str2double(get(hObject,'String')) returns contents of trial as a double


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



function N1_Callback(hObject, eventdata, handles)
% hObject    handle to N1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of N1 as text
%        str2double(get(hObject,'String')) returns contents of N1 as a double


% --- Executes during object creation, after setting all properties.
function N1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function N2_Callback(hObject, eventdata, handles)
% hObject    handle to N2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of N2 as text
%        str2double(get(hObject,'String')) returns contents of N2 as a double


% --- Executes during object creation, after setting all properties.
function N2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fn,path] = uigetfile('*','Select the first image in the set.');
set(handles.fn_path,'String',path(1:(end-1)));
ind = regexp(fn, '\.');
if (isempty(ind))
    img_ext = '';
else
    img_ext = fn((ind+1):end);
    fn(ind:end) = [];
end
ind = regexp(fn, '[0-9]$');
if (isempty(ind))
    fn1 = fn;
    N1 = '';
    N2 = '';
    set(handles.N2,'String',N2);
else
    fn1 = fn(1:(ind-2));
    N1 = fn(ind:end);
end
set(handles.trial,'String',fn1);
set(handles.N1,'String',N1);
set(handles.img_ext,'String',img_ext);


function img_ext_Callback(hObject, eventdata, handles)
% hObject    handle to img_ext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of img_ext as text
%        str2double(get(hObject,'String')) returns contents of img_ext as a double


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



function extract_dir1_Callback(hObject, eventdata, handles)
% hObject    handle to extract_dir1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of extract_dir1 as text
%        str2double(get(hObject,'String')) returns contents of extract_dir1 as a double


% --- Executes during object creation, after setting all properties.
function extract_dir1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to extract_dir1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rgb_data_Callback(hObject, eventdata, handles)
% hObject    handle to rgb_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rgb_data as text
%        str2double(get(hObject,'String')) returns contents of rgb_data as a double


% --- Executes during object creation, after setting all properties.
function rgb_data_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rgb_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sz_Callback(hObject, eventdata, handles)
% hObject    handle to sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sz as text
%        str2double(get(hObject,'String')) returns contents of sz as a double


% --- Executes during object creation, after setting all properties.
function sz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rmosaic.
function rmosaic_Callback(hObject, eventdata, handles)
% hObject    handle to rmosaic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pc = get(handles.pc,'Value');
fn_path = get(handles.fn_path,'String');
trial = get(handles.trial,'String');
N1 = get(handles.N1,'String');
N2 = get(handles.N2,'String');

if (~isempty(N1))
    N1 = str2double(N1);
    N2 = str2double(N2);
    for i = N1:N2
        if (pc == 1)
            fn = [fn_path '\' trial '_' num2str(i) '.img'];
        elseif (pc == 0)
            fn = [fn_path '/' trial '_' num2str(i) '.img'];
        end
        if (exist(fn,'file') == 2)
            rough_mosaic(fn,1);
        end
    end
else
    if (pc == 1)
        fn = [fn_path '\' trial '.img'];
    elseif (pc == 0)
        fn = [fn_path '/' trial '.img'];
    end
    if (exist(fn,'file') == 2)
        rough_mosaic(fn,1);
    end
end
