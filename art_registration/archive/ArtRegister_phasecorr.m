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

% Last Modified by GUIDE v2.5 17-Sep-2013 19:32:11

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

% UIWAIT makes ArtRegister wait for user response (see UIRESUME)
% uiwait(handles.figure1);
%axes(handles.axes1)
%imshow([])
%axes(handles.axes2)
%imshow([])
%axes(handles.axes4)
%imshow([])
%axes(handles.axes5)
%imshow([])
%set(handles.running,'String','');
%set(handles.stop,'Value',0);

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

function ir_fn_Callback(hObject, eventdata, handles)
% hObject    handle to ir_fn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function ir_fn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ir_fn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

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
ir_fn = get(handles.ir_fn,'String');
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
filteron = get(handles.filteron,'Value');
addsharp = get(handles.addsharp,'Value');
debugmode = get(handles.debugmode,'Value');
usegui = 1;
mxd = str2double(get(handles.mxd,'String'));
art_register(rgb_fn,sc,ir_fn,j01,j02,repcnt_max,max_shift,minpts,Npoly,Ntran,memory_limited,morepts,filteron,addsharp,debugmode,usegui,mxd)

if (get(handles.shift_ref,'Value') == 1)
    sc = 1;
    tmp = get(handles.ir_fn,'String');
    rgb_fn = [tmp(1:(end-4)) '_IR.tif'];
end
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
    art_register(rgb_fn,sc,ir_fn2,j01,j02,repcnt_max,max_shift,minpts,Npoly,Ntran,memory_limited,morepts,filteron,addsharp,debugmode,usegui,mxd)
end

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
    art_register(rgb_fn,sc,ir_fn3,j01,j02,repcnt_max,max_shift,minpts,Npoly,Ntran,memory_limited,morepts,filteron,addsharp,debugmode,usegui,mxd)
end

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

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (get(handles.bsets,'Value') == 1)
    [fn,path] = uigetfile('*.img','Select a template IMG file.');
elseif (get(handles.bsets,'Value') == 0)
    [fn,path] = uigetfile('*.tif','Select a template TIF file.');
end
set(handles.ir_fn,'String',[path fn]);

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

function ir_fn2_Callback(hObject, eventdata, handles)
% hObject    handle to ir_fn2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ir_fn2 as text
%        str2double(get(hObject,'String')) returns contents of ir_fn2 as a double

% --- Executes during object creation, after setting all properties.
function ir_fn2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ir_fn2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (get(handles.bsets,'Value') == 1)
    [fn,path] = uigetfile('*.img','Select a template IMG file.');
elseif (get(handles.bsets,'Value') == 0)
    [fn,path] = uigetfile('*.tif','Select a template TIF file.');
end
set(handles.ir_fn2,'String',[path fn]);

function ir_fn3_Callback(hObject, eventdata, handles)
% hObject    handle to ir_fn3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ir_fn3 as text
%        str2double(get(hObject,'String')) returns contents of ir_fn3 as a double

% --- Executes during object creation, after setting all properties.
function ir_fn3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ir_fn3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (get(handles.bsets,'Value') == 1)
    [fn,path] = uigetfile('*.img','Select a template IMG file.');
elseif (get(handles.bsets,'Value') == 0)
    [fn,path] = uigetfile('*.tif','Select a template TIF file.');
end
set(handles.ir_fn3,'String',[path fn]);

% --- Executes on button press in filteron.
function filteron_Callback(hObject, eventdata, handles)
% hObject    handle to filteron (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of filteron

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
if (get(handles.bsets,'Value') == 0)
    ir_extractX(get(handles.ir_fn,'String'),1,str2double(get(handles.sc1,'String')))
elseif (get(handles.bsets,'Value') == 1)
    ir_extract(get(handles.ir_fn,'String'),str2double(get(handles.extract_dir1,'String')),1,str2double(get(handles.sc1,'String')))
end

% --- Executes on button press in pushbutton_extract2.
function pushbutton_extract2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_extract2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (get(handles.bsets,'Value') == 0)
    ir_extractX(get(handles.ir_fn2,'String'),1,str2double(get(handles.sc2,'String')))
elseif (get(handles.bsets,'Value') == 1)
    ir_extract(get(handles.ir_fn2,'String'),str2double(get(handles.extract_dir2,'String')),1,str2double(get(handles.sc2,'String')))
end

% --- Executes on button press in pushbutton_extract3.
function pushbutton_extract3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_extract3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (get(handles.bsets,'Value') == 0)
    ir_extractX(get(handles.ir_fn3,'String'),1,str2double(get(handles.sc3,'String')))
elseif (get(handles.bsets,'Value') == 1)
    ir_extract(get(handles.ir_fn3,'String'),str2double(get(handles.extract_dir3,'String')),1,str2double(get(handles.sc3,'String')))
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

function extract_dir2_Callback(hObject, eventdata, handles)
% hObject    handle to extract_dir2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of extract_dir2 as text
%        str2double(get(hObject,'String')) returns contents of extract_dir2 as a double

% --- Executes during object creation, after setting all properties.
function extract_dir2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to extract_dir2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function extract_dir3_Callback(hObject, eventdata, handles)
% hObject    handle to extract_dir3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of extract_dir3 as text
%        str2double(get(hObject,'String')) returns contents of extract_dir3 as a double

% --- Executes during object creation, after setting all properties.
function extract_dir3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to extract_dir3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

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

function sc2_Callback(hObject, eventdata, handles)
% hObject    handle to sc2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sc2 as text
%        str2double(get(hObject,'String')) returns contents of sc2 as a double

% --- Executes during object creation, after setting all properties.
function sc2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sc2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sc3_Callback(hObject, eventdata, handles)
% hObject    handle to sc3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sc3 as text
%        str2double(get(hObject,'String')) returns contents of sc3 as a double

% --- Executes during object creation, after setting all properties.
function sc3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sc3 (see GCBO)
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

% --- Executes on button press in addsharp2.
function addsharp2_Callback(hObject, eventdata, handles)
% hObject    handle to addsharp2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of addsharp2

% --- Executes on button press in addsharp3.
function addsharp3_Callback(hObject, eventdata, handles)
% hObject    handle to addsharp3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of addsharp3

% --- Executes on button press in shift_ref.
function shift_ref_Callback(hObject, eventdata, handles)
% hObject    handle to shift_ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of shift_ref

% --- Executes on button press in bsets.
function bsets_Callback(hObject, eventdata, handles)
% hObject    handle to bsets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bsets

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
