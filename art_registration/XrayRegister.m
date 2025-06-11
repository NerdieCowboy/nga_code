function varargout = XrayRegister(varargin)
% XRAYREGISTER M-file for XrayRegister.fig
%      XRAYREGISTER, by itself, creates a new XRAYREGISTER or raises the existing
%      singleton*.
%
%      H = XRAYREGISTER returns the handle to a new XRAYREGISTER or the handle to
%      the existing singleton*.
%
%      XRAYREGISTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in XRAYREGISTER.M with the given input arguments.
%
%      XRAYREGISTER('Property','Value',...) creates a new XRAYREGISTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before XrayRegister_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to XrayRegister_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help XrayRegister

% Last Modified by GUIDE v2.5 05-Oct-2013 19:41:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @XrayRegister_OpeningFcn, ...
                   'gui_OutputFcn',  @XrayRegister_OutputFcn, ...
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

% --- Executes just before XrayRegister is made visible.
function XrayRegister_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to XrayRegister (see VARARGIN)

% Choose default command line output for XrayRegister
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes XrayRegister wait for user response (see UIRESUME)
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
function varargout = XrayRegister_OutputFcn(hObject, eventdata, handles) 
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
usegui = 3;
feather_width = str2double(get(handles.feather_width,'String'));
mxd = str2double(get(handles.mxd,'String'));
fn_path = get(handles.fn_path,'String');
trial = get(handles.trial,'String');
N1 = str2double(get(handles.N1,'String'));
N2 = str2double(get(handles.N2,'String'));
for j = N1:N2
    ir_fn = [fn_path trial '_' num2str(j) '.tif'];
    if (exist(ir_fn,'file') == 2)
        art_register(rgb_fn,sc,ir_fn,j01,j02,repcnt_max,max_shift,minpts,Npoly,Ntran,memory_limited,morepts,useresults,filteron,addsharp,debugmode,usegui,feather_width,mxd)
    end
    ir_last = [fn_path trial '_' num2str(j) '_IR.tif'];
    ir_next = [fn_path trial '_' num2str(j+1) '_IR.tif'];
    copyfile(ir_last,ir_next);
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

fn_path = get(handles.fn_path,'String');
trial = get(handles.trial,'String');
N1 = str2double(get(handles.N1,'String'));
N2 = str2double(get(handles.N2,'String'));
for j = N1:N2
    ir_fn = [fn_path trial '_' num2str(j) '.tif'];
    if (exist(ir_fn,'file') == 2)
        ir_extractX(ir_fn,3,str2double(get(handles.sc1,'String')))
    end
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


% --- Executes during object creation, after setting all properties.
function feather_width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to feather_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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


% --- Executes on button press in pushbutton_xray.
function pushbutton_xray_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_xray (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fn,path] = uigetfile('*','Select the first X-ray film in the set.');
set(handles.fn_path,'String',path);
fn((end-3):end) = [];
ind = regexp(fn, '[0-9]');
set(handles.trial,'String',fn(1:(ind-2)));
set(handles.N1,'String',fn(ind:end));


% --- Executes on button press in addsharp.
function add_sharp_Callback(hObject, eventdata, handles)
% hObject    handle to addsharp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of addsharp



function pc_Callback(hObject, eventdata, handles)
% hObject    handle to pc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pc as text
%        str2double(get(hObject,'String')) returns contents of pc as a double


% --- Executes during object creation, after setting all properties.
function pc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
