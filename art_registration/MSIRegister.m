function varargout = MSIRegister(varargin)
% MSIREGISTER MATLAB code for MSIRegister.fig
%      MSIREGISTER, by itself, creates a new MSIREGISTER or raises the existing
%      singleton*.
%
%      H = MSIREGISTER returns the handle to a new MSIREGISTER or the handle to
%      the existing singleton*.
%
%      MSIREGISTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MSIREGISTER.M with the given input arguments.
%
%      MSIREGISTER('Property','Value',...) creates a new MSIREGISTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MSIRegister_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MSIRegister_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MSIRegister

% Last Modified by GUIDE v2.5 16-Dec-2012 18:26:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MSIRegister_OpeningFcn, ...
                   'gui_OutputFcn',  @MSIRegister_OutputFcn, ...
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

% --- Executes just before MSIRegister is made visible.
function MSIRegister_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MSIRegister (see VARARGIN)

% Choose default command line output for MSIRegister
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MSIRegister wait for user response (see UIRESUME)
% uiwait(handles.figure1);
%axes(handles.axes1)
%imshow([])
set(handles.running,'Visible','on');
set(handles.running,'String','');
%ncores = getenv('NUMBER_OF_PROCESSORS');
ncores = feature('NumCores');
set(handles.ncores,'String',num2str(ncores));
set(handles.running,'String','');
warning('off','all')

% --- Outputs from this function are returned to the command line.
function varargout = MSIRegister_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.running,'Visible','on');
set(handles.running,'String','Running ...');
drawnow;
pause(0.1)

fn_path = get(handles.fn_path,'String');
fn = get(handles.fn,'String');
cube_fn = [fn_path fn];
j01 = str2double(get(handles.j01,'String'));
j02 = str2double(get(handles.j02,'String'));
repcnt_max = str2double(get(handles.repcnt_max,'String'));
max_shift = str2double(get(handles.max_shift,'String'));
minpts = str2double(get(handles.minpts,'String'));
morepts = get(handles.morepts,'Value');
filteron = get(handles.filteron,'Value');
usegui = 1;

cube_register(cube_fn,j01,j02,repcnt_max,minpts,max_shift,filteron,morepts,usegui);

%axes(handles.axes1)
%imshow([])
set(handles.running,'Visible','off');
set(handles.running,'String','');
set(handles.text17,'String','');
set(handles.text20,'String','');
set(handles.text18,'String','');
set(handles.text19,'String','');
drawnow;
pause(0.1)

function j01_Callback(hObject, eventdata, handles)
% hObject    handle to j01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of j01 as text
%        str2double(get(hObject,'String')) returns contents of j01 as a double

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

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fn,fn_path] = uigetfile('*','Select the first cube to be registered.');
set(handles.fn_path,'String',fn_path);
set(handles.fn,'String',fn);

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

% --- Executes on button press in morepts.
function morepts_Callback(hObject, eventdata, handles)
% hObject    handle to morepts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of morepts

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

function fn_Callback(hObject, eventdata, handles)
% hObject    handle to fn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fn as text
%        str2double(get(hObject,'String')) returns contents of fn as a double

% --- Executes during object creation, after setting all properties.
function fn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fn (see GCBO)
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
