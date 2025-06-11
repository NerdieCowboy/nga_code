function varargout = FlatFieldCube(varargin)
% FLATFIELDCUBE MATLAB code for FlatFieldCube.fig
%      FLATFIELDCUBE, by itself, creates a new FLATFIELDCUBE or raises the existing
%      singleton*.
%
%      H = FLATFIELDCUBE returns the handle to a new FLATFIELDCUBE or the handle to
%      the existing singleton*.
%
%      FLATFIELDCUBE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FLATFIELDCUBE.M with the given input arguments.
%
%      FLATFIELDCUBE('Property','Value',...) creates a new FLATFIELDCUBE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FlatFieldCube_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FlatFieldCube_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FlatFieldCube

% Last Modified by GUIDE v2.5 12-Mar-2013 14:04:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FlatFieldCube_OpeningFcn, ...
                   'gui_OutputFcn',  @FlatFieldCube_OutputFcn, ...
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

% --- Executes just before FlatFieldCube is made visible.
function FlatFieldCube_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FlatFieldCube (see VARARGIN)

% Choose default command line output for FlatFieldCube
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FlatFieldCube wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = FlatFieldCube_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function ir_fn_Callback(hObject, eventdata, handles)
% hObject    handle to ir_fn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ir_fn as text
%        str2double(get(hObject,'String')) returns contents of ir_fn as a double


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

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('*.*','Select the cube');
set(handles.ir_fn,'String',[PathName FileName]);
set(handles.PathName,'String',PathName);

function w_fn_Callback(hObject, eventdata, handles)
% hObject    handle to w_fn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of w_fn as text
%        str2double(get(hObject,'String')) returns contents of w_fn as a double

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

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PathName = get(handles.PathName,'String');
[FileName,PathName] = uigetfile([PathName '*.*'],'Select the white cube');
set(handles.w_fn,'String',[PathName FileName]);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text2,'String','Running ...');
set(handles.axes1,'Visible','off');
drawnow;

usegui = 1;
ir_fn = get(handles.ir_fn,'String');
w_fn = get(handles.w_fn,'String');
cube_rotate_ff(ir_fn,w_fn,usegui);
