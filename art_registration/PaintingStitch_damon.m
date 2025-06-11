function varargout = PaintingStitch(varargin)
% PAINTINGSTITCH M-file for PaintingStitch.fig
%      PAINTINGSTITCH, by itself, creates a new PAINTINGSTITCH or raises the existing
%      singleton*.
%
%      H = PAINTINGSTITCH returns the handle to a new PAINTINGSTITCH or the handle to
%      the existing singleton*.
%
%      PAINTINGSTITCH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PAINTINGSTITCH.M with the given input arguments.
%
%      PAINTINGSTITCH('Property','Value',...) creates a new PAINTINGSTITCH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PaintingStitch_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PaintingStitch_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PaintingStitch

% Last Modified by GUIDE v2.5 17-Feb-2013 23:36:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PaintingStitch_OpeningFcn, ...
                   'gui_OutputFcn',  @PaintingStitch_OutputFcn, ...
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

% --- Executes just before PaintingStitch is made visible.
function PaintingStitch_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PaintingStitch (see VARARGIN)

% Choose default command line output for PaintingStitch
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PaintingStitch wait for user response (see UIRESUME)
% uiwait(handles.figure1);
ncores = feature('NumCores');
set(handles.ncores,'String',num2str(ncores));
warning('off','all')

pc = 1;
set(handles.text16,'String',num2str(pc));
if (pc == 1)
    set(handles.text2,'String','\');
    set(handles.text3,'String','\');
elseif (pc == 0)
    set(handles.text2,'String','/');
    set(handles.text3,'String','/');
end

% --- Outputs from this function are returned to the command line.
function varargout = PaintingStitch_OutputFcn(hObject, eventdata, handles) 
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

function fn_path_Callback(hObject, eventdata, handles)
% hObject    handle to fn_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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

function set_name_Callback(hObject, eventdata, handles)
% hObject    handle to set_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function set_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_name (see GCBO)
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

function N_Callback(hObject, eventdata, handles)
% hObject    handle to N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function N_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fn_out_Callback(hObject, eventdata, handles)
% hObject    handle to fn_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function fn_out_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fn_out (see GCBO)
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

% --- Executes on button press in quick.
function quick_Callback(hObject, eventdata, handles)
% hObject    handle to quick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function mm_Callback(hObject, eventdata, handles)
% hObject    handle to mm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function mm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pixels_Callback(hObject, eventdata, handles)
% hObject    handle to pixels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function pixels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pixels (see GCBO)
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
set(handles.text15,'String','Running ...');
drawnow;

fn_path = get(handles.fn_path,'String');
set_name = get(handles.set_name,'String');
trial = get(handles.trial,'String');
fn_out = get(handles.fn_out,'String');
mmpixel = str2double(get(handles.mm,'String'))/str2double(get(handles.pixels,'String'));
N = str2double(get(handles.N,'String'));
feather_width = str2double(get(handles.feather_width,'String'));
pc = str2double(get(handles.text16,'String'));
ncores = str2double(get(handles.ncores,'String'));
turnonoutput = get(handles.turnonoutput,'Value');
quick = get(handles.quick,'Value');
usegui = 1;
mosaic_full(fn_path,set_name,trial,fn_out,mmpixel,N,feather_width,pc,ncores,turnonoutput,quick,usegui);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%{
if (get(handles.radiobutton1,'Value') == 1)
    pc = 1;
elseif (get(handles.radiobutton2,'Value') == 1)
    pc = 0;
end
%}
pc = str2double(get(handles.text16,'String'));
[fn,path] = uigetfile('*.tif','Select the last TIF file in the set.');

if (pc == 1)
    tmp = regexp(path,'\\');
elseif (pc == 0)
    tmp = regexp(path,'/');
end
set(handles.fn_path,'String',path(1:(tmp(end-1)-1)));
set(handles.set_name,'String',path((tmp(end-1)+1):(tmp(end)-1)));

ind = regexp(fn, '[0-9]');
set(handles.trial,'String',fn(1:(ind-1)));
set(handles.N,'String',fn(ind:(end-4)));

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
