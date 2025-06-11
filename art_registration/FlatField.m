function varargout = FlatField(varargin)
% FLATFIELD MATLAB code for FlatField.fig
%      FLATFIELD, by itself, creates a new FLATFIELD or raises the existing
%      singleton*.
%
%      H = FLATFIELD returns the handle to a new FLATFIELD or the handle to
%      the existing singleton*.
%
%      FLATFIELD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FLATFIELD.M with the given input arguments.
%
%      FLATFIELD('Property','Value',...) creates a new FLATFIELD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FlatField_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FlatField_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FlatField

% Last Modified by GUIDE v2.5 29-Jul-2014 13:59:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FlatField_OpeningFcn, ...
                   'gui_OutputFcn',  @FlatField_OutputFcn, ...
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


% --- Executes just before FlatField is made visible.
function FlatField_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FlatField (see VARARGIN)

% Choose default command line output for FlatField
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FlatField wait for user response (see UIRESUME)
% uiwait(handles.figure1);

axes(handles.axes1)
imshow([])

% --- Outputs from this function are returned to the command line.
function varargout = FlatField_OutputFcn(hObject, eventdata, handles) 
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

fn_path = get(handles.path1,'String');
pc = ~ismac();
if (pc == 1)
    fn_path2 = strcat(fn_path,'out/');
elseif (pc == 0)
    fn_path2 = strcat(fn_path,'out\');
end
B = get(handles.B,'String');
W = get(handles.W,'String');
N = str2double(get(handles.N,'String'));
trial = get(handles.fn,'String');

if (exist(fn_path2,'dir') == 7)
    rmdir(fn_path2,'s');
end
mkdir(fn_path2);

fn = strcat(fn_path,B);
black = double(imread(fn,'tif'));
fn = strcat(fn_path,W);
white = double(imread(fn,'tif'));

mn = str2double(get(handles.mn,'String'));
mx = str2double(get(handles.mx,'String'));
for i = 1:N
    if (N < 10)
        fn = [fn_path trial num2str(sprintf(['%01.0f'],i)) '.tif'];
        fn2 = [fn_path2 trial num2str(sprintf(['%01.0f'],i)) '.tif'];
    elseif (N < 100)
        fn = [fn_path trial num2str(sprintf(['%02.0f'],i)) '.tif'];
        fn2 = [fn_path2 trial num2str(sprintf(['%02.0f'],i)) '.tif'];
    elseif (N < 1000)
        fn = [fn_path trial num2str(sprintf(['%03.0f'],i)) '.tif'];
        fn2 = [fn_path2 trial num2str(sprintf(['%03.0f'],i)) '.tif'];
    elseif (N < 10000)
        fn = [fn_path trial num2str(sprintf(['%04.0f'],i)) '.tif'];
        fn2 = [fn_path2 trial num2str(sprintf(['%04.0f'],i)) '.tif'];
    end
    if (exist(fn,'file') == 2)
        img = double(imread(fn,'tif'));
        ff = (img - black)./(white - black);
        clear img
        ff = (ff - mn)/(mx - mn);
        ff = (2^16-1)*ff;
        ff = uint16(ff);

        imwrite(ff,fn2,'tif','Compression','None');
        clear ff
    end
end
axes(handles.axes1)
imshow([])
set(handles.mn,'String','');
set(handles.mx,'String','');
pause(0.1)

% --- Executes on button press in sel_file.
function sel_file_Callback(hObject, eventdata, handles)
% hObject    handle to sel_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fn,path] = uigetfile({'*.tif';'*.TIF'},'Select the last TIF file in the set.');

set(handles.path1,'String',path);

ind = regexp(fn, '[0-9]');
set(handles.fn,'String',fn(1:(ind-1)));
set(handles.N,'String',fn(ind:(end-4)));

function path1_Callback(hObject, eventdata, handles)
% hObject    handle to path1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function path1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to path1 (see GCBO)
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

% --- Executes on button press in sel_W.
function sel_W_Callback(hObject, eventdata, handles)
% hObject    handle to sel_W (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fn,path] = uigetfile({'*.tif';'*.TIF'},'Select the white card image.',get(handles.path1,'String'));
set(handles.path2,'String',path);
set(handles.W,'String',fn);

function path2_Callback(hObject, eventdata, handles)
% hObject    handle to path2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function path2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to path2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function W_Callback(hObject, eventdata, handles)
% hObject    handle to W (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function W_CreateFcn(hObject, eventdata, handles)
% hObject    handle to W (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in sel_B.
function sel_B_Callback(hObject, eventdata, handles)
% hObject    handle to sel_B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fn,path] = uigetfile({'*.tif';'*.TIF'},'Select the dark card image.',get(handles.path1,'String'));
set(handles.path3,'String',path);
set(handles.B,'String',fn);

function path3_Callback(hObject, eventdata, handles)
% hObject    handle to path3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function path3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to path3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function B_Callback(hObject, eventdata, handles)
% hObject    handle to B (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function B_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B (see GCBO)
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

function mn_Callback(hObject, eventdata, handles)
% hObject    handle to mn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function mn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mx_Callback(hObject, eventdata, handles)
% hObject    handle to mx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function mx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in preview.
function preview_Callback(hObject, eventdata, handles)
% hObject    handle to preview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fn_path = get(handles.path1,'String');
fn_path2 = strcat(fn_path,'out/');
B = get(handles.B,'String');
W = get(handles.W,'String');
N = str2double(get(handles.N,'String'));
trial = get(handles.fn,'String');

fn = strcat(fn_path,B);
black = double(imread(fn,'tif'));
fn = strcat(fn_path,W);
white = double(imread(fn,'tif'));

% find global max and min
mx = -2^16-1;
mn = 2^16-1;
for i = 1:N
    if (N < 10)
        fn = [fn_path trial num2str(sprintf(['%01.0f'],i)) '.tif'];
    elseif (N < 100)
        fn = [fn_path trial num2str(sprintf(['%02.0f'],i)) '.tif'];
    elseif (N < 1000)
        fn = [fn_path trial num2str(sprintf(['%03.0f'],i)) '.tif'];
    elseif (N < 10000)
        fn = [fn_path trial num2str(sprintf(['%04.0f'],i)) '.tif'];
    end
    if (exist(fn,'file') == 2)
        img = double(imread(fn,'tif'));
        ff = (img - black)./(white - black);
        clear img
        mx0 = max(ff(:))
        if (mx0 > mx)
            mx = mx0;
        end
        mn0 = min(ff(:));
        if (mn0 < mn)
            mn = mn0;
        end
        clear ff
    end
end

% compute overall histogram
ghist = double(zeros(1,1000));
x = linspace(mn,mx,1000);
for i = 1:N
    if (N < 10)
        fn = [fn_path trial num2str(sprintf(['%01.0f'],i)) '.tif'];
        fn2 = [fn_path2 trial num2str(sprintf(['%01.0f'],i)) '.tif'];
    elseif (N < 100)
        fn = [fn_path trial num2str(sprintf(['%02.0f'],i)) '.tif'];
        fn2 = [fn_path2 trial num2str(sprintf(['%02.0f'],i)) '.tif'];
    elseif (N < 1000)
        fn = [fn_path trial num2str(sprintf(['%03.0f'],i)) '.tif'];
        fn2 = [fn_path2 trial num2str(sprintf(['%03.0f'],i)) '.tif'];
    elseif (N < 10000)
        fn = [fn_path trial num2str(sprintf(['%04.0f'],i)) '.tif'];
        fn2 = [fn_path2 trial num2str(sprintf(['%04.0f'],i)) '.tif'];
    end
    if (exist(fn,'file') == 2)
        img = double(imread(fn,'tif'));
        ff = (img - black)./(white - black);
        clear img
        tmp = hist(ff(:),x);
        ghist = ghist + tmp;
        clear ff tmp
    end
end
axes(handles.axes1)
bar(x,ghist)
set(handles.mn,'String',num2str(mn));
set(handles.mx,'String',num2str(mx));


% --- Executes on button press in create_img.
function create_img_Callback(hObject, eventdata, handles)
% hObject    handle to create_img (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fn_path = get(handles.path1,'String');
B = get(handles.B,'String');
W = get(handles.W,'String');
N = str2double(get(handles.N,'String'));
trial = get(handles.fn,'String');

fn = strcat(fn_path,B);
black = double(imread(fn,'tif'));
fn = strcat(fn_path,W);
white = double(imread(fn,'tif'));

p0 = 0;
for i = 1:N
    if (N < 10)
        fn = [fn_path trial num2str(sprintf(['%01.0f'],i)) '.tif'];
    elseif (N < 100)
        fn = [fn_path trial num2str(sprintf(['%02.0f'],i)) '.tif'];
    elseif (N < 1000)
        fn = [fn_path trial num2str(sprintf(['%03.0f'],i)) '.tif'];
    elseif (N < 10000)
        fn = [fn_path trial num2str(sprintf(['%04.0f'],i)) '.tif'];
    end
    if (exist(fn,'file') == 2)
        img = imread(fn,'tif');
        [m0,n0] = size(img);
        p0 = p0 + 1;
    end
end

header = uint16(zeros(1,256));
header(2) = m0;
header(3) = n0;
header(22) = p0;

fn_full = [fn_path trial '1.img'];
fid = fopen(fn_full,'w');
fwrite(fid,header,'uint16');
clear header

mn = str2double(get(handles.mn,'String'));
mx = str2double(get(handles.mx,'String'));
for i = 1:N
    if (N < 10)
        fn = [fn_path trial num2str(sprintf(['%01.0f'],i)) '.tif'];
    elseif (N < 100)
        fn = [fn_path trial num2str(sprintf(['%02.0f'],i)) '.tif'];
    elseif (N < 1000)
        fn = [fn_path trial num2str(sprintf(['%03.0f'],i)) '.tif'];
    elseif (N < 10000)
        fn = [fn_path trial num2str(sprintf(['%04.0f'],i)) '.tif'];
    end
    if (exist(fn,'file') == 2)
        img = double(imread(fn,'tif'));
        ff = (img - black)./(white - black);
        clear img
        ff = (ff - mn)/(mx - mn);
        ff = (2^16-1)*ff;
        ff = uint16(ff);
        ff = ff(m0:-1:1,:);

        fwrite(fid,ff(:),'uint16');
        clear ff
    end
end
axes(handles.axes1)
imshow([])
set(handles.mn,'String','');
set(handles.mx,'String','');
pause(1)
fclose(fid);
