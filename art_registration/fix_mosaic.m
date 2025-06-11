function varargout = fix_mosaic(varargin)
% FIX_MOSAIC M-file for fix_mosaic.fig
%      FIX_MOSAIC, by itself, creates a new FIX_MOSAIC or raises the existing
%      singleton*.
%
%      H = FIX_MOSAIC returns the handle to a new FIX_MOSAIC or the handle to
%      the existing singleton*.
%
%      FIX_MOSAIC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIX_MOSAIC.M with the given input arguments.
%
%      FIX_MOSAIC('Property','Value',...) creates a new FIX_MOSAIC or raises the
%      existing singleton*.  Starting from the left, property value pairs
%      are
%      applied to the GUI before fix_mosaic_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fix_mosaic_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fix_mosaic

% Last Modified by GUIDE v2.5 21-Feb-2012 20:04:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fix_mosaic_OpeningFcn, ...
                   'gui_OutputFcn',  @fix_mosaic_OutputFcn, ...
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


% --- Executes just before fix_mosaic is made visible.
function fix_mosaic_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fix_mosaic (see VARARGIN)

% Choose default command line output for fix_mosaic
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fix_mosaic wait for user response (see UIRESUME)
% uiwait(handles.figure1);
img = [];
axes(handles.axes1)
imshow(img)

set(handles.pushbutton1,'enable','on');     %Start
set(handles.pushbutton2,'enable','off');    %Stop
set(handles.pushbutton3,'enable','off');    %Prev
set(handles.pushbutton4,'enable','off');    %Next
set(handles.slider1,'enable','off');        %vertical
set(handles.slider2,'enable','off');        %horizontal

% --- Outputs from this function are returned to the command line.
function varargout = fix_mosaic_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.pushbutton1,'enable','off');   %Start
set(handles.pushbutton2,'enable','on');    %Stop
set(handles.pushbutton3,'enable','on');    %Prev
set(handles.pushbutton4,'enable','on');    %Next
set(handles.slider1,'enable','on');        %vertical
set(handles.slider2,'enable','on');        %horizontal

tmp = regexp(get(handles.edit1,'String'),'\\');
pc = 1;
if (isempty(tmp))
    pc = 0;
end
if (pc == 1)
    filename = [get(handles.edit1,'String') '\' get(handles.edit2,'String') '\offset_values.csv'];
elseif (pc == 0)
    filename = [get(handles.edit1,'String') '/' get(handles.edit2,'String') '/offset_values.csv'];
end
rel_reg = double(csvread(filename));
mr = size(rel_reg,1);

i = 1;
I1 = rel_reg(i,1);
I2 = rel_reg(i+1,1);
dy = rel_reg(i+1,2);
dx = rel_reg(i+1,3);
set(handles.edit5,'String',num2str(i))

if (pc == 1)
    filename = [get(handles.edit1,'String') '\' get(handles.edit2,'String') '\' get(handles.edit2,'String') '_' num2str(sprintf('%03.0f',I1)) '.tif'];
elseif (pc == 0)
    filename = [get(handles.edit1,'String') '/' get(handles.edit2,'String') '/' get(handles.edit2,'String') '_' num2str(sprintf('%03d',I1)) '.tif'];
end
img1 = double(imread(filename,'tif'));
mx = max(img1(:));
img1 = img1/mx;
low_high = stretchlim(img1,[0.01 0.99]);
img1 = 255*imadjust(img1,low_high,[0 1],1);
[m01 n01] = size(img1);

if (pc == 1)
    filename = [get(handles.edit1,'String') '\' get(handles.edit2,'String') '\' get(handles.edit2,'String') '_' num2str(sprintf('%03.0f',I2)) '.tif'];
elseif (pc == 0)
    filename = [get(handles.edit1,'String') '/' get(handles.edit2,'String') '/' get(handles.edit2,'String') '_' num2str(sprintf('%03d',I2)) '.tif'];
end
img2 = double(imread(filename,'tif'));
%mx = max(img2(:));
img2 = img2/mx;
%low_high = stretchlim(img2,[0.01 0.99]);
img2 = 255*imadjust(img2,low_high,[0 1],1);
[m02 n02] = size(img2);
test0 = zeros((m01 + 2*m02 - 2),(n01 + 2*n01 - 2),3);
test0(m02:(m02+m01-1),n02:(n02+n01-1),3) = img1;
test0((m02+dy):(m02+dy+m02-1),(n02+dx):(n02+dx+n02-1),1) = img2;
clear img1 img2

axes(handles.axes1)
imshow(uint8(test0))

set(handles.slider1,'Min',-m01+1);
set(handles.slider1,'Max',m01-1);
set(handles.slider1,'Value',0);
set(handles.slider1,'SliderStep',[1/(2*m01-2) 10/(2*m01-2)]);
set(handles.slider2,'Min',-n01+1);
set(handles.slider2,'Max',n01-1);
set(handles.slider2,'Value',0);
set(handles.slider2,'SliderStep',[1/(2*n01-2) 10/(2*n01-2)]);    


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.pushbutton1,'enable','on');     %Start
set(handles.pushbutton2,'enable','off');    %Stop
set(handles.pushbutton3,'enable','off');    %Prev
set(handles.pushbutton4,'enable','off');    %Next
set(handles.slider1,'enable','off');        %vertical
set(handles.slider2,'enable','off');        %horizontal


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tmp = regexp(get(handles.edit1,'String'),'\\');
pc = 1;
if (isempty(tmp))
    pc = 0;
end
if (pc == 1)
    filename = [get(handles.edit1,'String') '\' get(handles.edit2,'String') '\offset_values.csv'];
elseif (pc == 0)
    filename = [get(handles.edit1,'String') '/' get(handles.edit2,'String') '/offset_values.csv'];
end
rel_reg = double(csvread(filename));
mr = size(rel_reg,1);

i = str2num(get(handles.edit5,'String'));
i = i - 1;
if (i > 1)
    set(handles.edit5,'String',num2str(i))
else
    i = 1;
    set(handles.edit5,'String',num2str(i))
end
I1 = rel_reg(i,1);
I2 = rel_reg(i+1,1);
dy = rel_reg(i+1,2);
dx = rel_reg(i+1,3);

if (pc == 1)
    filename = [get(handles.edit1,'String') '\' get(handles.edit2,'String') '\' get(handles.edit2,'String') '_' num2str(sprintf('%03.0f',I1)) '.tif'];
elseif (pc == 0)
    filename = [get(handles.edit1,'String') '/' get(handles.edit2,'String') '/' get(handles.edit2,'String') '_' num2str(sprintf('%03d',I1)) '.tif'];
end
img1 = double(imread(filename,'tif'));
mx = max(img1(:));
img1 = img1/mx;
low_high = stretchlim(img1,[0.01 0.99]);
img1 = 255*imadjust(img1,low_high,[0 1],1);
[m01 n01] = size(img1);

if (pc == 1)
    filename = [get(handles.edit1,'String') '\' get(handles.edit2,'String') '\' get(handles.edit2,'String') '_' num2str(sprintf('%03.0f',I2)) '.tif'];
elseif (pc == 0)
    filename = [get(handles.edit1,'String') '/' get(handles.edit2,'String') '/' get(handles.edit2,'String') '_' num2str(sprintf('%03d',I2)) '.tif'];
end
img2 = double(imread(filename,'tif'));
%mx = max(img2(:));
img2 = img2/mx;
%low_high = stretchlim(img2,[0.01 0.99]);
img2 = 255*imadjust(img2,low_high,[0 1],1);
[m02 n02] = size(img2);
test0 = zeros((m01 + 2*m02 - 2),(n01 + 2*n01 - 2),3);
test0(m02:(m02+m01-1),n02:(n02+n01-1),3) = img1;
test0((m02+dy):(m02+dy+m02-1),(n02+dx):(n02+dx+n02-1),1) = img2;
clear img1 img2

mx = max(test0(:));
mn = min(test0(:));
test0 = (test0 - mn)/(mx - mn)*255;
axes(handles.axes1)
imshow(uint8(test0))

set(handles.slider1,'Min',-m01+1);
set(handles.slider1,'Max',m01-1);
set(handles.slider1,'Value',0);
set(handles.slider1,'SliderStep',[1/(2*m01-2) 10/(2*m01-2)]);
set(handles.slider2,'Min',-n01+1);
set(handles.slider2,'Max',n01-1);
set(handles.slider2,'Value',0);
set(handles.slider2,'SliderStep',[1/(2*n01-2) 10/(2*n01-2)]);  


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tmp = regexp(get(handles.edit1,'String'),'\\');
pc = 1;
if (isempty(tmp))
    pc = 0;
end
if (pc == 1)
    filename = [get(handles.edit1,'String') '\' get(handles.edit2,'String') '\offset_values.csv'];
elseif (pc == 0)
    filename = [get(handles.edit1,'String') '/' get(handles.edit2,'String') '/offset_values.csv'];
end
rel_reg = double(csvread(filename));
mr = size(rel_reg,1);

i = str2num(get(handles.edit5,'String'));

dy = rel_reg(i+1,2);
dx = rel_reg(i+1,3);
deltay = -round(get(handles.slider1,'Value'));
rel_reg(i+1,2) = dy + deltay;
deltax = round(get(handles.slider2,'Value'));
rel_reg(i+1,3) = dx + deltax;
if (pc == 1)
    filename = [get(handles.edit1,'String') '\' get(handles.edit2,'String') '\offset_values.csv'];
elseif (pc == 0)
    filename = [get(handles.edit1,'String') '/' get(handles.edit2,'String') '/offset_values.csv'];
end
csvwrite(filename,rel_reg)

i = i + 1;
if (i < mr)
    set(handles.edit5,'String',num2str(i))
else
    i = i - 1;
    set(handles.edit5,'String',num2str(i))
end
I1 = rel_reg(i,1);
I2 = rel_reg(i+1,1);
dy = rel_reg(i+1,2);
dx = rel_reg(i+1,3);

if (pc == 1)
    filename = [get(handles.edit1,'String') '\' get(handles.edit2,'String') '\' get(handles.edit2,'String') '_' num2str(sprintf('%03.0f',I1)) '.tif'];
elseif (pc == 0)
    filename = [get(handles.edit1,'String') '/' get(handles.edit2,'String') '/' get(handles.edit2,'String') '_' num2str(sprintf('%03d',I1)) '.tif'];
end
img1 = double(imread(filename,'tif'));
mx = max(img1(:));
img1 = img1/mx;
low_high = stretchlim(img1,[0.01 0.99]);
img1 = 255*imadjust(img1,low_high,[0 1],1);
[m01 n01] = size(img1);

if (pc == 1)
    filename = [get(handles.edit1,'String') '\' get(handles.edit2,'String') '\' get(handles.edit2,'String') '_' num2str(sprintf('%03.0f',I2)) '.tif'];
elseif (pc == 0)
    filename = [get(handles.edit1,'String') '/' get(handles.edit2,'String') '/' get(handles.edit2,'String') '_' num2str(sprintf('%03d',I2)) '.tif'];
end
img2 = double(imread(filename,'tif'));
%mx = max(img2(:));
img2 = img2/mx;
%low_high = stretchlim(img2,[0.01 0.99]);
img2 = 255*imadjust(img2,low_high,[0 1],1);
[m02 n02] = size(img2);
test0 = zeros((m01 + 2*m02 - 2),(n01 + 2*n01 - 2),3);
test0(m02:(m02+m01-1),n02:(n02+n01-1),3) = img1;
test0((m02+dy):(m02+dy+m02-1),(n02+dx):(n02+dx+n02-1),1) = img2;
clear img1 img2

mx = max(test0(:));
mn = min(test0(:));
test0 = (test0 - mn)/(mx - mn)*255;
axes(handles.axes1)
imshow(uint8(test0))

set(handles.slider1,'Min',-m01+1);
set(handles.slider1,'Max',m01-1);
set(handles.slider1,'Value',0);
set(handles.slider1,'SliderStep',[1/(2*m01-2) 10/(2*m01-2)]);
set(handles.slider2,'Min',-n01+1);
set(handles.slider2,'Max',n01-1);
set(handles.slider2,'Value',0);
set(handles.slider2,'SliderStep',[1/(2*n01-2) 10/(2*n01-2)]); 


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

tmp = regexp(get(handles.edit1,'String'),'\\');
pc = 1;
if (isempty(tmp))
    pc = 0;
end
if (pc == 1)
    filename = [get(handles.edit1,'String') '\' get(handles.edit2,'String') '\offset_values.csv'];
elseif (pc == 0)
    filename = [get(handles.edit1,'String') '/' get(handles.edit2,'String') '/offset_values.csv'];
end
rel_reg = double(csvread(filename));

i = str2num(get(handles.edit5,'String'));
I1 = rel_reg(i,1);
I2 = rel_reg(i+1,1);
dy = rel_reg(i+1,2);
dx = rel_reg(i+1,3);
deltay = -round(get(handles.slider1,'Value'));
deltax = round(get(handles.slider2,'Value'));

if (pc == 1)
    filename = [get(handles.edit1,'String') '\' get(handles.edit2,'String') '\' get(handles.edit2,'String') '_' num2str(sprintf('%03.0f',I1)) '.tif'];
elseif (pc == 0)
    filename = [get(handles.edit1,'String') '/' get(handles.edit2,'String') '/' get(handles.edit2,'String') '_' num2str(sprintf('%03d',I1)) '.tif'];
end
img1 = double(imread(filename,'tif'));
mx = max(img1(:));
img1 = img1/mx;
low_high = stretchlim(img1,[0.01 0.99]);
img1 = 255*imadjust(img1,low_high,[0 1],1);
[m01 n01] = size(img1);
     
if (pc == 1)
    filename = [get(handles.edit1,'String') '\' get(handles.edit2,'String') '\' get(handles.edit2,'String') '_' num2str(sprintf('%03.0f',I2)) '.tif'];
elseif (pc == 0)
    filename = [get(handles.edit1,'String') '/' get(handles.edit2,'String') '/' get(handles.edit2,'String') '_' num2str(sprintf('%03d',I2)) '.tif'];
end
img2 = double(imread(filename,'tif'));
%mx = max(img2(:));
img2 = img2/mx;
%low_high = stretchlim(img2,[0.01 0.99]);
img2 = 255*imadjust(img2,low_high,[0 1],1);
[m02 n02] = size(img2);
test0 = zeros((m01 + 2*m02 - 2),(n01 + 2*n01 - 2),3);
test0(m02:(m02+m01-1),n02:(n02+n01-1),3) = img1;
test0((m02+dy+deltay):(m02+dy+deltay+m02-1),(n02+dx+deltax):(n02+dx+deltax+n02-1),1) = img2;
clear img1 img2

mx = max(test0(:));
mn = min(test0(:));
test0 = (test0 - mn)/(mx - mn)*255;
axes(handles.axes1)
imshow(uint8(test0))


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

tmp = regexp(get(handles.edit1,'String'),'\\');
pc = 1;
if (isempty(tmp))
    pc = 0;
end
if (pc == 1)
    filename = [get(handles.edit1,'String') '\' get(handles.edit2,'String') '\offset_values.csv'];
elseif (pc == 0)
    filename = [get(handles.edit1,'String') '/' get(handles.edit2,'String') '/offset_values.csv'];
end
rel_reg = double(csvread(filename));

i = str2num(get(handles.edit5,'String'));
I1 = rel_reg(i,1);
I2 = rel_reg(i+1,1);
dy = rel_reg(i+1,2);
dx = rel_reg(i+1,3);
deltay = -round(get(handles.slider1,'Value'));
deltax = round(get(handles.slider2,'Value'));

if (pc == 1)
    filename = [get(handles.edit1,'String') '\' get(handles.edit2,'String') '\' get(handles.edit2,'String') '_' num2str(sprintf('%03.0f',I1)) '.tif'];
elseif (pc == 0)
    filename = [get(handles.edit1,'String') '/' get(handles.edit2,'String') '/' get(handles.edit2,'String') '_' num2str(sprintf('%03d',I1)) '.tif'];
end
img1 = double(imread(filename,'tif'));
mx = max(img1(:));
img1 = img1/mx;
low_high = stretchlim(img1,[0.01 0.99]);
img1 = 255*imadjust(img1,low_high,[0 1],1);
[m01 n01] = size(img1);
     
if (pc == 1)
    filename = [get(handles.edit1,'String') '\' get(handles.edit2,'String') '\' get(handles.edit2,'String') '_' num2str(sprintf('%03.0f',I2)) '.tif'];
elseif (pc == 0)
    filename = [get(handles.edit1,'String') '/' get(handles.edit2,'String') '/' get(handles.edit2,'String') '_' num2str(sprintf('%03d',I2)) '.tif'];
end
img2 = double(imread(filename,'tif'));
%mx = max(img2(:));
img2 = img2/mx;
%low_high = stretchlim(img2,[0.01 0.99]);
img2 = 255*imadjust(img2,low_high,[0 1],1);
[m02 n02] = size(img2);
test0 = zeros((m01 + 2*m02 - 2),(n01 + 2*n01 - 2),3);
test0(m02:(m02+m01-1),n02:(n02+n01-1),3) = img1;
test0((m02+dy+deltay):(m02+dy+deltay+m02-1),(n02+dx+deltax):(n02+dx+deltax+n02-1),1) = img2;
clear img1 img2

mx = max(test0(:));
mn = min(test0(:));
test0 = (test0 - mn)/(mx - mn)*255;
axes(handles.axes1)
imshow(uint8(test0))


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
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
