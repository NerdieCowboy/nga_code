function varargout = HSIRegister(varargin)
% HSIREGISTER MATLAB code for HSIRegister.fig
%      HSIREGISTER, by itself, creates a new HSIREGISTER or raises the existing
%      singleton*.
%
%      H = HSIREGISTER returns the handle to a new HSIREGISTER or the handle to
%      the existing singleton*.
%
%      HSIREGISTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HSIREGISTER.M with the given input arguments.
%
%      HSIREGISTER('Property','Value',...) creates a new HSIREGISTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HSIRegister_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HSIRegister_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HSIRegister

% Last Modified by GUIDE v2.5 15-Feb-2014 17:38:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HSIRegister_OpeningFcn, ...
                   'gui_OutputFcn',  @HSIRegister_OutputFcn, ...
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


% --- Executes just before HSIRegister is made visible.
function HSIRegister_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HSIRegister (see VARARGIN)

% Choose default command line output for HSIRegister
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes HSIRegister wait for user response (see UIRESUME)
% uiwait(handles.figure1);
pc = 1;
set(handles.pc,'String',num2str(pc));
%axes(handles.axes1)
%imshow([])
%ncores = getenv('NUMBER_OF_PROCESSORS');
ncores = feature('NumCores');
set(handles.ncores,'String',num2str(ncores));
warning('off','all')

% --- Outputs from this function are returned to the command line.
function varargout = HSIRegister_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

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

% --- Executes on button press in pushbutton4.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fn,path] = uigetfile('*','Select the first cube in the set.');
set(handles.fn_path,'String',path);
ind = regexp(fn, '_[0-9]$');
set(handles.trial,'String',fn(1:(ind-1)));
set(handles.N1,'String',fn((ind+1):end));

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

function test_band_Callback(hObject, eventdata, handles)
% hObject    handle to test_band (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of test_band as text
%        str2double(get(hObject,'String')) returns contents of test_band as a double

% --- Executes during object creation, after setting all properties.
function test_band_CreateFcn(hObject, eventdata, handles)
% hObject    handle to test_band (see GCBO)
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

% --- Executes on button press in pushbutton4.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.running,'String','Running ...');
drawnow;
pause(0.1)

set(handles.stop,'Value',0);
rgb_fn = get(handles.rgb_fn,'String');
sc = str2double(get(handles.sc,'String'));
j01 = str2double(get(handles.j01,'String'));
j02 = str2double(get(handles.j02,'String'));
repcnt_max = str2double(get(handles.repcnt_max,'String'));
max_shift = str2double(get(handles.max_shift,'String'));
minpts = str2double(get(handles.minpts,'String'));
scan_mode = get(handles.scan_mode,'Value');
easel_mode = get(handles.easel_mode,'Value');
if (scan_mode == 1)
    Npoly = 3;
    Ntran = 3;
elseif (easel_mode == 1)
    Npoly = 0;
    Ntran = 0;
end
memory_limited = get(handles.memory_limited,'Value');
morepts = get(handles.morepts,'Value');
useresults = get(handles.useresults,'Value');
filteron = get(handles.filteron,'Value');
addsharp = get(handles.addsharp,'Value');
debugmode = get(handles.debugmode,'Value');
usegui = 2;
fn_path = get(handles.fn_path,'String');
trial = get(handles.trial,'String');
test_band = str2double(get(handles.test_band,'String'));
ir_fn = [fn_path trial '__' num2str(sprintf('%03.0f',test_band)) '.img'];
feather_width = str2double(get(handles.feather_width,'String'));
mxd = str2double(get(handles.mxd,'String'));
art_register(rgb_fn,sc,ir_fn,j01,j02,repcnt_max,max_shift,minpts,Npoly,Ntran,memory_limited,morepts,useresults,filteron,addsharp,debugmode,usegui,feather_width,mxd)
fn_blknum = [fn_path trial '__' num2str(sprintf('%03.0f',test_band)) '_blk.tif'];
blknum = uint16(imread(fn_blknum,'tif'));

fn = [fn_path trial '__' num2str(sprintf('%03.0f',test_band)) '_IR.tif'];
tmp = imread(fn,'tif');
[m2 n2] = size(tmp);
clear tmp

pc = str2double(get(handles.pc,'String'));
if (pc == 1)
    init_reg = [fn_path trial '__' num2str(sprintf('%03.0f',test_band)) '\init_reg.csv'];
    Xstar = [fn_path trial '__' num2str(sprintf('%03.0f',test_band)) '\Xstar.csv'];
elseif (pc == 0)
    init_reg = [fn_path trial '__' num2str(sprintf('%03.0f',test_band)) '/init_reg.csv'];
    Xstar = [fn_path trial '__' num2str(sprintf('%03.0f',test_band)) '/Xstar.csv'];
end
isOpen = matlabpool('size') > 0;
if (isOpen == 1)
    matlabpool close
end
N1 = str2double(get(handles.N1,'String'));
fnh = [fn_path trial '_' num2str(N1) '.hdr'];
fid = fopen(fnh);
while ~feof(fid)
    line = fgetl(fid);
    msk = isspace(line);
    line(msk==1) = '';
    [~,~,e] = regexp(line,'^bands=','match','start','end');
    if (e>0)
        p0 = str2double(line((e+1):end));
    end
end
fclose(fid);

set(handles.text17,'String','');
set(handles.text18,'String','');
set(handles.text19,'String','');
set(handles.running,'String','Running ...');
set(handles.text17,'Visible','off');
set(handles.text18,'Visible','off');
set(handles.text19,'Visible','off');
drawnow;
pause(0.1)

fn = [fn_path trial '_full.hdr'];
fid1 = fopen(fn,'w');
fid = fopen(fnh,'r');
while ~feof(fid)
    test = 0;
    line0 = fgetl(fid);
    line = line0;
    msk = isspace(line);
    line(msk==1) = '';
    [~,~,e] = regexp(line,'^samples=','match','start','end');
    if (e>0)
        test = 1;
        fprintf(fid1,'samples = %u\n',n2);
    end
    [~,~,e] = regexp(line,'^lines=','match','start','end');
    if (e>0)
        test = 1;
        fprintf(fid1,'lines = %u\n',m2);
    end
    [~,~,e] = regexp(line,'^bands=','match','start','end');
    if (e>0)
        test = 1;
        fprintf(fid1,'bands = %u\n',p0);
    end
    [~,~,e] = regexp(line,'^datatype=','match','start','end');
    if (e>0)
        test = 1;
        datatype0 = str2double(line((e+1):end));
        if ((datatype0 == 4) || (datatype0 == 5))
            fprintf(fid1,'data type = 4\n');
        else
            fprintf(fid1,'data type = 12\n');
        end
    end
    [~,~,e] = regexp(line,'^interleave=','match','start','end');
    if (e>0)
        test = 1;
        fprintf(fid1,'interleave = bsq\n');
    end
    [~,~,e] = regexp(line,'^byteorder=','match','start','end');
    if (e>0)
        test = 1;
        fprintf(fid1,'byte order = 0\n');
    end
    [~,~,e] = regexp(line,'^headeroffset=','match','start','end');
    if (e>0)
        test = 1;
        fprintf(fid1,'header offset = 0\n');
    end
    [~,~,e] = regexp(line,'^xstart=','match','start','end');
    if (e>0)
        test = 1;
    end
    [~,~,e] = regexp(line,'^ystart=','match','start','end');
    if (e>0)
        test = 1;
    end
    if (test == 0)
        fprintf(fid1,[line0 '\n']);
    end
end
fclose(fid);
fclose(fid1);

ncores = str2double(get(handles.ncores,'String'));
matlabpool(ncores)

fn = [fn_path trial '__' num2str(sprintf('%03.0f',test_band)) '_IR.tif'];
tmp = imread(fn,'tif');
fn = [fn_path trial '_full'];
if ((datatype0 == 4) || (datatype0 == 5))
    mxX = str2double(get(handles.mxX,'String'));
    mnX = str2double(get(handles.mnX,'String'));
    multibandwrite(((single(tmp)/(2^16-1))*(mxX-mnX))+mnX,fn,'bsq',[1,1,test_band],[m2,n2,p0],'machfmt','ieee-le')
else
    mxX = 2^16-1;
    mnX = 0;
    multibandwrite(uint16(tmp),fn,'bsq',[1,1,test_band],[m2,n2,p0],'machfmt','ieee-le')
end

stopcomm = get(handles.stop,'Value');
if (stopcomm == 1)
    set(handles.stop,'Value',0);
    error('Stop command issued')
end

parfor i = 1:p0
%for i = 1:p0
    if (i ~= test_band)
        fn = [fn_path trial '__' num2str(sprintf('%03.0f',i)) '.tif'];
        apply_transform(fn,blknum,m2,n2,Ntran,init_reg,Xstar,feather_width,0)
        fn = [fn_path trial '__' num2str(sprintf('%03.0f',i)) '1.tif'];
        tmp = imread(fn,'tif');
        fn_full = [fn_path trial '_full'];
        if ((datatype0 == 4) || (datatype0 == 5))
            multibandwrite(((single(tmp)/(2^16-1))*(mxX-mnX))+mnX,fn_full,'bsq',[1,1,i],[m2,n2,p0],'machfmt','ieee-le')
        else
            multibandwrite(uint16(tmp),fn_full,'bsq',[1,1,i],[m2,n2,p0],'machfmt','ieee-le')
        end
        %delete(fn)
        %folderName = [fn_path trial '__' num2str(sprintf('%03.0f',i))];
        %rmdir(folderName)
        %feval(@status_plot,p0);
    end
end
matlabpool close
clear tmp

set(handles.running,'String','');
drawnow;
pause(0.1)

% --- Executes on button press in prep.
function prep_Callback(hObject, eventdata, handles)
% hObject    handle to prep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.running,'String','Running ...');
set(handles.text17,'String','');
set(handles.text18,'String','');
set(handles.text19,'String','');
drawnow;
pause(0.1)

set(handles.stop,'Value',0);
fn_path = get(handles.fn_path,'String');
trial = get(handles.trial,'String');
N1 = str2double(get(handles.N1,'String'));
N2 = str2double(get(handles.N2,'String'));
test_band = str2double(get(handles.test_band,'String'));
ncores = str2double(get(handles.ncores,'String'));
pc = str2double(get(handles.pc,'String'));
fnh = [fn_path trial '_' num2str(N1) '.hdr'];
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
    [~,~,e] = regexp(line,'^datatype=','match','start','end');
    if (e>0)
        datatype0 = str2double(line((e+1):end));
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
        byteorder = str2double(line((e+1):end));
        if (byteorder == 0)
            byteorder = 'ieee-le';
        elseif (byteorder == 1)
            byteorder = 'ieee-be';
        end
    end
    [~,~,e] = regexp(line,'^headeroffset=','match','start','end');
    if (e>0)
        headeroffset = str2double(line((e+1):end));
    end
end
fclose(fid);

for i = 1:p0
    mkdir([fn_path trial '__' num2str(sprintf('%03.0f',i))]);
end
isOpen = matlabpool('size') > 0;
if (isOpen == 1)
    matlabpool close
end    
%matlabpool(ncores)

mxX = 0;
mnX = 2^16-1;
if ((datatype0 == 4) || (datatype0 == 5))
    for j = N1:N2
        fn = [fn_path trial '_' num2str(j)];
        if (exist(fn,'file') == 2)
            cube = multibandread(fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder);
            cube(isinf(cube)) = 0;
            cube(isnan(cube)) = 0;
            cube(cube<0) = 0;
            cube(cube>10) = 10;
            tmp = max(cube(:));
            if (tmp > mxX)
                mxX = tmp;
            end
            clear tmp
            tmp = min(cube(:));
            if (tmp < mnX)
                mnX = tmp;
            end
            clear tmp cube
        end
    end
    set(handles.mxX,'String',num2str(mxX));
    set(handles.mnX,'String',num2str(mnX));
end

for j = N1:N2
    fn = [fn_path trial '_' num2str(j)];
    if (exist(fn,'file') == 2)
        cube = multibandread(fn,[m0,n0,p0],datatype,headeroffset,interleave,byteorder);
        %parfor i = 1:p0
        for i = 1:p0
            if (pc == 1)
                if (datatype0 == 1)
                    imwrite(uint8(cube(:,:,i)),[fn_path trial '__' num2str(sprintf('%03.0f',i)) '\' trial '__' num2str(sprintf('%03.0f',i)) '_' num2str(sprintf('%03.0f',j)) '.tif'],'tif','Compression','None')
                elseif (datatype0 == 12)
                    imwrite(uint16(cube(:,:,i)),[fn_path trial '__' num2str(sprintf('%03.0f',i)) '\' trial '__' num2str(sprintf('%03.0f',i)) '_' num2str(sprintf('%03.0f',j)) '.tif'],'tif','Compression','None')
                elseif ((datatype0 == 4)  || (datatype0 == 5))
                    imwrite(uint16(((cube(:,:,i)-mnX)/(mxX-mnX))*(2^16-1)),[fn_path trial '__' num2str(sprintf('%03.0f',i)) '\' trial '__' num2str(sprintf('%03.0f',i)) '_' num2str(sprintf('%03.0f',j)) '.tif'],'tif','Compression','None')
                else
                    imwrite(uint16((2^14-1)*cube(:,:,i)),[fn_path trial '__' num2str(sprintf('%03.0f',i)) '\' trial '__' num2str(sprintf('%03.0f',i)) '_' num2str(sprintf('%03.0f',j)) '.tif'],'tif','Compression','None')
                end
            elseif (pc == 0)
                if (datatype0 == 1)
                    imwrite(uint8(cube(:,:,i)),[fn_path trial '__' num2str(sprintf('%03.0f',i)) '/' trial '__' num2str(sprintf('%03.0f',i)) '_' num2str(sprintf('%03.0f',j)) '.tif'],'tif','Compression','None')
                elseif (datatype0 == 12)
                    imwrite(uint16(cube(:,:,i)),[fn_path trial '__' num2str(sprintf('%03.0f',i)) '/' trial '__' num2str(sprintf('%03.0f',i)) '_' num2str(sprintf('%03.0f',j)) '.tif'],'tif','Compression','None')
                elseif ((datatype0 == 4) || (datatype0 == 5))
                    imwrite(uint16(((cube(:,:,i)-mnX)/(mxX-mnX))*(2^16-1)),[fn_path trial '__' num2str(sprintf('%03.0f',i)) '/' trial '__' num2str(sprintf('%03.0f',i)) '_' num2str(sprintf('%03.0f',j)) '.tif'],'tif','Compression','None')
                else
                    imwrite(uint16((2^14-1)*cube(:,:,i)),[fn_path trial '__' num2str(sprintf('%03.0f',i)) '/' trial '__' num2str(sprintf('%03.0f',i)) '_' num2str(sprintf('%03.0f',j)) '.tif'],'tif','Compression','None')
                end
            end
        end
        clear cube
    end
end
%matlabpool close
set(handles.running,'String','');
drawnow;
pause(0.1)

usegui = 2;
ir_fn1 = [fn_path trial '__' num2str(sprintf('%03.0f',test_band)) '.img'];
rough_mosaic(ir_fn1,usegui)

% --- Executes on button press in stop.
function stop_Callback(hObject, eventdata, handles)
% hObject    handle to stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of stop
set(handles.running,'String','');

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



function mxX_Callback(hObject, eventdata, handles)
% hObject    handle to mxX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mxX as text
%        str2double(get(hObject,'String')) returns contents of mxX as a double


% --- Executes during object creation, after setting all properties.
function mxX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mxX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mnX_Callback(hObject, eventdata, handles)
% hObject    handle to mnX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mnX as text
%        str2double(get(hObject,'String')) returns contents of mnX as a double


% --- Executes during object creation, after setting all properties.
function mnX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mnX (see GCBO)
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
