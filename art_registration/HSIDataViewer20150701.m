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

% Last Modified by GUIDE v2.5 08-Jun-2014 14:47:02

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

if (pc == 0)
    fnt = '/Applications/HSIDataViewer/application/HSIDataViewer.txt';
elseif (pc == 1)
    fnt = 'HSIDataViewer.txt';
end
if (exist(fnt,'file') == 2)
    fid = fopen(fnt);
    while ~feof(fid)
        line = fgetl(fid);
        msk = isspace(line);
        line(msk==1) = '';
        [~,~,e] = regexp(line,'^firstband=','match','start','end');
        if (e>0)
            p3 = str2double(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^lastband=','match','start','end');
        if (e>0)
            p2 = str2double(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^a0=','match','start','end');
        if (e>0)
            a0 = str2double(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^a1=','match','start','end');
        if (e>0)
            a1 = str2double(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^a2=','match','start','end');
        if (e>0)
            a2 = str2double(line((e+1):end));
        end
        [~,~,e] = regexp(line,'^a3=','match','start','end');
        if (e>0)
            a3 = str2double(line((e+1):end));
        end
    end
    fclose(fid);
    set(handles.p3,'String',num2str(p3));
    set(handles.p2,'String',num2str(p2));
    set(handles.a0,'String',num2str(a0));
    set(handles.a1,'String',num2str(a1));
    set(handles.a2,'String',num2str(a2));
    set(handles.a3,'String',num2str(a3));
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
end
trial_out = img_path;
if (read_only == 1)
    trial_out = [trial_out trial];
elseif (ff == 1) 
    trial_out = [trial_out 'ff_' trial];
elseif (ff == 0)
    trial_out = [trial_out trial];
end   
n0 = str2double(get(handles.n0,'string'));
m0 = str2double(get(handles.m0,'string'));
p0 = str2double(get(handles.p0,'string'));
datatype = get(handles.datatype,'string');
interleave = get(handles.interleave,'string');
byteorder = get(handles.byteorder,'string');
headeroffset = str2double(get(handles.headeroffset,'string'));

p3 = str2num(get(handles.p3,'String')); %#ok<ST2NM>
p2 = str2num(get(handles.p2,'String')); %#ok<ST2NM>
a0 = str2num(get(handles.a0,'String')); %#ok<ST2NM>
a1 = str2num(get(handles.a1,'String')); %#ok<ST2NM>
a2 = str2num(get(handles.a2,'String'))*(10^(-3)); %#ok<ST2NM>
a3 = str2num(get(handles.a3,'String'))*(10^(-6)); %#ok<ST2NM>
lambda3 = a0 + a1*p3 + a2*p3^2 + a3*p3^3;
lambda2 = a0 + a1*p2 + a2*p2^2 + a3*p2^3;
set(handles.lambda3,'String',num2str(lambda3));
set(handles.lambda2,'String',num2str(lambda2));
set(handles.slider1,'Min',p3); %#%#ok<MSNU> ok<ST2NM>
set(handles.slider1,'Max',p2);
ss = 1/(p2-p3+1);
set(handles.slider1,'SliderStep',[ss 20*ss]);

p1 = round(get(handles.slider1,'Value'));
set(handles.p1,'string',num2str(p1));
lambda1 = a0 + a1*p1 + a2*p1^2 + a3*p1^3;
set(handles.lambda1,'string',num2str(lambda1));

cube = multibandread(trial_out,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Band','Direct',p1});

val2 = get(handles.tog_img,'UserData');
axes(handles.axes1)
imagesc(cube)
colormap gray
if (val2 == 1)
    axis image
end
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

% --- Executes on button press in rot.
function rot_Callback(hObject, eventdata, handles)
% hObject    handle to rot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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

vnir_on = get(handles.VNIR,'Value');
xnir_on = get(handles.xNIR,'Value');
nir_on = get(handles.NIR,'Value');
read_only = get(handles.read_only,'Value');
if (read_only == 1)
    [fn,path] = uigetfile({'*.*'},'File Selector');
elseif (vnir_on == 1)
    [fn,path] = uigetfile({'*.cube';'*.img';'*.bsq';'*.bil',;'*.bip'},'File Selector');
elseif (xnir_on == 1)
    [fn,path] = uigetfile({'*.img'},'File Selector');
elseif (nir_on == 1)
    [fn,path] = uigetfile({'*.cube';'*.img';'*.bsq';'*.bil',;'*.bip'},'File Selector');
end
set(handles.img_path,'String',path);
ind = regexp(fn, '\.');
set(handles.img_ext,'String',fn((ind+1):end));
fn(ind:end) = [];
set(handles.trial,'String',fn);

% --- Executes on button press in file_info2.
function file_info2_Callback(hObject, eventdata, handles)
% hObject    handle to file_info2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tmp = get(handles.img_path,'String');
vnir_on = get(handles.VNIR,'Value');
xnir_on = get(handles.xNIR,'Value');
nir_on = get(handles.NIR,'Value');
if (vnir_on == 1)
    [fn,path] = uigetfile({'*.cube';'*.img';'*.bsq';'*.bil',;'*.bip'},'File Selector',tmp);
elseif (xnir_on == 1)
    [fn,path] = uigetfile({'*.img'},'File Selector',tmp);
elseif (nir_on == 1)
    [fn,path] = uigetfile({'*.cube';'*.img';'*.bsq';'*.bil',;'*.bip'},'File Selector',tmp);
end
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
vnir_on = get(handles.VNIR,'Value');
xnir_on = get(handles.xNIR,'Value');
nir_on = get(handles.NIR,'Value');
if (vnir_on == 1)
    [fn,path] = uigetfile('*.dark',tmp);
elseif (xnir_on == 1)
    [fn,path] = uigetfile({'*.img'},'File Selector',tmp);
elseif (nir_on == 1)
    [fn,path] = uigetfile({'*'},'File Selector',tmp);
end
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

deltan = 128;
vnir_on = get(handles.VNIR,'Value');
xnir_on = get(handles.xNIR,'Value');
nir_on = get(handles.NIR,'Value');
if (xnir_on == 1)
    %header_sz = 6553856;    %words
    %header_sz = 3277056;    %5*1024*640 + 256 words 
elseif (nir_on == 1)

elseif (vnir_on == 1)

end
read_only = get(handles.read_only,'Value');
p3 = str2num(get(handles.p3,'String')); %#ok<ST2NM>
p2 = str2num(get(handles.p2,'String')); %#ok<ST2NM>
a0 = str2num(get(handles.a0,'String')); %#ok<ST2NM>
a1 = str2num(get(handles.a1,'String')); %#ok<ST2NM>
a2 = str2num(get(handles.a2,'String'))*(10^(-3)); %#ok<ST2NM>
a3 = str2num(get(handles.a3,'String'))*(10^(-6)); %#ok<ST2NM>
img_path = get(handles.img_path,'String');
trial = get(handles.trial,'String');
img_ext = get(handles.img_ext,'String');
fn_full = [img_path trial '.' img_ext];
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
elseif ((b_on == 0) && (w_on == 1))
    ff = 1;
    div_only = 1;
else
    ff = 0;
end
set(handles.ff,'Value',ff);
if (read_only == 1)
    trial_out = img_path;
    trial_out = [trial_out trial];
    fnh = [trial_out '.hdr'];
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
    set(handles.interleave,'string',interleave);
elseif (ff == 1)
    if (xnir_on == 1)
        trial_out = [img_path 'ff_' trial];
        fnh = [trial_out '.hdr'];
        fidw = fopen(w_full,'r');
        header = fread(fidw,256,'uint16',0,'ieee-le');
        p0 = header(2);
        m0 = header(3);
        nw = header(22);
        clear header
        header = fread(fidw,5*p0*m0,'uint16',0,'ieee-le');
        clear header
        tmpw = zeros(p0,m0,nw,'single');
        if (div_only == 0)
            fidd = fopen(b_full,'r');
            header = fread(fidd,256 + 5*p0*m0,'uint16',0,'ieee-le');
            nd = header(22);
            clear header
            tmpd = zeros(p0,m0,nd,'single');
        end
        for n = 1:nw
            tmpw(:,:,n) = single(fread(fidw,[p0,m0],'uint16',0,'ieee-le'));
            if (div_only == 0)
                tmpd(:,:,n) = single(fread(fidd,[p0,m0],'uint16',0,'ieee-le'));
            end
            %if (~isempty(tmpw))
                %tmpw = single(tmpw);
            %end
        end
        fclose(fidw);
        w_tot = mean(tmpw,3);
        clear tmpw
        %w_tot = w_tot/nw;
        if (div_only == 0)
            fclose(fidd);
            %d_tot = d_tot/nd;
            d_tot = mean(tmpd,3);
            clear tmpd
        end

        fid = fopen(fn_full,'r');
        header = fread(fid,256 + 5*p0*m0,'uint16',0,'ieee-le');
        n0 = header(22);
        clear header
        cube = zeros(p0,m0,deltan,'single');
        n1 = 0;
        for n = 3600:n0
            n
            tmp = fread(fid,[p0,m0],'uint16',0,'ieee-le');
            n1 = n1 + 1;
            %if (~isempty(tmp))
                tmp = single(tmp);
                if (div_only == 0)
                    %tmp = ((tmp - d_tot)./(w_tot - d_tot))*(2^14-1);
                    tmp = (tmp - d_tot)./(w_tot - d_tot);
                elseif (div_only == 1)
                    %tmp = (tmp./w_tot)*(2^14-1);
                    tmp = tmp./w_tot;
                end
                %cube(:,:,n1) = uint16(tmp);
                cube(:,:,n1) = tmp;
                clear tmp
                if ((n1 == deltan) || (n == n0))
                    tmp = shiftdim(cube(:,:,1:n1),1);
                    tmp = tmp(:,:,p0:-1:1);
                    if (get(handles.bip_button,'value') == 1)
                        interleave1 = 'bip';
                    elseif (get(handles.bil_button,'value') == 1)
                        interleave1 = 'bil';
                    elseif (get(handles.bsq_button,'value') == 1)
                        interleave1 = 'bsq';
                    end
                    %multibandwrite(uint16(tmp),trial_out,interleave1,[1 (n-n1+1) 1],[m0 n0 p0],'machfmt','ieee-le')
                    multibandwrite(single(tmp),trial_out,interleave1,[1 (n-n1+1) 1],[m0 n0 p0],'machfmt','ieee-le')
                    clear tmp
                    n1 = 0;
                end
            %end
        end
        fclose(fid);
        clear cube w_tot
        if (div_only == 0)
            clear d_tot
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
        for i = 1:(p0-1)
            lambda1 = a0 + a1*i + a2*i^2 + a3*i^3;
            fprintf(fid1,'%f, ',lambda1);
        end
        fprintf(fid1,'%f}\n',a0 + a1*p0 + a2*p0^2 + a3*p0^3);
        fclose(fid1);
    elseif (nir_on == 1)
        
    elseif (vnir_on == 1)
        
    end
    if (ff == 0)
        datatype = 'uint16';
    else
        datatype = 'single';
    end
    if (get(handles.bip_button,'value') == 1)
        interleave1 = 'bip';
    elseif (get(handles.bil_button,'value') == 1)
        interleave1 = 'bil';
    elseif (get(handles.bsq_button,'value') == 1)
        interleave1 = 'bsq';
    end
    byteorder = 'ieee-le';
    headeroffset = 0;
    set(handles.interleave,'string',interleave1);
elseif (ff == 0)
    if (xnir_on == 1)
        trial_out = [img_path trial];
        fnh = [trial_out '.hdr'];
        fid = fopen(fn_full,'r');
        header = fread(fid,256,'uint16',0,'ieee-le');
        p0 = header(2);
        m0 = header(3);
        n0 = header(22);
        clear header
        header = fread(fid,5*p0*m0,'uint16',0,'ieee-le');
        clear header
        cube = zeros(p0,m0,deltan,'uint16');
        n1 = 0;
        for n = 1:n0
            tmp = fread(fid,[p0,m0],'uint16',0,'ieee-le');
            n1 = n1 + 1;
            %if (~isempty(tmp))
                cube(:,:,n1) = uint16(tmp);
                clear tmp

                if ((n1 == deltan) || (n == n0))
                    tmp = shiftdim(cube(:,:,1:n1),1);
                    tmp = tmp(:,:,p0:-1:1);
                    if (get(handles.bip_button,'value') == 1)
                        interleave1 = 'bip';
                    elseif (get(handles.bil_button,'value') == 1)
                        interleave1 = 'bil';
                    elseif (get(handles.bsq_button,'value') == 1)
                        interleave1 = 'bsq';
                    end
                    multibandwrite(uint16(tmp),trial_out,interleave1,[1 (n-n1+1) 1],[m0 n0 p0],'machfmt','ieee-le')
                    clear tmp
                    n1 = 0;
                end
            %end
        end
        fclose(fid);
        clear cube

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
        for i = 1:(p0-1)
            lambda1 = a0 + a1*i + a2*i^2 + a3*i^3;
            fprintf(fid1,'%f, ',lambda1);
        end
        fprintf(fid1,'%f}\n',a0 + a1*p0 + a2*p0^2 + a3*p0^3);
        fclose(fid1);
    elseif (nir_on == 1)
        
    elseif (vnir_on == 1)
        
    end
    if (ff == 0)
        datatype = 'uint16';
    else
        datatype = 'single';
    end
    if (get(handles.bip_button,'value') == 1)
        interleave1 = 'bip';
    elseif (get(handles.bil_button,'value') == 1)
        interleave1 = 'bil';
    elseif (get(handles.bsq_button,'value') == 1)
        interleave1 = 'bsq';
    end
    byteorder = 'ieee-le';
    headeroffset = 0;
    set(handles.interleave,'string',interleave1);
end
set(handles.n0,'string',num2str(n0));
set(handles.m0,'string',num2str(m0));
set(handles.p0,'string',num2str(p0));
set(handles.datatype,'string',datatype);
set(handles.byteorder,'string',byteorder);
set(handles.headeroffset,'string',headeroffset);

%set(handles.p3,'string',1);
%set(handles.p2,'string',p0);
%lambda3 = a0 + a1*1 + a2*1^2 + a3*1^3;
%lambda2 = a0 + a1*p0 + a2*p0^2 + a3*p0^3;
lambda3 = a0 + a1*p3 + a2*p3^2 + a3*p3^3;
lambda2 = a0 + a1*p2 + a2*p2^2 + a3*p2^3;
set(handles.lambda3,'String',num2str(lambda3));
set(handles.lambda2,'String',num2str(lambda2));
set(handles.slider1,'Enable','off');
set(handles.slider1,'Visible','on');
%set(handles.slider1,'Min',1); %#%#ok<MSNU> ok<ST2NM>
%set(handles.slider1,'Max',p0);
set(handles.slider1,'Min',p3); %#%#ok<MSNU> ok<ST2NM>
set(handles.slider1,'Max',p2);
%ss = 1/p0;
ss = 1/(p2-p3+1);
%set(handles.slider1,'Value',1);
set(handles.slider1,'Value',p3);
set(handles.slider1,'SliderStep',[ss 20*ss]);
set(handles.lambda1,'Visible','on');
set(handles.lambda1,'string',lambda3);
set(handles.p1,'Visible','off');
%set(handles.p1,'string',1);
set(handles.p1,'string',p3);
interleave = get(handles.interleave,'string');

axes(handles.axes1)
%tmp= multibandread(trial_out,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Band','Direct',1});
tmp= multibandread(trial_out,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Band','Direct',p3});

val2 = get(handles.tog_img,'UserData');
imagesc(tmp)
colormap gray
if (val2 == 1)
    axis image
end
axis off
clear tmp

%lambda = a0 + a1*(1:p0) + a2*((1:p0).^2) + a3*((1:p0).^3);
lambda = a0 + a1*(p3:p2) + a2*((p3:p2).^2) + a3*((p3:p2).^3);
axes(handles.axes2)
%plot(lambda,zeros(1,p0))
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

pathname = fileparts(which('HSIDataViewer.m'));
if (ismac == 0)
    fnt = [pathname '\HSIDataViewer.txt'];
elseif (ismac == 1)
    fnt = [pathname '/HSIDataViewer.txt'];
end
fid1 = fopen(fnt,'w');
fprintf(fid1,'first band = %u\n',p3);
fprintf(fid1,'last band = %u\n',p2);
fprintf(fid1,'a0 = %f\n',a0);
fprintf(fid1,'a1 = %f\n',a1);
fprintf(fid1,'a2 = %f\n',a2/(10^(-3)));
fprintf(fid1,'a3 = %f\n',a3/(10^(-6)));
fclose(fid1);

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

handles = guidata(HSIDataViewer);
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

    n0 = str2num(get(handles.n0,'String')); %#ok<ST2NM>
    m0 = str2num(get(handles.m0,'String')); %#ok<ST2NM>
    p0 = str2num(get(handles.p0,'String')); %#ok<ST2NM>
    datatype = get(handles.datatype,'string');
    interleave = get(handles.interleave,'string');
    byteorder = get(handles.byteorder,'string');
    headeroffset = str2double(get(handles.headeroffset,'string'));
    p2 = str2num(get(handles.p2,'String')); %#ok<ST2NM>
    p3 = str2num(get(handles.p3,'String')); %#ok<ST2NM>
    coordn = round(n0*x/pos(3));
    coordm = round(m0*y/pos(4));
    set(handles.x1,'String',num2str(coordn));
    set(handles.y1,'String',num2str(coordm));

    a0 = str2num(get(handles.a0,'String')); %#ok<ST2NM>
    a1 = str2num(get(handles.a1,'String')); %#ok<ST2NM>
    a2 = str2num(get(handles.a2,'String'))*(10^(-3)); %#ok<ST2NM>
    a3 = str2num(get(handles.a3,'String'))*(10^(-6)); %#ok<ST2NM>
    lambda3 = a0 + a1*p3 + a2*p3^2 + a3*p3^3;
    lambda2 = a0 + a1*p2 + a2*p2^2 + a3*p2^3;
    lambda = a0 + a1*(p3:p2) + a2*((p3:p2).^2) + a3*((p3:p2).^3);
    set(handles.lambda3,'string',num2str(lambda3));
    set(handles.lambda2,'string',num2str(lambda2));

    read_only = get(handles.read_only,'String');
    ff = get(handles.ff,'Value');
    trial_out = img_path;
    if (read_only == 1)
        trial_out = [trial_out trial];
    else
        if (ff == 1)
            trial_out = [trial_out 'ff_' trial];
        else
            trial_out = [trial_out trial];
        end
    end

    if (read_only == 1)
        cube = single(multibandread(trial_out,[m0,n0,p2-p3+1],datatype,headeroffset,interleave,byteorder,{'Row','Direct',coordm},{'Column','Direct',coordn},{'Band','Range',[p3,p2]}));
    else
        %cube = single(multibandread(trial_out,[m0,n0,p2-p3+1],datatype,headeroffset,interleave,byteorder,{'Row','Direct',coordm},{'Column','Direct',coordn},{'Band','Range',[p3,p2]}))/(2^14-1);
        cube = single(multibandread(trial_out,[m0,n0,p2-p3+1],datatype,headeroffset,interleave,byteorder,{'Row','Direct',coordm},{'Column','Direct',coordn},{'Band','Range',[p3,p2]}));
    end
    cube(cube<0) = 0;
    axes(handles.axes2)
    plot(lambda,cube(:));
    ylim([0 max(cube(:))]);
    xlabel('wavelength (nm)');
    ylabel('reflectance');
    xlim([lambda(1) lambda(end)]);
    clear cube
elseif ((x2>=1) && (x2<=pos2(3)) && (y2>=1) && (y2<=pos2(4)))
    set(handles.axes2,'Units','points');
    fig = figure(HSIDataViewer);
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
    
    fig = figure(HSIDataViewer);
    datacursormode on
    dcm_obj = datacursormode(fig);
    set(dcm_obj,'UpdateFcn',@myupdatefcn)
else
    set(handles.eid,'UserData',0);
    set(handles.eid,'String','Spatial pts');
    
    fig = figure(HSIDataViewer);
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
%        str2double(get(hObject,'String')) returns contents of n0 as a double


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
%        str2double(get(hObject,'String')) returns contents of datatype as a double


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
%        str2double(get(hObject,'String')) returns contents of byteorder as a double


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
%        str2double(get(hObject,'String')) returns contents of interleave as a double


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
%        str2double(get(hObject,'String')) returns contents of headeroffset as a double


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


% --- Executes on button press in tog_img.
function tog_img_Callback(hObject, eventdata, handles)
% hObject    handle to tog_img (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

read_only = get(handles.read_only,'Value');
w_path = get(handles.w_path,'String');
w_fn = get(handles.w_fn,'String');
w_ext = get(handles.w_ext,'String');
w_full = [w_path w_fn '.' w_ext];
if (exist(w_full,'file') == 2)
    w_on = 1;
else
    w_on = 0;
end
if (w_on == 1)
    ff = 1;
end
img_path = get(handles.img_path,'String');
trial = get(handles.trial,'String');
trial_out = img_path;
if (read_only == 1)
    trial_out = [trial_out trial];
elseif (ff == 1) 
    trial_out = [trial_out 'ff_' trial];
elseif (ff == 0)
    trial_out = [trial_out trial];
end
n0 = str2double(get(handles.n0,'string'));
m0 = str2double(get(handles.m0,'string'));
p0 = str2double(get(handles.p0,'string'));
datatype = get(handles.datatype,'string');
interleave = get(handles.interleave,'string');
byteorder = get(handles.byteorder,'string');
headeroffset = str2double(get(handles.headeroffset,'string'));

p1 = round(get(handles.slider1,'Value'));
val = get(handles.tog_img,'UserData');
if (val == 0)
    set(handles.tog_img,'UserData',1);
elseif (val == 1)
    set(handles.tog_img,'UserData',0);
end
val = get(handles.tog_img,'UserData');

if (val == 1)
    set(handles.eid,'enable','off');
    cube = multibandread(trial_out,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Band','Direct',p1});
    
    axes(handles.axes1)
    imagesc(cube)
    colormap gray
    axis image
    axis off
    clear cube
    
    fig = figure(HSIDataViewer);
    datacursormode on
    dcm_obj = datacursormode(fig);
    set(dcm_obj,'UpdateFcn',@myupdatefcn)
else
    set(handles.eid,'enable','on');
    cube = multibandread(trial_out,[m0,n0,p0],datatype,headeroffset,interleave,byteorder,{'Band','Direct',p1});

    val2 = get(handles.tog_img,'UserData');
    axes(handles.axes1)
    imagesc(cube)
    colormap gray
    if (val2 == 1)
        axis image
    end
    axis off
    clear cube
    
    fig = figure(HSIDataViewer);
    datacursormode off
end
