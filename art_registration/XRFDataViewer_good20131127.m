function varargout = XRFDataViewer(varargin)
% XRFDATAVIEWER MATLAB code for XRFDataViewer.fig
%      XRFDATAVIEWER, by itself, creates a new XRFDATAVIEWER or raises the existing
%      singleton*.
%
%      H = XRFDATAVIEWER returns the handle to a new XRFDATAVIEWER or the handle to
%      the existing singleton*.
%
%      XRFDATAVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in XRFDATAVIEWER.M with the given input arguments.
%
%      XRFDATAVIEWER('Property','Value',...) creates a new XRFDATAVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before XRFDataViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to XRFDataViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help XRFDataViewer

% Last Modified by GUIDE v2.5 14-Nov-2013 11:58:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @XRFDataViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @XRFDataViewer_OutputFcn, ...
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

% --- Executes just before XRFDataViewer is made visible.
function XRFDataViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to XRFDataViewer (see VARARGIN)

% Choose default command line output for XRFDataViewer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes XRFDataViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = XRFDataViewer_OutputFcn(hObject, eventdata, handles) 
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

p3 = str2num(get(handles.p3,'String')); %#ok<ST2NM>
p2 = str2num(get(handles.p2,'String')); %#ok<ST2NM>
set(handles.slider1,'Min',p3); %#%#ok<MSNU> ok<ST2NM>
set(handles.slider1,'Max',p2);
ss = 1/(p2-p3+1);
set(handles.slider1,'SliderStep',[ss 40*ss]);

slope = str2num(get(handles.slope,'String')); %#ok<ST2NM>
yint = str2num(get(handles.yint,'String')); %#ok<ST2NM>
p1 = round(get(handles.slider1,'Value'));
set(handles.p1,'string',num2str(p1));
energy1 = p1*slope + yint;
set(handles.energy1,'string',num2str(energy1));
e3 = p3*slope + yint;
e2 = p2*slope + yint;
set(handles.e3,'string',num2str(e3));
set(handles.e2,'string',num2str(e2));

n0 = str2num(get(handles.n0,'String')); %#ok<ST2NM>
m0 = str2num(get(handles.m0,'String')); %#ok<ST2NM>
%p0 = str2num(get(handles.p0,'String')); %#ok<ST2NM>

img_path = get(handles.img_path,'String');
trial = get(handles.trial,'String');
trial_out = [img_path trial];
cube = single(multibandread(trial_out,[m0,n0,p2-p3+1],'single',0,'bsq','ieee-le',{'Band','Direct',p1-p3+1}));

axes(handles.axes1)
imagesc(single(cube))
colormap gray
%axis image
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

function n0_Callback(hObject, eventdata, handles)
% hObject    handle to p0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p0 as text
%        str2double(get(hObject,'String')) returns contents of p0 as a double

% --- Executes during object creation, after setting all properties.
function n0_CreateFcn(hObject, eventdata, handles)
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

% Hints: get(hObject,'String') returns contents of m0 as text
%        str2double(get(hObject,'String')) returns contents of m0 as a double

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

function p0_Callback(hObject, eventdata, handles)
% hObject    handle to p0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p0 as text
%        str2double(get(hObject,'String')) returns contents of p0 as a double

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

function header_sz_Callback(hObject, eventdata, handles)
% hObject    handle to header_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of header_sz as text
%        str2double(get(hObject,'String')) returns contents of header_sz as a double

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

% --- Executes on button press in file_info.
function file_info_Callback(hObject, eventdata, handles)
% hObject    handle to file_info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

read_only = get(handles.read_only,'Value');
if (read_only == 0)
    [fn,path] = uigetfile('*.bin','Select the XRF files.');
else
    [fn,path] = uigetfile('*','Select the XRF files.');
end
set(handles.img_path,'String',path);
ind = regexp(fn, '\.');
set(handles.img_ext,'String',fn((ind+1):end));
fn(ind:end) = [];
set(handles.trial,'String',fn);

function img_path_Callback(hObject, eventdata, handles)
% hObject    handle to img_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of img_path as text
%        str2double(get(hObject,'String')) returns contents of img_path as a double

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

% --- Executes on button press in run_button.
function run_button_Callback(hObject, eventdata, handles)
% hObject    handle to run_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

p3 = str2num(get(handles.p3,'String')); %#ok<ST2NM>
p2 = str2num(get(handles.p2,'String')); %#ok<ST2NM>
set(handles.slider1,'Enable','off');
set(handles.slider1,'Visible','on');
set(handles.slider1,'Min',p3); %#%#ok<MSNU> ok<ST2NM>
set(handles.slider1,'Max',p2);
ss = 1/(p2-p3+1);
set(handles.slider1,'Value',p3);
set(handles.slider1,'SliderStep',[ss 40*ss]);
set(handles.slider1,'Enable','on');
set(handles.energy1,'Visible','on');
set(handles.energy1,'string',num2str(1));
set(handles.p1,'Visible','off');
set(handles.p1,'string',num2str(1));

img_path = get(handles.img_path,'String');
trial = get(handles.trial,'String');
img_ext = get(handles.img_ext,'String');

csv_path = get(handles.csv_path,'String');
csv_fn = get(handles.csv_fn,'String');
csv_ext = get(handles.csv_ext,'String');

n0 = str2num(get(handles.n0,'String')); %#ok<ST2NM>
m0 = str2num(get(handles.m0,'String')); %#ok<ST2NM>
p0 = str2num(get(handles.p0,'String')); %#ok<ST2NM>
pixel_sz = str2num(get(handles.pixel_sz,'String')); %#ok<ST2NM>
slope = str2num(get(handles.slope,'String')); %#ok<ST2NM>
yint = str2num(get(handles.yint,'String')); %#ok<ST2NM>
p2 = str2num(get(handles.p2,'String')); %#ok<ST2NM>
p3 = str2num(get(handles.p3,'String')); %#ok<ST2NM>
read_only = get(handles.read_only,'Value');
e3 = p3*slope + yint;
e2 = p2*slope + yint;
set(handles.e3,'string',num2str(e3));
set(handles.e2,'string',num2str(e2));

test = 2;

%p0 = p2 - p3 + 1;

if (read_only == 1)
    trial_out = [img_path trial];
    
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
    end
    fclose(fid);
    set(handles.n0,'string',num2str(n0));
    set(handles.m0,'string',num2str(m0));
    set(handles.p0,'string',num2str(p0));
    
    %cube2= single(multibandread(trial_out,[m0,n0,p0],'single',0,'bsq','ieee-le',{'Band','Direct',p3}));
elseif (read_only == 0)
    fncsv = [csv_path csv_fn '.' csv_ext];
    if (exist(fncsv,'file') == 2)
        tmp = csvread(fncsv,1,0);
        x = single(tmp(:,1));
        mnx = min(x);
        y = single(tmp(:,2));
        clear tmp
        mny = min(y);
        x = x - mnx;
        x = x/pixel_sz + 1.5;
        mxx = ceil(max(x)+0.5);
        y = y - mny;
        y = y/pixel_sz + 1.5;
        mxy = ceil(max(y)+0.5);
    end
    
    fn = [img_path '\' trial '.bin'];
    fid = fopen(fn,'r');
    dir = 0;
    cnt = 1;
    cube = single(zeros(m0,n0,(p2-p3+1)));
    if (exist(fncsv,'file') == 2)
        X = single(zeros(m0,n0));
        Y = single(zeros(m0,n0));
        if (test == 1)
            m = 1;
            for n = 1:n0
                X(m,n) = x(cnt) + 0.5;
                Y(m,n) = y(cnt);
                cnt = cnt + 1;
            end
            dir = 1;
            for m = 2:m0
                for n = 1:(n0+1)
                    if (dir == 0)
                        if (n ~= (n0+1))
                            X(m,n) = x(cnt) + 0.5;
                            Y(m,n) = y(cnt);
                        end
                    elseif (dir == 1)
                        if (n ~= (n0+1))
                            n1 = n0 - (n - 1);
                            X(m,n1) = x(cnt) - 0.5;
                            Y(m,n1) = y(cnt);
                        end
                    end
                    cnt = cnt + 1;
                end
                if (dir == 0)
                    dir = 1;
                elseif (dir == 1)
                    dir = 0;
                end
            end
        elseif (test == 2)
            for m = 1:m0
                for n = 1:(n0+1)
                    if (dir == 0)
                        if (n ~= (n0+1))
                            X(m,n) = x(cnt) + 0.5;
                            Y(m,n) = y(cnt);
                        end
                    elseif (dir == 1)
                        if (n ~= (n0+1))
                            n1 = n0 - (n - 1);
                            X(m,n1) = x(cnt) - 0.5;
                            Y(m,n1) = y(cnt);
                        end
                    end
                    cnt = cnt + 1;
                end
                if (dir == 0)
                    dir = 1;
                elseif (dir == 1)
                    dir = 0;
                end
            end
        end
        X = X(:,2:(n0-1));
        Y = Y(:,2:(n0-1));
    end
    
    for m = 1:m0
        for n = 1:n0
            tmp = single(fread(fid,p0,'uint32',0,'ieee-le'));
            tmp = tmp(p3:p2);

            if (dir == 0)
                if (n > 2)
                    cube(m,n,:) = tmp - tmp_prev;
                end
                if (n == 3)
                    cube(m,2,:) = cube(m,3,:);
                end
                %cube(m,n,:) = tmp;
            elseif (dir == 1)
                n1 = n0 - (n - 1);
                if (n1 < (n0 - 1))
                    cube(m,n1,:) = tmp - tmp_prev;
                end
                if (n1 == (n0 - 2))
                    cube(m,(n0-1),:) = cube(m,(n0-2),:);
                end
                %cube(m,n1,:) = tmp;
            end
            tmp_prev = tmp;
            clear tmp
        end
        if (dir == 0)
            dir = 1;
        elseif (dir == 1)
            dir = 0;
        end
    end
    cube = cube(:,2:(n0-1),:);
    n0 = n0 - 2;
    
    if (exist(fncsv,'file') == 2)
        [N,M] = meshgrid(1:mxx,1:mxy);
        [m0,n0] = size(N);
        %p0 = p2 - p3 + 1;
        cube2 = single(zeros(m0,n0,(p2-p3+1)));
        isOpen = matlabpool('size') > 0;
        if (isOpen == 1)
            matlabpool close
        end
        ncores = feature('numCores');
        matlabpool(ncores)
        parfor p = 1:(p2-p3+1)
            tmp = double(cube(:,:,p));
            F = TriScatteredInterp(double(X(:)),double(Y(:)),double(tmp(:)),'natural');
            cube2(:,:,p) = F(double(N),double(M));
            %tmp1 = F(double(N),double(M));
            %cube2(:,:,p) = tmp1(:,n0:-1:1);
        end
        clear tmp tmp1
        isOpen = matlabpool('size') > 0;
        if (isOpen == 1)
            matlabpool close
        end
        num_cols = n0;
        num_rows = m0;
        num_bands = p2-p3+1;
    end
        
    set(handles.n0,'string',num2str(n0));
    set(handles.m0,'string',num2str(m0));
    set(handles.p0,'string',num2str(p2-p3+1));
    
    trial_out = [img_path trial];
    fnh_out = [trial_out '.hdr'];
    
    if (exist(fncsv,'file') == 2)
        multibandwrite(single(cube2),trial_out,'bsq','machfmt','ieee-le')
    else
        multibandwrite(single(cube),trial_out,'bsq','machfmt','ieee-le')
    end

    fid1 = fopen(fnh_out,'w');
    fprintf(fid1,'ENVI\n');
    fprintf(fid1,'description = {}\n');
    fprintf(fid1,'samples = %u\n',n0);
    fprintf(fid1,'lines = %u\n',m0);
    fprintf(fid1,'bands = %u\n',p2-p3+1);
    fprintf(fid1,'header offset = 0\n');
    fprintf(fid1,'file type = ENVI Standard\n');
    fprintf(fid1,'data type = 4\n');
    fprintf(fid1,'interleave = bsq\n');
    fprintf(fid1,'byte order = 0\n');
    fprintf(fid1,'wavelength = {');
    for i = p3:(p2-1)
        lambda = i*slope+yint;
        fprintf(fid1,'%f, ',lambda);
    end
    fprintf(fid1,'%f}\n', p2*slope+yint);
    fclose(fid1);
end

axes(handles.axes1)
p1 = str2num(get(handles.p1,'string')); %#ok<ST2NM>
energy1 = p1*slope+yint;
set(handles.energy1,'string',num2str(energy1));
tmp= single(multibandread(trial_out,[m0,n0,p2-p3+1],'single',0,'bsq','ieee-le',{'Band','Direct',p1-p3+1}));
imagesc(single(tmp))
colormap gray
%axis image
axis off
clear tmp cube2

lambda = (p3:p2)*slope+yint;
axes(handles.axes2)
plot(lambda,zeros(1,(p2-p3+1)))
xlabel('energy');
ylabel('fluorescence');
xlim([lambda(1) lambda(end)]);
ylim([0 1]);

set(handles.axes1,'Visible','on');
set(handles.axes2,'Visible','on');

set(handles.axes1,'Units','pixels');
set(handles.axes2,'Units','points');

set(gcf,'WindowButtonDownFcn',@mouseMove);
set(gcf,'Units','pixels');
datacursormode off

function p1_Callback(hObject, eventdata, handles)
% hObject    handle to p1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p1 as text
%        str2double(get(hObject,'String')) returns contents of p1 as a double


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



function img_num_Callback(hObject, eventdata, handles)
% hObject    handle to img_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of img_num as text
%        str2double(get(hObject,'String')) returns contents of img_num as a double


% --- Executes during object creation, after setting all properties.
function img_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to img_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mouseMove (hobject, eventdata)

handles = guidata(XRFDataViewer);
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
    p2 = str2num(get(handles.p2,'String')); %#ok<ST2NM>
    p3 = str2num(get(handles.p3,'String')); %#ok<ST2NM>
    logon = get(handles.logon,'Value');
    coordn = round(n0*x/pos(3));
    coordm = round(m0*y/pos(4));
    set(handles.x1,'String',num2str(coordn));
    set(handles.y1,'String',num2str(coordm));

    slope = str2num(get(handles.slope,'String')); %#ok<ST2NM>
    yint = str2num(get(handles.yint,'String')); %#ok<ST2NM>
    lambda = (p3:p2)*slope+yint;
    e3 = p3*slope + yint;
    e2 = p2*slope + yint;
    set(handles.e3,'string',num2str(e3));
    set(handles.e2,'string',num2str(e2));


    trial_out = [img_path trial];
    cube = single(multibandread(trial_out,[m0,n0,p2-p3+1],'single',0,'bsq','ieee-le',{'Row','Range',[y,y]},{'Column','Range',[x,x]},{'Band','Range',[1,p2-p3+1]}));
    axes(handles.axes2)
    if (logon == 1)
        plot(lambda,log(cube(:)));
        ylim([0 max(log(cube(:)))]);
    else
        plot(lambda,cube(:));
        ylim([0 max(cube(:))]);
    end
    xlabel('energy (KeV)');
    ylabel('fluorescence');
    xlim([lambda(1) lambda(end)]);
elseif ((x2>=1) && (x2<=pos2(3)) && (y2>=1) && (y2<=pos2(4)))
    set(handles.axes2,'Units','points');
    fig = figure(XRFDataViewer);
    %datacursormode on
    %dcm_obj = datacursormode(fig);
    %set(dcm_obj,'UpdateFcn',@myupdatefcn)
%else
    %datacursormode(XRFDataViewer)
end

function txt = myupdatefcn(empt,event_obj)
% Customizes text of data tips

ptpos = get(event_obj,'Position');
logon = get(handles.logon,'Value');
if (logon == 1)
    txt = {['Energy: ',num2str(ptpos(1)),' KeV'],...
           ['Counts: ',num2str(round(exp(ptpos(2))))]};
else
    txt = {['Energy: ',num2str(ptpos(1)),' KeV'],...
           ['Counts: ',num2str(round(ptpos(2)))]};
end

function slope_Callback(hObject, eventdata, handles)
% hObject    handle to slope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of slope as text
%        str2double(get(hObject,'String')) returns contents of slope as a double


% --- Executes during object creation, after setting all properties.
function slope_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yint_Callback(hObject, eventdata, handles)
% hObject    handle to yint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yint as text
%        str2double(get(hObject,'String')) returns contents of yint as a double


% --- Executes during object creation, after setting all properties.
function yint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yint (see GCBO)
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

% Hints: get(hObject,'String') returns contents of p2 as text
%        str2double(get(hObject,'String')) returns contents of p2 as a double


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

% Hints: get(hObject,'String') returns contents of p3 as text
%        str2double(get(hObject,'String')) returns contents of p3 as a double


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


% --- Executes on button press in logon.
function logon_Callback(hObject, eventdata, handles)
% hObject    handle to logon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of logon


% --- Executes on button press in read_only.
function read_only_Callback(hObject, eventdata, handles)
% hObject    handle to read_only (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of read_only


% --- Executes on button press in csv_info.
function csv_info_Callback(hObject, eventdata, handles)
% hObject    handle to csv_info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fn,path] = uigetfile('*.csv','Select the CSV files.');
set(handles.csv_path,'String',path);
ind = regexp(fn, '\.');
set(handles.csv_ext,'String',fn((ind+1):end));
fn(ind:end) = [];
set(handles.csv_fn,'String',fn);


function csv_path_Callback(hObject, eventdata, handles)
% hObject    handle to csv_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of csv_path as text
%        str2double(get(hObject,'String')) returns contents of csv_path as a double


% --- Executes during object creation, after setting all properties.
function csv_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to csv_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function csv_fn_Callback(hObject, eventdata, handles)
% hObject    handle to csv_fn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of csv_fn as text
%        str2double(get(hObject,'String')) returns contents of csv_fn as a double


% --- Executes during object creation, after setting all properties.
function csv_fn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to csv_fn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function csv_ext_Callback(hObject, eventdata, handles)
% hObject    handle to csv_ext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of csv_ext as text
%        str2double(get(hObject,'String')) returns contents of csv_ext as a double


% --- Executes during object creation, after setting all properties.
function csv_ext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to csv_ext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pixel_sz_Callback(hObject, eventdata, handles)
% hObject    handle to pixel_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pixel_sz as text
%        str2double(get(hObject,'String')) returns contents of pixel_sz as a double


% --- Executes during object creation, after setting all properties.
function pixel_sz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pixel_sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in eid.
function eid_Callback(hObject, eventdata, handles)
% hObject    handle to eid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of eid

eid_state = get(handles.eid,'UserData');
if (eid_state == 0)
    set(handles.eid,'UserData',1);
    set(handles.eid,'String','Energy');
    
    fig = figure(XRFDataViewer);
    datacursormode on
    dcm_obj = datacursormode(fig);
    set(dcm_obj,'UpdateFcn',@myupdatefcn)
else
    set(handles.eid,'UserData',0);
    set(handles.eid,'String','Spatial pts');
    
    fig = figure(XRFDataViewer);
    datacursormode off
end
