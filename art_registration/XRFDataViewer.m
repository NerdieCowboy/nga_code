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

% Last Modified by GUIDE v2.5 07-Jul-2015 15:35:23

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

pc = ~ismac;
set(handles.pc,'String',num2str(pc));
ncores = feature('NumCores');

% read in previous calibration coefficients
if (pc == 0)
    fnc = '/Applications/XRFMapping/application/XRFMapping.txt';
elseif (pc == 1)
    fnc = 'XRFMapping.txt';
end
fid = fopen(fnc);
while ~feof(fid)
    line = fgetl(fid);
    msk = isspace(line);
    line(msk==1) = '';
    [~,~,e] = regexp(line,'^a1=','match','start','end');
    if (e>0)
        slope = str2double(line((e+1):end));
    end
    [~,~,e] = regexp(line,'^a0=','match','start','end');
    if (e>0)
        yint = str2double(line((e+1):end));
    end
end
fclose(fid);
set(handles.slope,'String',num2str(slope));
set(handles.yint,'String',num2str(yint));

set(handles.ncores,'String',num2str(ncores));
set(handles.axes1,'Units','pixels');
set(handles.axes2,'Units','pixels');

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
slope = str2num(get(handles.slope,'String')); %#ok<ST2NM>
yint = str2num(get(handles.yint,'String')); %#ok<ST2NM>
set(handles.slider1,'Min',p3); %#%#ok<MSNU> ok<ST2NM>
set(handles.slider1,'Max',p2);
ss = 1/(p2-p3+1);
set(handles.slider1,'SliderStep',[ss 40*ss]);

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

axes(handles.axes1)
cube = single(multibandread(trial_out,[m0,n0,p2-p3+1],'single',0,'bsq','ieee-le',{'Band','Direct',p1-p3+1}));

b0 = get(handles.axes1,'Position');
m1 = b0(4);
n1 = b0(3);
cube = imresize(cube,[m1 n1]);
imagesc(cube)
colormap gray
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
    set(handles.img_path,'String',path);
    ind = regexp(fn, '\.');
    set(handles.img_ext,'String',fn((ind+1):end));
    fn(ind:end) = [];
    set(handles.trial,'String',fn);
else
    [fn,path] = uigetfile('*','Select the XRF files.');
    set(handles.img_path,'String',path);
    set(handles.trial,'String',fn);
end

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

ncores = str2num(get(handles.ncores,'String'));
pc = str2num(get(handles.pc,'String'));
p3 = str2num(get(handles.p3,'String'));
p2 = str2num(get(handles.p2,'String'));
slope = str2num(get(handles.slope,'String'));
yint = str2num(get(handles.yint,'String'));
set(handles.slider1,'Enable','off');
set(handles.slider1,'Visible','on');
set(handles.slider1,'Min',p3);
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

n0 = str2num(get(handles.n0,'String'));
m0 = str2num(get(handles.m0,'String'));
p0 = str2num(get(handles.p0,'String'));
pixel_sz = str2num(get(handles.pixel_sz,'String'));
p2 = str2num(get(handles.p2,'String'));
p3 = str2num(get(handles.p3,'String'));
read_only = get(handles.read_only,'Value');
e3 = p3*slope + yint;
e2 = p2*slope + yint;
set(handles.e3,'string',num2str(e3));
set(handles.e2,'string',num2str(e2));

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
    p3 = 1;
    set(handles.p3,'string',num2str(p3));
    p2 = p0;
    set(handles.p2,'string',num2str(p2));
    e3 = p3*slope + yint;
    e2 = p2*slope + yint;
    set(handles.e3,'string',num2str(e3));
    set(handles.e2,'string',num2str(e2));
elseif (read_only == 0)
    if (ncores > 1)
        isOpen = matlabpool('size') > 0;
        if (isOpen == 1)
            matlabpool close
        end
        matlabpool(ncores)
    end
    
    m = zeros(m0*n0,1);
    n = zeros(m0*n0,1);
    fncsv = [csv_path csv_fn '.' csv_ext];
    if (exist(fncsv,'file') == 2)
        tmp = csvread(fncsv,1,0);
        x = single(tmp(:,1));
        y = single(tmp(:,2));
        clear tmp
        x = (x-min(x(:)))/pixel_sz + 1;
        y = (y-min(y(:)))/pixel_sz + 1;

        MN = size(x,1);
        deltam = abs(y(2:MN) - y(1:(MN-1)));
        deltam = [0;deltam];
        cnt = 0;
        pnt = 1;
        tmpm = zeros(n0+1,1);
        tmpn = zeros(n0+1,1);
        dir = 0;
        for i = 1:MN
            if (deltam(i) <= 0.5)
                cnt = cnt + 1;
            else
                if (cnt > n0)
                    if (dir == 0)
                        m(pnt:(pnt+n0-1)) = tmpm(1:n0);
                        n(pnt:(pnt+n0-1)) = tmpn(1:n0);
                    elseif (dir == 1)
                        m(pnt:(pnt+n0-1)) = tmpm(2:(n0+1));
                        n(pnt:(pnt+n0-1)) = tmpn(2:(n0+1));
                    end
                else
                    m(pnt:(pnt+n0-1)) = tmpm(1:n0);
                    n(pnt:(pnt+n0-1)) = tmpn(1:n0);
                end
                pnt = pnt + n0;
                cnt = 1;
                if (dir == 0)
                    dir = 1;
                elseif (dir == 1)
                    dir = 0;
                end
            end
            tmpm(cnt) = y(i);
            tmpn(cnt) = x(i);
        end
        if (cnt > n0)
            if (dir == 0)
                m(pnt:(pnt+n0-1)) = tmpm(1:n0);
                n(pnt:(pnt+n0-1)) = tmpn(1:n0);
            elseif (dir == 1)
                m(pnt:(pnt+n0-1)) = tmpm(2:(n0+1));
                n(pnt:(pnt+n0-1)) = tmpn(2:(n0+1));
            end
        else
            m(pnt:(pnt+n0-1)) = tmpm(1:n0);
            n(pnt:(pnt+n0-1)) = tmpn(1:n0);
        end
        clear x y
        x = n;
        y = m;
        clear m n
    end
    
    fn = [img_path '\' trial '.bin'];
    fid = fopen(fn,'r');
    dir = 0;
    cnt = 1;
    cube = single(zeros(m0,n0,(p2-p3+1)));
    if (exist(fncsv,'file') == 2)
        X = single(zeros(m0,n0));
        Y = single(zeros(m0,n0));
        for m = 1:m0
            for n = 1:n0
                if (dir == 0)
                    X(m,n) = x(cnt) + 0.5;
                    Y(m,n) = y(cnt);
                elseif (dir == 1)
                    n1 = n0 - (n - 1);
                    X(m,n1) = x(cnt) - 0.5;
                    Y(m,n1) = y(cnt);
                end
                cnt = cnt + 1;
            end
            if (dir == 0)
                dir = 1;
            elseif (dir == 1)
                dir = 0;
            end
        end
        
        X = X(:,2:(n0-1));
        Y = Y(:,2:(n0-1));
        mnx = ceil(min(X(:)));
        mxx = floor(max(X(:)));
        mny = ceil(min(Y(:)));
        mxy = floor(max(Y(:)));
        set(handles.mnx,'String',num2str(mnx));
        set(handles.mxx,'String',num2str(mxx));
        set(handles.mny,'String',num2str(mny));
        set(handles.mxy,'String',num2str(mxy));
    end
    
    for m = 1:m0
        for n = 1:n0
            tmp = single(fread(fid,p0,'uint32',0,'ieee-le'));
            tmp = tmp(p3:p2);
            if (dir == 0)
                if (n > 1)
                    cube(m,n,:) = tmp - tmp_prev;
                    %cube(m,n,:) = tmp;
                end
            elseif (dir == 1)
                n1 = n0 - (n - 1);
                if (n1 < n0)
                    cube(m,n1,:) = tmp - tmp_prev;
                    %cube(m,n1,:) = tmp;
                end
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
        [N,M] = meshgrid(mnx:mxx,mny:mxy);
        [m0,n0] = size(N);
        cube2 = single(zeros(m0,n0,(p2-p3+1)));
        parfor p = 1:(p2-p3+1)
            tmp = double(cube(:,:,p));
            F = TriScatteredInterp(double(X(:)),double(Y(:)),double(tmp(:)),'natural');
            cube2(:,:,p) = F(double(N),double(M));
        end
        clear tmp tmp1
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
    
    if (ncores > 1)
        isOpen = matlabpool('size') > 0;
        if (isOpen == 1)
            matlabpool close
        end
    end
end

p1 = str2num(get(handles.p1,'string'));
energy1 = p1*slope+yint;
set(handles.energy1,'string',num2str(energy1));

if (pc == 0)
    fn = '/Applications/XRFMapping/application/XRFMapping.txt';
elseif (pc == 1)
    fn = 'XRFMapping.txt';
end
if (pc == 0)
    fid = fopen(fn,'w');
elseif (pc == 1)
    fid = fopen(fn,'w+t','n');
end
fprintf(fid,'a1 = %s\n',num2str(slope));
fprintf(fid,'a0 = %s\n',num2str(yint));
fclose(fid);

axes(handles.axes1)
tmp= single(multibandread(trial_out,[m0,n0,p2-p3+1],'single',0,'bsq','ieee-le',{'Band','Direct',p1-p3+1}));

b0 = get(handles.axes1,'Position');
sc1 = b0(3)/n0;
sc2 = b0(4)/m0;
if (sc2 < sc1)
    sc = sc2;
else
    sc = sc1;
end
n1 = floor(n0*sc);
m1 = floor(m0*sc);
tmp = imresize(tmp,[m1 n1]);
b0(3) = n1;
b0(4) = m1;
set(handles.axes1,'Position',b0)

imagesc(tmp)
colormap gray
axis off
clear tmp cube2

lambda = (p3:p2)*slope+yint;
axes(handles.axes2)
plot(lambda,zeros(1,(p2-p3+1)))
xlabel('energy (KeV)');
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

    n0 = str2num(get(handles.n0,'String'));
    m0 = str2num(get(handles.m0,'String'));
    p0 = str2num(get(handles.p0,'String'));
    p2 = str2num(get(handles.p2,'String'));
    p3 = str2num(get(handles.p3,'String'));
    logon = get(handles.logon,'Value');
    coordn = round(n0*x/pos(3));
    coordm = round(m0*y/pos(4));
    set(handles.x1,'String',num2str(coordn));
    set(handles.y1,'String',num2str(coordm));

    slope = str2num(get(handles.slope,'String'));
    yint = str2num(get(handles.yint,'String'));
    lambda = (p3:p2)*slope+yint;
    e3 = p3*slope + yint;
    e2 = p2*slope + yint;
    set(handles.e3,'string',num2str(e3));
    set(handles.e2,'string',num2str(e2));

    trial_out = [img_path trial];
    cube = single(multibandread(trial_out,[m0,n0,p2-p3+1],'single',0,'bsq','ieee-le',{'Row','Direct',coordm},{'Column','Direct',coordn},{'Band','Range',[1,p2-p3+1]}));
    
    cube(cube<0) = 0;
    cube(isnan(cube)) = 0;
    cube(isinf(cube)) = 0;
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
txt = {['Energy: ',num2str(ptpos(1)),' KeV'],...
       ['Counts: ',num2str(round(exp(ptpos(2))))]};
   
function txt = myupdatefcn2(empt,event_obj)
% Customizes text of data tips

ptpos = get(event_obj,'Position');
txt = {['Energy: ',num2str(ptpos(1)),' KeV'],...
       ['Counts: ',num2str(round(ptpos(2)))]};

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

logon = get(handles.logon,'Value');
eid_state = get(handles.eid,'UserData');
if (eid_state == 0)
    set(handles.eid,'UserData',1);
    set(handles.eid,'String','Energy');
    
    fig = figure(XRFDataViewer);
    datacursormode on
    dcm_obj = datacursormode(fig);
    if (logon == 1)
        set(dcm_obj,'UpdateFcn',@myupdatefcn)
    else
        set(dcm_obj,'UpdateFcn',@myupdatefcn2)
    end
else
    set(handles.eid,'UserData',0);
    set(handles.eid,'String','Spatial pts');
    
    fig = figure(XRFDataViewer);
    datacursormode off
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



function mnx_Callback(hObject, eventdata, handles)
% hObject    handle to mnx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mnx as text
%        str2double(get(hObject,'String')) returns contents of mnx as a double


% --- Executes during object creation, after setting all properties.
function mnx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mnx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mxx_Callback(hObject, eventdata, handles)
% hObject    handle to mxx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mxx as text
%        str2double(get(hObject,'String')) returns contents of mxx as a double


% --- Executes during object creation, after setting all properties.
function mxx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mxx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mny_Callback(hObject, eventdata, handles)
% hObject    handle to mny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mny as text
%        str2double(get(hObject,'String')) returns contents of mny as a double


% --- Executes during object creation, after setting all properties.
function mny_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mxy_Callback(hObject, eventdata, handles)
% hObject    handle to mxy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mxy as text
%        str2double(get(hObject,'String')) returns contents of mxy as a double


% --- Executes during object creation, after setting all properties.
function mxy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mxy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
