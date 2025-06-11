function varargout = XRFMapping(varargin)
% XRFMAPPING MATLAB code for XRFMapping.fig
%      XRFMAPPING, by itself, creates a new XRFMAPPING or raises the existing
%      singleton*.
%
%      H = XRFMAPPING returns the handle to a new XRFMAPPING or the handle to
%      the existing singleton*.
%
%      XRFMAPPING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in XRFMAPPING.M with the given input arguments.
%
%      XRFMAPPING('Property','Value',...) creates a new XRFMAPPING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before XRFMapping_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to XRFMapping_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help XRFMapping

% Last Modified by GUIDE v2.5 24-Jun-2015 17:27:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @XRFMapping_OpeningFcn, ...
                   'gui_OutputFcn',  @XRFMapping_OutputFcn, ...
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


% --- Executes just before XRFMapping is made visible.
function XRFMapping_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to XRFMapping (see VARARGIN)

% Choose default command line output for XRFMapping
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes XRFMapping wait for user response (see UIRESUME)
% uiwait(handles.figure1);
pc = ~ismac;
set(handles.pc,'String',num2str(pc));
ncores = feature('NumCores');
set(handles.ncores,'String',num2str(ncores));

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
        a1 = str2double(line((e+1):end));
    end
    [~,~,e] = regexp(line,'^a0=','match','start','end');
    if (e>0)
        a0 = str2double(line((e+1):end));
    end
    [~,~,e] = regexp(line,'^elements=','match','start','end');
    if (e>0)
        Zlist = line((e+1):end);
    end
end
fclose(fid);
set(handles.a1,'String',num2str(a1));
set(handles.a0,'String',num2str(a0));

element_list = {'H'; 'He'; 'Li'; 'Be'; 'B'; 'C'; 'N'; 'O'; 'F'; ...
    'Ne'; 'Na'; 'Mg'; 'Al'; 'Si'; 'P'; 'S'; 'Cl'; 'Ar'; 'K'; ...
    'Ca'; 'Sc'; 'Ti'; 'V'; 'Cr'; 'Mn'; 'Fe'; 'Co'; 'Ni'; 'Cu'; ...
    'Zn'; 'Ga'; 'Ge'; 'As'; 'Se'; 'Br'; 'Kr'; 'Rb'; 'Sr'; 'Y'; ...
    'Zr'; 'Nb'; 'Mo'; 'Tc'; 'Ru'; 'Rh'; 'Pd'; 'Ag'; 'Cd'; 'In'; ...
    'Sn'; 'Sb'; 'Te'; 'I'; 'Xe'; 'Cs'; 'Ba'; 'La'; 'Ce'; 'Pr'; ...
    'Nd'; 'Pm'; 'Sm'; 'Eu'; 'Gd'; 'Tb'; 'Dy'; 'Ho'; 'Er'; 'Tm'; ...
    'Yb'; 'Lu'; 'Hf'; 'Ta'; 'W'; 'Re'; 'Os'; 'Ir'; 'Pt'; 'Au'; ...
    'Hg'; 'Tl'; 'Pb'; 'Bi'; 'Po'; 'At'; 'Rn'; 'Fr'; 'Ra'; 'Ac'; ...
    'Th'; 'Pa'; 'U'; 'Np'; 'Pu'; 'Am'; 'Cm'; 'Bk'; 'Cf'; 'Es'; ...
    'Fm'};
N = [];
tmp = regexp(Zlist,'\,','split');
for i = 1:length(tmp)
    tmp1 = strcmp(tmp(i),element_list);
    N = [N find(tmp1==1)];
end
N1 = [];
if (sum(N==12))
    set(handles.Mg0,'Value',1);
    N1 = [N1 12];
end
if (sum(N==13))
    set(handles.Al0,'Value',1);
    N1 = [N1 13];
end
if (sum(N==15))
    set(handles.P0,'Value',1);
    N1 = [N1 15];
end
if (sum(N==16))
    set(handles.S0,'Value',1);
    N1 = [N1 16];
end
if (sum(N==17))
    set(handles.Cl0,'Value',1);
    N1 = [N1 17];
end
if (sum(N==19))
    set(handles.K0,'Value',1);
    N1 = [N1 19];
end
if (sum(N==20))
    set(handles.Ca0,'Value',1);
    N1 = [N1 20];
end
if (sum(N==22))
    set(handles.Ti0,'Value',1);
    N1 = [N1 22];
end
if (sum(N==24))
    set(handles.Cr0,'Value',1);
    N1 = [N1 24];
end
if (sum(N==25))
    set(handles.Mn0,'Value',1);
    N1 = [N1 25];
end
if (sum(N==26))
    set(handles.Fe0,'Value',1);
    N1 = [N1 26];
end
if (sum(N==27))
    set(handles.Co0,'Value',1);
    N1 = [N1 27];
end
if (sum(N==28))
    set(handles.Ni0,'Value',1);
    N1 = [N1 28];
end
if (sum(N==29))
    set(handles.Cu0,'Value',1);
    N1 = [N1 29];
end
if (sum(N==30))
    set(handles.Zn0,'Value',1);
    N1 = [N1 30];
end
if (sum(N==33))
    set(handles.As0,'Value',1);
    N1 = [N1 33];
end
if (sum(N==34))
    set(handles.Se0,'Value',1);
    N1 = [N1 34];
end
if (sum(N==35))
    set(handles.Br0,'Value',1);
    N1 = [N1 35];
end
if (sum(N==47))
    set(handles.Ag0,'Value',1);
    N1 = [N1 47];
end
if (sum(N==48))
    set(handles.Cd0,'Value',1);
    N1 = [N1 48];
end
if (sum(N==50))
    set(handles.Sn0,'Value',1);
    N1 = [N1 50];
end
if (sum(N==51))
    set(handles.Sb0,'Value',1);
    N1 = [N1 51];
end
if (sum(N==56))
    set(handles.Ba0,'Value',1);
    N1 = [N1 56];
end
if (sum(N==79))
    set(handles.Au0,'Value',1);
    N1 = [N1 79];
end
if (sum(N==80))
    set(handles.Hg0,'Value',1);
    N1 = [N1 80];
end
if (sum(N==82))
    set(handles.Pb0,'Value',1);
    N1 = [N1 82];
end
[~,ia,~] = intersect(N,N1);
for i = length(ia):-1:1
    N(ia(i)) = [];
end
if (length(N) >= 1)
    Zother = '';
    Zother = [Zother char(element_list(N(1)))];
    for i = 2:length(N)
        Zother = [Zother ',' char(element_list(N(i)))];
    end
    set(handles.Zother,'String',Zother);
end

set(handles.running_status,'String','Idle');
drawnow;
pause(1)

% --- Outputs from this function are returned to the command line.
function varargout = XRFMapping_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in file_info.
function file_info_Callback(hObject, eventdata, handles)
% hObject    handle to file_info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fn,path] = uigetfile('*','Select the XRF cube.');
set(handles.fnpath,'String',path);
set(handles.trial,'String',fn);

function fnpath_Callback(hObject, eventdata, handles)
% hObject    handle to fnpath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fnpath as text
%        str2double(get(hObject,'String')) returns contents of fnpath as a double


% --- Executes during object creation, after setting all properties.
function fnpath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fnpath (see GCBO)
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


% --- Executes on button press in run_button.
function run_button_Callback(hObject, eventdata, handles)
% hObject    handle to run_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.running_status,'String','Running ...');
drawnow;
pause(1)

pc = str2double(get(handles.pc,'String'));
fnpath = get(handles.fnpath,'String');
trial = get(handles.trial,'String');
if (pc == 0)
    xrf_table_fn = '/Applications/XRFMapping/application/xrf_table.csv';
elseif (pc == 1)
    xrf_table_fn = 'xrf_table.csv';
end
minE = str2double(get(handles.minE,'String'));
maxE = str2double(get(handles.maxE,'String'));
ncores = str2double(get(handles.ncores,'String'));
a1 = str2double(get(handles.a1,'String'));
a0 = str2double(get(handles.a0,'String'));
update_coeffs = get(handles.update_coeffs,'Value');

element_list = {'H'; 'He'; 'Li'; 'Be'; 'B'; 'C'; 'N'; 'O'; 'F'; ...
    'Ne'; 'Na'; 'Mg'; 'Al'; 'Si'; 'P'; 'S'; 'Cl'; 'Ar'; 'K'; ...
    'Ca'; 'Sc'; 'Ti'; 'V'; 'Cr'; 'Mn'; 'Fe'; 'Co'; 'Ni'; 'Cu'; ...
    'Zn'; 'Ga'; 'Ge'; 'As'; 'Se'; 'Br'; 'Kr'; 'Rb'; 'Sr'; 'Y'; ...
    'Zr'; 'Nb'; 'Mo'; 'Tc'; 'Ru'; 'Rh'; 'Pd'; 'Ag'; 'Cd'; 'In'; ...
    'Sn'; 'Sb'; 'Te'; 'I'; 'Xe'; 'Cs'; 'Ba'; 'La'; 'Ce'; 'Pr'; ...
    'Nd'; 'Pm'; 'Sm'; 'Eu'; 'Gd'; 'Tb'; 'Dy'; 'Ho'; 'Er'; 'Tm'; ...
    'Yb'; 'Lu'; 'Hf'; 'Ta'; 'W'; 'Re'; 'Os'; 'Ir'; 'Pt'; 'Au'; ...
    'Hg'; 'Tl'; 'Pb'; 'Bi'; 'Po'; 'At'; 'Rn'; 'Fr'; 'Ra'; 'Ac'; ...
    'Th'; 'Pa'; 'U'; 'Np'; 'Pu'; 'Am'; 'Cm'; 'Bk'; 'Cf'; 'Es'; ...
    'Fm'};
N = [];
Mg0 = get(handles.Mg0,'Value');
Al0 = get(handles.Al0,'Value');
P0 = get(handles.P0,'Value');
S0 = get(handles.S0,'Value');
Cl0 = get(handles.Cl0,'Value');
K0 = get(handles.K0,'Value');
Ca0 = get(handles.Ca0,'Value');
Ti0 = get(handles.Ti0,'Value');
Cr0 = get(handles.Cr0,'Value');
Mn0 = get(handles.Mn0,'Value');
Fe0 = get(handles.Fe0,'Value');
Co0 = get(handles.Co0,'Value');
Ni0 = get(handles.Ni0,'Value');
Cu0 = get(handles.Cu0,'Value');
Zn0 = get(handles.Zn0,'Value');
As0 = get(handles.As0,'Value');
Se0 = get(handles.Se0,'Value');
Br0 = get(handles.Br0,'Value');
Ag0 = get(handles.Ag0,'Value');
Cd0 = get(handles.Cd0,'Value');
Sn0 = get(handles.Sn0,'Value');
Sb0 = get(handles.Sb0,'Value');
Ba0 = get(handles.Ba0,'Value');
Au0 = get(handles.Au0,'Value');
Hg0 = get(handles.Hg0,'Value');
Pb0 = get(handles.Pb0,'Value');
Zother = get(handles.Zother,'String');
if (Mg0 == 1)
    N = [N 12];
end
if (Al0 == 1)
    N = [N 13];
end
if (P0 == 1)
    N = [N 15];
end
if (S0 == 1)
    N = [N 16];
end
if (Cl0 == 1)
    N = [N 17];
end
if (K0 == 1)
    N = [N 19];
end
if (Ca0 == 1)
    N = [N 20];
end
if (Ti0 == 1)
    N = [N 22];
end
if (Cr0 == 1)
    N = [N 24];
end
if (Mn0 == 1)
    N = [N 25];
end
if (Fe0 == 1)
    N = [N 26];
end
if (Co0 == 1)
    N = [N 27];
end
if (Ni0 == 1)
    N = [N 28];
end
if (Cu0 == 1)
    N = [N 29];
end
if (Zn0 == 1)
    N = [N 30];
end
if (As0 == 1)
    N = [N 33];
end
if (Se0 == 1)
    N = [N 34];
end
if (Br0 == 1)
    N = [N 35];
end
if (Ag0 == 1)
    N = [N 47];
end
if (Cd0 == 1)
    N = [N 48];
end
if (Sn0 == 1)
    N = [N 50];
end
if (Sb0 == 1)
    N = [N 51];
end
if (Ba0 == 1)
    N = [N 56];
end
if (Au0 == 1)
    N = [N 79];
end
if (Hg0 == 1)
    N = [N 80];
end
if (Pb0 == 1)
    N = [N 82];
end
msk = isspace(Zother);
Zother(msk==1) = '';
tmp = regexp(Zother,'\,','split');
for i = 1:length(tmp)
    tmp1 = strcmp(tmp(i),element_list);
    N = [N find(tmp1==1)];
end

m1 = str2double(get(handles.m1,'String'));
n1 = str2double(get(handles.n1,'String'));
m2 = get(handles.m2,'String');
if (strcmpi(m2,'end'))
    m2 = 0;
else
    m2 = str2double(m2);
end
n2 = get(handles.n2,'String');
if (strcmpi(n2,'end'))
    n2 = 0;
else
    n2 = str2double(n2);
end

[map1,mapp,A1,fwhm1,c1,r1,element_list,LN,c,xrf_table,xray_label,a1,a0,r2] = xrf_classification1(fnpath,trial,xrf_table_fn,m1,m2,n1,n2,N,minE,maxE,a1,a0,update_coeffs,ncores);
if (update_coeffs == 1)
    set(handles.a1,'String',num2str(a1));
    set(handles.a0,'String',num2str(a0));
    set(handles.r2,'String',num2str(r2));
end
save([fnpath trial '.mat'],'map1','mapp','A1','fwhm1','c1','r1','element_list','N','LN','c','xrf_table','xray_label','a1','a0')
pc = str2double(get(handles.pc,'String'));
if (pc == 0)
    fn = '/Applications/XRFMapping/application/XRFMapping.txt';
elseif (pc == 1)
    fn = 'XRFMapping.txt';
end
if (pc == 0)
    fid = fopen(fn,'w+t','n');
elseif (pc == 1)
    fid = fopen(fn,'w+t','n');
end
fprintf(fid,'a1 = %s\n',num2str(a1));
fprintf(fid,'a0 = %s\n',num2str(a0));
fprintf(fid,'elements = %s,',char(element_list(N(1))));
for i = (2:(length(N)-1))
    fprintf(fid,'%s,',char(element_list(N(i))));
end
fprintf(fid,'%s\n',char(element_list(N(end))));
fclose(fid);

set(handles.running_status,'String','Idle');
drawnow;
pause(1)

% --- Executes on button press in Mg0.
function Mg0_Callback(hObject, eventdata, handles)
% hObject    handle to Mg0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Mg0


% --- Executes on button press in Al0.
function Al0_Callback(hObject, eventdata, handles)
% hObject    handle to Al0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Al0


% --- Executes on button press in P0.
function P0_Callback(hObject, eventdata, handles)
% hObject    handle to P0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of P0


% --- Executes on button press in S0.
function S0_Callback(hObject, eventdata, handles)
% hObject    handle to S0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of S0


% --- Executes on button press in Cl0.
function Cl0_Callback(hObject, eventdata, handles)
% hObject    handle to Cl0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Cl0


% --- Executes on button press in K0.
function K0_Callback(hObject, eventdata, handles)
% hObject    handle to K0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of K0


% --- Executes on button press in Ca0.
function Ca0_Callback(hObject, eventdata, handles)
% hObject    handle to Ca0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Ca0


% --- Executes on button press in Ti0.
function Ti0_Callback(hObject, eventdata, handles)
% hObject    handle to Ti0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Ti0


% --- Executes on button press in Cr0.
function Cr0_Callback(hObject, eventdata, handles)
% hObject    handle to Cr0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Cr0


% --- Executes on button press in Mn0.
function Mn0_Callback(hObject, eventdata, handles)
% hObject    handle to Mn0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Mn0


% --- Executes on button press in Fe0.
function Fe0_Callback(hObject, eventdata, handles)
% hObject    handle to Fe0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Fe0


% --- Executes on button press in Co0.
function Co0_Callback(hObject, eventdata, handles)
% hObject    handle to Co0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Co0


% --- Executes on button press in Ni0.
function Ni0_Callback(hObject, eventdata, handles)
% hObject    handle to Ni0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Ni0


% --- Executes on button press in Cu0.
function Cu0_Callback(hObject, eventdata, handles)
% hObject    handle to Cu0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Cu0


% --- Executes on button press in Zn0.
function Zn0_Callback(hObject, eventdata, handles)
% hObject    handle to Zn0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Zn0


% --- Executes on button press in As0.
function As0_Callback(hObject, eventdata, handles)
% hObject    handle to As0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of As0


% --- Executes on button press in Se0.
function Se0_Callback(hObject, eventdata, handles)
% hObject    handle to Se0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Se0


% --- Executes on button press in Br0.
function Br0_Callback(hObject, eventdata, handles)
% hObject    handle to Br0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Br0


% --- Executes on button press in Ag0.
function Ag0_Callback(hObject, eventdata, handles)
% hObject    handle to Ag0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Ag0


% --- Executes on button press in Cd0.
function Cd0_Callback(hObject, eventdata, handles)
% hObject    handle to Cd0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Cd0


% --- Executes on button press in Sn0.
function Sn0_Callback(hObject, eventdata, handles)
% hObject    handle to Sn0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Sn0


% --- Executes on button press in Sb0.
function Sb0_Callback(hObject, eventdata, handles)
% hObject    handle to Sb0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Sb0


% --- Executes on button press in Ba0.
function Ba0_Callback(hObject, eventdata, handles)
% hObject    handle to Ba0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Ba0


% --- Executes on button press in Au0.
function Au0_Callback(hObject, eventdata, handles)
% hObject    handle to Au0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Au0


% --- Executes on button press in Hg0.
function Hg0_Callback(hObject, eventdata, handles)
% hObject    handle to Hg0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Hg0


% --- Executes on button press in Pb0.
function Pb0_Callback(hObject, eventdata, handles)
% hObject    handle to Pb0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Pb0



function Zother_Callback(hObject, eventdata, handles)
% hObject    handle to Zother (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Zother as text
%        str2double(get(hObject,'String')) returns contents of Zother as a double


% --- Executes during object creation, after setting all properties.
function Zother_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Zother (see GCBO)
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

% Hints: get(hObject,'String') returns contents of a1 as text
%        str2double(get(hObject,'String')) returns contents of a1 as a double


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

% Hints: get(hObject,'String') returns contents of a0 as text
%        str2double(get(hObject,'String')) returns contents of a0 as a double


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


% --- Executes on button press in update_coeffs.
function update_coeffs_Callback(hObject, eventdata, handles)
% hObject    handle to update_coeffs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of update_coeffs



function minE_Callback(hObject, eventdata, handles)
% hObject    handle to minE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minE as text
%        str2double(get(hObject,'String')) returns contents of minE as a double


% --- Executes during object creation, after setting all properties.
function minE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxE_Callback(hObject, eventdata, handles)
% hObject    handle to maxE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxE as text
%        str2double(get(hObject,'String')) returns contents of maxE as a double


% --- Executes during object creation, after setting all properties.
function maxE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxE (see GCBO)
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



function n1_Callback(hObject, eventdata, handles)
% hObject    handle to n1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n1 as text
%        str2double(get(hObject,'String')) returns contents of n1 as a double


% --- Executes during object creation, after setting all properties.
function n1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function m1_Callback(hObject, eventdata, handles)
% hObject    handle to m1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of m1 as text
%        str2double(get(hObject,'String')) returns contents of m1 as a double


% --- Executes during object creation, after setting all properties.
function m1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to m1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function n2_Callback(hObject, eventdata, handles)
% hObject    handle to n2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n2 as text
%        str2double(get(hObject,'String')) returns contents of n2 as a double


% --- Executes during object creation, after setting all properties.
function n2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function m2_Callback(hObject, eventdata, handles)
% hObject    handle to m2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of m2 as text
%        str2double(get(hObject,'String')) returns contents of m2 as a double


% --- Executes during object creation, after setting all properties.
function m2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to m2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function r2_Callback(hObject, eventdata, handles)
% hObject    handle to r2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of r2 as text
%        str2double(get(hObject,'String')) returns contents of r2 as a double


% --- Executes during object creation, after setting all properties.
function r2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function running_status_Callback(hObject, eventdata, handles)
% hObject    handle to running_status (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of running_status as text
%        str2double(get(hObject,'String')) returns contents of running_status as a double


% --- Executes during object creation, after setting all properties.
function running_status_CreateFcn(hObject, eventdata, handles)
% hObject    handle to running_status (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
