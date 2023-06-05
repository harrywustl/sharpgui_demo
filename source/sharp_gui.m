% Sharp GUI ver0.4
% Usage
% need 'plotpkmarkers.m', 'sharp_gui.m', and 'sharp_gui.fig' in same folder
% use MATLAB import data tool, to load 'SharpElectrodeRecording'
% After loading data into the base workspace, in the box next to Data Name, 
% type in the data set name, such as 'data_P18', then click 'Import Data', 
% which loads into current program workspace. 
% After setting parameters, click 'Find Peaks', which plots data, peaks and other statistics.

% All Right Reserved
% 
% Release Date: Apr 20, 2016
% Created by:   Wenjie Wu
% Email:        wuwenjie@wustl.edu
% Organzation:  Washington University in St. Louis, School of Medicine
% Description:  GUI for processing sharp electrode data
%               improved on picking out lone peaks during resting

%%
function varargout = sharp_gui(varargin)
% SHARP_GUI MATLAB code for sharp_gui.fig
%      SHARP_GUI, by itself, creates a new SHARP_GUI or raises the existing
%      singleton*.
%
%      H = SHARP_GUI returns the handle to a new SHARP_GUI or the handle to
%      the existing singleton*.
%
%      SHARP_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHARP_GUI.M with the given input arguments.
%
%      SHARP_GUI('Property','Value',...) creates a new SHARP_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sharp_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sharp_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sharp_gui

% Last Modified by GUIDE v2.5 20-Apr-2016 17:42:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sharp_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @sharp_gui_OutputFcn, ...
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


% --- Executes just before sharp_gui is made visible.
function sharp_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sharp_gui (see VARARGIN)

% WENJIE: Add my own to 'handles' structure
handles.myfontsize = 10;
handles.data = [0 0;1 1;2 2]; % initialize with arbitrary data
handles.pkVals = []; % value of peaks
handles.pkLocs = []; % locations of peaks
handles.pkProm = []; % prominance of peaks
handles.pkW = []; % width of peaks
handles.minPkW = 20; % minimum peak width, initialize
handles.minPkP = 7; % minimum peak prominence, initialize
handles.minIntvl = 4000; % minimum resting interval, (ms), initialize
handles.minBurst = 2800; % minimum burst, (ms) , initialize
handles.restIntvl = 0; % initialize with arbitrary number
handles.burstDur = 0; % initialize with arbitrary number
handles.burstStart = 0; % indices of the start of burst
handles.burstEnd = 0; % indices of the end of burst
handles.wsVarName = ''; % varible name in the workspace
handles.sampleFreq = 1000; % sampling frequency

plot(handles.dataPlot,0,0); % plot a dummy as dataPlot

fz = handles.myfontsize;
xlabel('time (ms)', 'FontSize',fz);
ylabel('voltage (mV)','FontSize',fz);
grid on;
grid minor;
disp('sharp_gui_opening!');

% Choose default command line output for sharp_gui
handles.output = hObject;
guidata(hObject, handles); % update handles structure

% UIWAIT makes sharp_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = sharp_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% WENJIE: Executes on button press in uiImportData.
function uiImportData_Callback(hObject, eventdata, handles)
% hObject    handle to uiImportData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% WENJIE: import base workspace data into gui workspace
tmpStr = handles.wsVarName;

w = evalin('base','who'); % return all base varibles' name
isInBase = strcmp(tmpStr,w); % is the requested data in base workspace?

if sum(isInBase) ~= 0 % if at least one is in base workspace, get data from it
    handles.data = evalin('base',tmpStr); 
    plot(handles.dataPlot, handles.data(:,end));
    disp('NOTE: zoom out if you see nothing after import.');
else
    disp('ERROR: your requested data has not been loaded in base workspace.');
    plot(0,0); % plot a dummy
end

% WENJIE: important!
% the 'NextPlot' property of the axes has been set to 'replacechildren',
% so that the updated plot will not change zoom.

guidata(hObject,handles); % update gui handles

% WENJIE: Executes on button press in uiFindPks.
function uiFindPks_Callback(hObject, eventdata, handles)
% hObject    handle to uiFindPks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% WENJIE: get parameters from handles
minW = handles.minPkW; % min peak width
minP = handles.minPkP; % min peak prominence
minIntvl = handles.minIntvl; % min resting interval to be recognized
minBurst = handles.minBurst; % min burst to be recognized

% WENJIE: find peaks 
% If you specify neither 'x' nor 'Fs', then returned 'locs' is a vector of integer indices.
% If you specify a location vector 'x', then returned 'locs' contains the values of x at the peak indices.
% I didn't specify 'x' as [pks,locs] = findpeaks(data,x)
[pks,locs,pkW,pkProm] = findpeaks(handles.data(:,end),...
    'MinPeakWidth',minW,'MinPeakProminence',minP);

dx = diff(locs); % interval between to indices of peaks

% 1st round criteria, roughly find all burst packets
% burstStart and burstEnd stores absolute indices of the data
dxIdx = find(dx > minIntvl); % find idx of which is a candidate of resting interval
burstStart = [locs(1); locs(dxIdx+1)]; % start indices of data for all burst packets, idx+1 moves to the next peak
burstEnd = [locs(dxIdx); locs(end)]; % end indices of data for all burst packets
burstDur = burstEnd - burstStart; % duration for all bursts

% 2nd round criteria, pick out a few lone peaks within resting interval
% burstStart and burstEnd stores absolute indices of the data
notBurstIdx = find(burstDur < minBurst); % indices that is too short to be consider as burst
burstStart(notBurstIdx) = []; % delete indices that is too short to be consider as burst
burstEnd(notBurstIdx) = []; % delete indices that is too short to be consider as burst
burstDur = burstEnd - burstStart; % re-calculate real burst duration

% put it in handles
handles.burstStart = burstStart;
handles.burstEnd = burstEnd;
handles.burstDur = burstDur; 

% index of locs, where burst starts and ends
locStart = []; % index of locs where burst starts, initialize
for i = 1:size(burstStart)
    locStart = [locStart; find(locs == burstStart(i))];
end
locEnd = []; % index of locs where burst ends, initialize
for i = 1:size(burstEnd)
    locEnd = [locEnd; find(locs == burstEnd(i))];
end
% locStart % For DEBUG
% locEnd % For DEBUG

% pickout peak location that is not within burst packets
truePkLocs = []; % locations of peaks
truePkVals = []; % value of peaks
truePkProm = []; % prominance of peaks
truePkW = [];
% size(locStart)
if size(locStart) == size (locEnd)
    for i = 1:size(locStart)
        truePkLocs = [truePkLocs; locs(locStart(i):locEnd(i))];
        truePkVals = [truePkVals; pks(locStart(i):locEnd(i))];
        truePkProm = [truePkProm; pkProm(locStart(i):locEnd(i))];
        truePkW = [truePkW; pkW(locStart(i):locEnd(i))];
    end    
else
    disp('ERROR: locStart and locEnd have different dimesion');
end

% truePkLocs % For DEBUG
handles.pkLocs = truePkLocs; % locations of peaks
handles.pkVals = truePkVals; % value of peaks
handles.pkProm = truePkProm; % prominance of peaks
handles.pkW = truePkW; % width of peaks

% re-calculate real resting interval
restIntvl = burstStart(2:end) - burstEnd(1:end-1); %  start2-end1, start3-end2, ...
handles.restIntvl = restIntvl; % put it in handles

% WENJIE: plot, and mark peaks
% now it doesn't have x - the real time stamp
y = handles.data(:,2);

hAxe = handles.dataPlot; % retrieve dataPlot handle

plot(hAxe,y,'Tag','Signal'); % plot 
hAxes = ancestor(hAxe,'Axes'); 
color = get(hAxe,'Color');
hAxe = line(locs,y(locs),'Parent',hAxes, ...
     'Marker','o','LineStyle','none','Color',[1.0,0.4,0.0],'tag','Peak');
if coder.target('MATLAB')
    plotpkmarkers(hAxe,y(locs)); % plot offset peaks, need plotpkmarks function 
end %

hAxe = line(burstStart,y(burstStart),'Parent',hAxes, ...
     'Marker','o','LineStyle','none','Color','b','tag','Peak');
plotpkmarkers(hAxe,y(burstStart)); % plot offset peaks, need plotpkmarks function

hAxe = line(burstEnd,y(burstEnd),'Parent',hAxes, ...
     'Marker','o','LineStyle','none','Color','g','tag','Peak');
plotpkmarkers(hAxe,y(burstEnd)); % plot offset peaks, need plotpkmarks function

handles.dataPlot = hAxes;

% WENJIE: data stat calculation
riMean = mean(restIntvl); % mean of resting interval
riStd = std(restIntvl); % std of resting interval
bdMean = mean(burstDur); % mean of burst duration
bdStd = std(burstDur); % std of burst duration
pksMean = mean(pks); % mean of peaks
pksStd = std(pks); % std of peaks

allStat = [pksMean, pksStd; riMean, riStd; bdMean, bdStd];

% WENJIE: update data stat table
set(handles.uiRestIntvl, 'Data', restIntvl);
set(handles.uiBurstDur, 'Data', burstDur);
set(handles.uiAllStat, 'Data', allStat);

guidata(hObject,handles); % update gui handles

% WENJIE: Executes on button press in uiFFT.
function uiFFT_Callback(hObject, eventdata, handles)
% hObject    handle to uiFFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
y = handles.data(:,2); % data values

Fs = handles.sampleFreq; % sampling freuency

xLim = handles.dataPlot.XLim; % xlimit of dataPlot
xStart = uint32(xLim(1)) + 1; % plus 1 to prevent 0 index, uint can be used as index 
xEnd = uint32(xLim(end)) + 1; % plus 1 to prevent 0 index, uint can be used as index 
% check within the bound of data y, size(y,1) returns size of rows
if xEnd > size(y,1)
    xEnd = size(y,1);
end
% if the length from xStart to xEnd is odd
if mod((xEnd - xStart + 1),2) ~= 0
    y = y(xStart:xEnd - 1); % make it a even number of length
else
    y = y(xStart:xEnd);
end
L = size(y,1); % len of signal, must be even, returns size of rows

f = Fs*(0:(L/2))/L; % construct freq vector

% fft, abs gets magnitude of complex fft
% fft scaling, see http://www.mathworks.com/matlabcentral/answers/15770-scaling-the-fft-and-the-ifft
fftY = abs(fft(y)/Fs); % scale fft with 1/Fs(i.e. dt), Parseval's theorem
fftYpside = fftY(1:(L/2+1)); % pick positive freq of FFT
fftYpside(2:end-1) = 2 * fftYpside(2:end-1); % multiply by two expect DC(at 0 freq)
figure(1) % open a new figure
plot(f,fftYpside);
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('f (Hz)')
ylabel('|Y(f)|')
grid on;
% disp('FFT'); % for debug

% --- Executes on button press in uiSaveBursts.
function uiSaveBursts_Callback(hObject, eventdata, handles)
% hObject    handle to uiSaveBursts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% WENJIE: As Apr 19, 2016, MATLAB for Mac does not support saving excel file. 
% CSV file is saved instead
% Chartacters does not work well in csv, therefore, 0s are used as
% seperaters. 
bStart = handles.burstStart; % get burstStart
bEnd = handles.burstEnd; % get bursrEnd
sz = size(bEnd);
bInfo = zeros(sz(1),4); % make 4 column to match allInfo
bInfo(:,1) = bStart; % store bStart at first col
bInfo(:,2) = bEnd; % store bStart at 2nd col

% bInfo % for DEBUG

pkLocs = handles.pkLocs; % locations of peaks
pkVals = handles.pkVals; % value of peaks
pkProm = handles.pkProm; % prominance of peaks
pkW = handles.pkW; % width of peaks

wsVarName = handles.wsVarName; % get working variable name
filename = strcat(wsVarName,'_analysis.xlsx'); % save the file as its variable name

% allInfo
% [burstStart, burstEnd, 0,        0]
% [0,          0,        0,        0]
% [pkLocs,     pkVals,   pkProm, pkW]
allInfo = [bInfo;0,0,0,0; pkLocs, pkVals, pkProm, pkW];
xlswrite(filename,allInfo,1,'E2');
disp('As Apr 19, 2016, MATLAB for Mac does not support saving excel file. CSV file is saved instead');
% below do not work due to Mac version MATLAB's issue with  
% A = ['PeakTime','PeakVoltage','PeakProminence','PeakWidth'; pkLocs, pkVals, pkProm, pkW;]
% sheetTitle = ['PeakTime','PeakVoltage','PeakProminence','PeakWidth'];
% xlswrite(filename,sheetTitle,1,'E1'); % 


%% Callbacks
function uiWsVarName_Callback(hObject, eventdata, handles)
% hObject    handle to uiWsVarName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.wsVarName = get(hObject,'String');% get varName in the workspace
guidata(hObject,handles); % update gui handles

function uiMinPkWidth_Callback(hObject, eventdata, handles)
% hObject    handle to uiMinPkWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.minPkW = str2double(get(hObject,'String')); % get min peak width
guidata(hObject,handles); % update gui handles

function uiMinPkProm_Callback(hObject, eventdata, handles)
% hObject    handle to uiMinPkProm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.minPkP = str2double(get(hObject,'String')); % get min peak width
guidata(hObject,handles); % update gui handles

function uiMinIntvl_Callback(hObject, eventdata, handles)
% hObject    handle to uiMinIntvl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.minIntvl = str2double(get(hObject,'String')); % get min resting interval
guidata(hObject,handles); % update gui handles

function uiMinBurst_Callback(hObject, eventdata, handles)
% hObject    handle to uiMinBurst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of uiMinBurst as text
%        str2double(get(hObject,'String')) returns contents of uiMinBurst as a double
handles.minBurst = str2double(get(hObject,'String')); % get min resting interval
guidata(hObject,handles); % update gui handles

%% CreateFcn
% --- Executes during object creation, after setting all properties.
function uiWsVarName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uiWsVarName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function uiMinPkWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uiMinPkWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function uiMinPkProm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uiMinPkProm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function uiMinIntvl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uiMinIntvl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function dataPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: place code in OpeningFcn to populate dataPlot

% --- Executes during object creation, after setting all properties.
function uiMinBurst_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uiMinBurst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
