function varargout = mohsenhrv0(varargin)
% MOHSENHRV0 MATLAB code for mohsenhrv0.fig
%      MOHSENHRV0, by itself, creates a new MOHSENHRV0 or raises the existing
%      singleton*.
%
%      H = MOHSENHRV0 returns the handle to a new MOHSENHRV0 or the handle to
%      the existing singleton*.
%
%      MOHSENHRV0('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MOHSENHRV0.M with the given input arguments.
%
%      MOHSENHRV0('Property','Value',...) creates a new MOHSENHRV0 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mohsenhrv0_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mohsenhrv0_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mohsenhrv0

% Last Modified by GUIDE v2.5 30-Sep-2017 22:18:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mohsenhrv0_OpeningFcn, ...
                   'gui_OutputFcn',  @mohsenhrv0_OutputFcn, ...
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


% --- Executes just before mohsenhrv0 is made visible.
function mohsenhrv0_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mohsenhrv0 (see VARARGIN)

% Choose default command line output for mohsenhrv0
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mohsenhrv0 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mohsenhrv0_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
chkdf=[get(handles.checkbox9,'Value') get(handles.checkbox11,'Value') ...
    get(handles.checkbox10,'Value') get(handles.checkbox12,'Value')...
    get(handles.checkbox13,'Value') get(handles.checkbox14,'Value')...
    get(handles.checkbox32,'Value')...
    get(handles.checkbox17,'Value') get(handles.checkbox15,'Value')...
    get(handles.checkbox16,'Value') get(handles.checkbox33,'Value')...
    get(handles.checkbox18,'Value') get(handles.checkbox19,'Value')...
    get(handles.checkbox34,'Value') get(handles.checkbox35,'Value')...
    get(handles.checkbox24,'Value') get(handles.checkbox20,'Value')...
    get(handles.checkbox21,'Value') get(handles.checkbox22,'Value')...
    get(handles.checkbox23,'Value') ...
    get(handles.checkbox25,'Value') get(handles.checkbox26,'Value')...
    get(handles.checkbox27,'Value')];
blah = get(handles.slcrr,'SelectedObject');
whatever=get(blah,'String');
antyp=get(handles.slcan,'SelectedObject');
antypchs=get(antyp,'String');
if antypchs(1)=='E'
    set(handles.quantiles,'Enable','Off');
    q=1;
elseif antypchs(1)=='Q'
    set(handles.quantiles,'Enable','On');
    q=str2num(get(handles.quantiles,'String'));
elseif antypchs(1)=='S'
    set(handles.quantiles,'Enable','Off');
    q=0;
end
global fullpathnameRR 
if whatever(20)=='u'
% set(handles.edit3,'Enable','On');

RES=importdata(fullpathnameRR);
set(handles.text6,'String',['Loaded ' char(hex2dec('2713'))]);
sl=find(fullpathnameRR=='/' | fullpathnameRR=='\');
if q==1
set(handles.edit6,'String',fullpathnameRR(sl(end)+1:end-4));
elseif q>1
    set(handles.edit6,'String',[fullpathnameRR(sl(end)+1:end-4) 'Q']);
elseif q==0
    set(handles.edit6,'String',[fullpathnameRR(sl(end)+1:end-4) 'C']);
end
% assignin('base','q',q); % to pass to workspace
% assignin('base','RES',RES); % to pass to workspace

% if isfield(RES.CNT.rate,'ECG')
% fs=RES.CNT.rate.ECG;
% else
%     set(handles.edit3,'Enable','On');
%     h = msgbox('Enter sampling frequency and Run again! ');
%     fs=str2num(get(handles.edit3,'String'));
% end

% xECG=RES.CNT.ECG;
if ~exist('RES')
    RES=Res;
end
fgr=fieldnames(RES.CNT.rate);
fs=eval(['RES.CNT.rate' '.' fgr{1,1}])
set(handles.edit3,'String',fs);
% RR_tot_ind=round((RES.HRV.Data.T_RRs{1,1}-RES.CNT.Offset)*fs);
RR_tot_ind=round((RES.HRV.Data.T_RR-RES.CNT.Offset)*fs);

RR_tot=(RR_tot_ind(2:end)-RR_tot_ind(1:end-1))./fs;
RR_tot_ind(1)=[];
rj=find(RR_tot>1.7 | RR_tot<0.5);
RR_tot(rj)=[];
RR_tot_ind(rj)=[];

RR_tot_time=RR_tot_ind./fs;
RRts=spline(RR_tot_ind./fs,RR_tot,1/fs:0.25:RR_tot_ind(end)/fs); % 4 Hz interpolation
RRts_time=1/fs:0.25:RR_tot_ind(end)/fs;
elseif whatever(20)=='i'
    RR_tot_time=importdata(fullpathnameRR);
    set(handles.text6,'String',['Loaded ' char(hex2dec('2713'))]);
    sl=find(fullpathnameRR=='/' | fullpathnameRR=='\');
    set(handles.edit6,'String',fullpathnameRR(sl(end)+1:end-4));
    RR_tot=RR_tot_time(2:end)-RR_tot_time(1:end-1);
    RR_tot_time(1)=[];
    RRts=spline(RR_tot_time,RR_tot,0:0.25:RR_tot_time(end)); % 4 Hz interpolation
    RRts_time=0:0.25:RR_tot_time(end);
elseif whatever(20)=='e'
    RR_tot_ind=importdata(fullpathnameRR);
    set(handles.text6,'String',['Loaded ' char(hex2dec('2713'))]);
    sl=find(fullpathnameRR=='/' | fullpathnameRR=='\');
    set(handles.edit6,'String',fullpathnameRR(sl(end)+1:end-4));
    fs=str2num(get(handles.edit3,'String'));
    RR_tot_time=RR_tot_ind./fs;
    RR_tot=RR_tot_time(2:end)-RR_tot_time(1:end-1);
    RR_tot_time(1)=[];
    RRts=spline(RR_tot_time,RR_tot,0:0.25:RR_tot_time(end)); % 4 Hz interpolation
    RRts_time=0:0.25:RR_tot_time(end);
end
% RRts_time(end) 
global fullpathnameMRK
laur=get(handles.checkbox39,'Value');
newmrk=get(handles.checkbox40,'Value');
if (laur==0 && newmrk==0)
scrs=importdata(fullpathnameMRK);
end
if (laur==1 && newmrk==0)
    fid = fopen(fullpathnameMRK);
    C=textscan(fid,'%s');
    CC=[C{1,1}];
    CC2=CC(2:3:end);
    for i=1:length(CC2)
        scrs(i,1)=str2num(CC2{i,1});
    end
    assignin('base','testscrs',scrs);
end
if (laur==0 && newmrk==1)
mann=importdata(fullpathnameMRK);
scrs=mann.stages;
end

if (laur==1 && newmrk==1)
    h = msgbox('Please select only one or no checkbox for the sleep marker');
end
assignin('base','laur',laur);
assignin('base','newmrk',newmrk);
% assignin('base','scrs',scrs);  % to pass to workspace
set(handles.text7,'String',['Loaded ' char(hex2dec('2713'))]);

uw=str2num(get(handles.undisturbedmin,'String'));
aw=str2num(get(handles.analbin,'String'));
epl=str2num(get(handles.edit4,'String')); % epoch length in sec

[~,c]=size(scrs);
if c==3
mrkr=scrs(:,2);
elseif c==1
    mrkr=scrs;
end
strt=0;
%%% from here
if q>1
    bg=find(mrkr>0 & mrkr<7);
    strt=30*(bg(1)-1);
    mrkr(1:bg(1)-1)=[];
end
mdiv=floor(length(mrkr)/q);
assignin('base','mrkr',mrkr)
assignin('base','mdiv',mdiv)
    p=0;
for qi=1:q
    bg=(qi-1)*mdiv+1;
    ed=qi*mdiv;
    assignin('base','bg',bg)
    assignin('base','ed',ed)
% t=find((mrkr(2:end)-mrkr(1:end-1))~=0);
t=find((mrkr(bg+1:ed)-mrkr(bg:ed-1))~=0);
% t=[0 t' length(mrkr)];
t=bg(1)-1+[0 t' mdiv];
bnd=zeros(length(t)-1,2);
for i=1:length(t)-1
    bnd(i,:)=[t(i)+1 t(i+1)];
end
bmn=(bnd(:,2)-bnd(:,1)+1).*epl./60; %epl=30;
bndstg=mrkr(bnd(:,1));
mintof=find(bmn>(uw+aw)); % on pinchun request, was: find(bmn>(uw+aw) & bndstg~=7);
analepochs_mrk=bndstg(mintof);
analepochs_dur=bmn(mintof);
ContEpoch=floor((analepochs_dur-uw)./aw); %column 5
analepochs_bnd=bnd(mintof,:);
ll=length(analepochs_dur);
% HRmean=zeros(ll,1);HRstd=zeros(ll,1); RRmean=zeros(ll,1); 
% SDNN=zeros(ll,1); SDSD=zeros(ll,1); RMSSD=zeros(ll,1);
% NN50=zeros(ll,1); pNN50=zeros(ll,1);
% 
%  VLF=zeros(ll,1); LF=zeros(ll,1);
%  HF=zeros(ll,1);  Total=zeros(ll,1); 
%  LFovHF=zeros(ll,1);  LFnu=zeros(ll,1); 
%  HFnu=zeros(ll,1);  VLFperc=zeros(ll,1); 
%  LFperc=zeros(ll,1); HFperc=zeros(ll,1); 
%  LFpk=zeros(ll,1);  HFpk=zeros(ll,1); 
%  MPF=zeros(ll,1); 
%  SD1=zeros(ll,1);  SD2=zeros(ll,1); 
%  SD1ovSD2=zeros(ll,1); 
assignin('base','analepochs_dur',analepochs_dur)
assignin('base','bnd',bnd);
assignin('base','epl',epl)

if ~isempty(analepochs_dur)
%     SID=cell(length(analepochs_dur),1);
%     start_time=STD;
%     stage=STD;
for j=1:length(analepochs_dur)
    SID{j+p,1}=fullpathnameRR(sl(end)+1:end-4);
    t1=(analepochs_bnd(j,1)-1)*epl+uw*60;
    start_time{j+p,1}=[num2str(floor((t1+strt)/60)) ':' num2str((t1+strt)-60*floor((t1+strt)/60))]; %col3
    t3=analepochs_bnd(j,2)*epl;
    t2=t1+60*aw*floor((t3-t1)/(aw*60));
%     s1=t1*4;
%     s2=t2*4;
    x=RR_tot(find(RR_tot_time>t1+strt & RR_tot_time<=t2+strt));
    x=reshape(x,1,length(x));
    xts=RRts(find(RRts_time>t1+strt & RRts_time<=t2+strt)); 
    dxx=1000*abs(x(2:end)-x(1:end-1));
%     length(xts)/(4*60)
%     figure;subplot(2,1,1);plot(x);subplot(2,1,2);plot(xts);
    if analepochs_mrk(j)==1
        stage{j+p,1}='N1'; % col4
    elseif analepochs_mrk(j)==2
        stage{j+p,1}='N2';
    elseif analepochs_mrk(j)==3
        stage{j+p,1}='N3';
    elseif analepochs_mrk(j)==5
        stage{j+p,1}='REM';
    elseif analepochs_mrk(j)==0
        stage{j+p,1}='Wake';
    elseif analepochs_mrk(j)==-1
        stage{j+p,1}='NoStage';
    elseif analepochs_mrk(j)==7
        stage{j+p,1}='NoStage';
    end
%     sbjID{j,1}=filenameRR(1:end-4);
    x(find(x>2.5))=[]; x(find(x<0.35))=[];
    HRmean(j+p,1)=mean(60./x);
    RRmean(j+p,1)=mean(x*1000);
    SDNN(j+p,1)=std(x*1000);
    SDSD(j+p,1)=std(dxx);
    RMSSD(j+p,1)=sqrt(sum(dxx.^2)./(length(dxx)-1));
    NN50(j+p,1)=length(find(dxx>50));
    pNN50(j+p,1)=length(find(dxx>50))/length(dxx)*100;
    [px,f]=pyulear(detrend(xts*1000),16,0:1/1000:1,4);
%     figure;plot(f,px)
    VLF(j+p,1)=sum(px(find(f>=0.0033 & f<0.04)))/1000;
    LF(j+p,1)=sum(px(find(f>=0.04 & f<0.15)))/1000;
    HF(j+p,1)=sum(px(find(f>=0.15 & f<=0.4)))/1000;
    Total(j+p,1)=sum(px)/1000;
    LFovHF(j+p,1)=LF(j+p,1)/HF(j+p,1);
    LFnu(j+p,1)=LF(j+p,1)/(LF(j+p,1)+HF(j+p,1));
    HFnu(j+p,1)=HF(j+p,1)/(LF(j+p,1)+HF(j+p,1));
    VLFperc(j+p,1)=100*VLF(j+p,1)/(VLF(j+p,1)+LF(j+p,1)+HF(j+p,1));
    LFperc(j+p,1)=100*LF(j+p,1)/(VLF(j+p,1)+LF(j+p,1)+HF(j+p,1));
    HFperc(j+p,1)=100*HF(j+p,1)/(VLF(j+p,1)+LF(j+p,1)+HF(j+p,1));
    [~,kk]=max(px(find(f>=0.04 & f<0.15)));
    LFpk(j+p,1)=0.04+f(kk);
    [~,kk]=max(px(find(f>=0.15 & f<=0.4)));
    HFpk(j+p,1)=0.15+f(kk);
    MPF(j+p,1)=sum(f.*px')/sum(px);
    
    poincar=[1000*x(2:end)' 1000*x(1:end-1)' ];
    if ~isempty(poincar)
    sd=eig(cov(poincar));
    SD1(j+p,1)=sqrt(sd(1));
    SD2(j+p,1)=sqrt(sd(2));
    SD1ovSD2(j+p,1)=sqrt(sd(1))/sqrt(sd(2));
    else
        SD1(j+p,1)=0;
        SD2(j+p,1)=0;
        SD1ovSD2(j+p,1)=0;
    end
    Epoch(j+p,1)=j;
    quantile(j+p,1)=qi;
    ContEpochs(j+p,1)=ContEpoch(j);
end
% Epoch=p+[1:j]';
p=p+ll;
% global output
% output=table(SID,start_time,Epoch,stage,ContEpoch,RRmean,HRmean,SDNN,SDSD,...
%     RMSSD,NN50,pNN50,VLF,LF,HF,Total,LFovHF,LFnu,HFnu,VLFperc,LFperc,HFperc,LFpk,...
%     HFpk,MPF,SD1,SD2,SD1ovSD2);
% unch=find(chkdf==0);
% output(:,5+unch)=[];
% h = msgbox('Analysis finished. You can export the output now!');
% % assignin('base','output',output);
else
    h = msgbox('Sorry! There is no output by the selected parameters in Step 2. If you want to try 2-4 min analysis, you will need an upgrade to the GUI');
end
end
%for q=0 case
if q==0
    tmpscoring=zeros(1,length(mrkr)*epl*fs);
    for i=1:length(scrs)
        tmpscoring((i-1)*epl*fs+1:i*epl*fs)=scrs(i);
    end
    fsample=fs;
    SleepCycles = [];
    flag    = -1;    nrem1   = 1;
    
    while flag < 0
        nrem1   = nrem1 - 1 + find(ismember(tmpscoring(nrem1:end),[2 3 4 12 13 14]),1,'first');
        rem1    = nrem1 - 1 + find(ismember(tmpscoring(nrem1:end),[5 15]),1,'first');
        
        if isempty(nrem1)
            error('No 15 min NREM period existing');
        elseif (rem1 - nrem1) / (fsample * 60) > 15
            SleepCycles = [SleepCycles; nrem1, rem1];
            flag = 1;
        else
            nrem1 = rem1;
        end
    end
    cycType = 1;
    currPos = rem1 - 1 + find(ismember(tmpscoring(rem1:end),[2 3 4 12 13 14]),1,'first');
    
    while ~isempty(currPos)
        switch cycType
            case 1 % look for next 15 min NREM period
                nextPos = find(ismember(tmpscoring(currPos:end),[5 15]),1,'first');
                if ~isempty(nextPos)
                    nextPos = nextPos + currPos - 1;
                    if (nextPos - currPos) / (fsample * 60) > 15
                        SleepCycles = [SleepCycles; currPos, NaN];
                        cycType = 2;
                        currPos = nextPos;
                    else
                        currPos = nextPos - 1 + find(ismember(tmpscoring(1,nextPos:end),[2 3 4 12 13 14]),1,'first');
                    end
                else
                    currPos = nextPos;
                end
            case 2 % look for next 5 min REM period
                nextPos = find(ismember(tmpscoring(1,currPos:end),[2 3 4 12 13 14]),1,'first');
                if ~isempty(nextPos)
                    nextPos = nextPos + currPos - 1;
                    if (nextPos - currPos) / (fsample * 60) > 5
                        SleepCycles(end,2) = currPos;
                        cycType = 1;
                        currPos = nextPos;
                    else
                        currPos = nextPos - 1 + find(ismember(tmpscoring(1,nextPos:end),[5 15]),1,'first');
                    end
                else
                    currPos = nextPos;
                end
        end
    end
    switch cycType
        case 1 % Was trying to identify NREM but could not find suceeding REM sleep
            lastcyc = SleepCycles(end,2);
            lastpos = lastcyc - 1 + find(ismember(tmpscoring(1,lastcyc:end),[2 3 4 12 13 14]),1,'first');
            if (lastpos - lastcyc) / (fsample * 60) < 5
                SleepCycles(end,2) = NaN;
            else
                SleepCycles = [SleepCycles; lastpos, NaN];
            end
        case 2
            lastcyc = SleepCycles(end,1);
            lastpos = lastcyc - 1 + find(ismember(tmpscoring(1,lastcyc:end),[5 15]),1,'first');
            if (lastpos - lastcyc) / (fsample * 60) < 5
                SleepCycles(end,:) = [];
            else
                SleepCycles(end,2) = lastpos;
            end
    end
    
    numNREMCyc                  = size(SleepCycles,1) - mod(numel(SleepCycles(~isnan(SleepCycles))),2);
    NREMCycles   = SleepCycles(1:numNREMCyc,:);
    
    numREMCyc                   = size(SleepCycles,1) - 1;
    REMCycles    = [SleepCycles(1:numREMCyc,2), SleepCycles(2:numREMCyc+1,1)];
    remnrems=[NREMCycles;REMCycles];
    [~,k]=sort(remnrems(:,1));
    remnrems=remnrems(k,:);
    remnrems(:,2)=remnrems(:,2)-1;
    bgs=1+(remnrems(:,1)-1)./(30*256);
    eds=(remnrems(:,2))./(30*256);
    
%     mdiv=floor(length(mrkr)/q);
    pc=0;
    for ci=1:size(remnrems,1)
        t=find((mrkr(bgs(ci)+1:eds(ci))-mrkr(bgs(ci):eds(ci)-1))~=0);
        % t=[0 t' length(mrkr)];
        t=bgs(ci)-1+[0 t' eds(ci)-bgs(ci)+1];
        bnd=zeros(length(t)-1,2);
        for i=1:length(t)-1
            bnd(i,:)=[t(i)+1 t(i+1)];
        end
        bmn=(bnd(:,2)-bnd(:,1)+1).*epl./60; %epl=30;
        bndstg=mrkr(bnd(:,1));
        mintof=find(bmn>(uw+aw) & bndstg~=7);
        analepochs_mrk=bndstg(mintof);
        analepochs_dur=bmn(mintof);
        ContEpoch=floor((analepochs_dur-uw)./aw); %column 5
        analepochs_bnd=bnd(mintof,:);
        ll=length(analepochs_dur);
        
        if ~isempty(analepochs_dur)
           
            for j=1:length(analepochs_dur)
                SID{j+pc,1}=fullpathnameRR(sl(end)+1:end-4);
                t1=(analepochs_bnd(j,1)-1)*epl+uw*60;
                start_time{j+pc,1}=[num2str(floor((t1+strt)/60)) ':' num2str((t1+strt)-60*floor((t1+strt)/60))]; %col3
                t3=analepochs_bnd(j,2)*epl;
                t2=t1+60*aw*floor((t3-t1)/(aw*60));
                %     s1=t1*4;
                %     s2=t2*4;
                x=RR_tot(find(RR_tot_time>t1+strt & RR_tot_time<=t2+strt));
                x=reshape(x,1,length(x));
                xts=RRts(find(RRts_time>t1+strt & RRts_time<=t2+strt));
                dxx=1000*abs(x(2:end)-x(1:end-1));
                
                if analepochs_mrk(j)==1
                    stage{j+pc,1}='N1'; % col4
                elseif analepochs_mrk(j)==2
                    stage{j+pc,1}='N2';
                elseif analepochs_mrk(j)==3
                    stage{j+pc,1}='N3';
                elseif analepochs_mrk(j)==5
                    stage{j+pc,1}='REM';
                elseif analepochs_mrk(j)==0
                    stage{j+pc,1}='Wake';
                elseif analepochs_mrk(j)==-1
                    stage{j+p,1}='NoStage';
                elseif analepochs_mrk(j)==7
                    stage{j+p,1}='NoStage';
                end
                %     sbjID{j,1}=filenameRR(1:end-4);
                x(find(x>2.5))=[]; x(find(x<0.35))=[];
                HRmean(j+pc,1)=mean(60./x);
                RRmean(j+pc,1)=mean(x*1000);
                SDNN(j+pc,1)=std(x*1000);
                SDSD(j+pc,1)=std(dxx);
                RMSSD(j+pc,1)=sqrt(sum(dxx.^2)./(length(dxx)-1));
                NN50(j+pc,1)=length(find(dxx>50));
                pNN50(j+pc,1)=length(find(dxx>50))/length(dxx)*100;
                [px,f]=pyulear(detrend(xts*1000),16,0:1/1000:1,4);
                %     figure;plot(f,px)
                VLF(j+pc,1)=sum(px(find(f>=0.0033 & f<0.04)))/1000;
                LF(j+pc,1)=sum(px(find(f>=0.04 & f<0.15)))/1000;
                HF(j+pc,1)=sum(px(find(f>=0.15 & f<=0.4)))/1000;
                Total(j+pc,1)=sum(px)/1000;
                LFovHF(j+pc,1)=LF(j+pc,1)/HF(j+pc,1);
                LFnu(j+pc,1)=LF(j+pc,1)/(LF(j+pc,1)+HF(j+pc,1));
                HFnu(j+pc,1)=HF(j+pc,1)/(LF(j+pc,1)+HF(j+pc,1));
                VLFperc(j+pc,1)=100*VLF(j+pc,1)/(VLF(j+pc,1)+LF(j+pc,1)+HF(j+pc,1));
                LFperc(j+pc,1)=100*LF(j+pc,1)/(VLF(j+pc,1)+LF(j+pc,1)+HF(j+pc,1));
                HFperc(j+pc,1)=100*HF(j+pc,1)/(VLF(j+pc,1)+LF(j+pc,1)+HF(j+pc,1));
                [~,kk]=max(px(find(f>=0.04 & f<0.15)));
                LFpk(j+pc,1)=0.04+f(kk);
                [~,kk]=max(px(find(f>=0.15 & f<=0.4)));
                HFpk(j+pc,1)=0.15+f(kk);
                MPF(j+pc,1)=sum(f.*px')/sum(px);
                
                poincar=[1000*x(2:end)' 1000*x(1:end-1)' ];
                if ~isempty(poincar)
                    sd=eig(cov(poincar));
                    SD1(j+pc,1)=sqrt(sd(1));
                    SD2(j+pc,1)=sqrt(sd(2));
                    SD1ovSD2(j+pc,1)=sqrt(sd(1))/sqrt(sd(2));
                else
                    SD1(j+pc,1)=0;
                    SD2(j+pc,1)=0;
                    SD1ovSD2(j+pc,1)=0;
                end
                Epoch(j+pc,1)=j;
                if analepochs_mrk(j)==5
                cycles{j+pc,1}=['REM' num2str(ceil(ci/2))];
                else
                    cycles{j+pc,1}=['NREM' num2str(ceil(ci/2))];
                end
                ContEpochs(j+pc,1)=ContEpoch(j);
            end
            % Epoch=p+[1:j]';
            pc=pc+ll;
            % global output
            % output=table(SID,start_time,Epoch,stage,ContEpoch,RRmean,HRmean,SDNN,SDSD,...
            %     RMSSD,NN50,pNN50,VLF,LF,HF,Total,LFovHF,LFnu,HFnu,VLFperc,LFperc,HFperc,LFpk,...
            %     HFpk,MPF,SD1,SD2,SD1ovSD2);
            % unch=find(chkdf==0);
            % output(:,5+unch)=[];
            % h = msgbox('Analysis finished. You can export the output now!');
            % % assignin('base','output',output);
%         else
%             h = msgbox('Sorry! There is no output by the selected parameters in Step 2. If you want to try 2-4 min analysis, you will need an upgrade to the GUI');
        end
    end


end

%
%making output
global output
if q~=0
output=table(SID,quantile,start_time,Epoch,stage,ContEpochs,RRmean,HRmean,SDNN,SDSD,...
    RMSSD,NN50,pNN50,VLF,LF,HF,Total,LFovHF,LFnu,HFnu,VLFperc,LFperc,HFperc,LFpk,...
    HFpk,MPF,SD1,SD2,SD1ovSD2);
else
    output=table(SID,cycles,start_time,Epoch,stage,ContEpochs,RRmean,HRmean,SDNN,SDSD,...
    RMSSD,NN50,pNN50,VLF,LF,HF,Total,LFovHF,LFnu,HFnu,VLFperc,LFperc,HFperc,LFpk,...
    HFpk,MPF,SD1,SD2,SD1ovSD2);
end
unch=find(chkdf==0);
output(:,6+unch)=[];
h = msgbox('Analysis finished. You can export the output now!');








function undisturbedmin_Callback(hObject, eventdata, handles)
% hObject    handle to undisturbedmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of undisturbedmin as text
%        str2double(get(hObject,'String')) returns contents of undisturbedmin as a double



% --- Executes during object creation, after setting all properties.
function undisturbedmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to undisturbedmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function analbin_Callback(hObject, eventdata, handles)
% hObject    handle to analbin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of analbin as text
%        str2double(get(hObject,'String')) returns contents of analbin as a double


% --- Executes during object creation, after setting all properties.
function analbin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to analbin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
crd=date;
global filenameRR pathnameRR
global fullpathnameRR

% if datenum(crd)<737884
[filenameRR,pathnameRR] = uigetfile({'*.mat';'*.txt';'*.*'},'File Selector'); 
% filenameRR
fullpathnameRR=strcat(pathnameRR,filenameRR); 
% text= fileread(fullpathname);
set(handles.text6,'String',fullpathnameRR); % showing full path name
% else
%     set(handles.text6,'BackgroundColor','red');
%     set(handles.text6,'ForegroundColor','yellow');
%     wm=[84 104 101 32 115 111 102 116 119 97 114 101 32 105 115 32 101 ...
%         120 112 105 114 101 100 46 32 65 115 107 32 102 111 114 32 97 ...
%         110 111 116 104 101 114 32 118 101 114 115 105 111 110];
%     set(handles.text6,'String',char(wm));
%     set(handles.pushbutton1,'String',char([58 40]));
%     
% end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filenameMRK,pathnameMRK] = uigetfile({'*.txt';'*.vmrk';'*.mat';'*.*'},'File Selector');
global fullpathnameMRK
fullpathnameMRK=strcat(pathnameMRK,filenameMRK);
% text= fileread(fullpathname);
set(handles.text7,'String',fullpathnameMRK); % showing full path name



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in slcrr.
function slcrr_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in slcrr 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
blah = get(handles.slcrr,'SelectedObject');
whatever=get(blah,'String');
if whatever(20)~='e'
set(handles.edit3,'Enable','Off');
% set(handles.edit3,'Enable','On');
end
if whatever(20)=='e'
 set(handles.edit3,'Enable','On');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9


% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox10


% --- Executes on button press in checkbox11.
function checkbox11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox11


% --- Executes on button press in checkbox12.
function checkbox12_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox12


% --- Executes on button press in checkbox13.
function checkbox13_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox13


% --- Executes on button press in checkbox14.
function checkbox14_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox14


% --- Executes on button press in checkbox15.
function checkbox15_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox15


% --- Executes on button press in checkbox16.
function checkbox16_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox16


% --- Executes on button press in checkbox17.
function checkbox17_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox17


% --- Executes on button press in checkbox18.
function checkbox18_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox18


% --- Executes on button press in checkbox19.
function checkbox19_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox19


% --- Executes on button press in checkbox20.
function checkbox20_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox20


% --- Executes on button press in checkbox21.
function checkbox21_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox21


% --- Executes on button press in checkbox22.
function checkbox22_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox22


% --- Executes on button press in checkbox23.
function checkbox23_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox23


% --- Executes on button press in checkbox24.
function checkbox24_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox24


% --- Executes on button press in checkbox25.
function checkbox25_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox25


% --- Executes on button press in checkbox26.
function checkbox26_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox26


% --- Executes on button press in checkbox27.
function checkbox27_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox27


% --- Executes on button press in checkbox28.
function checkbox28_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox28


% --- Executes on button press in checkbox29.
function checkbox29_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox29


% --- Executes on button press in checkbox30.
function checkbox30_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox30


% --- Executes on button press in checkbox31.
function checkbox31_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox31


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3



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


% --- Executes on key press with focus on edit5 and none of its controls.
function edit5_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox33.
function checkbox33_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox33


% --- Executes on button press in checkbox34.
function checkbox34_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox34


% --- Executes on button press in checkbox35.
function checkbox35_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox35


% --- Executes on button press in checkbox32.
function checkbox32_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox32
%%%



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
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
csvok=get(handles.checkbox37,'Value');
matok=get(handles.checkbox36,'Value');
txtok=get(handles.checkbox38,'Value');
global output pathnameRR

outname=get(handles.edit6,'String');
if csvok==1
writetable(output,[pathnameRR outname '_table.csv'],'Delimiter',',','QuoteStrings',true);
end
if txtok==1
writetable(output,[pathnameRR outname '_table.txt'],'Delimiter',' ');
end
if matok==1
%     assignin('base','output',output);
    save([pathnameRR outname '_table'],'output');
end


% --- Executes on button press in checkbox36.
function checkbox36_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox36


% --- Executes on button press in checkbox37.
function checkbox37_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox37


% --- Executes on button press in checkbox38.
function checkbox38_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox38


% --- Executes during object creation, after setting all properties.
function text14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in checkbox39.
function checkbox39_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox39


% --- Executes on button press in checkbox40.
function checkbox40_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox40


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5


% --- Executes on button press in radiobutton8.
function radiobutton8_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton8


% --- Executes on button press in radiobutton7.
function radiobutton7_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton7


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton6



function quantiles_Callback(hObject, eventdata, handles)
% hObject    handle to quantiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of quantiles as text
%        str2double(get(hObject,'String')) returns contents of quantiles as a double


% --- Executes during object creation, after setting all properties.
function quantiles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to quantiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in slcan.
function slcan_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in slcan 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
antyp=get(handles.slcan,'SelectedObject');
antypchs=get(antyp,'String');
if antypchs(1)=='E'
    set(handles.quantiles,'Enable','Off');
    q=1;
elseif antypchs(1)=='Q'
    set(handles.quantiles,'Enable','On');
    q=str2num(get(handles.quantiles,'String'));
elseif antypchs(1)=='S'
    set(handles.quantiles,'Enable','Off');
    q=0;
end
