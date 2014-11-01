clear  %all
close all
NumbIntervals=221; % 441 --> 100Hz
Np=NumbIntervals;
t2_text=[];
AVGDLS=[]; STDDLS=[]; n17=0; StatDLS=[];
n=0;
tic
LL=[];
freq=40000;  %44100;   %42600; 96000; 50000; 48000;         % Sampling frequncy ############################
fn=freq/2;               % Nyquist Frequency
k105=1;
%************** Preferense **************
k105=1; % k105=1 - DLS and Press, k105=0 - DLS only
%##########################################################

% DLS transfer
'DLS files';
Exp_PathS=struct('name','xxx');
%Exp_PathS(1).name='F:\BP2012\SKsystem6\LabTest\';
%Exp_PathS(1).name='F:\BP2012\SKsystem6\CurrentTest\';
%Exp_PathS(1).name='D:\DataElfior\DeltranTransducer\Ymode_Msens-FingL2tip_DLSphalange2-OcclSensor\Test08_260614_DLSstation\HomeTest_260614\';
Exp_PathS(1).name='..\data';
%Exp_PathS(1).name='D:\BP2013a\DLSCurrentTest\';
%Exp_PathS(1).name='D:\DialysisTestStand160912\April 03i09 test1\';
PathI=1;   %Select the work dir path!!!

while ~exist(Exp_PathS(PathI).name) && PathI<=length(Exp_PathS)
    PathI=PathI+1;
end
Exp_Path=Exp_PathS(PathI).name;


cd(Exp_Path);
%####################################################
[f_name,f_path]=uigetfile('*.*');   % first and second group of files, initial input
cd(f_path)
fn_len=length(f_name);
data_path=f_path;
File_name=f_name;
if 0
    File_name12='DLS';
    
    test_files2=dir(fullfile(f_path,[File_name12,'*.txt']));
end
if 1
    %File_name12='Leo';
    File_name12='Meas';
    %File_name12='Meir';
    test_files2=dir(fullfile(f_path,[File_name12,'*.mat']));
end
%###############################################################################

DbgFlag=0;



f_name22=[];


if exist(f_path)
    LoopFlag_test=1;
    
    
    
    d_struc1 = size(test_files2(:,1));
    d_struc=d_struc1(1);
    
    %f_date1=test_files(1).date;
    %date_len=length(f_date1);
    if k105==1
        for k20=1:d_struc
            if File_name==test_files2(k20).name
                k21=k20;
                break;
            end
        end
        
    end
    
    for k=k21:d_struc    % d_struc - auto-number of files, 10 - manual number of files
        n=n+1;
        
        
        
        
        % Input 2 (ASCII) DLS files from device
        %*******************************************
        f_name22=[f_name22;test_files2(k).name];
        MeasFileName2=[data_path,test_files2(k).name];
        NumName22=f_name22(n,4:6);
        f_name22(n,:)
        
        if 1  % ASCII data and structure binary data
            
            
            %G= dlmread(MeasFileName2); %DLS file
            
            if 0
                G1=load(MeasFileName2); % DLS2 file (Fs2=44100Hz)
                c= struct2cell(G1);
                data = cell2mat(c);
                DLS2=data.signals.values;
            end
            %if 1
            G1=load(MeasFileName2); % PPG&DLS file (Fs=100000Hz)
            c= struct2cell(G1);
            %data = struct2cell(G1);
            data = cell2mat(c);
            clear G1 c
            
        end
        %************** PPD & SP *************
        if 1
            L = length(data);
            T=L/freq;
            TimeData=(1/freq)*(0:1:L-1);
            if 0
                fcDLS0=5;[b5,a5]=butter(2,fcDLS0/(fs/2), 'low');
                data(:,7)=filtfilt(b5,a5,data(:,7));
                data(:,8)=filtfilt(b5,a5,data(:,8));
            end
            figure;
            plotyy(TimeData,data(:,4),TimeData,data(:,5)); % SP and ECG
            %plot(TimeData,data(:,4:5)); % SP and ECG
            %plotyy(TimeData,data(:,4),TimeData,[data(:,7),data(:,8)]); %SP and PPG
            %plotyy(TimeData,data(:,5),TimeData,[data(:,7),data(:,8)]); %ECG and PPG
            SP2=data(:,4);
            ECG2=data(:,5);
        end
        dls3sens=[];
        dlsHsens=[];
        for t2=1:3
            if 0
                if t2==1
                    %S=(data1)';
                    %clear data1
                    S=(data.dls1)';
                    clear data.dls1
                elseif t2==2
                    %S=(data2)';
                    %clear data2
                    S=(data.dls2)';
                    clear data.dls2
                elseif t2==3
                    %S=(data3)';
                    %clear data3
                    S=(data.dls3)';
                    clear data.dls3
                end
            end
            %S=data(t2,:); LL=0; % 50kHz
            S=data(:,t2); LL=0;
            t2_text=[t2_text,num2str(t2)];
            %##########################################################
            
            
            clear MeasFileName2
            
            %clear G
            
            %---------------------- Processing -------------------------
            Signal=S;  %-mean(S); % reduces the mean
            L = length(Signal); % signal length
            Signal2=Signal';
            %-----------------------------------  End of FILE READING SESSION--------------------------------
            clear S
            
            
            
            %T111=mod(L,44408.16325);
            %freq=44408.16325;
            T=L/freq;% Measure Time
            
            
            %    Calculated paramters
            Fr=100; % Fr=10-500Hz   ##################################################
            NumbIntervals=Fr*T ; % Total number of pulsewave points
            
            % -----------------   End of USers Constants
            
            %     PULSEWAVE CALUCLATIONS and PLOT ###############################
            clear DOSE
            DOSE=fix(L/NumbIntervals);  % number of points of each  pulsewave interval
            
            if 1
                %FCH=[1,1000;1000,4000;5000,10000;10000,20000];  % 4 paramters for bandpasss filters intervals [1,200;1,1000;1000,4000;10000,20000];
                %FCH=[5000,10000];
                %FCH=[0.001,1000];
                FCH=[10000,fn-1];
                SIGNALafterFilter=[];
                for i=1:1
                    q=FCH(i,:);[b1S,a1S]=butter(2,q/fn, 'bandpass');
                    SIGNALafterFilter(i,:)=filter(b1S,a1S,Signal);      % Filter Paramters for the bandPass calulations
                end
            end
            %  ------ Transforms the signal to measurement points
            if 0
                SIGNALafterFilter=Signal2; % 1 - non filtered DLS signal
            end
            
            % calculates the STD for each interval
            VH=[]; PULSEwave=[]; VELOC=[];VELOC1=[];PULSEwave1=[]; VC1=[]; T12=[]; HEST0=[]; akm=1; VH1=[];
            for k=1:NumbIntervals
                
                OS=Signal(1+(k-1)*DOSE:k*DOSE); % Original signal before filtration
                for i=1:1
                    interv1=1+(k-1)*DOSE:k*DOSE;
                    VH(i,:)=SIGNALafterFilter(i,interv1);  % STD for each filter + removement of mean
                    %PULSEwave1(i,k)=sum(std((VH(i,:))));
                    PULSEwave1(i,k)=std((VH(i,:)));
                    if i==1 && akm==1
                        VH1(i,:)=Signal2(i,interv1);
                        %VH1(i,:)=Signal2(interv1,i)';
                        HEST0(i,k,:)=wfbmesti((VH1(i,:)));
                        
                    end
                end  % loop by i
                %                 T12(k)=std(OS); %std(OS);
                %                 VELOC1(k)=PULSEwave1(4,k)/PULSEwave1(2,k); % /PULSEwave1(2,k);
                %                 VELOC22(k)=PULSEwave1(4,k);
            end
            clear Signal Signal2 SIGNALafterFilter
            Len00=length(PULSEwave1);
            %             T13=T12';
            %             Vc0=PULSEwave1';
            %             VC1=VELOC1';
            %             Re0=VC1./T13;
            
            PULSEwave=PULSEwave1;
            if akm==1
                HEST=HEST0;
            end
            Len0=length(PULSEwave);
            PULSEwave00=PULSEwave';
            TimeDLS=(1/Fr)*(0:1:Len0-1);
            %*************************************************
            fs  = Fr;
            
            %dls0=[]; dls00=[]; % with LPF
            dls0=PULSEwave1; dls00=HEST; % without LPF
            if 1 % final LPF for DLSPulse
                %CutOffF(1,:)=[0.001,5]/(fs/2);  %[0.001,4]/(fs/2);
                %[b,a]=butter(2,CutOffF(1,:)); % order=2
                fcDLS0=5;[b,a]=butter(2,fcDLS0/(fs/2), 'low');
                for km1=1:1
                    %dls0=[dls0;filter(b,a,PULSEwave(km1,:))];
                    dls0(km1,:)=filtfilt(b,a,PULSEwave(km1,:));
                    dls0(km1,:)=dls0(km1,:)+abs(min(dls0(km1,:)));
                end
                if akm==1
                    %CutOffF1(1,:)=[0.001,5]/(fs/2); %[0.001,3]/(fs/2);
                    %[b1,a1]=butter(2,CutOffF1(1,:)); % order=2
                    fcDLS0=5;[b1,a1]=butter(2,fcDLS0/(fs/2), 'low');
                    if 1
                        for km1=1:3
                            dls00(1,:,km1)=filter(b1,a1,HEST(1,:,km1));
                            %dls00=[dls00;filtfilt(b1,a1,HEST(1,:,km1))];
                            if km1==3 && (any(isnan(dls00(1,:,km1)))==1 || any(isinf(dls00(1,:,km1)))==1)
                                akm=0;
                            end
                            dls00(1,:,km1)=dls00(1,:,km1)+abs(min(dls00(1,:,km1)));
                        end
                    end
                    if 0
                        for km1=1:3
                            dls00(1,:,km1)=HEST(1,:,km1);
                        end
                    end
                end
            end
            if 0
                T14=(filter(b,a,T12))';
                %Re2=filter(b,a,Re1);
                Re2=(slidingavg(Re1',20))'; % (Re1',20)
            end
            
            %dls50=PULSEwave(1:4,:)';
            %dls0(1,:)=PULSEwave(1,:);
            dls=dls0(1,:); %dls2=dls0(2,:);
            Ldls=1:length(dls); % dls=dls0(4,:);
            figure; plot(TimeDLS,dls0(1,:)); pause(1)
            %             figure(2); plot(Ldls,dls0(2,:)); pause(1)
            %             figure(3); plot(Ldls,dls0(3,:)); pause(1)
            %             figure(4); plot(Ldls,dls0(4,:)); pause(1)
            %             title(['BPF3, DLSsens # ',t2_text],'FontSize',12,'FontWeight','bold');
            %             figure(5); plot(Ldls,dls0(1,:),Ldls,dls0(2,:),Ldls,dls0(3,:),Ldls,dls0(4,:));
            %pause(1)
            if 0   %akm==1
                figure; plot(TimeDLS,dls00(1,:,1),TimeDLS,dls00(1,:,2),TimeDLS,dls00(1,:,3));
                title(['Hurst, DLSsens # ',t2_text],'FontSize',12,'FontWeight','bold');
            end
            pause(1)
            LL=(1:length(dls));
            dls3sens=[dls3sens,dls0(1,:)'];
            dlsHsens=[dlsHsens,dls00(1,:,1)'];
            %****************************************
            
            WorkTime1 = toc;
            WorkTime = round(toc/60);
            disp('     ***');
            TextStr = sprintf('     Work TimeRound = %d (min)', WorkTime);
            TextStr1 = sprintf('     Work Time = %d (sec)', WorkTime1);
            disp(TextStr);
            disp('     ***');
            disp(TextStr1);
            disp('     ******');
            a=12;
        end
        %pause(1)
        %PPdec=decimate(data(:,4),2500);
        %PP=-62.35242+149.45759.*PPdec;
        %pause(1)

        clear data
        Avg3DLS=median(dls3sens);
        std3DLS=std(dls3sens);
        AvgH3DLS=median(dlsHsens);
        stdH3DLS=std(dlsHsens);
        StatDLS=[StatDLS;[Avg3DLS(1),std3DLS(1),Avg3DLS(2),std3DLS(2),Avg3DLS(3),std3DLS(3), ...
            AvgH3DLS(1),stdH3DLS(1),AvgH3DLS(2),stdH3DLS(2),AvgH3DLS(3),stdH3DLS(3)]];
        AVGDLS=[AVGDLS;[Avg3DLS,AvgH3DLS]];
        STDDLS=[STDDLS;[std3DLS,stdH3DLS]];
        t2_text1=[t2_text(1),'-',t2_text(end)]; n17=n17+1; n17txt=num2str(n17);
        figure; plot(TimeDLS,dls3sens(:,1),TimeDLS,dls3sens(:,2),TimeDLS,dls3sens(:,3));
        title(['meas # ',n17txt,', DLS, DLSsens # ',t2_text1],'FontSize',12,'FontWeight','bold'); %hold on
        %figure(8); plot(Ldls',dls3sens(:,1),Ldls',dls3sens(:,2));
        %a=12;
        SP1=downsample(SP2,400);
        ECG1=downsample(ECG2,400);
        figure; 
        plotyy(TimeDLS,[dls3sens(:,1),dls3sens(:,2),dls3sens(:,3)],TimeDLS,[SP1,ECG1]);
        %title(['meas # ',n17txt,', DLS, DLSsens # ',t2_text1],'FontSize',12,'FontWeight','bold'); %hold on
        %figure(9); plot(TimeDLS,dlsHsens(:,1),TimeDLS,dlsHsens(:,2),TimeDLS,dlsHsens(:,3));
        %figure(9); plot(Ldls',dlsHsens(:,1),Ldls',dlsHsens(:,2));
        %pause(1)
        a=12;
        if 0
        AVGout=AvgH3DLS(1); STDout=stdH3DLS(1);  % 2 - "Hurst" 2nd probe
        AVGtxt=num2str(AVGout); %STDtxt=num2str(STDout1);
        STD1=STDout/AVGout*100; STD1txt=num2str(STD1);
        n17=n17+1; n17txt=num2str(n17);
        figure;
        %grid off
        errorbar(n17,AVGout,STDout,'--ob','MarkerSize',12)
        title(['meas # ',n17txt,',  Mean&Std(dlsH1) = ',AVGtxt,' & ',STD1txt,' %'],'FontSize',12,'FontWeight','bold')
        %grid on
        hold on
        figure; plot(TimeDLS,dls0(1,:),TimeDLS,PPdec); pause(1)
        a=12;
        pause(1)
        end
        break
    end
    
else
    LoopFlag_test=0;
end
data_path_result=[data_path,'Result','\'];
save([data_path_result,'StatDLS.dat'],'StatDLS','-ASCII');

if 0
    
    data_path_result=[data_path,'Result','\'];
    save([data_path_result,'PressDiastoleHarm.dat'],'BPdiastHarm','-ASCII');
    save([data_path_result,'PressSystoleHarm.dat'],'BPsystHarm','-ASCII');
    save([data_path_result,'Rq.dat'],'Rq0','-ASCII');
    save([data_path_result,'PressDiastoleEnv.dat'],'BPdiastEnv','-ASCII');
    save([data_path_result,'PressSystoleEnv.dat'],'BPsystEnv','-ASCII');
    %save([data_path_result,'Rq.dat'],'Rq0','-ASCII');
end
%------------------------------------------------------------------------------------------------
WorkTime1 = toc;
WorkTime = round(toc/60);
disp('     ***');
TextStr = sprintf('     Work TimeRound = %d (min)', WorkTime);
TextStr1 = sprintf('     Work Time = %d (sec)', WorkTime1);
disp(TextStr);
disp('     ***');
disp(TextStr1);
disp('     ******');
