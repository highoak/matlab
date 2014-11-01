clear all
close all

data_struct1 = load('..\data\Meas0001rawDLS.mat');
fs = 40000; %sample freq
fn = fs/2; %Nyquist freq
data = cell2mat(struct2cell(data_struct1));

DLS1_FL3 = data(:,1);
%DLS2_FOREHEAD = data(:,2);
%DLS3_FR3 = data(:,3);

SP_FL2 = data(:,4);
ECG = data(:,5);
%CROSS = data(:,6);
%PPG0 = data(:,7);
%PPG1 = data(:,8);
N = length(data);
time = (1/fs) * (0:1:(N-1));

clear data

Wn = [10000,fn-1] / fn;
[b,a] = butter(2, Wn, 'bandpass');
DLS1_bandpassed = filter(b,a,DLS1_FL3);


Fr=100; % Fr=10-500Hz 
NumbIntervals=Fr*N/fs; % Total number of pulsewave points


DOSE = fix(N / NumbIntervals);
for k = 1:NumbIntervals
    interval_indexes=1+(k-1)*DOSE:k*DOSE;
    HEST(k,:)=wfbmesti(DLS1_FL3(interval_indexes));
    PULSEwave(k) = std(DLS1_bandpassed(interval_indexes));
end
time_dls = (1/Fr)*(0:1:(length(PULSEwave)-1));


Wn = 5/(Fr/2);
[b,a]=butter(2,Wn,'low');
dls0=filtfilt(b,a,PULSEwave);
dls0=dls0 + abs(min(dls0));

for km1=1:3
    dls00(:,km1)=filter(b,a,HEST(:,km1));
    dls00(:,km1)=dls00(:,km1)+abs(min(dls00(:,km1)));
end

figure;
plotyy(time, ECG, time_dls, dls00);
