clear all
close all

data_struct1 = load('..\data\Meas0001rawDLS.mat');
fs = 40000;
data = cell2mat(struct2cell(data_struct1));

DLS1_FL3 = data(:,1);
%DLS2_FOREHEAD = data(:,2);
%DLS3_FR3 = data(:,3);

%SP_FL2 = data(:,4);
%ECG = data(:,5);
%CROSS = data(:,6);
%PPG0 = data(:,7);
%PPG1 = data(:,8);

time = (1/fs) * (0:1:(length(DLS1_FL3)-1));
Wn = 5/(fs/2);
[a,b] = butter(2, Wn, 'low');
DLS1_filtered = filtfilt(b,a,DLS1_FL3);
figure;
plot(time, DLS1_filtered);