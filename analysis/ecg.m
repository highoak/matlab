close all
clear all
disp('started')
data_struct1 = load('..\data\Meas0001rawDLS.mat');
fs = 40000; %sample freq
fn = fs/2; %Nyquist freq
data = cell2mat(struct2cell(data_struct1));
N = length(data);
t = (1/fs) * (0:1:(N-1));

ECG_data = data(:,5);

[~,locs_Rwave] = findpeaks(ECG_data,'MinPeakHeight',3,...
                                    'MinPeakDistance',0.5*fs);
figure
hold on
plot(t,ECG_data);
plot(t(locs_Rwave),ECG_data(locs_Rwave),'rv','MarkerFaceColor','r');
grid on
title('Thresholding Peaks in Signal')
xlabel('Time'); ylabel('Amplitude')
legend('ECG signal','R-wave');

disp('finished');