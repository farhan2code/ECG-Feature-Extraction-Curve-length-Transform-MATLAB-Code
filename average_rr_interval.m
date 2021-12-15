clear all
load ('100mml2.mat')              %Data from MIT-BIH Database
ecg = (val)/200;
fs = 360;                         %sampling frequency of 360 Hz
Ts = 1/fs;
t = (0:length(ecg) - 1)/fs;       %time variable

ax(1) = subplot(3,1,1);
plot(t,ecg,'linewidth',1.5);
xlabel('Time (in seconds)');
ylabel('Amplitude');
title('Raw ECG Signal');
grid


%%%%%%%%%%%%%    Bandpass filter (5-40 Hz)     %%%%%%%%%%%%%%%%%%%%%%%%%%%%

fl=5;                                                                      % cuttoff low frequency to get rid of baseline wander
fh=40;                                                                     % cuttoff frequency to discard high frequency noise
Wn=[fl fh]*2/fs;                                                           % cutt off based on fs                                                                    
[a,b] = butter(3,Wn);                                                      % bandpass filtering
ecg_h = filtfilt(a,b,ecg);
ecg_h = ecg_h/ max( abs(ecg_h));

ax(2) = subplot(3,1,2);
plot(t,ecg_h,'linewidth',1.5);
xlabel('Time (in seconds)');
ylabel('Amplitude');
title('Band Pass Filtered');
grid


%%%%%%%%%%%%%%%%%  Absolute CLT (ACLT) signal  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(ecg);
w = 18;
for i = 2+w:N
    for k = i-w:i
        dy = ecg_h(k)-ecg_h(k-1);
        l(k+w-i+1) = abs(Ts.^2 + abs(4.*(dy^2)));
    end
    aclt(i) = sum(l);
end

aclt = circshift(aclt, [0 -fix(w/2)]);


ax(3) = subplot(3,1,3);
plot(t,aclt,'linewidth',1.5);
xlabel('Time (in seconds)');
ylabel('Amplitude');
title('Absolute CLT signal');
grid

linkaxes(ax,'x');



%%%%%%%%%%%%%%%%%   Peak Detection of ACLT Signal   %%%%%%%%%%%%%%%%%%%%%%%

maxtab = [];
mintab = [];

aclt = aclt(:);

mn = Inf; mx = -Inf;
mnpos = NaN; mxpos = NaN;
mnpos_s = NaN; mxpos_s = NaN;    

lookformax = 1;
delta = 0.5;
x = (1:length(ecg));

for i=1:length(aclt)
  this = aclt(i);
  if this > mx, mx = this; mxpos = t(i); mxpos_s = x(i); end
  if this < mn, mn = this; mnpos = t(i); mnpos_s = x(i); end
  
  if lookformax
    if this < mx-delta
      maxtab = [maxtab ; mxpos_s mxpos];
      mn = this; mnpos = t(i); mnpos_x = x(i);
      lookformax = 0;
    end  
  else
    if this > mn+delta
      mintab = [mintab ; mnpos_s mnpos];
      mx = this; mxpos = t(i); mxpos_s = x(i);
      lookformax = 1;
    end
  end
end

qrs_time = [];
qrs_sample = [];

qrs_sample = maxtab(:,1);    %this vector contains the samples where there is a QRS Complex peak
qrs_time = maxtab(:,2);      %this vector contains the times where there is a QRS Complex peak




%%%%%%%%%%%%%   Peak detection of Raw ECG Signal (Fine tuning to find correct R peaks)   %%%%%%%%%%%%%%%%%%%%%%%%


maxtab1 = [];
mintab1 = [];

ecg = ecg(:);

mn = Inf; mx = -Inf;
mnpos = NaN; mxpos = NaN;
mnpos_s = NaN; mxpos_s = NaN; 

lookformax = 1;
delta1 = 0.5;


for n=1:length(qrs_sample)
    k = qrs_sample(n);
    
    for i = k-5:k+5
        this = ecg(i);
        if this > mx, mx = this; mxpos = t(i); mxpos_s = x(i); end
        if this < mn, mn = this; mnpos = t(i); mnpos_s = x(i); end
        
        if lookformax
            if this < mx-delta
                maxtab1 = [maxtab1 ; mxpos_s mxpos];
                mn = this; mnpos = t(i); mnpos_x = x(i);
                lookformax = 0;
            end
        else
            if this > mn+delta
                mintab1 = [mintab1 ; mnpos_s mnpos];
                mx = this; mxpos = t(i); mxpos_s = x(i);
                lookformax = 1;
            end
        end
    end
end


qrs_time_raw = [];
qrs_sample_raw = [];
r_peak = [];

qrs_sample_raw = maxtab1(:,1);    %this vector contains the samples where there is a QRS Complex peak in RAW ECG signal
qrs_time_raw = maxtab1(:,2);      %this vector contains the times where there is a QRS Complex peak in RAW ECG signal

for n = 1:length(qrs_sample_raw)
    i = qrs_sample_raw(n);
    r_peak(n) =  ecg(i);
end


figure(2)

plot(t,ecg);
ylim([-2.4 2])
xlabel('Time (in seconds)');
ylabel('Amplitude');
title('Detected R Peaks');
grid
hold on;

plot(qrs_time_raw, r_peak, 'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',8);

hold off;

legend({'Raw ECG Signal','R peaks'}, 'Location','southeast');
lgd = legend;
lgd.FontSize = 14;




%%%%%%%%%%%%%% Calculation of R-R interval and Heart Rate    %%%%%%%%%%%%%%
% 
% 
% beat_count = length(qrs_time);
% N = length(ecg);
% duration_in_seconds = N/fs; 
% duration_in_minutes = duration_in_seconds/60;
% BPM = beat_count/duration_in_minutes;
% 
% disp('Heart Rate (in BPM)');
% disp(BPM);
% 
% 
% 
% p = length(qrs_time);
% avg_rr_interval = (qrs_time(p) - qrs_time(1))/beat_count;
% disp('Average RR interval (in s)');
% disp(avg_rr_interval);

 