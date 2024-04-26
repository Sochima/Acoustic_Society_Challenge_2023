%Author: Sochima Omenkeukwu
% This script gets data from 3 hydrophones and uses the data to:
% 1. Plot the spectrogram of one of the hydrophes
% 2. Estimate diver's breath rate
% 3. Estimate time diver is closest to the middle hydrophone
% 4. Estimates diver's altitude
% 5. Estimate diver's swimming speed
%
%% Task 1a
% Reading in audio data for each mic
[N_Mic,Fs1] = audioread("../Audio_Signal/SCP23_ Hyd_N.wav");
[O_Mic,Fs2] = audioread("../Audio_Signal/SCP23_Hyd_O.wav");
[P_Mic,Fs3] = audioread("../Audio_Signal/SCP23_Hyd_P.wav");

%Speed of sound in water
sound_speed = 1520;

%getting the frame length for the audio data
framelen = 2048;

%window to be applied to data
window = hamming(framelen);
%overlap for frames
overlap = 0;
%zero padding 
zp = 0;
nfft = zp*framelen;

%Getting the length of data for each mic
n_num = length(N_Mic);
o_num = length(O_Mic);
p_num = length(P_Mic);

%getting the number of frames for each mic
n_frames1 = round(n_num/framelen);
n_frames2 = round(o_num/framelen);
n_frames3 = round(p_num/framelen);

%generating the time vector for each mic
tt1 = (1:n_num)*1/Fs1;
tt2 = (1:o_num)*1/Fs2;
tt3 = (1:p_num)*1/Fs3;

%Calculating the time frame in seconds for each mic
t_frame1 = framelen/Fs1;
t_frame2 = framelen/Fs2;
t_frame3 = framelen/Fs3;


%Calculating
f_res1 = Fs1/framelen;
f_res2 = Fs2/framelen;
f_res3 = Fs3/framelen;

figure()
plot(tt1,N_Mic)
xlabel("Time(s)")
ylabel("Amplitude")
title("Time Domain Plot of N Mic")
% 
% 
figure()
plot(tt2,O_Mic)
xlabel("Time(s)")
ylabel("Amplitude")
title("Time Domain Plot of O Mic")

figure()
plot(tt3,P_Mic)
xlabel("Time(s)")
ylabel("Amplitude")
title("Time Domain Plot of P Mic")


%stft

%% 
% Performing  STFT on the data of the 3 mics with windowing and overlap

ST_N = stft(N_Mic,'Window',window','OverlapLength',overlap);
ST_O = spectrogram(O_Mic,window,overlap);
ST_P = stft(P_Mic,'Window',window','OverlapLength',overlap);

%Getting magnitude and bin location from stft data for each mic
ST_N_dB = 20*log10(abs(ST_N));
[nBins,nFr] = size(ST_N);
f_bins_loc = (0:nBins-1)*f_res2;
frame_center1 = (t_frame1/2) + (1:nFr);

ST_O_dB = 20*log10(abs(ST_O));
[nBins2,nFr2] = size(ST_O);
f_bins_loc2 = (0:nBins2-1)*f_res2;
frame_center2 = (t_frame2/2) + (1:nFr2);

ST_P_dB = 20*log10(abs(ST_P));
[nBins3,nFr3] = size(ST_P);
f_bins_loc3 = (0:nBins3-1)*f_res3;
frame_center3 = (t_frame3/2) + (1:nFr3);




%Spectrogram Plot of O(middle) mic

imagesc(frame_center2,f_bins_loc2,ST_O_dB)
xlabel('Time (s)')
ylabel('Freq (Hz)')
set(gca,'YDir','normal')
title('O Mic Spectrogram')






%% Task 1b & 2a
%Estimating the time at each mic and the breath rate

%Calculating the toltal energy for each mic
total_energy_O = zeros(1,nFr2);
for n = 1:nFr2
    start = (n*framelen) - (framelen-1);
    stop = n*framelen;
    frame_data = O_Mic(start:stop);
    tot_e = sum(frame_data.^2);

    total_energy_O(n) = tot_e; 
end

total_energy_N = zeros(1,nFr);
for n = 1:nFr
    start = (n*framelen) - (framelen-1);
    stop = n*framelen;
    frame_data = N_Mic(start:stop);
    tot_e = sum(frame_data.^2);

    total_energy_N(n) = tot_e; 
end

total_energy_P = zeros(1,nFr3);
for n = 1:nFr3
    start = (n*framelen) - (framelen-1);
    stop = n*framelen;
    frame_data = P_Mic(start:stop);
    tot_e = sum(frame_data.^2);

    total_energy_P(n) = tot_e; 
end


%total energy plot
figure
% plot(frame_center1,total_energy_N)
% % hold on
plot(frame_center2,total_energy_O)
% plot(frame_center3,total_energy_P)
% hold off
xlabel("Time(s)")
ylabel("Amplitude")
title("Envelope of O Mic")
legend("O-Mic")


% 
% center_times = (1:nFr)+(t_frame2/2);

%Find peaks in total energy
[pks_O,locs_O,w_O,p_O] = findpeaks(total_energy_O);
[max_p_O,ind_O] = max(p_O);
[max_pk_O,ind_pk_O] = max(pks_O);
peak_time_O = frame_center2(ind_O);
time_at_mic_O = (locs_O(ind_O)*framelen)/Fs2;


[pks_N,locs_N,w_N,p_N] = findpeaks(total_energy_N);
[max_p_N,ind_N] = max(p_N);
peak_time_N = frame_center1(ind_N);
time_at_mic_N = (locs_N(ind_N)*framelen)/Fs1;

[pks_P,locs_P,w_P,p_P] = findpeaks(total_energy_P);
[max_p_P,ind_P] = max(p_P);
peak_time_P = frame_center1(ind_P);
time_at_mic_P = (locs_P(ind_P)*framelen)/Fs1;


%fft of total energy
te_fft = fft(total_energy_O);
Fs_new = Fs2/framelen;
freqs = (1:length(total_energy_O))*(Fs_new/nFr);


%fft plot
figure
plot(freqs,20*log10(abs(te_fft)))
xlabel("Frequency(Hz)")
ylabel("Magnitude(dB")
title("FFT plot of O-Mic total energy")

%pks of total energy
[pks,locs,w,p] = findpeaks(20*log10(abs(te_fft)));
[pk,ind] = max(pks);
in = locs(ind);
breath_rate = freqs(in);   %breath rate

% avg_p = mean(p);
% 
% 
% pks(ind)
% 
% locs(ind)
% 
% w(ind)
% 
% 



% 



%% Task 2b and 2c
%Estimating Diver speed and height
diver_speed = (2*14)/(time_at_mic_P-time_at_mic_N);
divers_NO = 14/(time_at_mic_O-time_at_mic_N);
divers_OP = 14/(time_at_mic_P-time_at_mic_O);
avg_speed = (divers_NO+divers_OP)/2;




% hyp = total_energy_O(ind_O);


hyp = sound_speed*time_at_mic_O;
height = sqrt((hyp^2)-(14^2));

[c,lag] = xcorr(total_energy_O,total_energy_P);
[max_c,ind_c] = max(c);
delay = lag(ind_c);



% height2=((time_at_mic_O-time_at_mic_N)*sound_speed);
% 
% 
% dist = total_energy_O(ind_O) - total_energy_N(ind_O);
% 
% height3 =  sqrt((dist^2)-(14^2));
% ang = angle(height3);
% height4 = 14/(tan(ang));
% 
% 

dist = max_pk_O - total_energy_N(ind_pk_O);

dist = max_pk_O;
height = sqrt((dist^2)-(14^2));

time = time_at_mic_O-time_at_mic_N;
dist = sound_speed*time;



pks_N(ind_O)
pks_O(ind_O)


time_idx = floor(204.08*Fs2);
[c,lags] = xcorr(O_Mic(time_idx-Fs2:time_idx+Fs2),P_Mic(time_idx-Fs2:time_idx+Fs2));
[max_c,ind_c] = max(c);

delay = lags(ind_c);
time_delay22 = delay/Fs2;

