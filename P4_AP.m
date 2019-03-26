%% AudioProcessing
recObj = audiorecorder;

while(true)
disp('Start speaking.');
recordblocking(recObj, 10);
disp('End of Recording.');


data = getaudiodata(recObj);
data_fft = fft(data);
subplot(2,1,2)
plot(abs(data_fft(:,1)));
audioInMono = mean(data,2);

t = (0:length(audioInMono)-1)/44100;
subplot(2,1,1)
plot(t,audioInMono)
ylabel('Amplitude') 
end


%% Supervised learning
asljgbkeabgr
erg
ewr
