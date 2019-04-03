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
ylabel('Magnitude')
xlabel('Frequency Values')

t = (0:length(audioInMono)-1)/44100;
subplot(2,1,1)
plot(t,audioInMono)

ylabel('Amplitude') 
ylabel('Amplitude')
xlabel('Time')
end


%% Supervised learning

GMMTestValues = importdata('bimodal_example.csv')
histfit(GMMTestValues)
GMModel = fitgmdist(GMMTestValues, 20)
gmPDF = @(x)reshape(pdf(GMModel,x(:)),size(x));

fsurf(gmPDF,[-10, 10])
