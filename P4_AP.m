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

GMMTestValues = importdata('bimodal_example.csv')
histfit(GMMTestValues, 20)
GMModel = fitgmdist(GMMTestValues, 103)
gmPDF = @(x)reshape(pdf(GMModel,x(:)),size(x)); 

fsurf(gmPDF,[-10 10])


>>>>>>> 0a60cb6638804e6c91b60ad68ea6a26c4b926e0a
