%% AudioProcessing
recObj = audiorecorder(44100, 8, 1);
 


disp('Start speaking.');
recordblocking(recObj, 1);
disp('End of Recording.');


data = getaudiodata(recObj, 'double');

% if data < data.length 
%    for (a[i] = 0; i*1323 < data; i++) {
%        testData = data(1323*i:1323*i+1323);
%        }
%    
% end
testArr = zeros(1323,1);
for i = 0:44100/1323
   if i < 33
    test = data(i*1323+1:i*1323+1323); 
    testArr = [testArr test];
   end
end
data_fft = fft(data);
subplot(2,1,2)
%plot(data_fft)
plot(abs(data_fft(:,1)));
audioInMono = mean(data,2);

t = (0:length(audioInMono)-1)/44100;
subplot(2,1,1)
%0p0lot(t,audioInMono)
ylabel('Amplitude')
fid=fopen('StressTestSheet.txt', 'w');
fmt = '%5d %5d %5d %5d\n';
fprintf(fid,fmt,data_fft);
fclose(fid);




%% Supervised learning

GMMTestValues = importdata('bimodal_example.csv')
histfit(GMMTestValues)
GMModel = fitgmdist(GMMTestValues, 20)
gmPDF = @(x)reshape(pdf(GMModel,x(:)),size(x));

fsurf(gmPDF,[-10, 10])
