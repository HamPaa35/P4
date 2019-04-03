clear
clc
recObj = audiorecorder;
name = "qwe";
commandID = 1;
ledPin = 'D7';
deltaT_blink = 0.5;
a = arduino('COM4','Uno');

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
%<<<<<<< HEAD
    for i=1:10
        a.writeDigitalPin(ledPin, 1);
        pause(deltaT_blink/2);
        a.writeDigitalPin(ledPin,0);
        pause(deltaT_blink/2);
    end
end
