clear
clc
recObj = audiorecorder;
name = "qwe";
commandID = 1;
ledPin = 'D7';
deltaT_blink = 0.5;
port = 'COM4';
board = 'UNO';
%a = arduino(port, board);
%for i=1:100
%    a.writeDigitalPin(ledPin, 1);
%    pause(deltaT_blink/2);
%    a.writeDigitalPin(ledPin,0);
%    pause(deltaT_blink/2);
%end

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
sendCommand(recObj,name,commandID,57)
end


%Supervised learning
%asljgbkeabgr
%erg
%ewr
%>>>>>>> 5fdb2cb7f66c4b5411ab0e64743041daccaa4bed
