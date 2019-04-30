%% Audioprocessing
h = audioProcessing();
%% Train Gmm-Hmm model
GmmModels = trainModels(h);
%% Eval model vs. data
[HelpTest,Fs] = audioread('Skrig_Rasmus.wav');
[skrigTest,Fs] = audioread('test_skrig.mp3');
recObj = audiorecorder(Fs, 8, 2);

disp('Start speaking.');
recordblocking(recObj, 2);
disp('End of Recording.');
 
liveData = getaudiodata(recObj, 'double');
Live = {audioProcessing(liveData, Fs)};

ht = {audioProcessing(HelpTest, Fs)};

st = {audioProcessing(skrigTest, Fs)};

[loglikFromAllModels, bestModel] = evalModels(Live, GmmModels);