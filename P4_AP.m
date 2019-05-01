%% Retrieval of training data
h = dataRetrieval();
%% Train Gmm-Hmm model
GmmModels = trainModels(h);
%% Eval model vs. data
addpath(genpath("./audio"))
[HelpTest,Fs] = audioread('Hjælp_GMM_refinement.wav');
[skrigTest,Fs] = audioread('Skrig_GMM_refinement.wav');
recObj = audiorecorder(8000, 16, 1);

disp('Start speaking.');
recordblocking(recObj, 2);
disp('End of Recording.');

liveData = getaudiodata(recObj, 'double');
Live = {audioProcessingNorm(liveData, 8000)};

ht = {audioProcessingNorm(HelpTest, 8000)};

st = {audioProcessingNorm(skrigTest, 8000)};

[loglikFromAllModels, bestModel] = evalModels(Live, GmmModels);
