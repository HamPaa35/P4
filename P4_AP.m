%% AudioProcessing
addpath(genpath("./audio"))
[Help1,Fs] = audioread('10_hjælp_Mathias.wav');
[Help2,Fs] = audioread('help-Glerup.wav');
[Help3,Fs] = audioread('Hjælp_Rasmus.wav');
[Help4,Fs] = audioread('10_helps_Female.wav');
[skrig1,Fs] = audioread('10_skrig_Mathias.wav');
[skrig2,Fs] = audioread('screech-Glerup.wav');
[skrig3,Fs] = audioread('Skrig_Rasmus.wav');
[skrig4,Fs] = audioread('10_screams_Female.wav');

% Not currently used
fs=8000;
framesize = 30/1000*fs;

% Do AP on help signals
h1 = audioProcessing(Help1, Fs);
h2 = audioProcessing(Help2, Fs);
h3 = audioProcessing(Help3, Fs);
h4 = audioProcessing(Help4, Fs);

AllH = {h1, h2, h3, h4};

% Do AP on scream signals
s1 = audioProcessing(skrig1, Fs);
s2 = audioProcessing(skrig2, Fs);
s3 = audioProcessing(skrig3, Fs);
s4 = audioProcessing(skrig4, Fs);

AllS = {s1, s2, s3, s4};

disp("done")
%% test zone
addpath(genpath("./audio"))
test = '10'
test2 = strcat(test, '_hjælp_Mathias.wav')
[Helptest, Fs] = audioread(test2);
sound(Helptest, 44200)
%% Supervised learning
GMMTestValues = importdata('bimodal_example.csv')
histfit(GMMTestValues)
GMModel = fitgmdist(GMMTestValues, 20)
gmPDF = @(x)reshape(pdf(GMModel,x(:)),size(x));

fsurf(gmPDF,[-10, 10])
%% Train Gmm-Hmm model
% Locates files
addpath(genpath("./matlab-hmm-master"))

% GMM and HMM implementation calls
Q = 6;      % state num
M = 3;      % mix num
[p_start, A, phi, loglik] = ChmmGmm(AllH, Q, M);
disp("hdone")
[tp_start, tA, tphi, tloglik] = ChmmGmm(AllS, Q, M);
disp("Trained")
model_1h = {p_start, A, phi};
model_2h = {tp_start, tA, tphi};
%% Eval model vs. data
[HelpTest,Fs] = audioread('Skrig_Rasmus.wav');
[skrigTest,Fs] = audioread('test_skrig.mp3');
recObj = audiorecorder(Fs, 8, 2);

disp('Start speaking.');
recordblocking(recObj, 2);
disp('End of Recording.');
% 
liveData = getaudiodata(recObj, 'double');
Live = {audioProcessing(liveData, Fs)}; 

ht = {audioProcessing(HelpTest, Fs)};

st = {audioProcessing(skrigTest, Fs)};
 
CalOfLoglik(ht, model_1h, model_2h)

function CalOfLoglik(evalData, model_l, model_2)
    obj_num = length(evalData);
    for r = 1:obj_num
        logp_xn_given_zn = Gmm_logp_xn_given_zn(evalData{r}, model_l{3});
        [~, ~, Loglik1] = LogForwardBackward(logp_xn_given_zn, model_l{1}, model_l{2});
        logp_xn_given_znt = Gmm_logp_xn_given_zn(evalData{r}, model_2{3});
        [~, ~, Loglik2] = LogForwardBackward(logp_xn_given_znt, model_2{1}, model_2{2});
    end
    disp(Loglik1)
    disp(Loglik2)
    if Loglik2 > Loglik1
        disp("2")
    end
    if Loglik1 > Loglik2
        disp("1")
    end
end

function CalOfLoglik2(evalData, model_1_phi, model_2_phi, model_1_p_start, model_2_p_start, model_1_A, model_2_A)
    obj_num = length(evalData);
    for r = 1:obj_num
        logp_xn_given_zn = Gmm_logp_xn_given_zn(evalData{r}, model_1_phi);
        [~, ~, Loglik1] = LogForwardBackward(logp_xn_given_zn, model_1_p_start, model_1_A);
        logp_xn_given_znt = Gmm_logp_xn_given_zn(evalData{r}, model_2_phi);
        [~, ~, Loglik2] = LogForwardBackward(logp_xn_given_znt, model_2_p_start, model_2_A);
    end
    if Loglik2 > Loglik1
        disp("2")
    end
    if Loglik1 > Loglik2
        disp("1")
    end
end

% Function for running all AP methods on a signal and returning them in a
% matrix.
function dataInMatrix = audioProcessing(soundToAnalyse, fs)
    % Converts Stereo signals into mono signals.
    if size(soundToAnalyse,2)>1
        soundToAnalyse= sum(soundToAnalyse, 2) / size(soundToAnalyse, 2);
    end
    hr = harmonicRatio(soundToAnalyse, fs);
    fTemp = pitch(soundToAnalyse, fs);
    hrSize = size(hr);
    fTempSize = size(fTemp);
    sizeDifference = hrSize(1) - fTempSize(1); 
    fzeros = zeros(sizeDifference, hrSize(2));
    f0 = [fTemp;fzeros];
    centroid = spectralCentroid(soundToAnalyse, fs);
    flux = spectralFlux(soundToAnalyse,fs);
    rolloffPoint = spectralRolloffPoint(soundToAnalyse,fs);
    flatness = spectralFlatness(soundToAnalyse,fs);
    dataInMatrix = [f0, hr,centroid,flux,rolloffPoint,flatness];
end