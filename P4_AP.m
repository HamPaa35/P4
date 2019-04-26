%% AudioProcessing
clc
clear

% retrieve names of individual files in folders
% Define a starting folder.
h = {1,3};
for k = 1:3
    start_path = fullfile(matlabroot, '\toolbox');
    if ~exist(start_path, 'dir')
        start_path = matlabroot;
    end
    % Ask user to confirm the folder, or change it.
    switch k
        case 1
            uiwait(msgbox('Pick the folder for screams in the following window.'));
        case 2
            uiwait(msgbox('Pick the help folder in the following window.'));
        case 3
            uiwait(msgbox('Pick the folder with falls in the following window.'));
        otherwise
            disp('Error')
    end
    topLevelFolder = uigetdir(start_path);
    if topLevelFolder == 0
        return;
    end
    fprintf('The top level folder is "%s".\n', topLevelFolder);
    % Retrieves the names of the audio files within the folder and its
    % subfolders
    names = fileRetrieve(topLevelFolder);

    % Do AP on the audio files
    h{1,k} = {1,length(names)};
    for i = 1:length(names)
        [Help,Fs] = audioread(names(i));
        h{1,k}{1,i} = audioProcessing(Help, Fs);
    end
end

disp("done")
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
