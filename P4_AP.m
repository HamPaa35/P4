%% AudioProcessing

% retrieve names of individual files in folders
% Define a starting folder.
h = {};
for k = 1:3
    start_path = fullfile("./audio");
    if ~exist(start_path, 'dir')
        start_path = matlabroot;
    end
    % Ask user to choose the folder.
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
    filePaths = fileRetrieve(topLevelFolder);

    % Do AP on retrieved audio files
    for i = 1:length(filePaths)
        [Help,Fs] = audioread(filePaths(i));
        h{k}{i} = audioProcessing(Help, Fs);
    end
end

disp("done")
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

CalOfLoglik(Live, GmmModels{1}, GmmModels{2})

%% Functions
function CalOfLoglik(evalData, model_l, model_2)
    obj_num = length(evalData);
    for r = 1:obj_num
        logp_xn_given_zn = Gmm_logp_xn_given_zn(evalData{r}, model_l.phi);
        [~, ~, Loglik1] = LogForwardBackward(logp_xn_given_zn, model_l.p_start, model_l.A);
        logp_xn_given_znt = Gmm_logp_xn_given_zn(evalData{r}, model_2.phi);
        [~, ~, Loglik2] = LogForwardBackward(logp_xn_given_znt, model_2.p_start, model_2.A);
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

function CalOfLoglik2(evalData, model_l, model_2)
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
%Should take the array of sound data as input
function modelMatrix = trainModels(audioData)
    addpath(genpath("./matlab-hmm-master"))
    dataSize = size(audioData);
    GmmModels = cell(dataSize(1), 1);
    for i =1:dataSize(2)
        % GMM and HMM implementation calls
        Q = 6;      % state num %Antal stavelser
        M = 3;      % mix num %Gausians per stavelser %Ikke mega vigtigt, men prøv lidt
        [p_start, A, phi, ~] = ChmmGmm(audioData{i}, Q, M);
        GmmModels(i, 1) = {Model(p_start, A, phi)};
        disp("current model done")
        disp(num2str(i))
    end
    disp("All models done")
    modelMatrix = GmmModels;
end

function [loglikMatrix, bestModelMatch] = evalModels(evalData, modelMatrix)
    obj_num = length(evalData);
    addpath(genpath("./matlab-hmm-master"))
    dataSize = size(modelMatrix);
    loglikOnModels = cell(dataSize(1), 2);
    for i =1:dataSize(1)
        for r = 1:obj_num
            logp_xn_given_zn = Gmm_logp_xn_given_zn(evalData{r}, modelMatrix{i}.phi);
            [~, ~, Loglik] = LogForwardBackward(logp_xn_given_zn, modelMatrix{i}.p_start, modelMatrix{i}.A);
        end
        loglikOnModels(i, 1) = {Loglik};
        loglikOnModels(i, 2) = {i};
        disp("eval current model done")
        disp(num2str(i))
    end
    disp("All models done")
    sortedLoglik = sortrows(loglikOnModels, [1]);
    sizeForOut = size(sortedLoglik);
    bestModelMatch = sortedLoglik{sizeForOut(1), sizeForOut(2)};
    loglikMatrix = sortedLoglik;

end