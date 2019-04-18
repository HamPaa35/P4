%% AudioProcessing
addpath(genpath("./audio"))
[Help1,Fs] = audioread('10_hjælp_Mathias.wav');
[Help2,Fs] = audioread('help-Glerup.wav');
[Help3,Fs] = audioread('Hjælp_Rasmus.wav');
[skrig1,Fs] = audioread('10_skrig_Mathias.wav');
[skrig2,Fs] = audioread('screech-Glerup.wav');
[skrig3,Fs] = audioread('Skrig_Rasmus.wav');

fs=8000;
framesize = 30/1000*fs;

recObj = audiorecorder(fs, 8, 1);

%%%%%%%%%%%mnjkn

% disp('Start speaking.');
% recordblocking(recObj, 1);
% disp('End of Recording.');

% data = getaudiodata(Help1, 'double');
[h1hr,h1centroid,h1flux,h1rolloffPoint,h1flatness] = audioProcessing(Help1, Fs);
h1 = [h1hr,h1centroid,h1flux,h1rolloffPoint,h1flatness];
[h2hr,h2centroid,h2flux,h2rolloffPoint,h2flatness] = audioProcessing(Help2, Fs);
h2 = [h2hr,h2centroid,h2flux,h2rolloffPoint,h2flatness];
[h3hr,h3centroid,h3flux,h3rolloffPoint,h3flatness] = audioProcessing(Help3, Fs);
h3 = [h3hr,h3centroid,h3flux,h3rolloffPoint,h3flatness];

AllH = {h1, h2, h3};

[s1hr,s1centroid,s1flux,s1rolloffPoint,s1flatness] = audioProcessing(skrig1, Fs);
s1 = [s1hr,s1centroid,s1flux,s1rolloffPoint,s1flatness];
[s2hr,s2centroid,s2flux,s2rolloffPoint,s2flatness] = audioProcessing(skrig2, Fs);
s2 = [s2hr,s2centroid,s2flux,s2rolloffPoint,s2flatness];
[s3hr,s3centroid,s3flux,s3rolloffPoint,s3flatness] = audioProcessing(skrig3, Fs);
s3 = [s3hr,s3centroid,s3flux,s3rolloffPoint,s3flatness];

AllS = {s1, s2, s3};

% s1 = audioProcessing(skrig1, Fs);
% s2 = audioProcessing(skrig2, Fs);
% s3 = audioProcessing(skrig3, Fs);

C2 = [hr,centroid,flux,rolloffPoint,flatness];

C = {f0,hr,centroid,flux,rolloffPoint,flatness,data_fft};

% if data < data.length
%    for (a[i] = 0; i*1323 < data; i++) {
%        testData = data(1323*i:1323*i+1323);
%        }
%
% end



testArr = zeros(framesize,1);
% for i = 0:fs/framesize
%    if i < 33
%     test = data(i*framesize+1:i*framesize+framesize);
%     testArr = [testArr test];
%    end
% end

% subplot(2,1,2)
% %plot(data_fft)
% plot(abs(data_fft(:,1)));
% % audioInMono = mean(data,2);
% 
% t = (0:length(audioInMono)-1)/fs;
% subplot(2,1,1)
% %0p0lot(t,audioInMono)
% ylabel('Amplitude')
% fclose(fid);
disp("done")
%% test zone
[h1f0,h1hr,h1centroid,h1flux,h1rolloffPoint,h1flatness] = audioProcessing(Help1, Fs);

%% Supervised learning
GMMTestValues = importdata('bimodal_example.csv')
histfit(GMMTestValues)
GMModel = fitgmdist(GMMTestValues, 20)
gmPDF = @(x)reshape(pdf(GMModel,x(:)),size(x));

fsurf(gmPDF,[-10, 10])


%% Train Gmm-Hmm model
% Generate Data
GMMTestValues = importdata('bimodal_example.csv');
addpath(genpath("./matlab-hmm-master"))
% for i1 = 1:2
%     X1 = mvnrnd([0,0], [0.5, 0.2; 0.2, 0.3]/5, 20);
%     X2 = mvnrnd([0,2], [0.3, -0.2; -0.2, 0.5]/5, 30);
%     X3 = mvnrnd([0,4], [0.5, 0; 0, 0.3]/5, 40);
%     X = [X1; X2; X3];
%     Data{i1} = X;
% end
% for i1 = 3:4
%     X1 = mvnrnd([2,0], [0.5, 0.2; 0.2, 0.3]/5, 20);
%     X2 = mvnrnd([2,2], [0.3, -0.2; -0.2, 0.5]/5, 30);
%     X3 = mvnrnd([2,4], [0.5, 0; 0, 0.3]/5, 30);
%     X = [X1; X2; X3];
%     Data{i1} = X;
% end
% Xall = cell2mat(Data');
% 
test = {GMMTestValues(1:50), GMMTestValues(51:100)};
% delAfTest = {10, 11};

Q = 6;      % state num
M = 3;      % mix num
[p_start, A, phi, loglik] = ChmmGmm(AllH, Q, M);
disp("hdone")
[tp_start, tA, tphi, tloglik] = ChmmGmm(AllS, Q, M);
disp("Trained")
model_1h = {p_start, A, phi};
model_2h = {tp_start, tA, tphi};
%% Eval model vs. data
[HelpTest,Fs] = audioread('test.wav');
[skrigTest,Fs] = audioread('test_skrig.mp3');
recObj = audiorecorder(Fs, 8, 2);

disp('Start speaking.');
recordblocking(recObj, 2);
disp('End of Recording.');

data = getaudiodata(recObj, 'double');

[lhr,lcentroid,lflux,lrolloffPoint,lflatness] = audioProcessing(data, Fs);
l = [lhr,lcentroid,lflux,lrolloffPoint,lflatness];
Live = {l};

[ht1hr,ht1centroid,ht1flux,ht1rolloffPoint,ht1flatness] = audioProcessing(HelpTest, Fs);
ht1 = [ht1hr,ht1centroid,ht1flux,ht1rolloffPoint,ht1flatness];
ht = {ht1};

[st1hr,st1centroid,st1flux,st1rolloffPoint,st1flatness] = audioProcessing(skrigTest, Fs);
st1 = [st1hr,st1centroid,st1flux,st1rolloffPoint,st1flatness];
st = {st1};
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

function [f0, hr,centroid,flux,rolloffPoint,flatness] = audioProcessing(soundToAnalyse, fs)
    hr = harmonicRatio(soundToAnalyse, fs);
    f0 = pitch(soundToAnalyse, fs);
    centroid = spectralCentroid(soundToAnalyse, fs);
    flux = spectralFlux(soundToAnalyse,fs);
    rolloffPoint = spectralRolloffPoint(soundToAnalyse,fs);
    flatness = spectralFlatness(soundToAnalyse,fs);
%     data_fft = fft(soundToAnalyse);
end