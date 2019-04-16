%% AudioProcessing
fs=8000;
framesize = 30/1000*fs;

recObj = audiorecorder(fs, 8, 1);



disp('Start speaking.');
recordblocking(recObj, 1);
disp('End of Recording.');

data = getaudiodata(recObj, 'double');

f0 = pitch(data, fs);
hr = harmonicRatio(data, fs);
centroid = spectralCentroid(data, fs);
flux = spectralFlux(data,fs);
rolloffPoint = spectralRolloffPoint(data,fs);
flatness = spectralFlatness(data,fs);
data_fft = fft(data);

C = {f0,hr,centroid,flux,rolloffPoint,flatness,data_fft};

% if data < data.length
%    for (a[i] = 0; i*1323 < data; i++) {
%        testData = data(1323*i:1323*i+1323);
%        }
%
% end



testArr = zeros(framesize,1);
for i = 0:fs/framesize
   if i < 33
    test = data(i*framesize+1:i*framesize+framesize);
    testArr = [testArr test];
   end
end

subplot(2,1,2)
%plot(data_fft)
plot(abs(data_fft(:,1)));
audioInMono = mean(data,2);

t = (0:length(audioInMono)-1)/fs;
subplot(2,1,1)
%0p0lot(t,audioInMono)
ylabel('Amplitude')
fclose(fid);




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
for i1 = 1:2
    X1 = mvnrnd([0,0], [0.5, 0.2; 0.2, 0.3]/5, 20);
    X2 = mvnrnd([0,2], [0.3, -0.2; -0.2, 0.5]/5, 30);
    X3 = mvnrnd([0,4], [0.5, 0; 0, 0.3]/5, 40);
    X = [X1; X2; X3];
    Data{i1} = X;
end
for i1 = 3:4
    X1 = mvnrnd([2,0], [0.5, 0.2; 0.2, 0.3]/5, 20);
    X2 = mvnrnd([2,2], [0.3, -0.2; -0.2, 0.5]/5, 30);
    X3 = mvnrnd([2,4], [0.5, 0; 0, 0.3]/5, 30);
    X = [X1; X2; X3];
    Data{i1} = X;
end
Xall = cell2mat(Data');
scatter(Xall(:,1), Xall(:,2), '.'); hold on

test = {GMMTestValues(1:50), GMMTestValues(51:100)};
delAfTest = {10, 11};

Q = 2;      % state num
M = 2;      % mix num
p = 2;      % feature dim
[p_start, A, phi, loglik] = ChmmGmm(test, Q, M);
[tp_start, tA, tphi, tloglik] = ChmmGmm(Data, Q, M);
%% Eval model vs. data
CalOfLoglik(Data, phi, tphi, p_start, tp_start, A, tA)

function CalOfLoglik(evalData, model_1_phi, model_2_phi, model_1_p_start, model_2_p_start, model_1_A, model_2_A)
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