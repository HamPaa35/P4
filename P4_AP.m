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
data_fft = fft(data);
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


%% test
% Generate Data
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

test = {GMMTestValues(1:50), GMMTestValues(51:100)}

Q = 3;      % state num
M = 2;      % mix num
p = 2;      % feature dim

% Train Gmm-Hmm model
[p_start, A, phi, loglik] = ChmmGmm(test, Q, M);
[tp_start, tA, tphi, tloglik] = ChmmGmm(Data, Q, M);
% [p_start, A, phi, loglik] = ChmmGmm(Data, Q, M, 'p_start0', p_start0, 'A0', A0, 'phi0', phi0, 'cov_type', 'diag', 'cov_thresh', 1e-1)

% Calculate p(X) & vertibi decode
logp_xn_given_zn = Gmm_logp_xn_given_zn(Data{1}, phi);
[~,~, loglik] = LogForwardBackward(logp_xn_given_zn, p_start, A);
path = LogViterbiDecode(logp_xn_given_zn, p_start, A);

p_start;
A;
phi;
loglik;
path';

color = {'r', 'g', 'k'};
for q = 1:Q
    for m = 1:M
        error_ellipse(phi.Sigma(:,:,m,q), phi.mu(:,m,q), 'style', color{q}); hold on
    end
end
