%% AudioProcessing
fs=8000;
framesize = 30/1000*fs;

recObj = audiorecorder(fs, 8, 1);
 


disp('Start speaking.');
recordblocking(recObj, 1);
disp('End of Recording.');


data = getaudiodata(recObj, 'double');

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

t = (0:length(audioInMono)-1)/44100;
subplot(2,1,1)
%0p0lot(t,audioInMono)
ylabel('Amplitude')
fid=fopen('StressTestSheet.txt', 'w');
fmt = '%5d %5d %5d %5d\n';
fprintf(fid,fmt,data_fft);
fclose(fid);




%% Supervised learning

GMMTestValues = importdata('bimodal_example.csv')
histfit(GMMTestValues)
GMModel = fitgmdist(GMMTestValues, 20)
gmPDF = @(x)reshape(pdf(GMModel,x(:)),size(x));

fsurf(gmPDF,[-10, 10])

%% test
% Generate Data
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

Q = 3;      % state num
M = 2;      % mix num
p = 2;      % feature dim

% Train Gmm-Hmm model
[p_start, A, phi, loglik] = ChmmGmm(Data, Q, M);
% [p_start, A, phi, loglik] = ChmmGmm(Data, Q, M, 'p_start0', p_start0, 'A0', A0, 'phi0', phi0, 'cov_type', 'diag', 'cov_thresh', 1e-1)

% Calculate p(X) & vertibi decode
logp_xn_given_zn = Gmm_logp_xn_given_zn(Data{1}, phi);
[~,~, loglik] = LogForwardBackward(logp_xn_given_zn, p_start, A);
path = LogViterbiDecode(logp_xn_given_zn, p_start, A);

p_start
A
phi
loglik
path'


color = {'r', 'g', 'k'};
for q = 1:Q
    for m = 1:M
        error_ellipse(phi.Sigma(:,:,m,q), phi.mu(:,m,q), 'style', color{q}); hold on
    end
end

