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


%% test
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

%% test af e
testOfEStep(test);

function testOfEStep(Data)
Q = 1;
M = 2;
p = size(Data{1},2);
tmp = rand(1,Q);
p_start = tmp / sum(tmp);
converge = 1 + 1e-4;
A = bsxfun(@rdivide, tmp, sum(tmp,2));
Xall = cell2mat(Data');
[prior_, mu_, Sigma_] = Gmm(Xall, M*Q, 'diag');
tmp = reshape(prior_,M,Q);
phi.B = bsxfun(@rdivide, tmp, sum(tmp, 1));
phi.mu = reshape(mu_,p,M,Q);
phi.Sigma = reshape(Sigma_,p,p,M,Q);
pre_ll = -inf;
obj_num = length(Data);
    % E STEP
    for r = 1:obj_num
        logp_xn_given_zn = Gmm_logp_xn_given_zn(Data{r}, phi);
        [LogGamma{r}, LogKsi{r}, Loglik{r}] = LogForwardBackward(logp_xn_given_zn, p_start, A);
        logp_xn_given_vn = Get_logp_xn_given_vn(Data{r}, phi);
        LogIta{r} = CalculateLogIta(logp_xn_given_vn, p_start, A, phi);
    end
% loglik
    loglik = 0;
    for r = 1:obj_num
        loglik = loglik + Loglik{r};
    end
    if (loglik-pre_ll<log(converge)) return;
    else pre_ll = loglik; end
end

function logp_xn_given_vn = Get_logp_xn_given_vn(X, phi)
    [N,p] = size(X);
    [M,Q] = size(phi.B);
    logp_xn_given_vn = zeros(N,M,Q);
    for q = 1:Q
        for m = 1:M
            x_minus_mu = bsxfun(@minus, X, phi.mu(:,m,q)');
            logp_xn_given_vn(:,m,q) = -0.5*p*log(2*pi) - 0.5*log(det(phi.Sigma(:,:,m,q))) - 0.5 * sum(x_minus_mu * inv(phi.Sigma(:,:,m,q)) .* x_minus_mu, 2);
        end
    end
end

function logita = CalculateLogIta(logp_xn_given_vn, p_start, A, phi)
    [N,M,Q] = size(logp_xn_given_vn);
    
    % reserve space
    logc = zeros(N,1);
    logalpha = zeros(N,M,Q);
    logbeta = zeros(N,M,Q);
    logita = zeros(N,M,Q);
    
    Tmp = bsxfun( @plus, log(phi.B) + reshape(logp_xn_given_vn(1,:,:),M,Q), log(p_start) );
    logc(1) = log( sum( sum( exp( Tmp - max(Tmp(:)) ) ) ) ) + max(Tmp(:));
    logalpha(1,:,:) = -logc(1) + Tmp;
    logbeta(N,:,:) = 0;
 
    % calculate c, alpha
    for n = 2:N
        T4 = zeros(M,Q,M,Q);    % dim 1,2: vn-1; dim 3,4: vn
        for q = 1:Q
            for m = 1:M
                T4(:,:,m,q) = logp_xn_given_vn(n,m,q) + log(phi.B) + bsxfun( @plus, reshape(logalpha(n-1,:,:),M,Q), log(A(:,q)') );
            end
        end
        tmp = exp( T4 - max(T4(:)) );
        logc(n) = log( sum(tmp(:)) ) + max(T4(:));
        
        for q = 1:Q
            for m = 1:M
                T2 = bsxfun( @plus, reshape(logalpha(n-1,:,:),M,Q), log(A(:,q)') );
                if isinf(max(T2(:)))
                    logalpha(n,m,q) = -inf;
                else
                    logalpha(n,m,q) = -logc(n) + logp_xn_given_vn(n,m,q) + log(phi.B(m,q)) + log( sum( sum( exp( T2 - max(T2(:)) ) ) ) ) + max(T2(:));
                end
            end
        end
    end
    
    for n = N-1:-1:1
        for q = 1:Q
            for m = 1:M
                T2 = bsxfun( @plus, reshape(logbeta(n+1,:,:),M,Q) + reshape(logp_xn_given_vn(n+1,:,:),M,Q) + log(phi.B), log(A(q,:)) );
                logbeta(n,m,q) = -logc(n+1) + log( sum( sum( exp( T2 - max(T2(:) ) ) ) ) ) + max(T2(:));
            end
        end
    end
    
    % calculate ita
    logita = logalpha + logbeta;
end