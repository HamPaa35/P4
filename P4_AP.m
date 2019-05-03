%% Retrieval of training data
% 1 = scream 2 = help 3 = fall
e = dataRetrieval();
%% Train Gmm-Hmm model
% for i = 7:11
%     Models{i} = trainModels(h,i);
%     disp("done")
%     disp(i)
% end
% AllGmm = trainModels(h,1);
% Models8 = trainModels(h,8);
% Models6 = trainModels(h,6);
% Models7 = trainModels(h,7);
% Models9 = trainModels(h,9);
% Models10 = trainModels(h,10);
% Models11 = trainModels(h,11);
%% Eval test
modelRow = 11;
OverallSize = size(e);
for j = 1:OverallSize(2)
    smallerSize = size(e{modelRow,j});
    for i = 1:smallerSize(2)
        [~, bestModel] = evalModels(e{modelRow, j}(1, i), Models11);
        z11(j,i) = bestModel;
    end
end
disp('du ser godt ud')

%% Optælling
OverallSize = size(z6);
res = 0;
for j = 1:OverallSize(1)
    for i = 1:OverallSize(2)
        if j == z11(j, i)
            res=res++1;
        end
    end
end
res11 = res;
%% Eval model vs. data
addpath(genpath("./audio"))
[HelpTest,Fs] = audioread('Hjælp_GMM_refinement.wav');
[skrigTest,Fs] = audioread('Skrig_GMM_refinement.wav');
recObj = audiorecorder(8000, 16, 1);

disp('Start speaking.');
recordblocking(recObj, 2);
disp('End of Recording.');

liveData = getaudiodata(recObj, 'double');
Live = {audioProcessing(liveData, 8000)};

ht = {audioProcessing(HelpTest, 8000)};

st = {audioProcessing(skrigTest, 8000)};

[loglikFromAllModels, bestModel] = evalModels(Live, GmmModels);
