%% Retrieval of training data
% 1 = scream 2 = help 3 = fall
h = dataRetrieval();
%% Train Gmm-Hmm model
AllAudioData = trainModels(h,1);
ModelTop5 = trainModels(h,2);
ModelTop4 = trainModels(h,3);
ModelTop3 = trainModels(h,4);
ModelTop2 = trainModels(h,5);
ModelButtom5 = trainModels(h,6);
ModelButtom4 = trainModels(h,7);
ModelButtom3 = trainModels(h,8);
ModelButtom2 = trainModels(h,9);
ModelButtom2 = trainModels(h,11);
UltraAudioData = trainModels(ultra,1);
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
recObj = audiorecorder(8000, 16, 1);
while true
    recordblocking(recObj, 2);
    liveData = getaudiodata(recObj, 'double');
    signalApmlitude = rms(liveData);
    if signalApmlitude>0.15
        disp("min reached")
        Live = {audioProcessing(liveData, 8000)};
        [loglikFromAllModels, bestModel] = evalModels(Live, AllAudioData);
        assignin("base", "bestModel", bestModel);
        disp(bestModel)
    end

end 
