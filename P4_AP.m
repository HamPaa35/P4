%% Retrieval of training data
% 1 = scream 2 = help 3 = fall
h = dataRetrieval();
%% Train Gmm-Hmm model
% AllAudioData = trainModels(h,1);
% ModelTop5 = trainModels(h,2);
% ModelTop4 = trainModels(h,3);
% ModelTop3 = trainModels(h,4);
% ModelTop2 = trainModels(h,5);
% ModelButtom5 = trainModels(h,6);
% ModelButtom4 = trainModels(h,7);
% ModelButtom3 = trainModels(h,8);
% ModelButtom2 = trainModels(h,9);
% ModelPitch = trainModels(h,11);
% UltraAudioData = trainModels(ultra,1);
%% Eval test
%allModels = {AllAudioData, ModelTop5, ModelTop4, ModelTop3, ModelTop2, ModelButtom5, ModelButtom4, ModelButtom3, ModelButtom2, ModelsHr, ModelPitch, ModelsRoll, ModelsFlux, ModelsCen, ModelsFlat}
% e = dataRetrieval();
OverallSize = size(e);
ModelSize = size(allModels);
modelRow = ModelSize(2);
precision = zeros(3,15);
z = zeros(3,110);
for k = 1:ModelSize(2)
    for j = 1:OverallSize(1, 2)
        smallerSize = size(e{modelRow,j});
        for i = 1:smallerSize(2)
            [~, bestModel] = evalModels(e{k, j}(1, i), allModels{1, k});
            z(j,i) = bestModel;
            OverallSize1 = size(z);
            for L = 1:OverallSize1(1)
                res = 0;
                for m = 1:OverallSize1(2)
                    if L == z(L, m)
                        res=res++1;
                    end
                    precision(L, k) = res;
                end
            end
        end
    end
end
disp('du ser godt ud')

%% Opt�lling
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
