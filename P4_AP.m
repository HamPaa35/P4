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
% UltraAudioData = trainModels(ultra,1);
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
clear a;
recObj = audiorecorder(8000, 16, 1);
disp("qwe");
while true
    disp("efefe")
    recordblocking(recObj, 2);
    liveData = getaudiodata(recObj, 'double');
    signalApmlitude = rms(liveData);
    disp(signalApmlitude)
    if signalApmlitude<0.15
        disp("min reached")
        Live = {audioProcessing(liveData, 8000)};
        [loglikFromAllModels, bestModel] = evalModels(Live, AllAudioData);
        disp(bestModel)


    end
       a = arduino('COM4', 'Uno');
       if bestModel == 1
       writeDigitalPin(a, 'D12', 1);
       end
       if bestModel == 2
           writeDigitalPin(a, 'D11', 1);
       end
       if bestModel == 3
       writeDigitalPin(a, 'D10', 1);
       end
       pause(1);
       clear a;
end
