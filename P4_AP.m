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
false = zeros(3,15);
z = zeros(3,110);
for k = 1:15
    sf = 0;
    hf = 0;
    ff = 0;
    for j = 1:OverallSize(2)
        smallerSize = size(e{modelRow,j});
        for i = 1:smallerSize(2)
            [~, bestModel] = evalModels(e{k, j}(1, i), allModels{1, k});
            z(j,i) = bestModel;
        end
    end
    rawData{1,k} = z;
end
        

%%
for k = 1:2
    sf = 0;
    hf = 0;
    ff = 0;
    for j = 1:OverallSize(2)
        smallerSize = size(e{modelRow,j});
        for i = 1:smallerSize(2)
            [~, bestModel] = evalModels(e{k, j}(1, i), allModels{1, k});
            z(j,i) = bestModel;
            OverallSize1 = size(z);
            for L = 1:OverallSize1(1)
                res = 0;
                for m = 1:OverallSize1(2)
                    if L == z(L, m)
                        res=res+1;
                    else
                        if z(L,m) == 1
                            sf = sf + 1;
                        elseif z(L,m) == 2
                            hf = hf + 1;
                        else
                            ff = ff + 1;
                        end
                    end 
                end
                precision(L, k) = res;
            end
            
        end
    end
    false(1, k) = sf;
    false(2, k) = hf;
    false(3, k) = ff;
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
clear a; %For the arduino
recObj = audiorecorder(8000, 16, 1);
while true
    disp("Recording...");
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
% Lights up LED depending on which type of emergency has been registered.
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
       disp("Done recording!");
end
