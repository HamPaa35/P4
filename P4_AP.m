%% Retrieval of training data
% 1 = scream 2 = help 3 = fall
% Function call that loads the training audio
h = dataRetrieval();
%% Train Gmm-Hmm model
% The code that calls the functions that train the relevant models on the
% sounds. Not in a loop so that the variables are saved if something goes
% wrong, and so that some of the models can be retrained easily.

% AllAudioData = trainModels(h,1);
% ModelTop5 = trainModels(h,2);
% ModelTop4 = trainModels(h,3);
% ModelTop3 = trainModels(h,4);
% ModelTop2 = trainModels(h,5);
% ModelButtom5 = trainModels(h,6);
% ModelButtom4 = trainModels(h,7);
% ModelButtom3 = trainModels(h,8);
% ModelButtom2 = trainModels(h,9);
% ModelHr = trainModels(h,10);
% ModelPitch = trainModels(h,11);
% ModelRoll = trainModels(h,12);
% ModelFlux = trainModels(h,13);
% ModelCen = trainModels(h,14);
% ModelFlat = trainModels(h,15);
% UltraAudioData = trainModels(ultra,1);
% allModels = {AllAudioData, ModelTop5, ModelTop4, ModelTop3, ModelTop2, 
% ModelButtom5, ModelButtom4, ModelButtom3, ModelButtom2, ModelsHr, 
% ModelPitch, ModelsRoll, ModelsFlux, ModelsCen, ModelsFlat}
%% Eval sound
% The sounds for the evaluation is loaded and analysed.
e = dataRetrieval();
evalAudioOverallSize = size(e);
ModelSize = size(allModels);
z = zeros(3,110);
% outer loop that goes through the audio data combinations
for k = 1:evalAudioOverallSize(1)
    % Loops through the data for each sound type
    for j = 1:evalAudioOverallSize(2)
        soundTypeSize = size(e{1,j});
        % Loops through all the audio recording data for each sound type
        for i = 1:soundTypeSize(2)
            % Evaluates the sound against the correct GMM model
            [~, bestModel] = evalModels(e{k, j}(1, i), allModels{1, k});
            % Saves the type of sound the system thinks the input is
            z(j,i) = bestModel;
        end
    end
    rawData{1,k} = z;
end    

%% Sort data
rawDataSize = size(rawData);
bestModelRawDataSize = size(rawData{1,1});
precision = zeros(3,15);
false = zeros(3,15);
% Loop through the cells that contain the best model raw data matricies.
for k = 1:rawDataSize(2)
    sf = 0;
    hf = 0;
    ff = 0;
    % Loops through the matrix containing the raw data
    for j = 1:bestModelRawDataSize(1)
        res = 0;
        %Loops through the rows in the above mentioned matrix
        for i = 1:bestModelRawDataSize(2)
            % If the row number and the bestModel number match the system
            % waas correct
            if j == rawData{1, k}(j,i)
                res=res+1;
            % If not classify what the system though it was.
            elseif rawData{1, k}(j,i) == 1
                sf = sf + 1;
            elseif rawData{1, k}(j,i) == 2
                hf = hf + 1;
            elseif rawData{1, k}(j,i) == 3
                ff = ff + 1;
            end
        end
        precision(j, k) = res;
    end
    false(1, k) = sf;
    false(2, k) = hf;
    false(3, k) = ff;
end
disp('du ser godt ud') 
%% Eval model vs. data
clear a; %For the arduino+
% The recorder that records audio for analysis
recObj = audiorecorder(8000, 16, 1);

while true
    disp("Recording...");
    % Record for 2 seconds
    recordblocking(recObj, 2);
    % Converts the recorded audio to audio data
    liveData = getaudiodata(recObj, 'double');
    % Measures the mean amplitude of the signal
    signalApmlitude = rms(liveData);
    disp(signalApmlitude)
    % If the audio is over the threashold, run analysis
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
