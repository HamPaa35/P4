%Should take the array of sound data as input
function modelMatrix = trainModels(audioData, rowToTrain)
    addpath(genpath("./matlab-hmm-master"))
    dataSize = size(audioData);
    GmmModelsSize = size(audioData(rowToTrain,:));
    GmmModels = cell(GmmModelsSize(2), 1);
    % Loop that trains models for each sound type in the audio array.
    parfor i = 1:dataSize(2)
        % GMM and HMM implementation calls
        Q = 30;      % state num, Settled on this through testing
        M = 1;      % mix num, Settled on this through testing
        [p_start, A, phi, ~] = ChmmGmm(audioData{rowToTrain,i}, Q, M);
        GmmModels(i, 1) = {Model(p_start, A, phi)};
        disp("current model done")
        disp(num2str(i))
    end
    disp("All models done")
    modelMatrix = GmmModels;
end