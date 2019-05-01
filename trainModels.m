%Should take the array of sound data as input
function modelMatrix = trainModels(audioData)
    addpath(genpath("./matlab-hmm-master"))
    dataSize = size(audioData);
    GmmModels = cell(dataSize(1), 1);
    for i = 1:dataSize(2)
        % GMM and HMM implementation calls
        Q = 6;      % state num %Antal stavelser
        M = 3;      % mix num %Gausians per stavelser %Ikke mega vigtigt, men prøv lidt
        [p_start, A, phi, ~] = ChmmGmm(audioData{i}, Q, M);
        GmmModels(i, 1) = {Model(p_start, A, phi)};
        disp("current model done")
        disp(num2str(i))
    end
    disp("All models done")
    modelMatrix = GmmModels;
end