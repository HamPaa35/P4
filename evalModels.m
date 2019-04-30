function [loglikMatrix, bestModelMatch] = evalModels(evalData, modelMatrix)
    obj_num = length(evalData);
    addpath(genpath("./matlab-hmm-master"))
    dataSize = size(modelMatrix);
    loglikOnModels = cell(dataSize(1), 2);
    for i =1:dataSize(1)
        for r = 1:obj_num
            logp_xn_given_zn = Gmm_logp_xn_given_zn(evalData{r}, modelMatrix{i}.phi);
            [~, ~, Loglik] = LogForwardBackward(logp_xn_given_zn, modelMatrix{i}.p_start, modelMatrix{i}.A);
        end
        loglikOnModels(i, 1) = {Loglik};
        loglikOnModels(i, 2) = {i};
        disp("eval current model done")
        disp(num2str(i))
    end
    disp("All models done")
    sortedLoglik = sortrows(loglikOnModels, [1]);
    sizeForOut = size(sortedLoglik);
    bestModelMatch = sortedLoglik{sizeForOut(1), sizeForOut(2)};
    loglikMatrix = sortedLoglik;

end