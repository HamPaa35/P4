% Should take signal to be evaluated and models to evaluate against as
% input
function [loglikMatrix, bestModelMatch] = evalModels(evalData, modelMatrix)
    obj_num = length(evalData);
    addpath(genpath("./matlab-hmm-master"))
    dataSize = size(modelMatrix);
    loglikOnModels = cell(dataSize(1), 2);
    % Loops though all models
    for i =1:dataSize(1)
        % Loops through all the data to be evaluated
        for r = 1:obj_num
            logp_xn_given_zn = Gmm_logp_xn_given_zn(evalData{r}, modelMatrix{i}.phi);
            [~, ~, Loglik] = LogForwardBackward(logp_xn_given_zn, modelMatrix{i}.p_start, modelMatrix{i}.A);
        end
        loglikOnModels(i, 1) = {Loglik};
        loglikOnModels(i, 2) = {i};
    end
    sortedLoglik = sortrows(loglikOnModels, [1]);
    sizeForOut = size(sortedLoglik);
    % returnes the model with best likelyhood
    bestModelMatch = sortedLoglik{sizeForOut(1), sizeForOut(2)};
    % returns the ranking of the data
    loglikMatrix = sortedLoglik;
end