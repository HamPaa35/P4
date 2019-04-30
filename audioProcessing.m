function res = audioProcessing()
    % Retrieve files and process them into a single variable allAudioRes.
    combNumb = 6;
    colStart = combNumb +1;
    soundTypes = 3;
    allAudioRes = cell(combNumb, soundTypes);
    for k = 1:soundTypes
        start_path = fullfile("./audio");
        if ~exist(start_path, 'dir')
            start_path = matlabroot;
        end
        % Ask user to choose the folder.
        switch k
            case 1
                uiwait(msgbox('Pick the folder for screams in the following window.'));
            case 2
                uiwait(msgbox('Pick the help folder in the following window.'));
            case 3
                uiwait(msgbox('Pick the folder with falls in the following window.'));
            otherwise
                disp('Error')
        end
        topLevelFolder = uigetdir(start_path);
        if topLevelFolder == 0
            return;
        end
        fprintf('The top level folder is "%s".\n', topLevelFolder);
        % Retrieves the names of the audio files within the folder and its
        % subfolders
        filePaths = fileRetrieve(topLevelFolder);

        % Do AP on retrieved audio files
        for i = 1:length(filePaths)
            [Help,Fs] = audioread(filePaths(i));
            allAudioRes{1,k}{i} = audioProcessingNorm(Help, Fs);
%             allAudioRes{1,k}{i} = audioProcessingRaw(Help, Fs);
        end
        % Split AP results into different combinations
        for j = 2:combNumb
            allAudioRes{j, k} = allAudioRes{1,k};
            for i = 1:length(allAudioRes{1,k})
%                 if j > 1 % Check if this can be removed!
                    for g = 1:j-1
                        allAudioRes{j, k}{i}(:, colStart-g) = [];
                    end
%                 end
            end
        end
    end
    res = allAudioRes;
    disp("done")
end

% Function for running all AP methods on a signal and returning them in a
% matrix.
function dataInMatrix = audioProcessingRaw(soundToAnalyse, fs)
    a = 0; % minimum value after mapping
    b = 1; % maximum value after mapping
    % Converts Stereo signals into mono signals.
    if size(soundToAnalyse,2)>1
        soundToAnalyse= sum(soundToAnalyse, 2) / size(soundToAnalyse, 2);
    end
    % Harmonic to noise ratio
    hr = harmonicRatio(soundToAnalyse, fs);
    % Pitch
    fTemp = pitch(soundToAnalyse, fs);
    hrSize = size(hr);
    fTempSize = size(fTemp);
    sizeDifference = hrSize(1) - fTempSize(1);
    fzeros = zeros(sizeDifference, hrSize(2));
    f0 = [fTemp;fzeros];
    % Spectral centroid
    centroid = spectralCentroid(soundToAnalyse, fs);
    %Flux
    flux = spectralFlux(soundToAnalyse,fs);
    % Roll off Point
    rolloffPoint = spectralRolloffPoint(soundToAnalyse,fs);
    % Spectral flatness
    flatness = spectralFlatness(soundToAnalyse,fs);
    % Output
    dataInMatrix = [f0, hr,centroid,flux,rolloffPoint,flatness];
end
% AP with normalization
function dataInMatrix = audioProcessingNorm(soundToAnalyse, fs)
    % Converts Stereo signals into mono signals.
    if size(soundToAnalyse,2)>1
        soundToAnalyse= sum(soundToAnalyse, 2) / size(soundToAnalyse, 2);
    end
    % Harmonic to noise ratio
    hr = harmonicRatio(soundToAnalyse, fs);
    % Pitch
    fTemp = pitch(soundToAnalyse, fs);
    hrSize = size(hr);
    fTempSize = size(fTemp);
    sizeDifference = hrSize(1) - fTempSize(1);
    fzeros = zeros(sizeDifference, hrSize(2));
    f0 = [fTemp;fzeros];
    f0 = f0/(fs/2);
    % Spectral centroid
    centroid = spectralCentroid(soundToAnalyse, fs);
    centroid = centroid/(fs/2);
    %Flux
%     flux = spectralFlux(soundToAnalyse, fs);
%     flux = flux*10;
    [s,cf,t] = melSpectrogram(soundToAnalyse,fs);
    flux = spectralFlux(s,cf);
    % Roll off Point
    rolloffPoint = spectralRolloffPoint(soundToAnalyse,fs);
    rolloffPoint = rolloffPoint/(fs/2);
    % Spectral flatness
    flatness = spectralFlatness(soundToAnalyse,fs);
    % Output
    dataInMatrix = [f0, hr,centroid,flux,rolloffPoint,flatness];
end
