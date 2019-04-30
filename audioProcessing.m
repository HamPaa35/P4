function res = audioProcessing()
    % Retrieve files and process them into a single variable allAudioRes.
    combNumb = 6;
    colStart = combNumb +1;
    soundTypes = 3;
    allAudioRes = cell(combNumb*2-1, soundTypes);
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
            allAudioRes{1,k}{i} = aPMethod(Help, Fs);
        end
        % Split AP results into different combinations
        for j = 2:combNumb
            allAudioRes{j, k} = allAudioRes{1,k};
            row = j+5;
            for i = 1:length(allAudioRes{1,k})
                allAudioRes{row, k}{i} = allAudioRes{1, k}{i}(:,j);
                    for g = 1:j-1
                        allAudioRes{j, k}{i}(:, colStart-g) = [];
                    end
            end
        end
        
    end
    res = allAudioRes;
    disp("done")
end