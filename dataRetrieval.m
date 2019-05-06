function res = dataRetrieval()
    % Retrieve files and process them into a single variable allAudioRes.
    combNumb = 6;
    soundTypes = 3;
    allAudioRes = cell((combNumb-2)*2+1+combNumb, soundTypes);
    titleString = "";
    for k = 1:soundTypes
        start_path = fullfile("./audio");
        if ~exist(start_path, 'dir')
            start_path = matlabroot;
        end
        % Ask user to choose the folder.
        switch k
            case 1
                uiwait(msgbox('Pick the folder for screams in the following window.'));
                titleString = "Select screams";
            case 2
                uiwait(msgbox('Pick the help folder in the following window.'));
                titleString = "Select calls for help";
            case 3
                uiwait(msgbox('Pick the folder with falls in the following window.'));
                titleString = "Select falls";
            otherwise
                disp('Error')
        end
        topLevelFolder = uigetdir(start_path, titleString);
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
            allAudioRes{1,k}{i} = audioProcessing(Help, Fs);
        end
        % Split AP results into different combinations
        % loops rows
        for j = 2:combNumb-1
            row = j+combNumb-2;
            allAudioRes{j, k} = allAudioRes{1,k};
            allAudioRes{row, k} = allAudioRes{1,k};
            % loops columns
            for i = 1:length(allAudioRes{1,k})
                % Loops columns for matrix
                for g = 0:j-2
                    allAudioRes{j, k}{i}(:, combNumb-g) = [];
                    allAudioRes{row, k}{i} = allAudioRes{row, k}{i}(:,2:end);
                end
            end
        end
        % adds the single values
        for j = 1:combNumb
            rowSingles = j+((combNumb-2)*2+1);
            for i = 1:length(allAudioRes{1,k})
                allAudioRes{rowSingles, k}{i} = allAudioRes{1, k}{i}(:,j);
            end
        end
        
    end
    res = allAudioRes;
    disp("done")
end