function fileNames = fileRetrieve(path)
    % Source: https://se.mathworks.com/matlabcentral/answers/uploaded_files/80696/recurse_subfolders_R2016b.m
    % Use with R2016b and later.
    % It's done differently with R2016a and earlier, with genpath().
    % Initialization steps:
    clc;    % Clear the command window.
    workspace;  % Make sure the workspace panel is showing.
    format long g;
    format compact;

    % Define a folder.
    topLevelFolder = path;
    if topLevelFolder == 0
        return;
    end
    fprintf('The top level folder is "%s".\n', topLevelFolder);

    % Specify the file pattern.
    % Get ALL files using the pattern *.*
    % Note the special file pattern.  It has /**/ in it if you want to get files in subfolders of the top level folder.
    % filePattern = sprintf('%s/**/*.m; %s/**/*.xml', topLevelFolder, topLevelFolder);
    % Finds all files in the topfolder.
    filePattern = sprintf('%s/**/*.wav', topLevelFolder);
    allFileInfo = dir(filePattern);

    % Throw out any folders.  We want files only, not folders.
    isFolder = [allFileInfo.isdir]; % Logical list of what item is a folder or not.
    % Now set those folder entries to null, essentially deleting/removing them from the list.
    allFileInfo(isFolder) = [];

    % Process all files in those folders.
    totalNumberOfFiles = length(allFileInfo);
    fullNames = strings(length(allFileInfo),1);
    % Now we have a list of all files, matching the pattern, in the top level folder and its subfolders.
    if totalNumberOfFiles >= 1
        for k = 1 : totalNumberOfFiles
            % Go through all those files and assign their full names to string array.
            thisFolder = allFileInfo(k).folder;
            thisBaseFileName = allFileInfo(k).name;
            fullFileName = fullfile(thisFolder, thisBaseFileName);
            fullNames(k) = fullFileName;
    % 		fprintf('     Processing file %d of %d : "%s".\n', k, totalNumberOfFiles, fullFileName);
        end
    else
        fprintf('     Folder %s has no files in it.\n', thisFolder);
    end
    fileNames = fullNames;
    fprintf('\nDone looking in all folders!\n');
end