% Function for running all AP methods on a signal and returning them in a
% matrix.
function dataInMatrix = audioProcessing(soundToAnalyse, fs)
    % Converts Stereo signals into mono signals.
    if size(soundToAnalyse,2)>1
        soundToAnalyse= sum(soundToAnalyse, 2) / size(soundToAnalyse, 2);
    end
    % Harmonic to noise ratio
    hr = harmonicRatio(soundToAnalyse, fs);
    % Pitch using Normalized Correlation Function
    fTemp = pitch(soundToAnalyse, fs, 'Range',[50, fs/2]);
    hrSize = size(hr);
    fTempSize = size(fTemp);
    sizeDifference = hrSize(1) - fTempSize(1);
    fzeros = zeros(sizeDifference, hrSize(2));
    f0 = [fTemp;fzeros];
    % Spectral centroid
    centroid = spectralCentroid(soundToAnalyse, fs);
    %Flux
    flux = spectralFlux(soundToAnalyse,fs);
    % Roll off Point (95%)
    rolloffPoint = spectralRolloffPoint(soundToAnalyse,fs);
    % Spectral flatness
    flatness = spectralFlatness(soundToAnalyse,fs);
    % Output
    dataInMatrix = [hr, f0, rolloffPoint, flux, centroid, flatness];
end

% ---- AP with normalization
% function dataInMatrix = audioProcessing(soundToAnalyse, fs)
%     % Converts Stereo signals into mono signals.
%     if size(soundToAnalyse,2)>1
%         soundToAnalyse= sum(soundToAnalyse, 2) / size(soundToAnalyse, 2);
%     end
%     % Harmonic to noise ratio
%     hr = harmonicRatio(soundToAnalyse, fs);
%     % Pitch
%     fTemp = pitch(soundToAnalyse, fs);
%     hrSize = size(hr);
%     fTempSize = size(fTemp);
%     sizeDifference = hrSize(1) - fTempSize(1);
%     fzeros = zeros(sizeDifference, hrSize(2));
%     f0 = [fTemp;fzeros];
%     f0 = f0/(fs/2);
%     % Spectral centroid
%     centroid = spectralCentroid(soundToAnalyse, fs);
%     centroid = centroid/(fs/2);
%     %Flux
% %     flux = spectralFlux(soundToAnalyse, fs);
% %     flux = flux*10;
%     [s,cf,t] = melSpectrogram(soundToAnalyse,fs);
%     flux = spectralFlux(s,cf);
%     % Roll off Point
%     rolloffPoint = spectralRolloffPoint(soundToAnalyse,fs);
%     rolloffPoint = rolloffPoint/(fs/2);
%     % Spectral flatness
%     flatness = spectralFlatness(soundToAnalyse,fs);
%     % Output
%     dataInMatrix = [f0, hr,centroid,flux,rolloffPoint,flatness];
% end
