function [templateMatches pkRatio_low corrsAll] = bz_temporalCompression(spkmat,template,avgTemplateTime,tcFactors)
% USAGE
%   [] = bz_temporalCompression(spkmat,template,tcFactors)
%
%
% INPUTS
%   
%       spkmat              -output of bz_SpktToSpkmat
%       template            -MxN matrix of mean FR template (M bins by N cells)
%       avgTemplateTime     -average duration (seconds) of template bins
%       tcFactors           -range of compression factors to examine
%       plotting            -true/false for plotting results (default: true)
%       status              -true/false for updating user w/ progress (default: true)
%
% OUTPUTS
%
%
%
% written by david tingley, 2018
plotting = true;
status = true;

nBins = size(template,1);
avgTimePerBin = avgTemplateTime./nBins;

% we have avgTimePerBin for the template(ratemap), and spkmat.dt for the spk
% times.  We can solve for the degree of 'warping' we must do such that the
% bin size for the template (seconds) is the same as the bin size for
% spkmat
warpingFactor = round(nBins ./ (avgTimePerBin ./ spkmat.dt)); % template should be thise long for avgTimePerBin==dt
warpingFactor = nBins;
% remove zero FR cells from template and spkmat
spkmat.data(:,sum(template)==0)=[];
template(:,sum(template)==0)=[];
for cell=1:size(template,2)
    % so lets warp the template to get bin sizes equivalent to
    % spkmat.dt, and Z-score all the things 
    template_norm(:,cell) = zscore(makeLength(template(:,cell),round(warpingFactor),'spline'));
    spkmat.dataZscored(:,cell) = zscore(spkmat.data(:,cell));
end
nCells = cell;

for factor = 1:length(tcFactors)
    %% compress spkmat data
    compressionFactor = round(size(template_norm,1)./tcFactors(factor));
    for cell=1:size(template,2)
       % lets warp the template by compression factor
        templateCompressed(:,cell) = makeLength(template_norm(:,cell),...
            compressionFactor);
    end
    
    %% get correlations for this compression factor
    for ts = 1:size(spkmat.dataZscored,1)-compressionFactor
        if sum(sum(spkmat.data(ts:ts+compressionFactor-1,:),2)) > 1 % exclude epochs w/ 0 spikes
           [corr2d(ts)] = corr(reshape(spkmat.dataZscored(ts:ts+compressionFactor-1,:),size(templateCompressed,1)*nCells,1),...
               templateCompressed(:),'rows','complete'); 
           if status
               if mod(ts,10000) == 0
                   disp(['done with ' num2str(ts) ' of ' num2str(size(spkmat.dataZscored,1))...
                       ' timestamps, for factor ' num2str(factor) ' of ' num2str(length(tcFactors))])
               end
           end
        else
            corr2d(ts) = NaN;
        end
    end
    
    % zscore and make original length
    templateMatches(factor,:) = [nanZscore(corr2d), zeros(1,compressionFactor+1)];
    corrsAll(factor,:) = [(corr2d), zeros(1,compressionFactor+1)];
    % Get ratio of peaks
%     [pks locs] = findpeaks(templateMatches(factor,:),'minpeakheight',2);
      [pks locs] = findpeaks(templateMatches(factor,:));
%     [maxtab mintab] = peakdet(templateMatches(factor,:),2);
    pkRatio_low(factor) = sum(pks>4)./length(pks);

%     pkRatio_hi(factor) = length(pks); sum(maxtab(:,2)>4)./length(maxtab);
    
    %% plot TCF and other stuff
    if plotting
        subplot(2,2,1)
        plot(tcFactors(1:factor),pkRatio_low(1:factor))
        title('Peak ratio (>4/total)')
        subplot(2,2,2)
        plot(templateMatches(factor,:))
%         plot(tcFactors(1:factor),pkRatio_hi(1:factor))
%         title('Peak ratio (>7/total)')
        subplot(2,2,3)
        format long g
        imagesc(spkmat.timestamps(1:1000),tcFactors(1:factor),templateMatches(:,1:1000))    
        title('example template matches')
        colorbar
        subplot(2,2,4)
%         plot(templateMatches(factor,:))
        histogram(corr2d)
        pause(.1)
    end
    disp(['done with CF: ' num2str(tcFactors(factor))])
    clear compressedData corr2d templateCompressed corr2d
end

clear template_norm 

