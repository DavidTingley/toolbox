function [wCorr] = weightedCorr(mat)



    wCorr = sum(sum(weightedCov(mat))./ sqrt(weightedVar(mat,1) * weightedVar(mat,2)));

% weighted mean
function [wM] = wMean(mat)
    for i = 1:size(mat,1)
        for j=1:size(mat,2)
           wM(i,j) =  mat(i,j) .* i ./ sum(mat(:));
        end
    end
    wM = sum(wM(:));
return

%weighted covariance
function [wCov] = weightedCov(mat)
    wM_x = wMean(mat);
    wM_y = wMean(mat');
    
    for i=1:size(mat,1)
        for j=1:size(mat,2)scr
            wCov(i,j) = mat(i,j) .* ((i - wM_x) .* (j - wM_y)) ./ sum(mat(:));
        end
    end
    wCov = sum(wCov(:));
return

function [wVar] = weightedVar(mat,axis)

    if axis == 1
        wM_x = wMean(mat);
    else
        wM_x = wMean(mat');
    end

    for i=1:size(mat,axis)
        for j = 1:size(mat,axis)
            wVar(i,j) = mat(i,j) .* (i - wM_x) .* (j - wM_x) ./ sum(mat(:));
        end
    end
    wVar = sum(wVar(:));
return