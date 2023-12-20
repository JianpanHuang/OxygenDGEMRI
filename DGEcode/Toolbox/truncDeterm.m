function [malInd, nelInd, medInd] = truncDeterm(sv)
% Adaptively determine the truncation thresold for MLSVD denoising based
% on Malinowskis, Nelson and Median criteria.
% INPUT:
%              svc: Singular values
%           nRatio: Noise ratio (referring to the last nRatio part of the 
%                   svc, which is used to fit the L-curve
%            delta: Difference threshold to determine the truncation index
% OUTPUT:
%           malInd: The truncation index suggested by Malinowskis criterion
%           nelInd: The truncation index suggested by Nelson criterion
%           medInd: The truncation index suggested by Median criterion
%
% REFERENCES:
%   Breitling J. et al., NMR Biomed 2019;32:e4133.
%
%   This function was written with reference to Breitling's code, more details
%   can be found here: https://github.com/jbreitling/CEST-AdaptiveDenoising
    
    d = sort(sv,'descend');
    dm = diag(d);
    %% Malinowskis criterion. Malinowski, E. R. et al., Anal Chem 1977;49:612-617.
    kInd = zeros(length(d)-1,1);  
    for mm = 1: length(d)-1
        kInd(mm) = (sum(diag(dm(mm:end, mm:end))./(length(d)-mm)))^0.5 / (length(d)-mm)^2;
    end
    [~, malInd] = min(kInd);
    
   
    %% Nelson criterion. Nelson, L.R. et al., J Educ Res Meas 2005;3:1-17.
    for nn = 1:length(d)-1
        xx = nn:length(d);
        dSub = d(nn:end);
        p = polyfit(xx', dSub, 1);
        dFit = p(1)*xx' + p(2);
        rsq(nn) = 1 - sum((dSub - dFit).^2)/sum((dSub - mean(dSub)).^2);
    end
    nelInd = find(rsq > 0.80, 1, 'first'); % Find index with R^2 > 0.80.
    
    
    %% Median criteria. Manjón, J.V. et al., Med Image Anal. 2015;22:35-47.
    beta = 1.29;
    dSub = d(sqrt(d) < 2*median(sqrt(d)));
    medInd = find(sqrt(d) >= beta*sqrt(median(dSub)),1, 'last');
end