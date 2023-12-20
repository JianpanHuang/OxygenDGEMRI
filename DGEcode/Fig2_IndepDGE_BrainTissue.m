clear all; close all; clc;
addpath(genpath('toolbox')); warning off

%% Initialize parameters
dataPath = 'Data/IndependentNormoxia'; % Change to independent_hyperoxia to see the hyperoxic results
roiName = 'BrainTissue';
baseNum = 30;
tResol = 15 ; % Time resolution, in sec
discaNum = 2;
cutNum = 0;
tBase = baseNum*tResol/60;
ifDenoi = 1; % Decide whether to denoise, 1 for denoise
ifRoiFilt = 0; % Decide whether to fliter low SNR pixels in ROI, 1 with ROI filtering, 0 without filtering
filtThre = 0.2; % Set the ratio factor for filtering, value within 0~1
if ifRoiFilt == 0
    filtThre = 0;
end

%% Load image data
load([dataPath,filesep,'Result.mat']);
img = Result.image;
img(:,:,1:discaNum*2)=[];
imgParen = img(:,:,1:2:end);
[xn, yn, tn] = size(imgParen);

%% Denoising
if ifDenoi == 1
    % [~, ~, sv] = mlsvd(imgParen); % Singular value
    % svn{1} = sv{1}/max(sv{1}); % Normalized singular value
    % svn{2} = sv{2}/max(sv{2}); 
    % svn{3} = sv{3}/max(sv{3}); 
    % % Determine the truncation indexes using Malinowskis, Nelson and Median criteria.
    % [malInd(1,1), nelInd(1,1), medInd(1,1)] = truncDeterm(svn{1});
    % [malInd(1,2), nelInd(1,2), medInd(1,2)] = truncDeterm(svn{2});
    % [malInd(1,3), nelInd(1,3), medInd(1,3)] = truncDeterm(svn{3});
    % % Denoise
    % [u, dgeFun] = mlsvd(imgParen, [medInd(1), medInd(2), medInd(3)]);
    [u, dgeFun] = mlsvd(imgParen, round([xn/2, yn/2, tn/3*(yn/tn)]));
    imgParen = lmlragen(u, dgeFun);
end

%% Draw ROI
[roiMask, roiNum] = drawRoi(imgParen(:,:,1), dataPath, roiName, 'gray', 'm', filtThre);
load([dataPath, filesep, 'csf.mat'],'mask');
roiMask(mask==1) = 0; % Remove CSF

%% Calculate DGE maps
imgBase = mean(imgParen(:,:,1:baseNum), 3);
for mm = 1:tn
    imgDge(:,:,mm)=(imgBase - imgParen(:,:,mm))./imgBase*100; % detaS = (Sbase-S)/Sbase*100
end

%% Calculate DGE curves
dgeCur = zeros(tn,roiNum);
for mm = 1:tn
    imgTemp = imgDge(:,:,mm);
	for nr = 1:roiNum
        roiTemp = roiMask(:,:,nr);
        dgeCur(mm,nr) = mean2(imgTemp(roiTemp==1)); 
	end
end

%% Fit
t0 = tBase;
IV = [1       0.1];
LB = [1e-5    0  ];
UB = [20      0.5];
reluFun = @(t) max(t,0);
sFun = @(p,t)  p(1)-p(1)*exp(-reluFun(t-t0)*p(2));
tMin = (tResol*(1:tn))'/60; % Time in minutes
tCut = tMin(1:end-cutNum);
for nr = 1:roiNum
    dgeRaw(:,nr) = dgeCur(1:end-cutNum, nr);
    [pFit(:,nr), rn(1,nr)] = lsqcurvefit(sFun, IV, tCut([1:baseNum,baseNum+10*60/tResol+1:end]), dgeRaw([1:baseNum,baseNum+10*60/tResol+1:end],nr), LB, UB);
    dgeFit(:,nr) = sFun(pFit(:,nr), tCut);
    rsq(1,nr) = 1 - rn(1,nr)/sum(dgeRaw(:,nr).^2);
    varia(1,nr) = var(dgeRaw(:,nr) - dgeFit(:,nr));
    snrc(1,nr) = 10*log10(pFit(1,nr)^2./varia(1,nr));
end

%% Show DGE curve results
tScal = [(1-baseNum)*tResol/60, max(tCut)-tBase];
sScal = [-0.5, 10];
figure, set(gcf,'Position',[300 300 1300 500]);
subplot(1,2,1),imagesc(imgParen(:,:,1)); axis image; colormap(gray); hold on;
contour(sum(roiMask,3),1,'m-','LineWidth',2); title('ROI'); set(gca,'Fontsize',18);
subplot(1,2,2)
for nr = 1:roiNum
    plot(tCut-tBase,dgeRaw(:,nr),'o',tCut-tBase,dgeFit(:,nr),'-','LineWidth',2); hold on
end 
hold off
xlabel('Time (min)','FontName', 'Arial');
ylabel('Normalized Intensity (%)','FontName', 'Arial');
title(['Independent DGE cureve of ',roiName]), axis([tScal, sScal]);
set(gca,'Fontsize',18);
% Show the fit parameters on figure
txtStr = ['Smax = ',num2str(pFit(1,:)),', uin = ',num2str(pFit(2,:)),', SNRC = ',num2str(snrc)];
txtUi = uicontrol('Style', 'text', 'Fontname', 'Arial', 'FontSize', 16,...
                  'Position', [748 400 400 40], 'String', txtStr);
% Save DGE curve(s)
saveas(gcf,[dataPath, filesep, roiName,'_DGEcurve.jpg']);
xlswrite([dataPath, filesep, roiName,'_DGEcurve.xls'],[tCut-tBase, dgeRaw, tCut-tBase, dgeFit]);
% % Save fitting parameters
% savePath = [dataPath, filesep,'BrainTissue_Smax_uin_R2_SNRC_',roiName,'.txt'];
% saveMat = [pFit(1,:)', pFit(2,:)', rsq', snrc'];
% saveTxt(savePath,saveMat);

%% Show DGE map results
sampleNum = 10;
rowNum = 6;
xScal = 3:53;
yScal = 18:81; 
cScal = [-2, 6];
answ = questdlg('Do you want to obtain DGE maps?', ...
    'Decision on DGE maps', ...
    'No','Yes','No');
switch answ
    case 'Yes'
        mask = sum(roiMask,3);
        mask_sub = mask;
        mask_sub(mask==1) = 0;
        mask_sub(mask==0) = 1e8;
        imgDge = imgDge.*repmat(mask,[1,1,tn]);
        imgDge = imgDge - repmat(mask_sub,[1,1,tn]);
        for nn = 1:floor(tn/sampleNum)
            imgDgeSample(:,:,nn) = mean(imgDge(:,:,(nn-1)*sampleNum+1:nn*sampleNum),3);
        end
        lastImg = mean(imgDge(:,:,(nn-1)*sampleNum+1:end),3) - mask_sub;
        imgDgeSampleBrain = cat(3,imgDgeSample(xScal,yScal,:),lastImg(xScal,yScal));
        firstImg = zeros(length(xScal),length(yScal)) - mask_sub(xScal,yScal) - 1e8;
        imgDgeSampleBrain = cat(3,firstImg,imgDgeSampleBrain);
        imgDgeSampleBrain2D = imshow3dimage(imgDgeSampleBrain, rowNum);
        res_num = rowNum - mod(size(imgDgeSampleBrain,3),rowNum);
        imgDgeSampleBrain2D(end-size(firstImg,1)+1:end,end-res_num*size(firstImg,2)+1:end)=-1e8;
        h = figure('numbertitle','off','name','DGE maps','color','white');
        imshow(imgDgeSampleBrain2D,[],'InitialMagnification','fit'); 
        mycolormap(1); colorbar; caxis(cScal); colorbar('FontName','Arial','FontSize',18, 'LineWidth', 2);
        title(['Independent DGE maps of ',roiName]);
%         set(gcf,'Position',[100 300 600 800]);
        export_fig([dataPath,filesep,roiName,'_DGEmap_Every',num2str(sampleNum)], '-tif', '-r200')
        anatoImg = imgParen(:,:,1).*mask ;
        anatoImgSave = anatoImg(xScal,yScal);
        imwrite(anatoImgSave/max(anatoImgSave(:)),[dataPath,filesep,roiName,'_Anatomy.tiff'],'tiff');
    case 'No'
        disp('DGE processing is done!')
end