clear all;clc; close all;
addpath(genpath('toolbox'))

%% Parameters
dataPath = 'Data/Interleaved';
roiName = 'BrainTissue';
tResol = 15;   % in s time resolution 
discaNum = 2;
baseO2Num = 30;
baseAirNum = 32;
cutNum = 0;
baseO2Time = 7.5; % minx
baseTransTime = 14; % min
postO2Time = 30;
baseTime = baseO2Time+baseTransTime+baseAirNum*tResol/60;
ifRoiFilt = 0; % Decide if to fliter low SNR pixels in ROI, 1 with ROI filtering, 0 without filtering
filtThre = 0.2; % Set the ratio factor for filtering, value within 0~1
if ifRoiFilt == 0
    filtThre = 0;
end
ifDenoi = 1;

%% Read image data
load([dataPath,filesep,'Result.mat']);
imgTmp = Result.image;
imgTmp(:,:,1:discaNum*2)=[];
imgParen = imgTmp(:,:,1:2:end);
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
[roiMask,roiNum] = drawRoi(imgParen(:,:,1),dataPath, roiName, 'gray', 'm', filtThre);
load([dataPath, filesep, 'CSF.mat'],'mask');
roiMask(mask==1) = 0; % Remove CSF

%% Calculate DGE maps
baseImgAir = mean(imgParen(:,:,(baseO2Time+baseTransTime)*4+1:(baseO2Time+baseTransTime)*4+baseO2Num),3);
baseImgO2 = mean(imgParen(:,:,1:baseO2Num),3);
for nd = 1:tn
    imgDge(:,:,nd)=(baseImgO2-imgParen(:,:,nd))./baseImgO2*100;
%     dgeImgCorr(:,:,nd)=(baseImgO2-imgParen(:,:,nd))./baseImgO2;
end
% dgeImgCorr(:,:,end-postO2Time*4+1:end) = 0 ; 
% for nd = tn-postO2Time*4+1:tn
%     dgeImgCorr(:,:,nd) = (baseImgO2-imgParen(:,:,nd))./baseImgO2;
% end

%% DGE point-wise calculation
dgeVal = zeros(tn,roiNum);
dgeValCorr = zeros(tn,roiNum);
for m = 1:tn
    img = imgDge(:,:,m);
%     imgCorr = dgeImgCorr(:,:,m);
	for n = 1:roiNum
        roiTmp = roiMask(:,:,n);
        dgeVal(m,n) = mean2(img(roiTmp>0)); 
%         dgeValCorr(m,n) = mean2(imgCorr(roiTmp>0)); 
	end
end
timePoint = (tResol*(1:tn))'/60; % sec to minute
tCorr = timePoint-baseTime;

%% Show DGE curve results
sScal = [-0.5, 10];
tScal = [min(tCorr) max(tCorr)];
figure, imagesc(imgParen(:,:,1)); axis image;colormap(gray)
hold on
contour(roiMask,1,'m-','LineWidth',2);
title('ROI')
set(gca,'Fontsize',15)
figure
set(0,'defaultfigurecolor','w') 
for n = 1:roiNum
    plot(tCorr,dgeVal(:,n),'b-o','LineWidth',2)
    hold on
end 
xlabel('Time (min)','Fontname', 'Times New Roman','FontSize',16);
ylabel('Normalized Intensity (%)','Fontname', 'Times New Roman','FontSize',16);
title(['Interleaved DGE curve of ',roiName])
axis([tScal sScal]);
% set(gcf,'Position',[300 300 1000 500]);
set(gca,'Fontsize',18)
saveas(gcf,[dataPath, filesep, roiName,'_DGEcurve.jpg']);
xlswrite( [dataPath,filesep,roiName,'_DGEcurve.xls'],[tCorr,dgeVal]);

%% Show DGE map results
sampleNum = 10;
rowNum = 6;
xScal = 3:53;
yScal = 18:81; 
cScal = [-3, 9];
answ = questdlg('Do you want to obtain DGE maps?', ...
    'Decision on DGE maps', ...
    'No','Yes','No');
switch answ
    case 'Yes'
        mask = sum(roiMask,3);
        mask_sub = mask;
        mask_sub(mask==1) = 0;
        mask_sub(mask==0) = 1e8;
        imgDge = imgDge.*repmat(sum(mask,3),[1,1,tn]);
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
        h = figure('numbertitle','off','name','DGEImg','color','white');
        imshow(imgDgeSampleBrain2D,[],'InitialMagnification','fit'); 
        mycolormap(1); colorbar; caxis(cScal); colorbar('FontName','Arial','FontSize',18, 'LineWidth', 2);
        title(['Interleaved DGE maps of ',roiName]);
%         set(gcf,'Position',[100 300 600 800]);
        export_fig([dataPath,filesep,roiName,'_DGEmap_Every',num2str(sampleNum)], '-tif', '-r200')
        anatoImg = imgParen(:,:,1).*mask;
        anatoImgSave = anatoImg(xScal,yScal);
        imwrite(anatoImgSave/max(anatoImgSave(:)),[dataPath,filesep,roiName,'_Anatomy','.tiff'],'tiff');
    case 'No'
        disp('DGE processing is done!')
end