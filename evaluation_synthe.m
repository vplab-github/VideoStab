% read images 
% groundtruth
folder='./data/frog/D/';
GT = 'E:\research\DATA\datasets\SegTrackv2\JPEGImages\frog\';
GT_dir = dir(GT);
% jittery
JT = './data/frog/D/frames/';
JT_dir = dir(JT);
% stabilized
% warp stabilizer
ST = './data/frog/D/warpSt/';
ST_dir = dir(ST);
numFrames=99;
for i=3:101
    GT_img(:,:,:,i-2) = imread([GT GT_dir(i).name]);
    JT_img(:,:,:,i-2) = imresize(imread([JT JT_dir(i).name]),[size(GT_img,1) size(GT_img,2)]);
    ST_img(:,:,:,i-2) = imresize(imread([ST ST_dir(i).name]),[size(GT_img,1) size(GT_img,2)]);
end
load([folder 'resultFrames.mat'])
STo_img=result(:,:,:,1:numFrames);
[h,w,c,~] = size(GT_img);
frames2=uint8(STo_img);
frames2=(ST_img);

% if exist(strcat(folder,'/homography.mat'),'file')==0
    FListO = CalcMotion(folder,GT_img,frames2);
numFrames=99;
%%cropping
for i=1:numFrames
% get scale component
a=FListO(1,1,i);b=FListO(1,2,i);c(i)=FListO(1,3,i);
d=FListO(2,1,i);e=FListO(2,2,i);f(i)=FListO(2,3,i);

 p = sqrt(a*a +b*b);  r = (a*e - b*d)/(p);
  scalex = p; scaley = r;
[u s v]=svd(FListO(1:2,1:2,i));
ratioeigen = s(1,1)/s(2,2);
end
% get average over frames
avgscale = (sum(scalex)/size(scalex,2)+ sum(scaley)/size(scaley,2))/2
avgdist = sum(ratioeigen)/size(scaley,2)
%% stability
    FListst = CalcMotion_st(frames2);
% %        FListun = CalcMotion_st(folder,frames1);
        grayI = rgb2gray(frames2(:,:,:,1));
% get fft
numFrames=99;
magf=abs(fft(reshape(FListst(2,3,1:numFrames),1,numFrames))).^2;
Etx = sum(magf(2:6))/sum(magf(2:ceil(numFrames/2)))
magc=abs(fft(reshape(FListst(1,3,1:numFrames),1,numFrames))).^2;
Ety = sum(magc(2:6))/sum(magc(2:ceil(numFrames/2)))
% magt=abs(fft(reshape(FListst(1,2,1:numFrames),1,numFrames))).^2;
% Es = sum(mags(2:6))/sum(mags(2:numFrames/2))

magt=abs(fft(atan2d(reshape(FListst(1,2,1:numFrames),1,numFrames),reshape(FListst(1,1,1:numFrames),1,numFrames)))).^2;
Et = sum(magt(2:6))/sum(magt(2:ceil(numFrames/2)))
% 
% magfftc=abs(fft(c)).^2;
% Etx = sum(magfftc(2:6))/sum(magfftc);
% magfftf=abs(fft(f)).^2;
% Ety = sum(magfftf(2:6))/sum(magfftf);
% magfftt=abs(fft(theta)).^2;
% Etheta = sum(magfftt(2:6))/sum(magfftt);
Stability_score= min([Etx Ety Et])
avgscale
avgdist
% get energy percentage
% save([folder 'scores.mat'],'avgscale','avgdist','Etx','Ety','Et','Stability_score')
% scores of input
% scores of output