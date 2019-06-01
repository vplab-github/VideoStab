%% get global homography
clear all;
clc;
addpath('./MotionEstimation/');
addpath('./MotionEstimation/RANSAC');
addpath('./videoprocessing');
folder = './data/D/A/';  
% folder = 
% imdir1=dir([folder 'frames/*.jpg']);
f='simple';f1='6';f2='6_ours';
% imdir2 = dir([folder 'warped/*.png']);%'F:\research\DATA\datasets\SegTrackv2\JPEGImages\frog\*.png');
folder1=['E:\research\results\Stabilization\bundled_paths\' f '\' f f1 '\'];
% folder1=['E:\research\results\Stabilization\YoutubeResults\' f '\'];
% folder1=['E:\research\results\Stabilization\warpStabilizerOutput\' f '\'];
% folder1 = './comparison/warpStabilizer/';
% folder1 = folder;
folder2 = folder;
%output
InputVideoName1 = ['6_ours.avi'];
%input
InputVideoName2=['6.avi'];
[numFrames1,frameRate1,fr2]=ReadVideoAVI(InputVideoName1,folder1);
[numFrames2,frameRate2,fr1]=ReadVideoAVI(InputVideoName2,folder1);

numFrames=min(min(numFrames1,numFrames2),99);
frame1=fr1(:,:,:,1:numFrames);
frame2 = fr2(:,:,:,1:numFrames);
%% resize to the input video size
frames2=zeros(360,640,3,numFrames);
for i=1:numFrames
    frames1(:,:,:,i) = imresize(frame1(:,:,:,i),[360 640]);
    frames2(:,:,:,i) = imresize(frame2(:,:,:,i),[360 640]);
end
frames2=uint8(frames2);

%estimate homography between neighboring frames
[h,w,c,~] = size(frames1(:,:,:,:));
% if exist(strcat(folder,'/homography.mat'),'file')==0
    FList = CalcMotion(folder,frames1,frames2);

%%cropping
for i=1:numFrames%numFrames
% get scale component
a=FList(1,1,i);b=FList(1,2,i);c=FList(1,3,i);
d=FList(2,1,i);e=FList(2,2,i);f=FList(2,3,i);
 p = sqrt(a*a + b*b);  r = (a*e - b*d)/(p);
  scalex(i) = p; scaley(i) = r;
[u s v]=svd(FList(1:2,1:2,i));
ratioeigen(i) = s(1,1)/s(2,2);
end
% get average over frames
avgscale = (sum(scalex)/size(scalex,2)+ sum(scaley)/size(scaley,2))/2
avgdist = sum(ratioeigen)/size(scalex,2)
%%distortion
% get affine part
% getfirst two eigen values and ratio

%% stability
% get the stabilized video feature tracks from surf features
    FListst = CalcMotion_st(frames2);
% %        FListun = CalcMotion_st(folder,frames1);
        grayI = rgb2gray(frames2(:,:,:,2));

    points = detectSURFFeatures(grayI);
coordinates=points.Location;
% for i=1:numFrames
% %     homocoord(:,1:2)=coordinates;
% %     homocoord(:,3) = ones(size(coordinates,1),1);
% %     new(:,:,i) = double(FListst(:,:,i)*homocoord');
% 
% a(i)=FListst(1,1,i);b(i)=FListst(1,2,i);c(i)=FListst(1,3,i);
% d=FListst(2,1,i);e=FListst(2,2,i);f(i)=FListst(2,3,i);
% theta(i)=atan2d(b(i),a(i));
% end
%  newcoord=permute(new,[1 3 2]);
%  traj = newcoord(1:2,:,i);
%  fft(traj(1,:))
% get fft
magf=abs(fft(reshape(FListst(2,3,1:numFrames),1,numFrames))).^2;
Etx = sum(magf(2:6))/sum(magf(2:numFrames/2))
magc=abs(fft(reshape(FListst(1,3,1:numFrames),1,numFrames))).^2;
Ety = sum(magc(2:6))/sum(magc(2:numFrames/2))
% mags=abs(fft(reshape(FListst(1,2,1:99),1,99))).^2;
% Es = sum(mags(2:6))/sum(mags(2:numFrames/2))

magt=abs(fft(atan2d(reshape(FListst(1,2,1:numFrames),1,numFrames),reshape(FListst(1,1,1:numFrames),1,numFrames)))).^2;
Et = sum(magt(2:6))/sum(magt(2:numFrames/2))
% 
% magfftc=abs(fft(c)).^2;
% Etx = sum(magfftc(2:6))/sum(magfftc);
% magfftf=abs(fft(f)).^2;
% Ety = sum(magfftf(2:6))/sum(magfftf);
% magfftt=abs(fft(theta)).^2;
% Etheta = sum(magfftt(2:6))/sum(magfftt);
Stability_score= min([Etx Ety Et])
% get energy percentage
[Etx Ety Et avgscale avgdist]
save([folder1 'scores.mat'],'avgscale','avgdist','Etx','Ety','Et','Stability_score')