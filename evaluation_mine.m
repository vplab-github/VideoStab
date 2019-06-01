%% get global homography
clear all;
clc;
addpath('./MotionEstimation/');
addpath('./MotionEstimation/RANSAC');
addpath('./videoprocessing');
folder = './data/D/A/';            
imdir1=dir([folder 'frames/*.ppm']);
% imdir2 = dir([folder 'warped/*.png']);%'F:\research\DATA\datasets\SegTrackv2\JPEGImages\frog\*.png');
% % % folder1 = './comparison/sample_videos/';
% % % folder2 = 'F:\research\DATA\STABILIZATION\input\with_moving_object\unstable\C\';
% % % InputVideoName1 = 'C/CYoutube.mp4';
% % % InputVideoName2='example7_input.avi';
% % % [numFrames,frameRate,fr2]=ReadVideoAVI(InputVideoName1,folder1);
% % % [numFrames,frameRate,fr1]=ReadVideoAVI(InputVideoName2,folder2);
numFrames = 99;
% % % frames1=fr1(:,:,:,29:numFrames+28);
% % % frame2 = fr2(:,:,:,1:numFrames);
for i=1:numFrames
    frames1(:,:,:,i)=imread([folder 'frames/' imdir1(i).name]);
%         frames2(:,:,:,i)=imread(['F:\research\DATA\datasets\SegTrackv2\JPEGImages\frog\' imdir2(i).name]);

%      frame2(:,:,:,i)=imread([folder 'warped/' imdir2(i).name]);
end
load([folder 'resultFrames.mat']);
frame2=result;
% resize to the input video size
frames2=zeros(size(frames1));
for i=1:size(frame2,4)
    frames2(:,:,:,i) = imresize(frame2(:,:,:,i),[size(frames1,1) size(frames1,2)]);
end
frames2=uint8(frames2);

%estimate homography between neighboring frames
[h,w,c,~] = size(frames1(:,:,:,:));
% if exist(strcat(folder,'/homography.mat'),'file')==0
    FList = CalcMotion(folder,frames1,frames2);

%%cropping
for i=1:numFrames
% get scale component
a=FList(1,1,i);b=FList(1,2,i);c(i)=FList(1,3,i);
d=FList(2,1,i);e=FList(2,2,i);f(i)=FList(2,3,i);

 p = sqrt(a*a +b*b);  r = (a*e - b*d)/(p);
  scalex = p; scaley = r;
[u s v]=svd(FList(1:2,1:2,i));
ratioeigen = s(1,1)/s(2,2);
end
% get average over frames
avgscale = (sum(scalex)/size(scalex,2)+ sum(scaley)/size(scaley,2))/2
avgdist = sum(ratioeigen)/size(scaley,2)
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
[avgscale avgdist Etx Ety Et]
% get energy percentage
save([folder 'scores.mat'],'avgscale','avgdist','Etx','Ety','Et','Stability_score')