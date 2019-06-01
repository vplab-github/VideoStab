%stitch frames together
% read videos to frames1, frames2, frames3
clear all;
addpath('./videoprocessing/')
InputVideoName1='lim2.avi';
InputVideoName2='limSteadyflow.avi';
InputVideoName3='output.avi';

% input
% folder1='E:\research\results\Stabilization\bundled_paths\running\run5\';
% folder2='E:\research\results\Stabilization\YoutubeResults\running\';%bundled_paths\running\run5\';%WarpStabilizerOutput\running\';
% folder3='E:\research\codes\MY CODES\Work3\data\running\run5\';
folder2 = 'E:\research\DATA\datasets\SegTrackv2\JPEGImages\frog\';
folder3='E:\research\results\Stabilization\warpStabilizerOutput\lim\'
% folder3= 'E:\research\codes\steadyflo_from_site\SinglePathStabilizerV1.0\data\frogD\';
folder4='./data/limitations/lim2/'
load([folder4 'resultFrames.mat'])
folderout = './resultppt/';
[numFrames1,frameRate1,frames1]=ReadVideoAVI(InputVideoName1,folder4);
% [numFrames2,frameRate2,frames2]=ReadVideoAVI(InputVideoName2,folder2);
[numFrames3,frameRate3,frames3]=ReadVideoAVI(InputVideoName2,folder3);
frameno = min([numFrames1 numFrames3 size(result,4)]);
for i=1:frameno
   fr = video_horizontal(frames1(:,:,:,i),imresize(frames3(:,:,:,i),[size(frames1,1) size(frames1,2)]),50);
%    fr1= video_horizontal(fr,imresize(frames3(:,:,:,i),[size(frames1,1) size(frames1,2)]),50);
   ComparisonIM(:,:,:,i) = video_horizontal(fr,imresize(result(:,:,:,i),[size(frames1,1) size(frames1,2)]),50);

clear fr
end
WriteVideoAVI('lim2steadyflow.avi',folderout,2,ComparisonIM);
