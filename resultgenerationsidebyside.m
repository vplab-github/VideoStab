folder ='E:\research\DATA\STABILIZATION\input\datasets\crowd\crowd7\test\';
InputVideoName1='input.avi';
InputVideoName2='7.avi';
InputVideoName3='ours.avi';
[numFrames1,frameRate1,frames1]=ReadVideoAVI(InputVideoName1,folder);
[numFrames2,frameRate2,frames2]=ReadVideoAVI(InputVideoName2,folder);
[numFrames3,frameRate3,frames3]=ReadVideoAVI(InputVideoName3,folder);

for i=1:200
     k = video_horizontal(frames1(:,:,:,i),imresize(frames2(:,:,:,i),[size(frames,1) size(frames,2)]),30);
    ComparisionIM(:,:,:,i) = video_horizontal(k,imresize(frames3(:,:,:,i),[size(frames,1) size(frames,2)]),30);
end
WriteVideoAVI('output.avi',folder,frameRate1,ComparisionIM);
