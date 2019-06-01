clc;clear all;close all;
%----------------------------------
rand('state',0);
%----------------------------Set path
ip_dir = './data/C/A/';
img_dir = dir([ip_dir 'frames/*.ppm']);%['E:\research\DATA\STABILIZATION\frames_unstable\C\*.jpg']);
op_dir  = './output/';
load([ip_dir 'validTrajFull_CA_uns.mat'])
%----------------------------------
addpath(ip_dir);
%----------------------------------------
%img_dir=sort_nat(img_dir);
img1 = imread([ip_dir 'frames/' img_dir(1).name]);
initstate = [0 0 size(img1,1) size(img1,2)];

close all;
%-------------------------------------------------
num = length(img_dir);% number of frames
%--------------------------------------------------------
x = initstate(1);% x axis at the Top left corner
y = initstate(2);
w = initstate(3);% width of the rectangle
h = initstate(4);% height of the rectangle
%--------------------------------------------------------
flag = 0;
for i = 1:num-1
    
img1 = imread([ip_dir 'frames/' img_dir(1).name]);

% % % % if size(img1,3) == 3
% % % %    img1 = rgb2gray(img1);
% % % % end

imgP=img_dir(i).name;
frames{i}=img1;    
    img2=img1;
    img = img_dir(i+1).name;%imread(img_dir(i).name);
%     img1=imread(img);
% % % %    if size(img1,3) == 3
% % % %    img1 = rgb2gray(img1);
% % % %    end
% % % % frames{i+1}=img1;
% [desc1,desc2,corresp1,corresp2]=sift_corresp(imgP,img);
% % % %             sdesc{i}=desc1;
% % % %             sdesc{i+1}=desc2;
% % % %           sqrd_dist=(corresp1-corresp2).^2;
% % % %           distances=sqrd_dist(:,1)+sqrd_dist(:,2);
% % % %           matches{i,1}=corresp1(distances<5000,:);
% % % %           matches{i,2}=corresp2(distances<5000,:);
           corresp1=reshape(validTraj{1}(1:2,1:40:end,i),size(validTraj{1}(1:2,1:40:end,i),2),2);
           corresp2=reshape(validTraj{1}(1:2,1:40:end,i+1),size(validTraj{1}(1:2,1:40:end,i+1),2),2);
%     %plotting only 10th frame
     if i>=2
      
        % Plot only the first frame
        if i==2
            
            plotFrames3D(initstate,img1,i-2);
        % For all other patches
        else
           % cd sift
            
           % mCost = matchSift(des1,des2);
            PlotFrameAndSift_mine(initstate,corresp1,corresp2,i,img1);
            flag = 0;
            %length(loc1)
            %cd ..
        end
        imgP=img;
        x_prev = initstate(1);
        z_prev = initstate(2);
        
        
     end
   
    
end  