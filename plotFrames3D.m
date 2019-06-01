function plotFrames3D(bbox,imgData,f_ind)
%     clc;clear all;close all;
%     %----------------------------Set path
%     ip_dir = '//media/Local Disk/Joy/data/biking1/';
%     img_dir = dir([ip_dir '*.jpg']);
%     op_dir  = './output/';
%     %----------------------------------
%     addpath(ip_dir);
%
%     num = length(img_dir);% number of frames
%     % For all frames of the video
%     i=1;
%     while i<=num
%         img = imread(img_dir(i).name);
%         xImage = [-0.5 0.5; -0.5 0.5];   %# The x data for the image corners
%         yImage = i*ones(2);             %# The y data for the image corners
%         zImage = [0.5 0.5; -0.5 -0.5];   %# The z data for the image corners
%         surf(xImage,yImage,zImage,...    %# Plot the surface
%             'CData',img,...
%             'FaceColor','texturemap');
%         hold on
%         i=i+5;
%     end

xImage = [bbox(1) bbox(1)+bbox(4);bbox(1) bbox(1)+bbox(4)];
yImage = f_ind*ones(2);
zImage = [bbox(2) bbox(2);bbox(2)+bbox(4) bbox(2)+bbox(4)];


surf(xImage,yImage,zImage,...    %# Plot the surface
    'CData',flipdim(imgData,1),...
    'FaceColor','texturemap');
hold on

end
 

