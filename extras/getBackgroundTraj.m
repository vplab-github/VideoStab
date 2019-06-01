%% get the trajectories belonging to the background
addpath(genpath('./'));
clear all;f='A';set = 'worm';
folder = ['./data/' set '/' f '/'];
fileLocation = [folder set '_' f 'Tracks243.dat'];
readTraj = ['traj_coord_' set f '_uns.mat'];
if exist(strcat(folder,readTraj),'file')==0
    readpoints(folder,fileLocation,readTraj);
end
load([folder readTraj]);
%% blockwise division
validTrajname = [folder 'validTrajFull_' set f '_uns.mat'];
% validTraj=cell(0);trajids=cell(0);excluded_ids=cell(0);
% save(validTrajname,'validTraj','trajids','excluded_ids');
diffe=10;i=1;
% %   fullTraj1([folder readTraj],7,noOfFrames,validTrajname,1);
% fullTraj_adaptive([folder readTraj],i,diffe,noOfFrames,validTrajname);
load(validTrajname);
%% get all the trajectories
trajectories = validTraj{1}(1:2,:,:);trajids = validTraj{1}(3,:,1);
pertraj = permute(trajectories,[2 3 1]);
%% accumulated motion vector: get the histogram separate and get the threshold
% deltraj=zeros(size(pertraj));
% accumtraj=zeros(size(pertraj,1),2);
for i=1:size(pertraj,1)
   deltraj(i,:,1)=diff(pertraj(i,:,1));
   deltraj(i,:,2)=diff(pertraj(i,:,2));
   accumtraj(i,1) = sum((deltraj(i,:,1)));
   accumtraj(i,2) = sum((deltraj(i,:,2)));  
%    h=fspecial('gaussian');
%    v=[accumtraj(i,1) accumtraj(i,2)];
%         %gaussian smoothing
%         ct=filter2(h,v);
%         difference=v-[ct(1) ct(2)];
%         dist(i)=norm(difference);
   test(i) = norm(accumtraj(i,:));
end

%% get the ids for the camera motion group
hist=histogram(test,50);
threshold = find(hist.Values==min(hist.Values),1);%hist.BinEdges(2);%
ids = find(test<threshold);
backtrajectories = trajectories(:,ids,:);
%% model in kendall's shape space
ppp=backtrajectories;%kshapespace;%

%% get representative trajectory of camera motion
mean_shape = mean_shape(permute(ppp,[3 1 2]));
% [mean_shape pc_shape std_shape pc_projection new_pt_out c scale]=tangent_pca_shape(permute(p,[3 1 2]));
figure;
plot3(mean_shape(:,1),mean_shape(:,2),[1:size(mean_shape,1)],'LineWidth',1)
hold on;
x=[1:size(mean_shape,1)];
y=mean_shape(:,1);
sm_mean_shape(:,1)= smooth(x,y,0.7,'rloess');
y=mean_shape(:,2);
sm_mean_shape(:,2)= smooth(x,y,0.7,'rloess');
plot3(sm_mean_shape(:,1),sm_mean_shape(:,2),[1:size(mean_shape,1)],'LineWidth',1);
% axis('equal');
pause(0.5);
difference = mean_shape-sm_mean_shape;
%% regression of mean shape
TrajCoordinates = permute(backtrajectories,[3 1 2]);
% smoothedTrajectories = zeros(size(TrajectoryCoordinates{1},1),size(TrajectoryCoordinates{1},2),size(TrajectoryCoordinates,1));
smoothedTrajectories = zeros(size(TrajCoordinates));
for i=1:size(TrajCoordinates,3)%size(TrajectoryCoordinates,1)%
%% take every trajectory to the pre_shape
%  TrajCoordinates(:,:,i)=TrajectoryCoordinates{i};
  [row,~] = find(TrajCoordinates(:,:,i)~=0);
  uniquerows = unique(row);
 
  nzrows=uniquerows(uniquerows<=127);
  [traj c scale]=pre_shape(TrajCoordinates(nzrows,:,i));
% %   x=[1:size(traj,1)];
% % y=traj(:,1);
% % p(:,1)= smooth(x,y,0.3,'rloess');
% % y=traj(:,2);
% % p(:,2)= smooth(x,y,0.3,'rloess');

   origp = TrajCoordinates(nzrows,:,i);
   
%% add the  difference
 p=traj-difference(nzrows,:);
 
%% convert it back to the original space
  pp  = p*scale+repmat(c,size(p,1),1);
  smoothedTrajectories(nzrows,:,i)=pp;
  
%% plot smoothed and original
% figure;  
% plot3(origp(:,1),origp(:,2),[1:size(origp,1)],'LineWidth',1);
% hold on
% plot3(pp(:,1),pp(:,2),[1:size(pp,1)],'LineWidth',1);
% % axis('equal');
% pause(0.5)
% % 
% close all;
% clear p;clear pp; clear traj;clear origp;
end
hold off

img_dir = dir([folder 'frames/']);
%% warp CPW
% % % % for k=1:168
% % % %     (diff-overlap)*(l-1)+k
% % % X1 = nonzeros(TrajCoordinates(k,1,:));
% % % Y1 = nonzeros(TrajCoordinates(k,2,:));
% % % smoothX = nonzeros(smoothedTrajectories(k,1,:));
% % % smoothY = nonzeros(smoothedTrajectories(k,2,:));
% % % %% asap warping
% % % frames(:,:,:,k)=imread([folder 'frames/' img_dir(k+2).name]);

% % % I1 = imresize(img_rgb,[360 640]);
% % % % I2 = imresize(frames(:,:,:,2),[360 640]);
% % % clear I1_features;clear I2_features;
% % % I1_features(:,1) = X1;
% % % I1_features(:,2) = Y1;
% % % I2_features(:,1) = smoothX;
% % % I2_features(:,2) = smoothY;
% % % % % % 
% % % [height,width,~] = size(I1);
% % % %3x3 mesh
% % % quadWidth = width/(2^5);
% % % quadHeight = height/(2^5);
% % % 
% % % % %4x4 mesh
% % % % quadWidth = width/(2^4);
% % % % quadHeight = height/(2^4);
% % % 
% % % lamda = 1; %mesh more rigid if larger value. [0.2~5]
% % % asap = AsSimilarAsPossibleWarping(height,width,quadWidth,quadHeight,lamda);
% % % asap.SetControlPts(I1_features,I2_features);%set matched features
% % % asap.Solve();            %solve Ax=b for as similar as possible
% % % homos = asap.CalcHomos();% calc local hommograph transform
% % % 
% % % gap = 50;
% % % I1 = asap.Warp(I1,gap);                     %warp source image to target image
% % % I1warpmesh = asap.destin.drawMesh(I1,gap);  %draw mesh on the warped source image
% % % imshow(I1warpmesh);
% % % 
% % % gap = 50;
% % % I1warp(:,:,:,k) = asap.Warp(I1,gap);                     
% % % imwrite(I1warp(:,:,:,k),[folder 'warped/warp' num2str(k) '.png']);
% % % imshow(I1warp(:,:,:,k));pause(0.5);
% % % %access local homography
% % % [h,w,~,~] = size(homos);
% % % for i=1:h-1
% % %     for j=1:w-1
% % %        H(:,:,k) = homos(i,j,:,:);
% % % %        fprintf('Quad=[%d %d]\n',i,j);
% % % %        H(:,:,k)
% % %     end
% % % end
% % % end
% % % postprocessing(I1warp,[folder 'result/'],1)

disp('rendering...');
numFrames=127;
for i=1:numFrames
  fprintf('%3d / %3d\n',i,numFrames);
  frames(:,:,:,i)=imread([folder 'frames/' img_dir(i+2).name]);
   X1 = nonzeros(TrajCoordinates(i,1,:));
Y1 = nonzeros(TrajCoordinates(i,2,:));
smoothX = nonzeros(smoothedTrajectories(i,1,:));
smoothY = nonzeros(smoothedTrajectories(i,2,:));
I1_features(:,1) = X1;
I1_features(:,2) = Y1;
I2_features(:,1) = smoothX;
I2_features(:,2) = smoothY;

  [H,~] = EstimateHomographyByRANSAC(I1_features',I2_features', 0.001);
updateSet=H;
  result(:,:,:,i) = HomographyWarp(double(frames(:,:,:,i)),inv(updateSet));
  ComparisionIM(:,:,:,i) = video_horizontal(frames(:,:,:,i),result(:,:,:,i),20);
  imshow(result(:,:,:,i));pause(0.5)
  imwrite(result(:,:,:,i),[folder 'warped/warp' num2str(i) '.png']);
  clear I1_features;clear I2_features;
end
result = double(result);
% postprocessing(result,[folder,'finalFrames/'])
WriteVideoAVI('output.avi',folder,0.5,result);
WriteVideoAVI('comparison.avi',folder,0.5,ComparisionIM);

disp('[DONE]');

