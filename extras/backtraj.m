%% get the trajectories belonging to the background
addpath(genpath('./'));
clear all;f='A';set = 'A';
folder = ['./data/' set '/' f '/'];
fileLocation = [folder set '_' f 'Tracks291.dat'];
readTraj = ['traj_coord_' set f '_uns.mat'];
if exist(strcat(folder,readTraj),'file')==0
    readpoints(folder,fileLocation,readTraj);
end
load([folder readTraj]);
%% blockwise division
validTrajname = [folder 'validTrajFull_' set f '_uns.mat'];
% validTraj=cell(0);trajids=cell(0);excluded_ids=cell(0);
% save(validTrajname,'validTraj','trajids','excluded_ids');
% diffe=10;i=1;
%   fullTraj1([folder readTraj],100,noOfFrames,validTrajname,1);
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
   accumtraj(i,1) = sum((deltraj(i,:,1)),2);
   accumtraj(i,2) = sum((deltraj(i,:,2)),2);  
   h=fspecial('gaussian');
   v=[accumtraj(i,1) accumtraj(i,2)];
        %gaussian smoothing
        ct=filter2(h,v);
        difference=v-[ct(1) ct(2)];
        dist1(i)=norm(difference);
        if dist1(i)<19
            back(i)=1;
        else
            back(i)=0;
        end
   test(i) = norm(accumtraj(i,:));
end
%% plot the trajectories in the first image
trajidsback = trajids(find(back==1));
 img=imread([folder 'frames/0001.jpg']);
   imshow(img);
hold on;
for i=1:size(trajidsback,2)
   coordinates = TrajectoryCoordinates{trajidsback(i)}(1,:);
      plot(coordinates(1),coordinates(2),'r*');
end
trajidsback = trajids(find(back==0));
for i=1:size(trajidsback,2)
   coordinates = TrajectoryCoordinates{trajidsback(i)}(1,:);
      plot(coordinates(1),coordinates(2),'b*');
end
hold off
hist=histogram(dist1,10);

threshold = hist.BinEdges(find(hist.Values==min(hist.Values),1));%hist.BinEdges(2);%
ncl1 = nnz(dist1<threshold);
ncl2 = nnz(dist1>=threshold);
if ncl1>ncl2
    ids = find(dist1<threshold);
else
    ids = find(dist1>=threshold);
end
% ids=trajidsback;
backtrajectories = trajectories(:,ids,:);
%% model in kendall's shape space
ppp=backtrajectories;%kshapespace;%
%  figure;hold on;
% X=permute(reshape(backtrajectories(:,1,:),2,size(backtrajectories,3)),[2 1]);%Ypre = pre_shape(X);
% figure; hold on;%plot(Xpre(:,1),Xpre(:,2),'LineWidth',1);
for i=1:size(ids,2)
Y=permute(reshape(backtrajectories(:,i,:),2,size(backtrajectories,3)),[2 1]);
Ypre(:,:,i) = pre_shape(Y);
plot3(Ypre(:,1,i),Ypre(:,2,i),[1:size(Ypre,1)],'LineWidth',1);hold on
end
% X=permute(reshape(backtrajectories(:,1,:),2,size(backtrajectories,3)),[2 1]);
% for i=2:size(ids,2)
%     Y=permute(reshape(backtrajectories(:,i,:),2,size(backtrajectories,3)),[2 1]);
%    [d(i-1) kk(:,:,i-1) tr(i-1)]= procrustes(X,Y);
% kk(:,:,i-1)=pre_shape(kk(:,:,i-1));
% if d(i-1)<0.001
% plot3(kk(:,1,i-1),kk(:,2,i-1),[1:size(backtrajectories,3)],'LineWidth',1)
% end
%  hold on
% end
hold off;
%% get representative trajectory of camera motion
mean_shape = mean_shape1(Ypre);%permute(ppp,[3 1 2]));
% [mean_shape pc_shape std_shape pc_projection new_pt_out c scale]=tangent_pca_shape(permute(p,[3 1 2]));
figure;
plot3(mean_shape(:,1),mean_shape(:,2),[1:size(mean_shape,1)],'LineWidth',4)
hold on;
x=[1:size(mean_shape,1)]';
y=mean_shape(:,1);
sm_mean_shape(:,1)= smooth(x,y,0.6,'rloess');
y=mean_shape(:,2);
sm_mean_shape(:,2)= smooth(x,y,0.6,'rloess');
plot3(sm_mean_shape(:,1),sm_mean_shape(:,2),[1:size(mean_shape,1)],'LineWidth',1);
X=mean_shape;

for i=1:size(ids,2)
%% get the procrustes transformation from the mean to the traj
% Y=permute(reshape(backtrajectories(:,i,:),2,size(backtrajectories,3)),[2 1]);
% Ypre(:,:,i) = pre_shape(Y);
Y=Ypre(:,:,i);
[d Z tr] = procrustes(Y,X);
Ysmooth(:,:,i) = tr.b * sm_mean_shape * tr.T + tr.c;
% figure;hold on
% plot3(X(:,1),X(:,2),[1:size(mean_shape,1)],'LineWidth',1);
% plot3(sm_mean_shape(:,1),sm_mean_shape(:,2),[1:size(mean_shape,1)],'LineWidth',1);
% plot3(Y(:,1),Y(:,2),[1:size(mean_shape,1)],'LineWidth',1);
% plot3(Ysmooth(:,1,i),Ysmooth(:,2,i),[1:size(mean_shape,1)],'LineWidth',1);
% 
% pause(0.5); close all;
end
%% convert the smooth traj to the trajectory only if the procrustes distance is < threshold
%% for every trajectory except for those with the procrustes distance is < threshold, do the above
Trajectoryinframe=zeros(size(TrajectoryCoordinates,1),size(mean_shape,1));
Ys= zeros(size(TrajectoryCoordinates{1},1),2,size(TrajectoryCoordinates,1));
for i=1:size(Fulltrajectories,3)
%     Fulltrajectories(:,:,i)=TrajectoryCoordinates{i};
Trajectory = Fulltrajectories(:,:,i);
ind = find(Trajectory(1:size(mean_shape,1),1)~=0);
[Y c scale]=pre_shape(Trajectory(ind,:));
X=mean_shape(ind,:);
[d Z tr] = procrustes(Y,X);
% Ys=zeros(size(TrajectoryCoordinates,1),2);
if d<0.01 & nnz(trajidsback==i)==0
    Trajectoryinframe(i,ind)=1;
p = tr.b * sm_mean_shape(ind,:) * tr.T + tr.c;
% come back from preshape
Ys(ind,:,i)  = p*scale+repmat(c,size(p,1),1);
% Ys(ind,:,i) = 
%% in frame which all trajectories are available
% figure;hold on
% plot3(mean_shape(ind,1),mean_shape(ind,2),[1:size(ind,1)],'LineWidth',1);
% plot3(sm_mean_shape(ind,1),sm_mean_shape(ind,2),[1:size(ind,1)],'LineWidth',1);
% plot3(Y(:,1),Y(:,2),[1:size(Y,1)],'LineWidth',1);
% plot3(p(:,1),p(:,2),[1:size(Y,1)],'LineWidth',1);
% pause(0.5); 
end
% close all;
end
% clear TrajectoryCoordinates;
close all;
img_dir = dir([folder 'frames/']);
%% CPW
clear I1warp;
for i=1:size(mean_shape,1)
%% get the trajectory ids in frame number i
i
   testcol = Trajectoryinframe(:,i);
   trajids = find(testcol==1);
   originaltraj = permute(reshape(Fulltrajectories(i,:,trajids),2,size(trajids,1)),[2 1]);
   smoothedtraj = permute(reshape(Ys(i,:,trajids),2,size(trajids,1)),[2 1]); 
   X1 = originaltraj(:,1);
Y1 = originaltraj(:,2);
smoothX = smoothedtraj(:,1);
smoothX(smoothX<0)=1;
smoothY = smoothedtraj(:,2);
smooth(smoothY<0)=1;

% % % %% asap warping
frames(:,:,:,i)=imread([folder 'frames/' img_dir(i+2).name]);
figure;  imshow(frames(:,:,:,i));hold on;
plot(X1,Y1,'r*');
plot(smoothX,smoothY,'b*');
pause(0.5);close all;
% % I1 = imresize(img_rgb,[360 640]);
I1 = frames(:,:,:,i);%imresize(frames(:,:,:,i),[360 640]);
clear I1_features;clear I2_features;
I1_features(:,1) = X1;
I1_features(:,2) = Y1;
I2_features(:,1) = smoothX;
I2_features(:,2) = smoothY;


[height,width,~] = size(I1);
%3x3 mesh
quadWidth = width/(2^3);
quadHeight = height/(2^3);

% % % % %4x4 mesh
% % % % quadWidth = width/(2^4);
% % % % quadHeight = height/(2^4);

lamda = 1; %mesh more rigid if larger value. [0.2~5]
asap = AsSimilarAsPossibleWarping(height,width,quadWidth,quadHeight,lamda);
asap.SetControlPts(I1_features,I2_features);%set matched features
asap.Solve();            %solve Ax=b for as similar as possible
homos = asap.CalcHomos();% calc local hommograph transform

gap = 50;
I1 = asap.Warp(I1,gap);                     %warp source image to target image
I1warpmesh = asap.destin.drawMesh(I1,gap);  %draw mesh on the warped source image
imshow(I1warpmesh);

gap = 50;
I1warp{i} = asap.Warp(I1,gap); %= I1;%                    
imwrite(I1warp{i},[folder 'warped/warp' num2str(i) '.png']);
imshow(I1warp{i});pause(0.5);
%access local homography
[h,w,~,~] = size(homos);
for k=1:h-1
    for j=1:w-1
       H(:,:,i) = homos(k,j,:,:);
%        fprintf('Quad=[%d %d]\n',i,j);
%        H(:,:,k)
    end
end
close all;
end
cc = postprocessing(I1warp);

for i=1:size(mean_shape,1)
    cutI1warp{i} = I1warp{i}(50:size(I1warp{i},1)-50,100:size(I1warp{i},2)-100,:);
    imwrite(cutI1warp{i},[folder 'warped/warp' num2str(i) '.png']);

end

result = double(cutI1warp);
% postprocessing(result,[folder,'finalFrames/'])
WriteVideoAVI('output.avi',folder,0.5,result);
WriteVideoAVI('comparison.avi',folder,0.5,ComparisionIM);

disp('[DONE]');


