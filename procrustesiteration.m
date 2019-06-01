addpath(genpath('./'));
clear all;close all;f='D';set = 'frog';
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
diffe=10;i=1;l=1;
%   fullTraj1([folder readTraj],100,noOfFrames,validTrajname,1);
% fullTraj_adaptive([folder readTraj],i,diffe,noOfFrames,validTrajname);
load(validTrajname);
%% get all the trajectories
trajectories = validTraj{1}(1:2,:,:);trajids = validTraj{1}(3,:,1);
pertraj = permute(trajectories,[2 3 1]);
%% accumulated motion vector: get the histogram separate and get the threshold

for i=1:size(pertraj,1)
    pptraj = pre_shape(reshape(pertraj(i,:,:),size(pertraj,2),size(pertraj,3)));
 back(i)=1;
end
for k=1:2
trajidsback = trajids(find(back==1));
 img=imread([folder 'frames/0001.ppm']);
   imshow(img);
hold on;
for i=1:size(trajidsback,2)
   coordinates = TrajectoryCoordinates{trajidsback(i)}(1,:);
      plot(coordinates(1),coordinates(2),'r*');
      TrajCoord(:,:,i)= TrajectoryCoordinates{trajidsback(i)};
end
trajidsfore = trajids(find(back==0));
for i=1:size(trajidsfore,2)
   coordinates = TrajectoryCoordinates{trajidsfore(i)}(1,:);
      plot(coordinates(1),coordinates(2),'b*');
end
hold off

%% frechet mean
for i=1:nnz(back)
    trajid = trajidsback(i);
validtrajid = find(validTraj{l}(3,:,1)==trajid);
% TrajCoord(:,:,i) = permute(reshape(validTraj{l}(1:2,validtrajid,:),2,size(validTraj{l},3)),[2 1]);

Y(:,:,i)=reshape(TrajCoord(:,:,i),size(TrajCoord,1),2);
end

hold off;
[mean_shape,tran]= mean_shape1(Y);%permute(ppp,[3 1 2]));
% [mean_shape pc_shape std_shape pc_projection new_pt_out c scale]=tangent_pca_shape(permute(p,[3 1 2]));
figure;
plot3(mean_shape(:,1),mean_shape(:,2),[1:size(mean_shape,1)],'LineWidth',4)
hold on;
x=[1:size(mean_shape,1)]';
y=mean_shape(:,1);
sm_mean_shape(:,1)= smooth(x,y,0.7,'rloess');
y=mean_shape(:,2);
sm_mean_shape(:,2)= smooth(x,y,0.7,'rloess');
plot3(sm_mean_shape(:,1),sm_mean_shape(:,2),[1:size(mean_shape,1)],'LineWidth',1);
legend('frechet mean','smoothed frechet mean')
X=mean_shape;
hold off
%% get the procrustes distance 
for i=1:nnz(back)
%     Fulltrajectories(:,:,i)=TrajectoryCoordinates{i};
Trajectory = Y(:,:,i);%TrajCoord(:,:,i);
ind = find(Trajectory(1:size(mean_shape,1),1)~=0);
[Yy c scale]=pre_shape(Trajectory(ind,:));
% X=mean_shape(ind,:);
[d(i) Z tr] = procrustes(Yy,X(ind,:));
p = sm_mean_shape(ind,:) * tran(:,:,i)'+tr.c;
% come back from preshape
  kkk= p*scale+repmat(c,size(p,1),1);
  StabOTDM(ind,:,i)=kkk;


end
back1=back;
back(d>(0.7)/k)=0;
% back=difference<0.1;
trajidsback = trajids(find(back1==1));

% for i=1:nnz(back)
%     trajid = trajidsback(i);
% validtrajid = find(validTraj{l}(3,:,1)==trajid);
% % TrajCoord(:,:,i) = permute(reshape(validTraj{l}(1:2,validtrajid,:),2,size(validTraj{l},3)),[2 1]);
% 
% Y(:,:,i)=reshape(TrajCoord(:,:,i),size(TrajCoord,1),2);
% end

end
