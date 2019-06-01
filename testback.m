%% get the trajectories belonging to the background
addpath(genpath('./'));
clear all;f='D';set = 'frog';
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
diffe=10;i=1;
%   fullTraj1([folder readTraj],50,noOfFrames,validTrajname,1);
% fullTraj_adaptive([folder readTraj],i,diffe,noOfFrames,validTrajname);
load(validTrajname);
%% get all the trajectories
trajectories = validTraj{1}(1:2,:,:);trajids = validTraj{1}(3,:,1);
pertraj = permute(trajectories,[2 3 1]);
%% accumulated motion vector: get the histogram separate and get the threshold
% deltraj=zeros(size(pertraj));
% accumtraj=zeros(size(pertraj,1),2);
for i=1:size(pertraj,1)
     pptraj = pre_shape(reshape(pertraj(i,:,:),size(pertraj,2),size(pertraj,3)));
   deltraj(i,:,1)=diff(pptraj(:,1));
   deltraj(i,:,2)=diff(pptraj(:,2));
   accumtraj(i,1) = sum(((deltraj(i,:,1))),2);
   accumtraj(i,2) = sum((deltraj(i,:,2)),2); 
   h=fspecial('gaussian');
   v=[accumtraj(i,1) accumtraj(i,2)];
%    v=v./norm(v);
        %gaussian smoothing
        ct(1)=filter2(h,v(1));
         ct(2)=filter2(h,v(2));
        difference=v-[ct(1) ct(2)];
        dist1(i)=norm(difference);
       
   test(i) = norm(accumtraj(i,:));
end
hist(dist1)

%% plot the trajectories in the first image
averageAcc = 1/size(accumtraj,1) * sum(accumtraj,1);
difference = sqrt(sum((accumtraj-repmat(averageAcc,size(accumtraj,1),1)).^2,2))
hist(difference,5)
back=difference<0.3;
trajidsback = trajids(find(back==1));
for kk=1:100
 img=imread([folder 'frames/000' num2str(kk) '.ppm']);
  h(:,:,:,kk)= imshow(img);
hold on;
for i=1:size(trajidsback,2)
   coordinates = TrajectoryCoordinates{trajidsback(i)}(1,:);
      plot(coordinates(1),coordinates(2),'r*');
end
trajidsfore = trajids(find(back==0));
for i=1:size(trajidsfore,2)
   coordinates = TrajectoryCoordinates{trajidsfore(i)}(1,:);
      plot(coordinates(1),coordinates(2),'b*');
end
hold off
end