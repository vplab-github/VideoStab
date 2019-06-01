addpath(genpath('./'));
clear all;f='A';set = 'car1';
folder = ['./data/' set '/' f '/'];
fileLocation = [folder f 'Tracks200.dat'];
InputVideoName1='1.avi';thresholdtraj=0.7;iterations=20;lambda1=0.5;
diffe=50;overlap1=0;
% fileLocation = [folder set '_' f 'Tracks291.dat'];
readTraj = ['traj_coord_' set f '_uns.mat'];
 if exist(strcat(folder,readTraj),'file')==0
    readpoints(folder,fileLocation,readTraj);
 end
 
load([folder readTraj]);
%% blockwise division
validTrajname = [folder 'validTrajFull_' set f '_uns.mat'];
validTraj=cell(0);trajids=cell(0);excluded_ids=cell(0);
save(validTrajname,'validTraj','trajids','excluded_ids');
i=1;
  fullTraj1([folder readTraj],diffe,noOfFrames,validTrajname,overlap1);
% fullTraj_adaptive([folder readTraj],i,diffe,noOfFrames,validTrajname);
load(validTrajname);
%% get all the trajectories
for l=1:size(validTraj,2)
    
trajectories = validTraj{l}(1:2,:,:);trajids = validTraj{l}(3,:,1);
pertraj = permute(trajectories,[2 3 1]);
% [numFrames,frameRate,frames]=ReadVideoAVI(InputVideoName1,folder);
jj=1;clear backcoord
for indexbak=1:2
clear TrajCoord;clear accumtraj;clear back;
    
for i=1:size(pertraj,1)
    pptraj = pre_shape(reshape(pertraj(i,:,:),size(pertraj,2),size(pertraj,3)));
     deltraj(:,1)=diff(pptraj(:,1));
   deltraj(:,2)=diff(pptraj(:,2));
   accumtraj(i,1) = sum(((deltraj(:,1))));
   accumtraj(i,2) = sum((deltraj(:,2))); 
   h=fspecial('gaussian');
   v=[accumtraj(i,1) accumtraj(i,2)];
%    v=v./norm(v);
        %gaussian smoothing
        ct(1)=filter2(h,v(1));
         ct(2)=filter2(h,v(2));
        difference=v-[ct(1) ct(2)];
        dist1(i)=norm(difference);
%        plot(1:size(deltraj,1),deltraj(:,1),'r');
%        hold on;
   test(i) = norm(accumtraj(i,:));
end
%% average of accumulated traj
averageAcc = 1/size(accumtraj,1) * sum(accumtraj,1);
difference = sqrt(sum((accumtraj-repmat(averageAcc,size(accumtraj,1),1)).^2,2));
[histo,centres]=hist(difference,10)

an=[histo(1) sum(histo(1:2)) sum(histo(1:3)) sum(histo(1:4)) sum(histo(1:4)) sum(histo(1:5)) sum(histo(1:6)) sum(histo(1:7)) sum(histo(1:8)) sum(histo(1:9)) sum(histo(1:10))]<thresholdtraj*size(difference,1);
thresholdinitial = centres(find(an==1,1,'last'));
if nnz(an==1)==0
    thresholdinitial=centres(1);
end
% earlierback=back;
back=difference<1;%thresholdinitial;
% back=~back;
trajidsback = trajids(back==1);
n=sprintf('%04d.ppm',(l-1)*(diffe-overlap1)+1);
 img=imread([folder 'frames/' n]);
%  img=frames(:,:,:,1);
   imshow(img);
hold on;
for i=1:size(trajidsback,2)
   coordinates = TrajectoryCoordinates{trajidsback(i)}(l*diffe-overlap1,:);
      plot(coordinates(1),coordinates(2),'r*');

end
trajidsfore = trajids(back==0);
% add previous foreground points also
backcoord=[]
for i=1:size(trajidsfore,2)
   coordinates = TrajectoryCoordinates{trajidsfore(i)}(l*diffe-overlap1,:);
            backcoord(jj,:)=coordinates;
            jj=jj+1;
end
kkkk=1;
for kkkk=1:size(backcoord,1)
  plot(backcoord(kkkk,1),backcoord(kkkk,2),'b*');
end
hold off
% pause(0.5)

for i=1:nnz(back==1)
%    i
trajid = trajidsback(i);
validtrajid = find(validTraj{l}(3,:,1)==trajid);
TrajCoord{l}(:,:,i) = permute(reshape(validTraj{l}(1:2,validtrajid,:),2,size(validTraj{l},3)),[2 1]);
end
for i=1:nnz(back==1)
Y(:,:,i)=TrajCoord{l}(:,:,i);%reshape(,size(OTDM,1),2);
Ypre(:,:,i) = pre_shape(Y(:,:,i));
end

% hold off;
[mean_shape,tran]= mean_shape1(Y);%permute(ppp,[3 1 2]));

%% get representative trajectory of camera motion
[mean_shape,tran]= mean_shape1(Y);%permute(ppp,[3 1 2]));
% [mean_shape pc_shape std_shape pc_projection new_pt_out c scale]=tangent_pca_shape(permute(p,[3 1 2]));
figure;
plot3([1:size(mean_shape,1)],mean_shape(:,1),mean_shape(:,2),'r','LineWidth',4)
hold on;
x=[1:size(mean_shape,1)]';
y=mean_shape(:,1);
sm_mean_shape(:,1)= smooth(x,y,0.6,'rloess');
y=mean_shape(:,2);
sm_mean_shape(:,2)= smooth(x,y,0.6,'rloess');
plot3([1:size(mean_shape,1)],sm_mean_shape(:,1),sm_mean_shape(:,2),'r','LineWidth',4);
legend('frechet mean','smoothed frechet mean')
X=mean_shape;
Ypretran=zeros(size(Ypre,1),size(Ypre,2));
for j=1:nnz(back==1)
    Ypretran = Ypretran+Ypre(:,:,j)*tran(:,:,j);
end
norm(mean_shape*nnz(back==1)-Ypretran)
hold off;
clear trajids;clear Stab;clear ppp;
for i=1:size(TrajCoord{l},3)
%     Fulltrajectories(:,:,i)=TrajectoryCoordinates{i};
Trajectory = TrajCoord{l}(:,:,i);
ind = find(Trajectory(1:size(mean_shape,1),1)~=0);
[Y c scale]=pre_shape(Trajectory(ind,:));
A=Y'*mean_shape;
  [L, D, M] = svd(A);
T = M * L';
pp = mean_shape(ind,:)* T;
p = sm_mean_shape(ind,:)* T;%tran(:,:,i);%+tr.c;
d(i)= norm(p-sm_mean_shape(ind,:));
% come back from preshape
  kkk= p*scale+repmat(c,size(p,1),1);
  ppp(:,:,i)= pp*scale+repmat(c,size(p,1),1);
  Stab(ind,:,i)=kkk;
  hold on 

end
clear pertraj;
pertraj = permute(Stab,[3 1 2]);
trajids=trajidsback;
% close all;

end
Traj{l}=TrajCoord{l};
mean_shape_final{l}=mean_shape;
sm_mean_shape_final{l}=sm_mean_shape;
StabOTDM{l}=Stab;
clear mean_shapel;clear sm_mean_shape;
end
clear TrajCoord;
TrajCoord=Traj;
clear Traj
clear mean_shape;clear sm_mean_shape;clear x;clear y;
% mean_shape=zeros(size(validTraj,2)*diffe,2);
l=1;
    mean_shape(1:diffe,1)=mean_shape_final{l}(1:end,1);
    mean_shape(1:l*diffe,2)=mean_shape_final{l}(:,2);
    sm_mean_shape(1:l*diffe,1)=sm_mean_shape_final{l}(:,1);
    sm_mean_shape(1:l*diffe,2)=sm_mean_shape_final{l}(:,2);
  difference=zeros(size(validTraj,2)-1,2);
for l=2:size(validTraj,2)
    start =(l-1)*(diffe)-(l-1)*overlap1+overlap1
    difference(l,:)= mean_shape(end,:)-mean_shape_final{l}(1,:)
    mean_shape(start+1:start+diffe-overlap1,1)=mean_shape_final{l}(overlap1+1:end,1)+difference(l,1);
    mean_shape(start+1:start+diffe-overlap1,2)=mean_shape_final{l}(overlap1+1:end,2)+difference(l,2);
    sm_mean_shape(start+1:start+diffe-overlap1,1)=sm_mean_shape_final{l}(overlap1+1:end,1);
    sm_mean_shape(start+1:start+diffe-overlap1,2)=sm_mean_shape_final{l}(overlap1+1:end,2);
end

figure;
plot3([1:size(mean_shape,1)],mean_shape(:,1),mean_shape(:,2),'r','LineWidth',4)
hold on;
clear sm_mean_shape;
x=[1:size(mean_shape,1)]';
y=mean_shape(:,1);
sm_mean_shape(:,1)= smooth(x,y,0.8,'rloess');
y=mean_shape(:,2);
sm_mean_shape(:,2)= smooth(x,y,0.8,'rloess');
plot3([1:size(mean_shape,1)],sm_mean_shape(:,1),sm_mean_shape(:,2),'g','LineWidth',4);
hold off
X=mean_shape;
for l=1:size(validTraj,2)-1
    start= (l-1)*(diffe-overlap1)
ind = start+1:start+diffe;%find(Trajectory(1:diffe,1)~=0);
for i=1:size(TrajCoord{l},3)
Trajectory = TrajCoord{l}(:,:,i);
[Y c scale]=pre_shape(Trajectory);
diffmatt = mean_shape(ind,:)-Y;
A=Y'*(mean_shape(ind,:));
  [L, D, M] = svd(A);
T = M * L';
pp = (mean_shape(ind,:))* T;
p = (sm_mean_shape(ind,:))*T;%-diffmatt;%tran(:,:,i);%+tr.c;

% come back from preshape
  ppp(:,:,i)= pp*scale+repmat(c,size(p,1),1);
  diffmatt = ppp(:,:,i)-Trajectory;
    kkk= p*scale+repmat(c,size(p,1),1);

  StabOTDM{l}(:,:,i)=kkk;
%   plot3([1:diffe],Trajectory(:,1),Trajectory(:,2),'r','LineWidth',1);
%   hold on
%   plot3([1:diffe],kkk(:,1),kkk(:,2),'b','LineWidth',1);
% %   

end
end
hold off;

if ~exist([folder 'warped/'],'dir')
    mkdir([folder 'warped/']);
end
close all;
img_dir = dir([folder 'frames/']);
%% CPW
clear I1warp;clear set;
i=1;l=1;