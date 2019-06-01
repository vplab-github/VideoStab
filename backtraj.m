addpath(genpath('./'));
clear all;f='D';set = 'frog';
folder = ['./data/' set '/' f '/'];
fileLocation = [folder f 'Tracks773.dat'];
InputVideoName1='1.avi';thresholdtraj=0.95;iterations=20;lambda1=0.5;
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
diffe=10;i=1;
  fullTraj1([folder readTraj],200,noOfFrames,validTrajname,1);
% fullTraj_adaptive([folder readTraj],i,diffe,noOfFrames,validTrajname);
load(validTrajname);
%% get all the trajectories
trajectories = validTraj{1}(1:2,:,:);trajids = validTraj{1}(3,:,1);
pertraj = permute(trajectories,[2 3 1]);
% [numFrames,frameRate,frames]=ReadVideoAVI(InputVideoName1,folder);
jj=1;
for indexbak=1:3
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
       
   test(i) = norm(accumtraj(i,:));
end
%% average of accumulated traj
averageAcc = 1/size(accumtraj,1) * sum(accumtraj,1);
difference = sqrt(sum((accumtraj-repmat(averageAcc,size(accumtraj,1),1)).^2,2));
[histo,centres]=hist(difference,10)
an=[histo(1) sum(histo(1:2)) sum(histo(1:3)) sum(histo(1:4)) sum(histo(1:4)) sum(histo(1:5)) sum(histo(1:6)) sum(histo(1:7)) sum(histo(1:8)) sum(histo(1:9)) sum(histo(1:10))]<thresholdtraj*size(difference,1);
thresholdinitial = centres(find(an==1,1,'last'));
% earlierback=back;
back=difference<thresholdinitial;
trajidsback = trajids(back==1);
 img=imread([folder 'frames/0001.ppm']);
%  img=frames(:,:,:,1);
   imshow(img);
hold on;
for i=1:size(trajidsback,2)
   coordinates = TrajectoryCoordinates{trajidsback(i)}(1,:);
      plot(coordinates(1),coordinates(2),'r*');

end
trajidsfore = trajids(back==0);
% add previous foreground points also
for i=1:size(trajidsfore,2)
   coordinates = TrajectoryCoordinates{trajidsfore(i)}(1,:);
%       plot(coordinates(1),coordinates(2),'b*');
            backcoord(jj,:)=coordinates;
            jj=jj+1;
end
kkkk=1;
for kkkk=1:size(backcoord,1)
  plot(backcoord(kkkk,1),backcoord(kkkk,2),'b*');
end
hold off
pause(0.5)
l=1;
for i=1:nnz(back==1)
%    i
trajid = trajidsback(i);
validtrajid = find(validTraj{l}(3,:,1)==trajid);
TrajCoord(:,:,i) = permute(reshape(validTraj{l}(1:2,validtrajid,:),2,size(validTraj{l},3)),[2 1]);
end
for i=1:nnz(back==1)
Y(:,:,i)=TrajCoord(:,:,i);%reshape(,size(OTDM,1),2);
Ypre(:,:,i) = pre_shape(Y(:,:,i));
%  plot(Y(:,1,i),Y(:,2,i),'LineWidth',1);hold on
% plot(Ypre(:,1,i),Ypre(:,2,i),'LineWidth',1);hold on
% plot3([1:size(Y,1)],Y(:,1,i),Y(:,2,i),'LineWidth',1);hold on
% plot3([1:size(Y,1)],Ypre(:,1,i),Ypre(:,2,i),'LineWidth',1);hold on
end

hold off;
%% get representative trajectory of camera motion
[mean_shape,tran]= mean_shape1(Y);%permute(ppp,[3 1 2]));
% [mean_shape pc_shape std_shape pc_projection new_pt_out c scale]=tangent_pca_shape(permute(p,[3 1 2]));
figure;
plot3([1:size(mean_shape,1)],mean_shape(:,1),mean_shape(:,2),'r','LineWidth',4)
hold on;
x=[1:size(mean_shape,1)]';
y=mean_shape(:,1);
sm_mean_shape(:,1)= smooth(x,y,0.8,'rloess');
y=mean_shape(:,2);
sm_mean_shape(:,2)= smooth(x,y,0.8,'rloess');
plot3([1:size(mean_shape,1)],sm_mean_shape(:,1),sm_mean_shape(:,2),'g','LineWidth',4);
legend('frechet mean','smoothed frechet mean')
X=mean_shape;
Ypretran=zeros(size(Ypre,1),size(Ypre,2));
for j=1:nnz(back==1)
    Ypretran = Ypretran+Ypre(:,:,j)*tran(:,:,j);
end
norm(mean_shape*nnz(back==1)-Ypretran)
% % % [mean_shape,tran]= mean_shape1(Y);%permute(ppp,[3 1 2]));
% % % % [mean_shape pc_shape std_shape pc_projection new_pt_out c scale]=tangent_pca_shape(permute(p,[3 1 2]));
% % % figure;
% % % plot3([1:size(mean_shape,1)],mean_shape(:,1),mean_shape(:,2),'r','LineWidth',4)
% % % hold on;
% % % x=[1:size(mean_shape,1)]';
% % % y=mean_shape(:,1);
% % % % get optimization function
% % % % y1 = mean_shape(:,2);
% % % % sm_mean_shape(:,1) = y;
% % % % sm_mean_shape(:,2) = y1;
% % % % for k=1:iterations
% % % %    sm_mean_shape(:,1) = CalcOneUpdate(y,sm_mean_shape(:,1),100,lambda1);
% % % %    sm_mean_shape(:,2) = CalcOneUpdate(y1,sm_mean_shape(:,2),100,lambda1);
% % % % %    optPathy(i,:) = CalcOneUpdate(M{l}(size(mag(i,:),2)+1:2*size(mag(i,:),2),i),optPathy(i,:),diffe,lambda1);
% % % % 
% % % % end
% % % sm_mean_shape(:,1)= smooth(x,y,0.7,'rloess');
% % % y=mean_shape(:,2);
% % % sm_mean_shape(:,2)= smooth(x,y,0.7,'rloess');
% % % plot3([1:size(mean_shape,1)],sm_mean_shape(:,1),sm_mean_shape(:,2),'g','LineWidth',4);
% % % legend('frechet mean','smoothed frechet mean')
% % % X=mean_shape;
hold off;clear trajids;clear StabOTDM;clear ppp;
for i=1:size(TrajCoord,3)
%     Fulltrajectories(:,:,i)=TrajectoryCoordinates{i};
Trajectory = TrajCoord(:,:,i);
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
  StabOTDM(ind,:,i)=kkk;
  plot3([1:size(mean_shape,1)],Trajectory(:,1),Trajectory(:,2),'r','LineWidth',1);
  hold on
  plot3([1:size(mean_shape,1)],kkk(:,1),kkk(:,2),'b','LineWidth',1);
  

end
clear pertraj;
pertraj = permute(StabOTDM,[3 1 2]);
trajids=trajidsback;
close all;
end
if ~exist([folder 'warped/'],'dir')
    mkdir([folder 'warped/']);
end
close all;
img_dir = dir([folder 'frames/']);
%% CPW
clear I1warp;clear set;
i=1;
for i=1:size(mean_shape,1)
%% get the trajectory ids in frame number i
i
   originaltraj = permute(reshape(TrajCoord(i,:,:),size(TrajCoord,2),size(TrajCoord,3)), [2 1]);
%     originaltraj = permute(reshape(ppp(i,:,:),size(ppp,2),size(ppp,3)), [2 1]);
   smoothedtraj = permute(reshape(StabOTDM(i,:,:),size(TrajCoord,2),size(TrajCoord,3)), [2 1]);
   X1 = round(originaltraj(:,1));
Y1 = round(originaltraj(:,2));
smoothX = round(smoothedtraj(:,1));
smoothX(smoothX<0)=1;
smoothY = round(smoothedtraj(:,2));
smooth(smoothY<0)=1;

% % % %% asap warping
frames(:,:,:,i)=imread([folder 'frames/' img_dir(i+2).name]);
% figure;  imshow(frames(:,:,:,i));hold on;
% plot(X1,Y1,'r*');
% plot(smoothX,smoothY,'b*');
%pause(0.5);
% close all;
% % I1 = imresize(img_rgb,[360 640]);
I1 = frames(:,:,:,i);%imresize(frames(:,:,:,i),[360 640]);
clear I1_features;clear I2_features;
I1_features(:,1) = X1;
I1_features(:,2) = Y1;
I2_features(:,1) = smoothX;
I2_features(:,2) = smoothY;
I1=frames(:,:,:,i);
imshow(I1);hold on
plot(I1_features(:,1),I1_features(:,2),'r*');
plot(I2_features(:,1),I2_features(:,2),'b*');
% pause(0.5)
if length(I1_features) < 20
error('not enough matched features');
return;
end
[height,width,~] = size(I1);
%3x3 mesh
quadWidth = width/(2^3);
quadHeight = height/(2^3);
% %4x4 mesh
% quadWidth = width/(2^4);
% quadHeight = height/(2^4);
lamda = 1; %mesh more rigid if larger value. [0.2~5]
asap = AsSimilarAsPossibleWarping(height,width,quadWidth,quadHeight,lamda);
asap.SetControlPts(I1_features,I2_features);%set matched features
asap.Solve();            %solve Ax=b for as similar as possible
homos = asap.CalcHomos();% calc local hommograph transform
gap = 100;
I1warp = asap.Warp(I1,gap);                     %warp source image to target image
I1warpmesh = asap.destin.drawMesh(I1warp,gap);  %draw mesh on the warped source image
% imshow(I1warpmesh);
gap = 100;
I1warp = asap.Warp(I1,gap);
h=imshow(I1warp);

% 
% h=imshow(I1);
set(h,'AlphaData',0.5); 
pause(0.1)
I1warp1{i}=I1warp;
                 
% imwrite(I1warp1{i},[folder 'warped/warp' num2str(i) '.png']);
% imshow(I1warp1{i});%pause(0.5);
%access local homography
[h,w,~,~] = size(homos);
for k=1:h-1
    for j=1:w-1
       H(:,:,i) = homos(k,j,:,:);
%        fprintf('Quad=[%d %d]\n',i,j);
%        H(:,:,k)
    end
end
% prevX1=smoothX;prevY1=smoothY;
% close all;
end
% cc = postprocessing(I1warp1);
s1=130;s2=130;
clear I2;clear ComparisionIM;
for i=1:size(mean_shape,1)
%     cutI1warp{i} = I1warp1{i}(20:size(I1warp1{i},1)-20,20:size(I1warp1{i},2)-20,:);
   
I2(:,:,:,i)=I1warp1{i}(s1:size(I1warp1{i},1)-s1,s2:size(I1warp1{i},2)-s2,:);
ki=sprintf('%04d',i);
%  imwrite(I2(:,:,:,i),[folder 'result/warp' ki '.png']);
resI2(:,:,:,i) = imresize(double(I2(:,:,:,i)),[size(frames,1) size(frames,2)]);
%  ComparisionIM(:,:,:,i) = video_horizontal(frames(:,:,:,i),resI2,20);
end

result = double(I2);
result = postprocessing(result);
save([folder 'resultFrames.mat'],'result')
for i=1:size(mean_shape,1)
    ComparisionIM(:,:,:,i) = video_horizontal(frames(:,:,:,i),imresize(result(:,:,:,i),[size(frames,1) size(frames,2)]),50);
end
WriteVideoAVI('output.avi',folder,0.5,result);
WriteVideoAVI('comparison.avi',folder,25,ComparisionIM);

disp('[DONE]');


