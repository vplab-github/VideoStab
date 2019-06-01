%% get the trajectories belonging to the background
addpath(genpath('./'));
f='D';set = 'frog';
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
  fullTraj1([folder readTraj],100,noOfFrames,validTrajname,1);
% fullTraj_adaptive([folder readTraj],i,diffe,noOfFrames,validTrajname);
load(validTrajname);
%% get all the trajectories
trajectories = validTraj{1}(1:2,:,:);trajids = validTraj{1}(3,:,1);
pertraj = permute(trajectories,[2 3 1]);
%% accumulated motion vector: get the histogram separate and get the threshold
%  deltraj=zeros(size(pertraj));
%  accumtraj=zeros(size(pertraj,1),2);
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
 img=imread([folder 'frames/0001.ppm']);
   imshow(img);
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
l=1;
for i=1:nnz(back==1)
%    i
trajid = trajidsback(i);
validtrajid = find(validTraj{l}(3,:,1)==trajid);
TrajCoord(:,:,i) = permute(reshape(validTraj{l}(1:2,validtrajid,:),2,size(validTraj{l},3)),[2 1]);
        tempx=TrajCoord(:,1,i);
        tempy=TrajCoord(:,2,i);
%           if nnz(diff(tempx)==0)==0
        %% ignore trajectories stagnant
        in = find(diff(tempx)==0);
        gradtx=diff(tempx);
        gradtx(in) = 0.0001;
        clear in;
         in = find(diff(tempy)==0);
        gradty=diff(tempy);
        gradty(in) = 0.0001;
%         tempx(in)=[];tempy(in)=[];
        x1(i,:) = zeros(1,size(tempx,1));
        y1(i,:) = zeros(1,size(tempx,1));
        x1(i,2:end) = gradtx;%diff(tempx);
        y1(i,2:end) = gradty;%diff(tempy); 
         if nnz(gradtx)~=0
        x1(i,1) = x1(i,2);
       end
       if nnz(gradtx)~=0
        y1(i,1) = y1(i,2);
       end
        mag(i,:) = sqrt(x1(i,:).^2+y1(i,:).^2)';
    cosi(i,:) = (x1(i,:))./mag(i,:);
    sini(i,:) = (y1(i,:))./mag(i,:);
     OTDM(:,1,i) = mag(i,:).*(sini(i,:));%x1(k,:)+300;%
     OTDM(:,2,i) = mag(i,:).*(cosi(i,:));%y1(k,:)+300;
end

for i=1:10:nnz(back==1)
Y(:,:,i)=reshape(TrajCoord(:,:,i),size(OTDM,1),2);
Ypre(:,:,i) = pre_shape(Y(:,:,i));
%  plot(Y(:,1,i),Y(:,2,i),'b','LineWidth',1);hold on
% plot(Ypre(:,1,i),Ypre(:,2,i),'LineWidth',1);hold on
plot3(Y(:,1,i),[1:size(Y,1)],Y(:,2,i),'b','LineWidth',1);hold on
% plot3([1:size(Y,1)],Ypre(:,1,i),Ypre(:,2,i),'b','LineWidth',1);hold on
end

hold off;
%% get representative trajectory of camera motion
[mean_shape,tran]= mean_shape1(Y);%permute(ppp,[3 1 2]));
% [mean_shape pc_shape std_shape pc_projection new_pt_out c scale]=tangent_pca_shape(permute(p,[3 1 2]));
figure;
plot3([1:size(mean_shape,1)],mean_shape(:,1),mean_shape(:,2),'r','LineWidth',2)
hold on;
x=[1:size(mean_shape,1)]';
y=mean_shape(:,1);
sm_mean_shape(:,1)= smooth(x,y,0.5,'rloess');
y=mean_shape(:,2);
sm_mean_shape(:,2)= smooth(x,y,0.5,'rloess');
plot3([1:size(mean_shape,1)],sm_mean_shape(:,1),sm_mean_shape(:,2),'g','LineWidth',3);
legend('frechet mean','smoothed frechet mean')
X=mean_shape;
optPathx=mean_shape(:,1);
optPathy=mean_shape(:,2);lambda1=0.5;
% for k=1:iterations
%    optPathx = CalcOneUpdate(mean_shape(:,1),optPathx,diffe,lambda1);
%    optPathy = CalcOneUpdate(mean_shape(:,2),optPathy,diffe,lambda1);
% end
% % % for i=1:size(OTDM,3)
% % %     i
% % % %% get the procrustes transformation from the mean to the traj
% % % % Y=Ypre(:,:,i);
% % % [d Z tr] = procrustes(Y(:,:,i),X);
% % % % Ysmooth(:,:,i) = tr.b * sm_mean_shape * tr.T + tr.c;
% % % % % figure;hold on
% % % % % plot3(X(:,1),X(:,2),[1:size(mean_shape,1)],'LineWidth',1);
% % % % % plot3(sm_mean_shape(:,1),sm_mean_shape(:,2),[1:size(mean_shape,1)],'LineWidth',1);
% % % plot3(Z(:,1),Z(:,2),[1:size(mean_shape,1)],'LineWidth',1);hold on
% % % % % plot3(Ysmooth(:,1,i),Ysmooth(:,2,i),[1:size(mean_shape,1)],'LineWidth',1);
% % % % % 
% % % % % pause(0.5);
% % % % % close all;
% % % end

Trajectoryinframe=zeros(size(OTDM,3),size(mean_shape,1));
Ys= zeros(size(OTDM,1),2,size(OTDM,3));
figure;hold on

for i=1:size(OTDM,3)
%     Fulltrajectories(:,:,i)=TrajectoryCoordinates{i};
Trajectory = TrajCoord(:,:,i);
ind = find(Trajectory(1:size(mean_shape,1),1)~=0);
[Y c scale]=pre_shape(Trajectory(ind,:));
% X=mean_shape(ind,:);
% [d Z tr] = procrustes(Y,X);
% if d<0.01 & nnz(trajidsback==i)==0
    Trajectoryinframe(i,ind)=1;
p = sm_mean_shape(ind,:)* tran(:,:,i);%+tr.c;
% come back from preshape
  kkk= p*scale+repmat(c,size(p,1),1);
  StabOTDM(ind,:,i)=kkk;
%% in frame which all trajectories are available
% plot3([1:size(ind,1)],mean_shape(ind,1),mean_shape(ind,2),'g','LineWidth',1);
% plot3([1:size(ind,1)],sm_mean_shape(ind,1),sm_mean_shape(ind,2),'g','LineWidth',4);
% plot3(Trajectory(ind,1),Trajectory(ind,2),[1:size(Y,1)],'LineWidth',1);
plot3([1:size(Y,1)],p(:,1),p(:,2),'b','LineWidth',1); hold on
% pause(0.5); 
% end
% close all;
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
%    testcol = Trajectoryinframe(:,i);
%    trajids = find(testcol==1);
%    originaldel = permute(reshape(OTDM(i,:,:),2,size(OTDM,3)),[2 1]);
%    originaldel = [0 0;originaldel]
%    smootheddel = permute(reshape(StabOTDM(i,:,:),2,size(OTDM,3)),[2 1]); 
   originaltraj = permute(reshape(TrajCoord(i,:,:),size(OTDM,2),size(OTDM,3)), [2 1]);
   smoothedtraj = permute(reshape(StabOTDM(i,:,:),size(OTDM,2),size(OTDM,3)), [2 1]);
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
gap = 0;
I1warp = asap.Warp(I1,gap);                     %warp source image to target image
I1warpmesh = asap.destin.drawMesh(I1warp,gap);  %draw mesh on the warped source image
imshow(I1warpmesh);
gap = 0;
I1warp = asap.Warp(I1,gap);
g=imshow(I1warp);
hold on
h=imshow(I1);
set(h,'AlphaData',0.5); 
I1warp1{i}=I1warp;
                 
imwrite(I1warp1{i},[folder 'warped/warp' num2str(i) '.png']);
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
close all;
end
% cc = postprocessing(I1warp1);
s=1;
clear I2;clear ComparisionIM;
for i=1:size(mean_shape,1)
%     cutI1warp{i} = I1warp1{i}(20:size(I1warp1{i},1)-20,20:size(I1warp1{i},2)-20,:);
   
I2(:,:,:,i)=I1warp1{i}(s:size(I1warp1{i},1)-s,s:size(I1warp1{i},2)-s,:);
ki=sprintf('%04d',i);
%  imwrite(I2(:,:,:,i),[folder 'result/warp' ki '.png']);
resI2 = imresize(double(I2(:,:,:,i)),[size(frames,1) size(frames,2)]);
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


