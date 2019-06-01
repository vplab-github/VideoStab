%% get valid trajectories with non homogeneous pixels and throughout the frames
% read continuous trajectories every 20 frames
% get the angles for each frame
function [] = fullTraj_adaptive(Trajname,in,diffe,numFrames,validTrajname,dim)
% for each set plot the 
%----------------------------------
rand('state',0);
overlap(1)=0;
in
diffe
% flag=0;
%----------------------------Set path
% addpath('E:\research\codes\MY CODES\Work2\Pixel_profile_plot\plot_points\')
load(Trajname);load(validTrajname);
k=0;
frames = zeros(numFrames,size(TrajectoryCoordinates,1));
% Organize trajectories based on frames
for i = 1: size(TrajectoryCoordinates,1)
    [row,~] = find(TrajectoryCoordinates{i}~=0);
    uniqrow = unique(row);
    frames(uniqrow,i)=1;
%     k = k+size(unique(row),1);
    
end
in=1;diffe=numFrames+3;k=1;
while in<numFrames
    in
% columns = 1:size(TrajectoryCoordinates,1);
% subframes = frames(in:min(in+diffe-1,numFrames),:);
%  [r,cols] = find(subframes==0);
%   
%     uniqcols = unique(cols);
%      ind = find(ismember(columns,uniqcols));
%     col=columns(setdiff(1:length(columns),ind));
%% get the total no of trajectory points in that frame
    % find the number of trajectories in the subframes
  flag=0;
  
while flag==0 & diffe-1-in>0
% k=size(validTraj,2)+1;
diffe=diffe-1;

columns = 1:size(TrajectoryCoordinates,1);
subframes = frames(in:min(diffe,numFrames),:);
 [r,cols] = find(subframes==0);
  
    uniqcols = unique(cols);
     ind = find(ismember(columns,uniqcols));
    col=columns(setdiff(1:length(columns),ind));
    sub=zeros(1,size(subframes,1));
% for i=1:size(subframes,1)
%     noOfTrajInFrame=nnz(subframes(i,:));
%     sub(i) = nnz(col)>3000;
    if nnz(col)>0.1*size(TrajectoryCoordinates,1)
     flag=1;
      break;
    end
% end
 
end
%         fullTraj_adaptive(Trajname,in,diffe,numFrames,validTrajname,overlap,dim);
%         fullTraj_adaptive(Trajname,in,diffe,numFrames,validTrajname,overlap,dim);
    
%     flag=1;
    ind = find(ismember(columns,uniqcols));
%     col=columns(setdiff(1:length(columns),ind));
    l = col;
    trajids{k}=l;
    excluded_ids{k} = uniqcols;
    if in==275
        disp('hai');
    end
    for j=1:size(l,2)
       
        traj(:,1:2) = TrajectoryCoordinates{l(j)};
        traj(:,3) = repmat(l(j),[size(traj,1) 1]);
      kk = traj(in:min(diffe,numFrames),:);
        validTraj{k}(:,j,:) = permute(kk,[2 3 1]);

    end
  
    k=k+1;
    
    overlap(k) = max(round((diffe-in+1)*0.3),1);
   
  
    
    in
    diffe
     in = diffe-overlap(k)+1;
    diffe=numFrames+1;
 end
    

% save the valid trajectory coordinates and trajectory ID in validTraj and the span of the frames  
 save(validTrajname,'validTraj','trajids','excluded_ids','overlap');
 clear all;
end