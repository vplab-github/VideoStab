%% get valid trajectories with non homogeneous pixels and throughout the frames
% read continuous trajectories every 20 frames
% get the angles for each frame
function [] = fullTraj1(Trajname,diff,numFrames,validTrajname,o)
% for each set plot the 
%----------------------------------
rand('state',0);
%----------------------------Set path
% addpath('E:\research\codes\MY CODES\Work2\Pixel_profile_plot\plot_points\')
load(Trajname);
k=0;
frames = zeros(numFrames,size(TrajectoryCoordinates,1));

% Organize trajectories based on frames
for i = 1: size(TrajectoryCoordinates,1)
    [row,~] = find(TrajectoryCoordinates{i}~=0);
    uniqrow = unique(row);
    frames(uniqrow,i)=1;
    k = k+size(unique(row),1);
    
end
k=0;
columns = 1:size(TrajectoryCoordinates,1);
% get the common trajectories in a set of frames
for i=1:diff-o:numFrames-mod(numFrames,diff)
    k=k+1;
%     
    subframes = frames(i:min(i+diff-1,numFrames),:);
   
    [r,cols] = find(subframes==0);
    uniqcols = unique(cols);
    ind = find(ismember(columns,uniqcols));
    col{k}=columns(setdiff(1:length(columns),ind));
    l = col{k};
    trajids{k}=l;
    excluded_ids{k} = uniqcols;
    for j=1:size(l,2)
       
        traj(:,1:2) = TrajectoryCoordinates{l(j)};
        traj(:,3) = repmat(l(j),[size(traj,1) 1]);  kk = traj(i:min(i+diff-1,numFrames),:);
        validTraj{k}(:,j,:) = permute(kk,[2 3 1]);

    
    end
    overlap(k)=o;
    %subframes(i:i+diff-1,col{k});
end
% save the valid trajectory coordinates and trajectory ID in validTraj and the span of the frames  
 save(validTrajname,'validTraj','trajids','excluded_ids','overlap');
 clear all;
end