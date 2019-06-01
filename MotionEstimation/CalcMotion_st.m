function FList = CalcMotion_st(input_frames)

[h,w,c,numFrames] = size(input_frames(:,:,:,:));
% numFrames=99;
FList = zeros(3,3,numFrames);
FList(:,:,1) = eye(3);

    for i=2:numFrames
         s = sprintf('compute homography %d / %d \n',i-1,numFrames);
         s
     FF = meshhom(input_frames(:,:,:,i),input_frames(:,:,:,i-1));
        FList(:,:,i) = FF(2,2,:,:);%EstimateHomography_SURF_RANSAC(input_frames(:,:,:,i),input_frames(:,:,:,i-1));
        FList(:,:,i)
    end
%     save(strcat(name,'\homography'),'FList');
end