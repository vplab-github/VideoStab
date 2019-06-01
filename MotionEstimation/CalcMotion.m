function FList = CalcMotion(name,input_frames,output_frames)

[h,w,c,numFrames] = size(input_frames(:,:,:,:));
% [h,w,c,numFrames] = size(output_frames(:,:,:,:));
FList = zeros(3,3,numFrames);
FList(:,:,1) = eye(3);

    fprintf('compute homography \n');
    for i=1:numFrames
         fprintf(' %3d / %d\n',i-1,numFrames);
         %fprintf('%03d\b\b\b',i-1);
%              FF = meshhom(output_frames(:,:,:,i),input_frames(:,:,:,i));
 FList(:,:,i) = EstimateHomography_SURF_RANSAC(output_frames(:,:,:,i),input_frames(:,:,:,i));
%         FList(:,:,i) = FF(2,2,:,:);%EstimateHomography_SURF_RANSAC(output_frames(:,:,:,i),input_frames(:,:,:,i));
    end
%     save(strcat(name,'\homography'),'FList');
end