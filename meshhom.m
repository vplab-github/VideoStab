function [homos]= meshhom(I1,I2)
% fprintf('detect surf features...');
    grayI1 = rgb2gray(I1);
    grayI2 = rgb2gray(I2);

      points1 = detectSURFFeatures(grayI1);
    points2 = detectSURFFeatures(grayI2);

    [f1, vpts1] = extractFeatures(grayI1, points1);
    [f2, vpts2] = extractFeatures(grayI2, points2);

    index_pairs = matchFeatures(f1, f2) ;
    matched_pts1 = vpts1(index_pairs(:, 1));
    matched_pts2 = vpts2(index_pairs(:, 2));

    [n,~] = size(matched_pts1);

    I1_features = zeros(n,2);
    I2_features = zeros(n,2);

    for i=1:n
        I1_features(i,:) = matched_pts1(i).Location;
        I2_features(i,:) = matched_pts2(i).Location;
    end
%  plot(I1_features(:,1),I1_features(:,2),'r*');
% plot(I2_features(:,1),I2_features(:,2),'b*');
% fprintf('[DONE]');

% if length(I1_features) < 20
%     error('not enough matched features');
%     return;
% end

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
end