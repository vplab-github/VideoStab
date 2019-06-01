function []=PlotFrameAndSift_mine(bbox,pts1,pts2,f_ind,imgData)
xImage = [bbox(1) bbox(1)+bbox(3);bbox(1) bbox(1)+bbox(3)];
yImage = f_ind*5*ones(2);
zImage = [bbox(2) bbox(2);bbox(2)+bbox(4) bbox(2)+bbox(4)];

%surf(xImage,yImage,zImage,...    %# Plot the surface
 %   'CData',flipdim(imgData,1),...
  %  'FaceColor','texturemap');
hold on

for i=1:size(pts1,1)
    if f_ind~=3
    line([pts1(i,2) pts2(i,2)],[yImage(1,1)-5 yImage(1,1)],[pts1(i,1)  pts2(i,1)]);
    else
            line([pts1(i,2) pts2(i,2)],[yImage(1,1)-15 yImage(1,1)],[pts1(i,1)  pts2(i,1)]);

    end
   % plot3(pts1(i,1),yImage(1,1)-10,pts1(i,2),'ro');
    %plot3(pts2(i,1),yImage(1,1),pts2(i,2),'ro');
end
end