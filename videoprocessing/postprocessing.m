function []=postprocessing(result,imgFolder_final,number)

addpath('N:\research\codes\Toolbox\','N:\research\codes\MY CODES\PerMIn\with_steadyflow\synthetic\result_final\test\');
% % % % img_dir=dir([imgFolder '*.png']);
% % % % cell_array=struct2cell(img_dir);
% % % % numFrames=size(cell_array,2);
% % % % C=cell(1,size(cell_array,2));
% % % % for i=1:numFrames
% % % % C{i}=cell_array{1,i};
% % % % end
% % % % [cs,index]=sort_nat(C);
% % % % numFrames=size(cs,2);
 im1=result{1}; %%imread([imgFolder cs{1}]);
numFrames=size(result,2);
% numFrames=121;
 new_col=ones(numFrames,1);  
new_row=ones(numFrames,1);  
cut_col=zeros(numFrames,size(im1,2));
cut_row=zeros(numFrames,size(im1,1));
%imwrite(im1,[imgFolder_final 'warp00',num2str(1),'.png']);
cut_final_max(1)=1;
cut_final_min(1)=size(im1,2);
cut_final_max_row(1)=1;
cut_final_min_row(1)=size(im1,1);

for i=1:numFrames-1
% % % %     %% read the current frame
% % % %     %% check if any black col extra to prev frame
% % % %     %% if yes, shift one col left,get one col from right previous frame  
i
%% chop off black part
  im2=result{i};%%imread([imgFolder cs{i}]);
%% every col
new_size=size(im2,2);
new_size_row=size(im2,1);

min_col = round(0.40*size(im2,2));
max_col = round(0.80*size(im2,2));
min_row = round(0.40*size(im2,1));
max_row = round(0.80*size(im2,1));
new_im=[];
new_im1=[];
close all;
for col=1:size(im2,2)
    array_of_col=im2(:,col,1);
    if nnz(array_of_col)<0.9*size(array_of_col,1)
       new_size=new_size-1;
       cut_col(i,col)=1;
    else
    new_im(:,new_col(i),1)=im2(:,col,1);
    new_im(:,new_col(i),2)=im2(:,col,2);
    new_im(:,new_col(i),3)=im2(:,col,3);
    new_col(i)=new_col(i)+1;
    end
end
if isempty(find(cut_col(i,1:min_col)))
cut_final_max(i)=1;
else
cut_final_max(i) = max(find(cut_col(i,1:min_col)));
end
if isempty(find(cut_col(i,max_col:size(im2,2))))
    cut_final_min(i)=size(im2,2);
else
cut_final_min(i) = max_col-1+min(find(cut_col(i,max_col:size(im2,2))));
end

for row=1:size(im2,1)
    array_of_row=im2(row,:,1);
    if nnz(array_of_row)<0.9*size(array_of_row,2)
       new_size_row=new_size_row-1;
       cut_row(i,row)=1;
    else
    new_im1(new_row(i),:,1)=new_im(row,:,1);
    new_im1(new_row(i),:,2)=new_im(row,:,2);
    new_im1(new_row(i),:,3)=new_im(row,:,3);
    new_row(i)=new_row(i)+1;
    end
    
end
if isempty(find(cut_row(i,1:min_row)))
cut_final_max_row(i)=1;
else
cut_final_max_row(i) = max(find(cut_row(i,1:min_row)));
end
if isempty(find(cut_row(i,max_row:size(im2,1))))
cut_final_min_row(i)=size(im2,1);
else
cut_final_min_row(i) = max_row-1+min(find(cut_row(i,max_row:size(im2,1))));
end
%imwrite(uint8(new_im1),[imgFolder_final 'final00',num2str(i),'.png']);
end
col_min=max(cut_final_max);
col_max=min(cut_final_min);

row_min=max(cut_final_max_row);
row_max=min(cut_final_min_row);
for i=1:numFrames
im1=result{i};
% imshow
new_im1=[];
new_im1=im1(row_min+1:row_max-1,col_min+1:col_max-1,:);
imwrite(uint8(new_im1),[imgFolder_final 'warp00',num2str(number),'.png']);
%  imshow(uint8(new_im1));
%  pause(0.3);
number=number+1;
end
end