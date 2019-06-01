function [kshape] = kendalls(backtraj)
for i=1:size(backtraj,2)
   kshape(:,:,i)= pre_shape(reshape(permute(backtraj(:,i,:),[3 1 2]),size(backtraj,3),size(backtraj,1)));
end
end