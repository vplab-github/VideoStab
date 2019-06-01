function resultPath = CalcOneSmoothing2(camPath,optPath,numFrames,lambda1,lambda2,optPathprev,overlap)
numFrames = size(optPath,2);
resultPath = zeros(1,numFrames);
L(1:overlap)=1;
L(overlap+1:numFrames)=0;
for i=1:numFrames
    den = numFrames^2*(overlap+lambda2*L(i))+lambda1*overlap*(numFrames-1);
    resultPath(i) = camPath(i)*((numFrames^2*overlap)/den);
    
    for j=1:i-1
        resultPath(i) = resultPath(i)+((lambda1*overlap)/den)*optPath(j);
    end
    for j=i+1:numFrames
        resultPath(i) = resultPath(i)+((lambda1*overlap)/den)*optPath(j);
    end
    if L(i)==1 & numFrames-overlap+i<=size(optPathprev,2)
        resultPath(i) = resultPath(i)+((lambda2*L(i)*numFrames^2)/den)*optPathprev(size(optPathprev,2)-overlap+i);
    end
end

end

