function resultPath = CalcOneSmoothing_wo2R(camPath,optPath,numFrames,lambda)
resultPath = zeros(1,size(optPath,2));
R=numFrames;
for i=1:numFrames
    alpha=(R)/(lambda*R-lambda+R);
    resultPath(i) = camPath(i)*(alpha);
    for j=1:i-1
        resultPath(i) = resultPath(i)+(lambda/(R)*alpha)*optPath(j);
    end
    for j=i+1:numFrames
        resultPath(i) = resultPath(i)+(lambda/(R)*alpha)*optPath(j);
    end
end
end

