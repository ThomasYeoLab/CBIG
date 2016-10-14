RandPhi = 2 * pi * rand(100000,1);
RandTheta = pi * rand(100000,1);
RandPts = [cos(RandPhi), sin(RandPhi).*cos(RandTheta), sin(RandPhi).*sin(RandTheta)];

tic
[tmp, tmp, TreeRoot] = kdtree(RandPts, []);

TestPhi = 2*pi*rand(100,1);
TestTheta = pi*rand(100,1);
TestPts = [cos(TestPhi), sin(TestPhi).*cos(TestTheta), sin(TestPhi).*sin(TestTheta)];


[ ClosestPts, DistA, TreeRoot ] = kdtree([], TestPts, TreeRoot);

toc

DistB = zeros(size(TestPts, 1),1);
for i = 1:size(TestPts, 1)
    
    D = sqrt(sum((RandPts - repmat(TestPts(i,:), [size(RandPts, 1) 1])).^2,2));
    
    DistB(i,1) = min(D);
    
end

DistA - DistB
    