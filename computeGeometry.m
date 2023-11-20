function [yP,yC] = computeGeometry(N,b)


yP = zeros(N+1,1);
yC = zeros(N,1);
d = b/N;
sum = -b/2;

for i=1:N
    if i == N
        yP(i,1) = sum;
        yP(i+1,1) = sum + d;
        yC(i,1) = yP(i,1) + d/2;
    else
        yP(i,1) = sum;
        yC(i,1) = yP(i,1) + d/2;
        sum = sum + d;
    end
end