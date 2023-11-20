function [yP,yC,deltaY] = computeGeometryCosine(N,b)

yP      = zeros(N+1,1);
yC      = zeros(N,1);
deltaY  = zeros(N,1);
beta    = linspace(0,180,N+1)';

for i=1:N
    if i == 1
        yP(i) =  -cosd(beta(i))*b/2;
        yP(i+1)   = -cosd(beta(i+1))*b/2;
        yC(i)     = (yP(i+1) + yP(i))/2;
        deltaY(i)   = yP(i+1) - yP(i);
    else
        yP(i+1)   = -cosd(beta(i+1))*b/2;
        yC(i)     = (yP(i+1) + yP(i))/2;
        deltaY(i)   = yP(i+1) - yP(i);
    end
end