function [yP, yC, deltaY,c,c12,theta,aoaE] = computeGeometryUniform(N,b,cR,cT,thetaT,aoa)


yP     = zeros(N+1,1);
yC     = zeros(N,1);
deltaY = zeros(N,1);
d      = b/N;
sum    = -b/2;
c12    = zeros(N,1);

for i=1:N
    if i == N
        yP(i,1)   = sum;
        yP(i+1,1) = sum + d;
        yC(i,1)   = yP(i,1) + d/2;
    else
        yP(i,1)   = sum;
        yC(i,1)   = yP(i,1) + d/2;
        sum       = sum + d;
    end
    deltaY (i)    = d;

end
theta = [thetaT:-thetaT/(N/2):0,thetaT/(N/2):thetaT/(N/2):thetaT];
aoaE  = theta + aoa; % Efective angle of attack
c     = [cT:(cR-cT)/(N/2):cR,cR-(cR-cT)/(N/2):-(cR-cT)/(N/2):cT];
for i = 1:N
    c12(i) = (c(i)+c(i+1))/2;
end
end

