function [coordsP,coordsC,deltaY,c,c12,theta,aoaE] = computeGeometryUniform(N,b,cR,cT,thetaT,aoa)


coordsP      = zeros(N+1,3);
coordsC      = zeros(N,3);
deltaY  = zeros(N,1);
d      = b/N;
sum    = -b/2;
c12    = zeros(N,1);

% Spanwise coordinates - Uniform Distribution
for i=1:N
    if i == N
        coordsP(i,2)   = sum;
        coordsP(i+1,2) = sum + d;
        coordsC(i,2)   = coordsP(i,2) + d/2;
    else
        coordsP(i,2)   = sum;
        coordsC(i,2)   = coordsP(i,2) + d/2;
        sum       = sum + d;
    end
    deltaY (i)    = d;
end

% Torsion and Efective Angle of Attack
if thetaT == 0
    theta = zeros(1,N+1);
else
    theta = [thetaT:-thetaT/(N/2):0,thetaT/(N/2):thetaT/(N/2):thetaT];
end
aoaE  = theta + aoa; % Efective angle of attack

% Chord length of each slice
c     = [cT:(cR-cT)/(N/2):cR,cR-(cR-cT)/(N/2):-(cR-cT)/(N/2):cT];
for i = 1:N
    c12(i) = (c(i)+c(i+1))/2;
end

end

