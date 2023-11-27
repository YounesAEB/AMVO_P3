function [coordsP,coordsC,deltaY,c,c12,theta,aoaE] = computeGeometryCosine(N,b,cR,cT,thetaT,aoa)

coordsP      = zeros(N+1,3);
coordsC      = zeros(N,3);
deltaY  = zeros(N,1);
beta    = linspace(0,180,N+1)';
c12     = zeros(N,1);

% Spanwise coordinates - Cosine Distribution
for i=1:N
    if i == 1
        coordsP(i,2)     =  -cosd(beta(i))*b/2;
        coordsP(i+1,2)   = -cosd(beta(i+1))*b/2;
        coordsC(i,2)     = (coordsP(i+1,2) + coordsP(i,2))/2;
        deltaY(i) = coordsP(i+1,2) - coordsP(i,2);
    else
        coordsP(i+1,2)   = -cosd(beta(i+1))*b/2;
        coordsC(i,2)     = (coordsP(i+1,2) + coordsP(i,2))/2;
        deltaY(i) = coordsP(i+1,2) - coordsP(i,2);
    end
end
% Torsion and Efective Angle of Attack
if thetaT == 0
    theta = zeros(1,N+1);
else
    theta = [thetaT:-thetaT/(N/2):0,thetaT/(N/2):thetaT/(N/2):thetaT];
end
aoaE  = theta + aoa; % Efective angle of attack

% Chord length of each slice -  NOT CORRECT
c     = [cT:(cR-cT)/(N/2):cR,cR-(cR-cT)/(N/2):-(cR-cT)/(N/2):cT];
for i = 1:N
    c12(i) = (c(i)+c(i+1))/2;
end

end