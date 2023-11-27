%--------------------------------------------------------------------------
%  /  POTENTIAL AERODYNAMICS - LIFTING LINE METHOD -AMVO 
%  /  Matlab code to assess the numerical solution via LLM - Part 1                                            
%  /  ESEIAAT_UPC                                           
%  /  MUEA - MQ1 - Younes Akhazzan - Joel Rajo - Pol Ruiz                         
%--------------------------------------------------------------------------
clc; clear; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% Data given by the excercise statement
b       = 6;     % Wingspan of the main wing
bh      = 2.2;   % Wingspan of the horizontal tail plane (HTP)
ba      = 2;     % Airleron Width of the semi-wing
cR      = 1.3;   % Root chord of the main wing
cT      = 0.7;   % Tip chord of the main wing
cRh     = 0.65;  % Root chord of HTP  
cTh     = 0.45;  % Tip chord of HTP 
lh      = 3;     % Main wing - HTP separation
thetaT  = 0;     % Twist at the tip of the main wing
thetaTh  = 0;    % Twist at the tip of the HTP
iw      = 0;     % Main wing incidence angle
it      =-2;     % HTP incidence angle
aoa     = 4;     % Angle of attack of the main wing central section
rho     = 1.225; % Air density
Uinf    = 1;   % Freestream Velocity field module
Qinf    = Uinf*[cosd(aoa);sind(aoa)]; % Freestream Velocity field
% Parabolic drag: Cd = Cd0 + K*Cl^2
Cd0     = 0.0075;  % Zero lift drag coefficient
K       = 0.0055;  % Drag coefficient constant 
% NACA 0010 Lift Coefficient: Cl = Clalpha*aoaE+Cl0
Clalpha = 0.117306319973439; % Lift coefficient slope with aoa
Cl0     = 0.000308895559508056; % Zero aoa lift coefficient

% Geometry definition
N       = 10; % Number of span slices main wing
M       = 10; % Number of span slices HTP
[MW.coordsP,MW.coordsC,MW.deltaY,MW.c,MW.c12,MW.theta,MW.aoaE] = computeGeometryUniform(N,b,cR,cT,thetaT,aoa+iw);
[HTP.coordsP,HTP.coordsC,HTP.deltaY,HTP.c,HTP.c12,HTP.theta,HTP.aoaE] = computeGeometryUniform(M,bh,cRh,cTh,thetaTh,aoa+it);
coordsP = [MW.coordsP;HTP.coordsP];
coordsP(N+1:end,1) = coordsP(N+1:end,1) + lh; % HTP displacement
coordsC = [MW.coordsC;HTP.coordsC];
coordsC(N+1:end,1) = coordsC(N+1:end,1) + lh; % HTP displacement
deltaY  = [MW.deltaY;HTP.deltaY];
c       = [MW.c';HTP.c'];
c12     = [MW.c12;HTP.c12];
theta   = [MW.theta';HTP.theta'];
aoaE    = [MW.aoaE';HTP.aoaE'];

% Variable definition
q       = zeros(N+M,1); % Vector of independent terms changed notation from "b" to "q"
A       = zeros(N+M,N+M); % Influence matrix
aoaInd  = zeros(N+M,1);

% System of equations resolution
for i= 1:N
    q(i,1) = 1/2*c12(i)*norm(Qinf)*(Cl0+Clalpha*(aoaE(i)));
    for j = 1:N
        if i==j
            v = computeHorseshoeSelf(coordsP,coordsC,i,j,aoa);
            A(i,i) = -1/2*Clalpha*c12(i)*v*[-sind(aoa),0,cosd(aoa)]' + 1;
        else
            v = computeHorseshoe(coordsP,coordsC,i,j,aoa);
            A(i,j) = -1/2*Clalpha*c12(i)*v*[-sind(aoa),0,cosd(aoa)]'; 
        end
    end
    for j = N+1:N+M
            v = computeHorseshoe(coordsP,coordsC,i,j+1,aoa);
            A(i,j) = -1/2*Clalpha*c12(i)*v*[-sind(aoa),0,cosd(aoa)]'; 
    end
end
for i= N+1:N+M
    q(i,1) = 1/2*c12(i)*norm(Qinf)*(Cl0+Clalpha*(aoaE(i)));
    for j = 1:N
            v = computeHorseshoe(coordsP,coordsC,i,j,aoa);
            A(i,j) = -1/2*Clalpha*c12(i)*v*[-sind(aoa),0,cosd(aoa)]'; 
    end
    for j = N+1:N+M
        if i==j
            v = computeHorseshoeSelf(coordsP,coordsC,i,j+1,aoa);
            A(i,i) = -1/2*Clalpha*c12(i)*v*[-sind(aoa),0,cosd(aoa)]' + 1;
        else
            v = computeHorseshoe(coordsP,coordsC,i,j+1,aoa);
            A(i,j) = -1/2*Clalpha*c12(i)*v*[-sind(aoa),0,cosd(aoa)]'; 
        end
    end
end
T = A\q;

% % Individual slice bidimensional lift coefficient
% Cl12   = 2*T./(c12*norm(Qinf));
% for i = 1:N
%     aoaInd(i,1) = (Cl12(i) - Cl0)/Clalpha - (aoaE(i+1)+aoaE(i))/2;
% end
% 
% % Total wing lift verification
% L = rho*norm(Qinf)*sum(T.*deltaY);
% 
% % Induced Drag calculation 
% S = 3.375;
% Dind = -rho*norm(Qinf)*sum(T.*deltaY.*aoaInd);
% 
% CDind = Dind/(0.5*rho*norm(Qinf)^2*S);
% 
% % Total Lift coefficient calculation
% CL = 2*sum(T.*deltaY/(norm(Qinf)*S));
% 
% % Parameter
% parameter = CL^2/(pi*(b^2/S)*CDind);