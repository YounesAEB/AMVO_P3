%--------------------------------------------------------------------------
%  /  POTENTIAL AERODYNAMICS - LIFTING LINE METHOD -AMVO 
%  /  Matlab code to assess the numerical solution of LLM for a wing.                                            
%  /  ESEIAAT_UPC                                           
%  /  MUEA - MQ1 - Younes Akhazzan - Joel Rajo - Pol Ruiz                         
%--------------------------------------------------------------------------
clc; clear; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% Data given by the excercise statement
b       = 4.5;   % Wingspan
cR      = 0.9;   % Root chord
cT      = 0.6;   % Tip chord
thetaT  = 0;     % Twist at the tip
aoa     = 0;     % Angle of attack
m       = 40;    % Mass in kilograms
rho     = 1.225; % Air density
g0      = 9.81;  % Gravity's acceleration
W       = m*g0;  % Wing's weight
Cl0     = 0.24;  % Zero degree lift coefficient
Clalpha = 6.7;   % Lift coefficient slope with angle of attack

% Desired Freestream velocity process
delta = 0.00001;
L     = 0; % Lift initialization
Uinf    = 10;   % Freestream Velocity field module
while abs(L-W)>delta
Qinf    = Uinf*[cosd(aoa);sind(aoa)]; % Freestream Velocity field

% Geometry definition
N = 10; % Number of span slices


%[coordsP,coordsC, deltaY,c,c12,theta,aoaE] = computeGeometryUniform(N,b,cR,cT,thetaT,aoa);
[coordsP, coordsC, deltaY,c,c12,theta,aoaE] = computeGeometryCosine(N,b,cR,cT,thetaT,aoa);

% Variable definition
q    = zeros(N,1); % Vector of independent terms changed notation from "b" to "q"
A    = zeros(N,N); % Influence matrix
aoaInd = zeros(N,1);

% System of equations resolution
for i= 1:N
    q(i,1) = 1/2*c12(i)*norm(Qinf)*(Cl0+Clalpha*(aoaE(i)));
    for j = 1:N
        if i==j
            v = computeHorseshoeSelf(coordsP,coordsC,i,j);
            A(i,i) = -1/2*Clalpha*c12(i)*v*[-sind(aoa),0,cosd(aoa)]' + 1;
        else
            v = computeHorseshoe(coordsP,coordsC,i,j);
            A(i,j) = -1/2*Clalpha*c12(i)*v*[-sind(aoa),0,cosd(aoa)]'; 
        end
    end
end
T = A\q;

% Individual slice bidimensional lift coefficient
Cl12   = 2*T./(c12*norm(Qinf));
for i = 1:N
    aoaInd(i,1) = (Cl12(i) - Cl0)/Clalpha - (aoaE(i+1)+aoaE(i))/2;
end
% Total wing lift verification
L = rho*norm(Qinf)*sum(T.*deltaY);

if abs(L-W)<delta
    msg = sprintf('Desired freestream velocity: Uinf=%i with L=W=%i \n', Uinf, L);
    disp(msg);
else 
    msg = sprintf('Inappropiate freestream velocity with L= %i and W=%i \n', L, W);
    disp(msg);   
    Uinf = (Uinf+(m*g0)/(rho*sum(T.*deltaY)))/2;  
end
end

% Induced Drag calculation 
Dind = -rho*norm(Qinf)*sum(T.*deltaY.*aoaInd);