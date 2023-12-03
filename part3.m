%--------------------------------------------------------------------------
%  /  POTENTIAL AERODYNAMICS - LIFTING LINE METHOD -AMVO 
%  /  Matlab code to assess the numerical solution via LLM - Part 3                                            
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
thetaT  = -2.5;  % Twist at the tip of the main wing
thetaTh = 0;     % Twist at the tip of the HTP
iw      = 0;     % Main wing incidence angle
it      =-2;     % HTP incidence angle
aoa     = 4;     % Angle of attack of the main wing central section
delta_l   = -10; % Elevator deflection angle
delta_r   = 10;  % Elevator deflection angle
rho     = 1.225; % Air density
Uinf    = 1;     % Freestream Velocity field module
Qinf    = Uinf*[cosd(aoa);sind(aoa)]; % Freestream Velocity field
% NACA 0015 Lift Coefficient: Cl = Clalpha*aoaE+Cl0+Cld*d
Clalpha_15 = 0.115491628925204; % Lift coefficient slope with aoa
Cl0_15     = 0.000272585908644561; % Zero aoa lift coefficient
Cld        = 0.0724135854767064; % Lift coefficient slope with flap deflection 
% NACA 0015 Momentum Coefficient:
Cm0_15    = 0; % Zero pitching moment about the aerodynamic center in symetric airfoils 
% NACA 0010 Lift Coefficient: Cl = Clalpha*aoaE+Cl0
Clalpha_10 = 0.117306319973439; % Lift coefficient slope with aoa
Cl0_10     = 0.000308895559508056; % Zero aoa lift coefficient
% NACA 0010 Momentum Coefficient:
Cm0_10   = 0; % Zero pitching moment about the aerodynamic center in symetric airfoils

% Geometry definition
N       = 512; % Number of span slices main wing
M       = 256; % Number of span slices HTP
[MW.coordsP,MW.coordsC,MW.deltaY,MW.c,MW.c12,MW.theta,MW.aoaE] = computeGeometryUniform(N,b,cR,cT,thetaT,aoa+iw);
[HTP.coordsP,HTP.coordsC,HTP.deltaY,HTP.c,HTP.c12,HTP.theta,HTP.aoaE] = computeGeometryUniform(M,bh,cRh,cTh,thetaTh,aoa+it);
coordsP = [MW.coordsP;HTP.coordsP];
coordsP(N+2:end,1) = coordsP(N+2:end,1) + lh; % HTP displacement
coordsP(N+2:end,3) = coordsP(N+2:end,3) - 0.05; % Zero angle interference correction
coordsC = [MW.coordsC;HTP.coordsC];
coordsC(N+1:end,1) = coordsC(N+1:end,1) + lh; % HTP displacement
coordsC(N+1:end,3) = coordsC(N+1:end,3) - 0.05; % Zero angle interference correction
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
    if (coordsC(i,2)<=-(b/2)+ba)
        q(i,1) = 1/2*c12(i)*norm(Qinf)*(Cl0_15+Clalpha_15*((aoaE(i)+aoaE(i+1))/2)+Cld*delta_l);
        for j = 1:N
            if i==j
                v = computeHorseshoeSelf(coordsP,coordsC,i,j,aoa);
                A(i,i) = -1/2*Clalpha_15*c12(i)*v*[-sind(aoa),0,cosd(aoa)]' + 1;
            else
                v = computeHorseshoe(coordsP,coordsC,i,j,aoa);
                A(i,j) = -1/2*Clalpha_15*c12(i)*v*[-sind(aoa),0,cosd(aoa)]'; 
            end
        end
        for j = N+1:N+M
                v = computeHorseshoe(coordsP,coordsC,i,j+1,aoa);
                A(i,j) = -1/2*Clalpha_10*c12(i)*v*[-sind(aoa),0,cosd(aoa)]'; 
        end
    elseif (coordsC(i,2)>=(b/2)-ba)
        q(i,1) = 1/2*c12(i)*norm(Qinf)*(Cl0_15+Clalpha_15*((aoaE(i)+aoaE(i+1))/2)+Cld*delta_r);
        for j = 1:N
            if i==j
                v = computeHorseshoeSelf(coordsP,coordsC,i,j,aoa);
                A(i,i) = -1/2*Clalpha_15*c12(i)*v*[-sind(aoa),0,cosd(aoa)]' + 1;
            else
                v = computeHorseshoe(coordsP,coordsC,i,j,aoa);
                A(i,j) = -1/2*Clalpha_15*c12(i)*v*[-sind(aoa),0,cosd(aoa)]'; 
            end
        end
        for j = N+1:N+M
                v = computeHorseshoe(coordsP,coordsC,i,j+1,aoa);
                A(i,j) = -1/2*Clalpha_10*c12(i)*v*[-sind(aoa),0,cosd(aoa)]'; 
        end
    else
        q(i,1) = 1/2*c12(i)*norm(Qinf)*(Cl0_10+Clalpha_10*((aoaE(i)+aoaE(i+1))/2));
        for j = 1:N
            if i==j
                v = computeHorseshoeSelf(coordsP,coordsC,i,j,aoa);
                A(i,i) = -1/2*Clalpha_10*c12(i)*v*[-sind(aoa),0,cosd(aoa)]' + 1;
            else
                v = computeHorseshoe(coordsP,coordsC,i,j,aoa);
                A(i,j) = -1/2*Clalpha_10*c12(i)*v*[-sind(aoa),0,cosd(aoa)]'; 
            end
        end
        for j = N+1:N+M
                v = computeHorseshoe(coordsP,coordsC,i,j+1,aoa);
                A(i,j) = -1/2*Clalpha_10*c12(i)*v*[-sind(aoa),0,cosd(aoa)]'; 
        end
    end    
end
for i= N+1:N+M
    q(i,1) = 1/2*c12(i)*norm(Qinf)*(Cl0_10+Clalpha_10*((aoaE(i+1)+aoaE(i+2))/2));
    for j = 1:N
        if abs(coordsC(i,2))>=((b/2)-ba) 
            v = computeHorseshoe(coordsP,coordsC,i,j,aoa);
            A(i,j) = -1/2*Clalpha_15*c12(i)*v*[-sind(aoa),0,cosd(aoa)]';
        else
            v = computeHorseshoe(coordsP,coordsC,i,j,aoa);
            A(i,j) = -1/2*Clalpha_10*c12(i)*v*[-sind(aoa),0,cosd(aoa)]';
        end
    end
    for j = N+1:N+M
        if i==j
            v = computeHorseshoeSelf(coordsP,coordsC,i,j+1,aoa);
            A(i,i) = -1/2*Clalpha_10*c12(i)*v*[-sind(aoa),0,cosd(aoa)]' + 1;
        else
            v = computeHorseshoe(coordsP,coordsC,i,j+1,aoa);
            A(i,j) = -1/2*Clalpha_10*c12(i)*v*[-sind(aoa),0,cosd(aoa)]'; 
        end
    end
end
T = A\q;

% Total Lift coefficient calculation
Sw = 2*(b/2*(cR+cT)/2);    % Main wing surface
Sh = 2*(bh/2*(cRh+cTh)/2); % HTP surface
CL = 2*sum(T.*deltaY/(norm(Qinf)*Sw));
% Total Lift 
L = rho*norm(Qinf)*sum(T.*deltaY);
% Individual slice bidimensional lift coefficient
Cl12   = 2*T./(c12*norm(Qinf));

% Individual slice induced angle of attack
for i = 1:N
    if (coordsC(i,2)<=-(b/2)+ba)
        aoaInd(i,1) = (Cl12(i) - Cl0_15 - Cld*delta_l)/Clalpha_15 - (aoaE(i+1)+aoaE(i))/2;
    elseif (coordsC(i,2)<=(b/2)-ba)
        aoaInd(i,1) = (Cl12(i) - Cl0_15 - Cld*delta_r)/Clalpha_15 - (aoaE(i+1)+aoaE(i))/2;
    else
        aoaInd(i,1) = (Cl12(i) - Cl0_10)/Clalpha_10 - (aoaE(i+1)+aoaE(i))/2;
    end
end
for i = N+1:N+M
    aoaInd(i,1) = (Cl12(i) - Cl0_10)/Clalpha_10 - (aoaE(i+2)+aoaE(i+1))/2;
end

% Individual slice bidimensional induced drag coefficient
Cdind   = -2*T.*aoaInd./(norm(Qinf).*c12);
% Induced Drag calculation 
Dind = -rho*norm(Qinf)*sum(T.*deltaY.*aoaInd);
CDind = Dind/(0.5*rho*norm(Qinf)^2*Sw);

 
% % Pitching moment coefficient 
% lambda = cT/cR; % Tip-to-Root chord ratio
% mac = 2/3*cR*(1+lambda+lambda^2)/(1+lambda); % Mean aerodynamic chord
% CM0 = Cm14 -2*sum(coordsC(:,1).*T.*deltaY)/(norm(Qinf)*Sw*mac);
% M0  = CM0*0.5*rho*norm(Qinf)^2*Sw*mac;
% msg =sprintf("Global CL=%i, CD=%i and CM0=%i",CL,CD,CM0);
% disp(msg);

% Plot of the lift coefficients per slice
figure
hold on
title("Spanwise distribution of the local coefficients of lift")
plot((2/b)*[-b/2;coordsC(1:N,2);b/2],[0;Cl12(1:N,1);0]);
plot((2/bh)*[-bh/2;coordsC(N+1:N+M,2);bh/2],[0;Cl12(N+1:N+M,1);0]);
xlabel("$2y/b$");
ylabel("Lift Coefficient $C_{l}$");
legend("Main Wing","Horizontal Tail Plane","Location","south");
xlim([-1,1]);
grid on;
grid minor;
box on;
axis padded
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize',13);
hold off;

% % Plot of the viscous drag coefficients per slice
% figure
% hold on
% title("Spanwise distribution of the local coefficients","of viscous drag")
% plot((2/b)*[-b/2;coordsC(1:N,2);b/2],[0;Cdv(1:N,1);0]);
% plot((2/bh)*[-bh/2;coordsC(N+1:N+M,2);bh/2],[0;Cdv(N+1:N+M,1);0]);
% xlabel("$2y/b$");
% ylabel("Viscous Drag Coefficient $C_{d_v}$");
% legend("Main Wing","Horizontal Tail Plane","Location","south");
% xlim([-1,1]);
% grid on;
% grid minor;
% box on;
% axis padded
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize',13);
% hold off;
% % Plot of the induced drag coefficients per slice
% figure
% hold on
% title("Spanwise distribution of the local coefficients","of induced drag")
% plot((2/b)*[-b/2;coordsC(1:N,2);b/2],[0;Cdv(1:N,1);0]);
% plot((2/bh)*[-bh/2;coordsC(N+1:N+M,2);bh/2],[0;Cdv(N+1:N+M,1);0]);
% xlabel("$2y/b$");
% ylabel("Induced Drag Coefficient $C_{d_{ind}}$");
% legend("Main Wing","Horizontal Tail Plane","Location","south");
% xlim([-1,1]);
% grid on;
% grid minor;
% box on;
% axis padded
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize',13);
% hold off;
% % Plot of the induced angle of attack per slice
% figure
% hold on
% title("Spanwise distribution of the local","induced angle of attack")
% plot((2/b)*[coordsC(1:N,2)],[aoaInd(1:N,1)]);
% plot((2/bh)*[coordsC(N+1:N+M,2)],[aoaInd(N+1:N+M,1)]);
% xlabel("$2y/b$");
% ylabel("Induced Angle of Attack $\alpha_{ind}$");
% legend("Main Wing","Horizontal Tail Plane","Location","south");
% xlim([-1,1]);
% grid on;
% grid minor;
% box on;
% axis padded
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize',13);
% hold off;