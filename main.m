b  = 4.5;  % Wingspan
cR = 0.9;  % Root chord
cT = 0.6;  % Tip chord
theta = 0; % Twist

m = 40; % Mass in kilograms

N = 2;

[yP, yC] = computeGeometry(N,b);