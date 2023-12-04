function v = computeHorseshoeSelf(coordsP,coordsC,i,j,aoa)
r1 = coordsC(i,:)-coordsP(j,:);
r2 = coordsC(i,:)-coordsP(j+1,:);
u = -[cos(aoa),0,sin(aoa)];

vInfA = computeSemiVortex(r1,u);
vInfB = computeSemiVortex(r2,u);

v = vInfA - vInfB;
end