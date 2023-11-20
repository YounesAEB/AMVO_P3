function v = computeHorseshoe(coordsP,coordsC,i,j)
r1 = coordsC(i,:)-coordsP(j,:);
r2 = coordsC(i,:)-coordsP(j+1,:);
u = [-1,0,0];

vInfA = computeSemiVortex(r1,u);
vAB = computeFiniteVortex(r1,r2);
vInfB = computeSemiVortex(r2,u);

v = vInfA + vAB - vInfB;
end