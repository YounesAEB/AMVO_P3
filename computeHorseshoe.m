function v = computeHorseshoe(coordsP,coordsC,i,j)
r1 = coordsC(i,:)-coordsP(j,:);
r2 = coordsC(i,:)-coordsP(j+1,:);
uA = [1,0,0];
uB = [-1,0,0];

vInfA = computeSemiVortex(r1,uA);
vAB = computeFiniteVortex(r1,r2);
vInfB = computeSemiVortex(r2,uB);

v = vInfA + vAB - vInfB;
end