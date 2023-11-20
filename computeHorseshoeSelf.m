function v = computeHorseshoeSelf(coordsP,coordsC,i,j)
r1 = coordsC(i,:)-coordsP(j,:);
r2 = coordsC(i,:)-coordsP(j+1,:);
u = [-1,0,0];
vInfA = computeSemiVortex(r1,u);
vInfB = computeSemiVortex(r2,u);

v = vInfA - vInfB;
end