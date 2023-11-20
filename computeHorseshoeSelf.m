function v = computeHorseshoeSelf(coordsP,coordsC,i,j)
r1 = coordsC(i,:)-coordsP(j,:);
r2 = coordsC(i,:)-coordsP(j+1,:);
uA = [-1,0,0];
uB = [-1,0,0];
vInfA = computeSemiVortex(r1,uA);
vInfB = computeSemiVortex(r2,uB);

v = vInfA - vInfB;
end