function v = computeFiniteVortex(r1,r2)
v = 1/(4*pi)*(norm(r1)+norm(r2))/...
    (norm(r1)*norm(r2)*(norm(r1)*norm(r2)+dot(r1,r2)))*cross(r1,r2);
end