function v = computeSemiVortex(r,u)
v = 1/(4*pi)*((1-dot(u,(r/norm(r))))/(norm(cross(u,r))^2))*(cross(u,r));

end