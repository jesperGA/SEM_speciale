function mesh = box_mesh_gen(X,Y,nelx,nely,N)
% mesh = box_mesh_gen(X,Y,nelx,nely,N)
% Input: X = [x0 x1] vector containing start and end x coordinates
%        Y = [y0 y1] vector containing start and end x coordinates
%        nelx and nely = number of elements in mesh
%        N = polynomial order of velocity
% 

[iglobP, xNP,yNP] = MeshBox_mod(2,2,2,2,study.n_GL,2);
[iglobV, xNV,yNV] = MeshBox_mod(2,2,2,2,n_GLL,1);

end