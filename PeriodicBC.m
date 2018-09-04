function PeriodicBC(nx, ny, nz, dumx, dumy, dumz, fname)
% Ran Ma
% 7/5/2018
% write periodic condition for Warp3d
% 7/6/2018
% modified to remove redundant constraint
%
%                          write output file
%           equation to get node ID from 3D coordinate: 
%               (a,b,c) => a + (b-1)*nx + (c-1) * nx * ny
%
fid = fopen(fname,'w');
% ============================= write surface =============================
fprintf(fid,'\nc     Now writing surface\n');
% x direction
for kk = 1:(nz-2)
    for jj = 1:(ny-2)
        node1 = nx*jj + nx*ny*kk + 1;                          % (1 ,jj,kk)
        node2 = node1 + nx - 1;                                % (nx,jj,kk)
        fprintf(fid,'  + %d 1.0 u - %d 1.0 u - %d 1.0 u = 0.0\n',node2,node1,dumx);
        fprintf(fid,'  + %d 1.0 v - %d 1.0 v - %d 1.0 v = 0.0\n',node2,node1,dumx);
        fprintf(fid,'  + %d 1.0 w - %d 1.0 w - %d 1.0 w = 0.0\n',node2,node1,dumx);
    end
end

% y direction
for kk = 1:(nz-2)
    for ii = 1:(nx-2)
        node1 = ii + nx*ny*kk + 1;                             % (ii,1 ,kk)
        node2 = node1 + nx*(ny-1);                             % (ii,ny,kk)
        fprintf(fid,'  + %d 1.0 u - %d 1.0 u - %d 1.0 u = 0.0\n',node2,node1,dumy);
        fprintf(fid,'  + %d 1.0 v - %d 1.0 v - %d 1.0 v = 0.0\n',node2,node1,dumy);
        fprintf(fid,'  + %d 1.0 w - %d 1.0 w - %d 1.0 w = 0.0\n',node2,node1,dumy);
    end
end

% z direction
for jj = 1:(ny-2)
    for ii = 1:(nx-2)
        node1 = ii + nx*jj + 1;                                % (ii,jj,1 )
        node2 = node1 + nx*ny*(nz-1);                          % (ii,jj,nz)
        fprintf(fid,'  + %d 1.0 u - %d 1.0 u - %d 1.0 u = 0.0\n',node2,node1,dumz);
        fprintf(fid,'  + %d 1.0 v - %d 1.0 v - %d 1.0 v = 0.0\n',node2,node1,dumz);
        fprintf(fid,'  + %d 1.0 w - %d 1.0 w - %d 1.0 w = 0.0\n',node2,node1,dumz);
    end
end

% ============================== write edge ===============================
fprintf(fid,'\nc     Now writing edge\n');
% edge x
for ii = 2:(nx-1)
    % edge pair 1
    node1 = ii;                                                % (ii,1 ,1 )
    node2 = ii + (ny-1)*nx;                                    % (ii,ny,1 )
    fprintf(fid,'  + %d 1.0 u - %d 1.0 u - %d 1.0 u = 0.0\n',node2,node1,dumy);
    fprintf(fid,'  + %d 1.0 v - %d 1.0 v - %d 1.0 v = 0.0\n',node2,node1,dumy);
    fprintf(fid,'  + %d 1.0 w - %d 1.0 w - %d 1.0 w = 0.0\n',node2,node1,dumy);
    % edge pair 2
    node1 = ii;                                                % (ii,1 ,1 )
    node2 = ii + nx*ny*nz - nx;                                % (ii,ny,nz)
    fprintf(fid,'  + %d 1.0 u - %d 1.0 u - %d 1.0 u - %d 1.0 u = 0.0\n',node2,node1,dumy,dumz);
    fprintf(fid,'  + %d 1.0 v - %d 1.0 v - %d 1.0 v - %d 1.0 v = 0.0\n',node2,node1,dumy,dumz);
    fprintf(fid,'  + %d 1.0 w - %d 1.0 w - %d 1.0 w - %d 1.0 w = 0.0\n',node2,node1,dumy,dumz);
    % edge pair 3
    node1 = ii;                                                % (ii,1 ,1 )
    node2 = ii + (nz-1)*nx*ny;                                 % (ii,1 ,nz)
    fprintf(fid,'  + %d 1.0 u - %d 1.0 u - %d 1.0 u = 0.0\n',node2,node1,dumz);
    fprintf(fid,'  + %d 1.0 v - %d 1.0 v - %d 1.0 v = 0.0\n',node2,node1,dumz);
    fprintf(fid,'  + %d 1.0 w - %d 1.0 w - %d 1.0 w = 0.0\n',node2,node1,dumz);
end

% edge y
for jj = 2:(ny-1)
    % edge pair 1
    node1 = 1 + (jj-1)*nx;                                     % ( 1,jj, 1)
    node2 = jj*nx;                                             % (nx,jj, 1)
    fprintf(fid,'  + %d 1.0 u - %d 1.0 u - %d 1.0 u = 0.0\n',node2,node1,dumx);
    fprintf(fid,'  + %d 1.0 v - %d 1.0 v - %d 1.0 v = 0.0\n',node2,node1,dumx);
    fprintf(fid,'  + %d 1.0 w - %d 1.0 w - %d 1.0 w = 0.0\n',node2,node1,dumx);
    % edge pair 2
    node1 = 1 + (jj-1)*nx;                                     % ( 1,jj, 1)
    node2 = jj*nx + (nz-1)*nx*ny;                              % (nx,jj,nz)
    fprintf(fid,'  + %d 1.0 u - %d 1.0 u - %d 1.0 u - %d 1.0 u = 0.0\n',node2,node1,dumx,dumz);
    fprintf(fid,'  + %d 1.0 v - %d 1.0 v - %d 1.0 v - %d 1.0 v = 0.0\n',node2,node1,dumx,dumz);
    fprintf(fid,'  + %d 1.0 w - %d 1.0 w - %d 1.0 w - %d 1.0 w = 0.0\n',node2,node1,dumx,dumz);
    % edge pair 3
    node1 = 1 + (jj-1)*nx;                                     % ( 1,jj, 1)
    node2 = 1 + (jj-1)*nx + (nz-1)*nx*ny;                      % ( 1,jj,nz)
    fprintf(fid,'  + %d 1.0 u - %d 1.0 u - %d 1.0 u = 0.0\n',node2,node1,dumz);
    fprintf(fid,'  + %d 1.0 v - %d 1.0 v - %d 1.0 v = 0.0\n',node2,node1,dumz);
    fprintf(fid,'  + %d 1.0 w - %d 1.0 w - %d 1.0 w = 0.0\n',node2,node1,dumz);
end

% edge x
for kk = 2:(nz-1)
    % edge pair 1
    node1 = 1 + (kk-1)*nx*ny;                                  % ( 1, 1,kk)
    node2 = kk*nx*ny - nx + 1;                                 % ( 1,ny,kk)
    fprintf(fid,'  + %d 1.0 u - %d 1.0 u - %d 1.0 u = 0.0\n',node2,node1,dumy);
    fprintf(fid,'  + %d 1.0 v - %d 1.0 v - %d 1.0 v = 0.0\n',node2,node1,dumy);
    fprintf(fid,'  + %d 1.0 w - %d 1.0 w - %d 1.0 w = 0.0\n',node2,node1,dumy);
    % edge pair 2
    node1 = 1 + (kk-1)*nx*ny;                                  % ( 1, 1,kk)
    node2 = kk*nx*ny;                                          % (nx,ny,kk)
    fprintf(fid,'  + %d 1.0 u - %d 1.0 u - %d 1.0 u - %d 1.0 u = 0.0\n',node2,node1,dumx,dumy);
    fprintf(fid,'  + %d 1.0 v - %d 1.0 v - %d 1.0 v - %d 1.0 v = 0.0\n',node2,node1,dumx,dumy);
    fprintf(fid,'  + %d 1.0 w - %d 1.0 w - %d 1.0 w - %d 1.0 w = 0.0\n',node2,node1,dumx,dumy);
    % edge pair 3
    node1 = 1 + (kk-1)*nx*ny;                                  % ( 1, 1,kk)
    node2 = nx + (kk-1)*nx*ny;                                 % (nx, 1,kk)
    fprintf(fid,'  + %d 1.0 u - %d 1.0 u - %d 1.0 u = 0.0\n',node2,node1,dumx);
    fprintf(fid,'  + %d 1.0 v - %d 1.0 v - %d 1.0 v = 0.0\n',node2,node1,dumx);
    fprintf(fid,'  + %d 1.0 w - %d 1.0 w - %d 1.0 w = 0.0\n',node2,node1,dumx);
end

% ============================= write corner ==============================
fprintf(fid,'\nc     Now writing corner\n');
%                      1  1  1    =>    1
%                     nx  1  1    =>    nx
%                     nx ny  1    =>    nx*ny
%                      1 ny  1    =>    (ny-1)*nx+1
%                      1  1 nz    =>    (nz-1)*nx*ny+1
%                     nx  1 nz    =>    (nz-1)*nx*ny+nx
%                     nx ny nz    =>    nx*ny*nz
%                      1 ny nz    =>    nx*ny*nz-nx+1
node1 = 1;                                                     % ( 1, 1, 1)
node2 = (nz-1)*nx*ny + 1;                                      % ( 1, 1,nz)
fprintf(fid,'  + %d 1.0 u - %d 1.0 u - %d 1.0 u = 0.0\n',node2,node1,dumz);
fprintf(fid,'  + %d 1.0 v - %d 1.0 v - %d 1.0 v = 0.0\n',node2,node1,dumz);
fprintf(fid,'  + %d 1.0 w - %d 1.0 w - %d 1.0 w = 0.0\n',node2,node1,dumz);

node1 = 1;                                                    % ( 1, 1, 1)
node2 = (nz-1)*nx*ny + nx;                                     % (nx, 1,nz)
fprintf(fid,'  + %d 1.0 u - %d 1.0 u - %d 1.0 u - %d 1.0 u = 0.0\n',node2,node1,dumx,dumz);
fprintf(fid,'  + %d 1.0 v - %d 1.0 v - %d 1.0 v - %d 1.0 v = 0.0\n',node2,node1,dumx,dumz);
fprintf(fid,'  + %d 1.0 w - %d 1.0 w - %d 1.0 w - %d 1.0 w = 0.0\n',node2,node1,dumx,dumz);

node1 = 1;                                         % ( 1, 1, 1)
node2 = nx*ny*nz-nx+1;                                         % ( 1,ny,nz)
fprintf(fid,'  + %d 1.0 u - %d 1.0 u - %d 1.0 u - %d 1.0 u = 0.0\n',node2,node1,dumy,dumz);
fprintf(fid,'  + %d 1.0 v - %d 1.0 v - %d 1.0 v - %d 1.0 v = 0.0\n',node2,node1,dumy,dumz);
fprintf(fid,'  + %d 1.0 w - %d 1.0 w - %d 1.0 w - %d 1.0 w = 0.0\n',node2,node1,dumy,dumz);

node1 = 1;                                                 % ( 1, 1, 1)
node2 = nx*ny*nz;                                              % (nx,ny,nz)
fprintf(fid,'  + %d 1.0 u - %d 1.0 u - %d 1.0 u - %d 1.0 u - %d 1.0 u = 0.0\n',node2,node1,dumx,dumy,dumz);
fprintf(fid,'  + %d 1.0 v - %d 1.0 v - %d 1.0 v - %d 1.0 v - %d 1.0 v = 0.0\n',node2,node1,dumx,dumy,dumz);
fprintf(fid,'  + %d 1.0 w - %d 1.0 w - %d 1.0 w - %d 1.0 w - %d 1.0 w = 0.0\n',node2,node1,dumx,dumy,dumz);

node1 = 1;                                                     % ( 1, 1, 1)
node2 = nx;                                                    % (nx, 1, 1)
fprintf(fid,'  + %d 1.0 u - %d 1.0 u - %d 1.0 u = 0.0\n',node2,node1,dumx);
fprintf(fid,'  + %d 1.0 v - %d 1.0 v - %d 1.0 v = 0.0\n',node2,node1,dumx);
fprintf(fid,'  + %d 1.0 w - %d 1.0 w - %d 1.0 w = 0.0\n',node2,node1,dumx);

node1 = 1;                                                     % ( 1, 1, 1)
node2 = (ny-1)*nx+1;                                           % ( 1,ny, 1)
fprintf(fid,'  + %d 1.0 u - %d 1.0 u - %d 1.0 u = 0.0\n',node2,node1,dumy);
fprintf(fid,'  + %d 1.0 v - %d 1.0 v - %d 1.0 v = 0.0\n',node2,node1,dumy);
fprintf(fid,'  + %d 1.0 w - %d 1.0 w - %d 1.0 w = 0.0\n',node2,node1,dumy);

node1 = 1;                                                    % ( 1, 1, 1)
node2 = nx*ny;                                                 % (nx,ny, 1)
fprintf(fid,'  + %d 1.0 u - %d 1.0 u - %d 1.0 u - %d 1.0 u = 0.0\n',node2,node1,dumx,dumy);
fprintf(fid,'  + %d 1.0 v - %d 1.0 v - %d 1.0 v - %d 1.0 v = 0.0\n',node2,node1,dumx,dumy);
fprintf(fid,'  + %d 1.0 w - %d 1.0 w - %d 1.0 w - %d 1.0 w = 0.0\n',node2,node1,dumx,dumy);

fclose(fid);
end