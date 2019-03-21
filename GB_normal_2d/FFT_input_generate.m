% Ran Ma
% 3/2/2019
%
% generate FFT input file with grain boundary thickness
%
clear
clc

% user input
GB_thickness = 2; % thickness of GB in element
GB_smoothness = 10; % Laplacian smoothing iterations
PF_iso_factor_GB = 10; % anisotropic factor
fname = 'RVE.tesr'; % neper input file
fname_ebsd = 'EBSD_RVE'; % EBSD output file
fname_gb = 'EBSD_GB'; % EBSD output file with GB
fname_PF = 'phase_0'; % FFT input file
cs = crystalSymmetry('cubic');
ss = specimenSymmetry('triclinic');

% read neper output file
neut_tesr_fscanf

% write EBSD file required by MTEX to extract grain boundary
dataAll = zeros(7,numel);
% dataAll(1,:)   % index
% dataAll(2,:)   % phase
% dataAll(3:4,:) % coordinate
% dataAll(5:7,:) % Euler angle
WriteEBSD

% call MTEX function to detect grain boundary
ebsd = loadEBSD_generic(fname_ebsd,'CS',cs,'SS',ss, 'ColumnNames', ...
{'Index' 'Phase' 'x' 'y' 'Euler1' 'Euler2' 'Euler3'}, 'Bunge');
[grains,ebsd.grainId] = calcGrains(ebsd);
grains = smooth(grains,GB_smoothness); % Laplacian smoothing
gB = grains.boundary;

% assign different material to grain boundary material
% gb.F - vertice to face
% gb.ebsdID - element besides interface
% gb.V - vertice coordinate
% gb.I_VF - incidence matrix, vertice x edge
% gb.A_V = I_VF * transpose(I_VF), extract adjacent vertice along GB
%          diag: sum of adjacent vertice
num_gb = length(gB.ebsdId);
gb_n = zeros(2,numel);
for ii = 1:num_gb
    jj = gB.ebsdId(ii,1);
    kk = gB.ebsdId(ii,2);
    if jj > 0 && kk > 0 % neither jj or kk is surrounding matrix
        % assign phase 1 Euler angle 0,0,0 to these 2 elements
        dataAll(2,jj) = 1;
        dataAll(5:7,jj) = 0;
        dataAll(2,kk) = 1;
        dataAll(5:7,kk) = 0;
        
        % grain boundary normal direction
        gb_n_local = [gB.direction.y(ii); -gB.direction.x(ii)];
        gb_n(:,jj) = gb_n_local;
        gb_n(:,kk) = gb_n_local;

        % element on one side of the GB
        % element(jj) => grid(aa,bb)
        bb = floor((jj-1)/numgrid(1)) + 1;
        aa = mod(jj-1,numgrid(1)) + 1;
        % grounp of GB element: grid(xx,yy)
        xx = max(aa-GB_thickness,1):min(aa+GB_thickness,numgrid(1));
        yy = max(bb-GB_thickness,1):min(bb+GB_thickness,numgrid(2));
        xx2 = transpose(xx) * ones(1,length(yy));
        yy2 = ones(length(xx),1) * yy;
        % grid(xx,yy) => element(temp)
        temp = (yy2-1)*numgrid(1)+xx2;
        temp = reshape(temp,[],1);
        % assign GB phase and Euler angle to GB
        dataAll(2,temp) = 1;
        dataAll(5:7,temp) = 0;
        % assign GB normal direction
        gb_n_local = [gB.direction.x(ii); gB.direction.y(ii)];
        gb_n(1,temp) = gb_n_local(1);
        gb_n(2,temp) = gb_n_local(2);
        
        % the other element on the other side of GB
        bb = floor((jj-1)/numgrid(1)) + 1;
        aa = mod(jj-1,numgrid(1)) + 1;
        xx = max(aa-GB_thickness,1):min(aa+GB_thickness,numgrid(1));
        yy = max(bb-GB_thickness,1):min(bb+GB_thickness,numgrid(2));
        xx2 = transpose(xx) * ones(1,length(yy));
        yy2 = ones(length(xx),1) * yy;
        temp = (yy2-1)*numgrid(1)+xx2;
        temp = reshape(temp,[],1);
        dataAll(2,temp) = 1;
        dataAll(5:7,temp) = 0;
        gb_n_local = [gB.direction.x(ii); gB.direction.y(ii)];
        gb_n(1,temp) = gb_n_local(1);
        gb_n(2,temp) = gb_n_local(2);
            
    end
end

% write new EBSD file for verification
fid = fopen(fname_gb,'wt');
fprintf(fid,'Index, Phase, X, Y, Euler1, Euler2, Euler3\n');
temp = dataAll(2,:);
dataAll(2,:) = 0;
fprintf(fid,'%4i\t%1i\t%10.8f\t%10.8f\t%5.2f\t%5.2f\t%5.2f\n',dataAll);
dataAll(2,:) = temp;
fclose(fid);
ebsd = loadEBSD_generic(fname_gb,'CS',cs,'SS',ss, 'ColumnNames', ...
{'Index' 'Phase' 'x' 'y' 'Euler1' 'Euler2' 'Euler3'}, 'Bunge');
plot(ebsd,ebsd.orientations.angle./degree);

% write FFT phase field input file
fid = fopen(fname_PF,'wt');
phase_field_0 = genCirc(numgrid);
PF_iso_factor = dataAll(2,:) * PF_iso_factor_GB;
elemtext = '%5.2f, %5.2f, %10.8f, %10.8f\n';
fprintf(fid,elemtext,[phase_field_0; PF_iso_factor; gb_n]);
fclose(fid);


function eta = genCirc(num_grid)
%
% circle in a square matrix
%
eta = zeros(num_grid(1),num_grid(2));
radius = 5;

for ii = 1 : num_grid(1)
    for jj = 1 : num_grid(2)
        dist = sqrt( ( ii - 49 ) ^ 2 + ( jj - 40 ) ^ 2 );
        if dist < radius
            eta(ii,jj) = 1;
        end
    end
end

eta = reshape(eta,1,[]);

end
