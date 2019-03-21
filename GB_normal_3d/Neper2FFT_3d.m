% Ran ma
% 3/8/2019
%
% extract grain boundary element and assign normal direction
%
% run the following bash script first
%
% #!/bin/bash
% NEPER="neper --rcfile none"
% if [ -a "gene_morp_2.tess" ]
% then
%   $NEPER -T -loadtess gene_morp_2.tess -format tess,tesr -tesrsize 64 -o gene_form_1
% else
%   $NEPER -T -n 20 -morpho gg -o gene_form_1 -format tess,tesr -tesrsize 64
% fi
% mv gene_form_1.tesr gene_form_2.tesr
% C="-datacellcol id -cameraangle 12 -imagesize 600:600"
% $NEPER -V gene_form_1.tess $C                   -print gene_form_1
% $NEPER -V gene_form_2.tesr $C -datarptedgerad 0 -print gene_form_2
% convert +append gene_form_?.png gene_form.png
% exit 0
%
clear
clc
fclose('all');

thick_half = 0.03; % GB thickness, relative to coordinate
fname_tess = 'gene_form_1.tess'; % Neper geometry output file
fname_tesr = 'gene_form_2.tesr'; % Neper RVE mesh output file
radius = 3; % inclusion radius, in element
PF_iso_factor = 3; % anisotropic factor

%% read tess file to get
% vertex, NodeOnEdge, NodeOnFace, EdgeOnFace, FaceOnRegion
fid = fopen(fname_tess,'r');

if(fid==-1)
    error(strcat('>>>Error: ',fname_tess,' does not exist ...'));
end

currLine = fgetl(fid);
if(~strcmp(currLine,'***tess'))
    error('>>> Error: This is not a tess file ...');
end

while (1)
    currLine = fgetl(fid);
    currLine = strtrim(currLine);
    if(strcmp(currLine,'***tess'))
        continue;
    elseif(strcmp(currLine,'**format'))
        currLine = fgetl(fid);
        if(~strcmp(currLine(end-2:end),'2.0'))
            error(strcat('>>>Error: unknown input file format ',currLine));
        end
    elseif(strcmp(currLine,'**general'))
        currLine = fgetl(fid);
        ndim = sscanf(currLine,'%f');
    elseif(strcmp(currLine,'**cell'))
        currLine = fgetl(fid);
        numCry = sscanf(currLine,'%f');
        for temp = 1:(2*numCry+9) % adjust according to tess file or rewrite
            currLine = fgetl(fid);
        end
    elseif(strcmp(currLine,'**vertex'))
        currLine = fgetl(fid);
        num_vertex = sscanf(currLine,'%f');
        temp = textscan(fid,'%f%f%f%f%f');
        vertex = cell2mat(temp);
        if(size(vertex,1) ~= num_vertex)
            error('>>> Error: Invalid number of vertex');
        end
    elseif(strcmp(currLine,'**edge'))
        currLine = fgetl(fid);
        num_edge = sscanf(currLine,'%f');
        temp = textscan(fid,'%f%f%f%f');
        NodeOnEdge = cell2mat(temp);
        if(size(NodeOnEdge,1) ~= num_edge)
            error('>>> Error: Invalid number of edge');
        end
    elseif(strcmp(currLine,'**face'))
        currLine = fgetl(fid);
        num_face = sscanf(currLine,'%f');
        NodeOnFace = zeros(num_face,10);
        EdgeOnFace = zeros(num_face,10);
        face_direction = zeros(num_face,4);
        for ii = 1:num_face
            
            % read vertex
            currLine = fgetl(fid);
            temp = sscanf(currLine,'%f');
            NodeOnFace(ii,1:length(temp)) = temp;
            
            % read edge
            currLine = fgetl(fid);
            temp = sscanf(currLine,'%f');
            EdgeOnFace(ii,1:length(temp)) = temp;
            
            % read normal direction
            currLine = fgetl(fid);
            temp = sscanf(currLine,'%f');
            face_direction(ii,1:4) = temp;
            
            % garbage
            fgetl(fid);
            
        end
    elseif(strcmp(currLine,'**polyhedron'))
        currLine = fgetl(fid);
        num_cell = sscanf(currLine,'%f');
        FaceOnRegion = zeros(num_cell,20);
        if(num_cell ~= numCry)
            error('>>> Error: polyhedron ~= crystal');
        end
        for ii = 1:num_cell
            currLine = fgetl(fid);
            temp = sscanf(currLine,'%f');
            FaceOnRegion(ii,1:length(temp)) = temp;
        end
    elseif(strcmp(currLine,'**domain'))
        warning(strcat('skipping domain part of',32,fname_tess));
        break;
    elseif(strcmp(currLine,'***end'))
        break;
    else
        fclose(fid);
        error(strcat('>>>Error: unknown command ',currLine));
    end
end

fclose(fid);

%% read tesr file to get RegionOnElement
fid = fopen(fname_tesr,'r');

if(fid==-1)
    error(strcat('>>>Error: ',fname_tesr,' does not exist ...'));
end

while (1)
    currLine = fgetl(fid);
    if(strcmp(currLine,'***tesr'))
        continue;
    elseif(strcmp(currLine,' **format'))
        currLine = fgetl(fid);
        if(strcmp(currLine(end-4:end),'ascii'))
            format = 1;
        elseif(strcmp(currLine(end-6:end),'binary8'))
            format = 2;
        elseif(strcmp(currLine(end-7:end),'binary16'))
            format = 3;
        elseif(strcmp(currLine(end-7:end),'binary32'))
            format = 4;
        else
            fclose(fid);
            error(strcat('>>>Error: unknown option ',currLine));
        end
    elseif(strcmp(currLine,' **general'))
        currLine = fgetl(fid);
        ndim = sscanf(currLine,'%f');
        currLine = fgetl(fid);
        numgrid = sscanf(currLine,'%f%f%f');
        numel = prod(numgrid);
        fgetl(fid);
    elseif(strcmp(currLine,' **cell'))
        currLine = fgetl(fid);
        numCry = sscanf(currLine,'%f');
        while(~strcmp(currLine,'  *ori'))
            currLine = fgetl(fid);
        end
        for temp = 1:(numCry+1)
            fgetl(fid);
        end
    elseif(strcmp(currLine,' **data'))
        if(format == 1)
            temp = textscan(fid,'%f');
            RegionOnElement = cell2mat(temp);
        elseif(format == 2)
            temp = fread(fid,numel,'*ubit8');
            RegionOnElement = double(temp);
            fgetl(fid);
        elseif(format == 3)
            temp = fread(fid,numel,'*ubit16');
            RegionOnElement = double(temp);
            fgetl(fid);
        elseif(format == 4)
            temp = fread(fid,numel,'*ubit32');
            RegionOnElement = double(temp);
            fgetl(fid);
        end
    elseif(strcmp(currLine,'***end'))
        break;
    else
        fclose(fid);
        error(strcat('>>>Error: unknown command ',currLine));
    end
end

fclose(fid);

%% remove domain surface from FaceOnRegion
% set FaceOnRegion to 0 for RVE surface
% criteria: all = 0 or all = 1
tol = 1.0e-8;
FaceOnDomain = zeros(num_face,1);
for ii = 1:num_face
    num_node_face = NodeOnFace(ii,2);
    node_list = NodeOnFace(ii,3:num_node_face+2);
    face_coord = sum( vertex(node_list,2:4) );
    FaceOnX = abs( face_coord(1) ) < tol || ...
        abs( face_coord(1)-num_node_face ) < tol;
    FaceOnY = abs( face_coord(2) ) < tol || ...
        abs( face_coord(2)-num_node_face ) < tol;
    FaceOnZ = abs( face_coord(3) ) < tol || ...
        abs( face_coord(3)-num_node_face ) < tol;
    if FaceOnX || FaceOnY || FaceOnZ
        FaceOnDomain(ii) = 1;
    end
end
FaceOnDomain = logical(FaceOnDomain);

for ii = 1:numCry
    num_face_region = FaceOnRegion(ii,2);
    for jj = 3:num_face_region+2
        currFace = abs( FaceOnRegion(ii,jj) );
        if FaceOnDomain(currFace)
            FaceOnRegion(ii,jj) = 0;
        end
    end
end

%% extract grain boundary element and assign normal direction
temp = 1; % current element ID
currNode = zeros(1,3); % element center coordinates
GB_n = zeros(numel,3); % GB normal direction
for kk = 1:numgrid(3) % loop over 3rd direction
    currNode(3) = (kk-0.5) / numgrid(3);
    for jj = 1:numgrid(2) % loop over 2nd direction
        currNode(2) = (jj-0.5) / numgrid(2);
        for ii = 1:numgrid(1) % loop over 1st direction
            currRegion = RegionOnElement(temp); % grain ID of this element
            currNode(1) = (ii-0.5) / numgrid(1);
            num_face_region = FaceOnRegion(currRegion,2);
            face_count = 0; % count if current element is triple junction
            % loop over all surface of this grain
            for ff = 1:num_face_region
                currFace = abs( FaceOnRegion(currRegion,ff+2) );
                if currFace == 0
                    % This is RVE boundary, skip it
                    continue
                end
                % distance between element center and grain boundary
                dist = abs( dot( face_direction(currFace,2:4), currNode ) ...
                    - face_direction(currFace,1) );
                if dist < thick_half
                    % This is a GB element, assign normal to it
                    GB_n(temp,:) = face_direction(currFace,2:4);
                    RegionOnElement(temp) = 0;
                    face_count = face_count + 1;
                end
            end
            % if more than 1 direction is assigned, it's a triple junction
            % triple junction is assumed to be isotropic
            if face_count > 1
                GB_n(temp,:) = 0;
            end
            temp = temp + 1;
        end
    end
end

%% output vtk file for visualization
fid = fopen('Polycrystal_RVE.vtk', 'w');
fprintf(fid, '# vtk DataFile Version 2.0\n');
fprintf(fid, 'VTK from Matlab\n');
fprintf(fid, 'ASCII\n');
precision = '2';
fprintf(fid, 'DATASET STRUCTURED_GRID\n');
fprintf(fid, 'DIMENSIONS %d %d %d\n', numgrid(1)+1, numgrid(2)+1, numgrid(3)+1);
num_node = (numgrid(1)+1)*( numgrid(2)+1)*( numgrid(3)+1);
fprintf(fid, ['POINTS ' num2str(num_node) ' float\n']);

x = 0:numgrid(1);
y = 0:numgrid(2);
z = 0:numgrid(3);
[xx,yy,zz] = ndgrid(x,y,z);
FFT_mesh = [reshape(xx,num_node,1),reshape(yy,num_node,1),reshape(zz,num_node,1)];
fprintf(fid, '%0.2f %0.2f %0.2f\n', transpose(FFT_mesh));

fprintf(fid, ['CELL_DATA ' num2str(numel)]);
fprintf(fid, '\nSCALARS mat_id integer\n');
fprintf(fid, 'LOOKUP_TABLE default\n');
fprintf(fid, '%3i', RegionOnElement);

fprintf(fid, '\nVECTORS GB_normal float\n');
% fprintf(fid, 'LOOKUP_TABLE default\n');
fprintf(fid, '%5.2f %5.2f %5.2f\n', transpose(GB_n));
fprintf(fid,'\n');
fclose(fid);

%% FFT input file
phase_field_0 = genCirc(numgrid, radius);
PF_iso_factor_all = PF_iso_factor * double(RegionOnElement<0.1);
fid = fopen('phase_0','wt');
elemtext = '%5.2f, %5.2f, %10.8f, %10.8f, %10.8f\n';
fprintf(fid,elemtext,[phase_field_0; PF_iso_factor_all'; GB_n']);
fclose(fid);

%% compute initial phase field with circular inclusion
function eta = genCirc(num_grid, radius)
%
% circle in a square matrix
%
c = zeros(1,3);
if mod(num_grid(1),2) % odd
    c(1) = 0.5 * ( num_grid(1) + 1 );
else % even
    c(1) = num_grid(1)*0.5 + 0.5;
end
if mod(num_grid(2),2) % odd
    c(2) = 0.5 * ( num_grid(2) + 1 );
else % even
    c(2) = num_grid(2)*0.5 + 0.5;
end
if mod(num_grid(3),2) % odd
    c(3) = 0.5 * ( num_grid(3) + 1 );
else % even
    c(3) = num_grid(3)*0.5 + 0.5;
end

[xx,yy,zz] = ndgrid(1:num_grid(1),1:num_grid(2),1:num_grid(3));
xx = xx - c(1);
yy = yy - c(2);
zz = zz - c(3);
dist = sqrt(xx.*xx + yy.*yy + zz.*zz);

eta = reshape(double(dist<radius),1,[]);

end
