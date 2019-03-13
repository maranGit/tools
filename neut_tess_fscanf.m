% Ran Ma
% 3/8/2019
% read Neper .tess file

if(~exist('fname','var'))
    fname = 'n10-id1.tess';
end

fid = fopen(fname,'r');

if(fid==-1)
    error(strcat('>>>Error: ',fname,' does not exist ...'));
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
        for temp = 1:(2*numCry+9)
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
        warning(strcat('skipping domain part of ',fname));
        break;
    elseif(strcmp(currLine,'***end'))
        break;
    else
        fclose(fid);
        error(strcat('>>>Error: unknown command ',currLine));
    end
end

fclose(fid);
