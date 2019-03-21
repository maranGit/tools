% Ran Ma
% 7/4/2018
% read Neper .tesr file

if(~exist('fname','var'))
    fname = 'n10-id1.tesr';
end

fid = fopen(fname,'r');

if(fid==-1)
    error(strcat('>>>Error: ',fname,' does not exist ...'));
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
        elseif(strcmp(currLine(end-6:end),'binary16'))
            format = 3;
        elseif(strcmp(currLine(end-6:end),'binary32'))
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