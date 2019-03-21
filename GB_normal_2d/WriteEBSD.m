% Ran Ma
% 3/2/2019
%
% generate EBSD file

if (~exist('numgrid','var'))
    numgrid = [101,101];
    numel = prod(numgrid);
end

if (~exist('ndim','var'))
    ndim = 2;
end

if (~exist('numCry','var'))
    numCry = 3;
end

if (~exist('RegionOnElement','var'))
    RegionOnElement = ones(1,numel);
    RegionOnElement(101*30:end) = 2;
    RegionOnElement(101*60:end) = 3;
end

if (~exist('fname_ebsd','var'))
    fname_ebsd = 'EBSD_RVE';
end

coordinate_1 = 0 : ( 1 / ( numgrid(1) - 1 ) ) : 1;
coordinate_2 = 0 : ( 1 / ( numgrid(2) - 1 ) ) : 1;

c = zeros(numel,2);

temp = 1;
for jj = 1:numgrid(2)
    for ii = 1:numgrid(1)
        c(temp,1) = coordinate_1(ii);
        c(temp,2) = coordinate_2(jj);
        temp = temp + 1;
    end
end

% random initial orientation
WARP3Dodf = uniformODF(cs);
oriWARP3D = calcOrientations(WARP3Dodf,numCry);
[phi1,Phi,phi2] = Euler(oriWARP3D);
EulAng = [(1:numCry)' phi1*180/pi Phi*180/pi phi2*180/pi];
EAngles = EulAng(:,2:4);
BigList = EAngles(RegionOnElement,1:3); % ix(:,nen1) is the grain ID #

dataAll = zeros(7,numel);
dataAll(1,:) = 1:numel;              % index
dataAll(2,:) = 0;                    % phase
dataAll(3:4,:) = transpose(c);       % coordinate
dataAll(5:7,:) = transpose(BigList); % Euler angle


fid = fopen(fname_ebsd,'wt');
fprintf(fid,'Index, Phase, X, Y, Euler1, Euler2, Euler3\n');
fprintf(fid,'%4i\t%1i\t%10.8f\t%10.8f\t%5.2f\t%5.2f\t%5.2f\n',dataAll);
fclose(fid);