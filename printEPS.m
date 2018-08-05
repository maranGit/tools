function printEPS(fig,fname,width,height,preview)
% Ran Ma
% 8/5/2018
% print figure in eps format
%
% Change the size of the figure
%
figuresize( fig, width, height, 'cm' );
%
% Save the current figure as an eps file and return
%
if ~preview
    print(fname, 'depsc');
    return
end
%
% Save the current figure as an eps file and add a TIFF preview.
%
print(fname, '-depsc', '-tiff', '-r300');
%
% Extract TIFF preview
%
source = strcat(fname,'.eps');
target = strcat(fname,'.tif');
f = fopen(source,'rb');
d = fread(f,'uint8');
fclose(f);
if ~isequal(d(1:4),[197;208;211;198])
    error 'This does not appear to be a vaild TIFF file.'
end
tiffStart = sum(d(21:24).*256.^(0:3)')+1;
tiffLength = sum(d(25:28).*256.^(0:3)');
f = fopen(target,'w');
fwrite(f,d(tiffStart:tiffStart-1+tiffLength),'uint8','b');
fclose(f);

end

function figuresize( fig , w , h , u )
%FIGURESIZE Set a figure to a specific size
%
% When saving a figure as a PDF, it is necessary to set the
% figure size appropriately. This function sets the "paper size"
% sufficient that the figure is saved with a tight bounding box.
% It will also set the figure size on screen correspondingly.
%
% figuresize(width,height) - sets the figure size in centimeters
% figuresize(width,height,units) - sets the figure size in <units>
%
% <units> can be any of the standard Matlab lengths.
%
% Will Robertson
% 28 April 2010
%
% Copyright and licence information appended.

p = inputParser;
p.addRequired('width', @(x) isnumeric(x) && all(size(x)==1) );
p.addRequired('height',@(x) isnumeric(x) && all(size(x)==1) );
p.addOptional('units','centimeters',...
  @(x) any(strcmpi(x,{'normalized','centimeters','cm','inches','in','points','pt'})) );

p.parse( w, h, u );
w = p.Results.width;
h = p.Results.height;
u = p.Results.units;

switch u
  case 'cm', u = 'centimeters';
  case 'in', u = 'inches';
  case 'pt', u = 'points';
end
  
p = 0.01;

set(fig,'Units',u);
screenpos = get(fig,'Position');

set(fig,...
  'Position',[screenpos(1:2) w h],...
  'PaperUnits',u,...
  'PaperPosition',[p*w p*h w h],...
  'PaperSize',[w*(1+2*p) h*(1+2*p)]);

end


% Copyright (c) 2010 Will Robertson, wspr 81 at gmail dot com
% All rights reserved.
%
% Distributed under the BSD licence in accordance with the wishes of the
% Matlab File Exchange. (Usually I'd pick the Apache License.)
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER ''AS IS'' AND ANY
% EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
% THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.