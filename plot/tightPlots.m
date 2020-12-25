function axs = tightPlots(varargin)
% axs(nrow,ncol) = tightPlots(nrow, ncol, width, AR, gap, xedge, yedge, units)

p = inputParser;
addOptional(p, 'nrow', 1, @(x) isnumeric(x) && isscalar(x));
addOptional(p, 'ncol', 1, @(x) isnumeric(x) && isscalar(x));
addOptional(p, 'width', 10, @(x) isnumeric(x) && isscalar(x));
addOptional(p, 'AR', [1, 1], @(x) isnumeric(x));
addOptional(p, 'gap', [1.6, 1.8], @(x) isnumeric(x));
addOptional(p, 'xedge', [1.3, 0.1], @(x) isnumeric(x));
addOptional(p, 'yedge', [1.0, 0.5], @(x) isnumeric(x));
addOptional(p, 'units', 'centimeters', @(x) ischar(x));
addParameter(p, 'preset', '', @(x) ischar(x));
parse(p, varargin{:});

nrow = p.Results.nrow;
ncol = p.Results.ncol;
figwidth = p.Results.width;

AR = p.Results.AR;
if numel(AR) > 1
    AR = AR(2) / AR(1);
end

gap = p.Results.gap;
if numel(gap) < 2
    gap = [gap, gap];
end

xedge = p.Results.xedge;
if numel(xedge) < 2
    xedge = [xedge, xedge];
else
    xedge = xedge(1:2);
end

yedge = p.Results.yedge;
if numel(yedge) < 2
    yedge = [yedge, yedge];
else
    yedge = yedge(1:2);
end

units = p.Results.units;

preset = p.Results.preset;
if strcmp(preset, 'compact')
    gap = [0.2, 0.2];
end

axw = (figwidth - sum(xedge) - (ncol - 1) * gap(1)) / ncol;
axh = axw * AR;

figheight = nrow * axh + (nrow - 1) * gap(2) + sum(yedge);

set(0, 'units', units);
screensize = get(0, 'screensize');
figSize = [screensize(3) / 2 - figwidth / 2, screensize(4) / 2 - figheight / 2, figwidth, figheight];

set(gcf, 'Units', units);
set(gcf, 'Position', figSize);
set(gcf, 'PaperUnits', units, 'PaperSize', [figwidth, figheight]);
set(gcf, 'PaperPositionMode', 'manual', 'PaperPosition', [0, 0, figwidth, figheight]);

py = figheight - yedge(2) - axh;
axs = zeros(nrow, ncol);
for irow = 1:nrow
    px = xedge(1);
    for jcol = 1:ncol
        axs(irow, jcol) = axes;
        set(gca, 'Units', units, 'Position', [px, py, axw, axh]);
        px = px + axw + gap(1);
        if strcmp(preset, 'compact')
            if irow < nrow
                set(gca, 'xticklabel', '');
            end
            if jcol > 1
                set(gca, 'yticklabel', '');
            end
        end
    end
    py = py - axh - gap(2);
end

end
