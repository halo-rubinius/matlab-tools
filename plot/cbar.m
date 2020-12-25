function hbar = cbar(axs, location, label)
% hbar = cbar(axs, location, label) create colobar

if ~exist('axs', 'var') || isempty(axs)
    axs = gca;
end
if ~exist('location', 'var') || isempty(location)
    location = 'eastoutside';
end
if ~ischar(location) || ~any(strcmp(location, {'eastoutside', 'southoutside'}))
    location = 'eastoutside';
end
if ~exist('label', 'var') || isempty(label)
    label = '';
end

fontsize = 10;


n = numel(axs);
pos = zeros(n, 4);
units = cell(n, 1);

for ii = 1:n
    units{ii} = get(axs(ii), 'units');
    set(axs(ii), 'units', 'centimeters');
    pos(ii, :) = get(axs(ii), 'position');
end

left = min(pos(:, 1));
right = max(pos(:, 1)+pos(:, 3));
bottom = min(pos(:, 2));
top = max(pos(:, 2)+pos(:, 4));

barwidth = 0.35;
if strcmp(location, 'eastoutside')
    [~, idx] = max(pos(:, 1)+pos(:, 3));
    ti = get(axs(idx), 'TightInset');
    p = [pos(idx, 1) + pos(idx, 3) + ti(3) + 0.1, bottom, ...
        barwidth, top - bottom];
else
    [~, idx] = min(pos(:, 2));
    ti = get(axs(idx), 'TightInset');
    p = [left, bottom - ti(2) - 0.8, right - left, barwidth];
end
hbar = colorbar(axs(idx), 'fontsize', fontsize, ...
    'units', 'centimeters', 'location', location, ...
    'position', p);

if strcmp(location, 'eastoutside')
    if ~isempty(label)
        ylabel(hbar, label, 'fontsize', fontsize);
    end
else
    if ~isempty(label)
        xlabel(hbar, label, 'fontsize', fontsize);
    end
end

for ii = 1:n
    set(axs(ii), 'units', units{ii});
end

end
