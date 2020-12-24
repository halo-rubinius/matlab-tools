function savefigure(filename, fig)
% savefigure(filename, fig)

if ~exist('fig', 'var')
    fig = gcf;
end

if ~isempty(findall(fig, 'Type', 'surface')) ...
        || ~isempty(findall(fig, 'Type', 'contour'))
    resolution = '-r600';
else
    resolution = '-r1200';
end

savefig(fig, filename); % .fig
print(fig, filename, '-dpng', resolution, '-opengl'); % .png
print(fig, filename, '-depsc', resolution, '-painters', '-loose'); % .eps
print(fig, filename, '-dmeta', resolution, '-painters'); % .emf

disp(['Figure saved: ', filename]);

end
