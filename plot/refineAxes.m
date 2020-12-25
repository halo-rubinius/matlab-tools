function refineAxes(axs)
% refineAxes(axs)

fontSize = 10;

if ~exist('axs', 'var') || isempty(axs)
    axs = gca;
end

set(axs, 'LineWidth', 1, 'TickDir', 'out' ...
    , 'FontName', 'Times', 'FontSize', fontSize ...
    , 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');

end
