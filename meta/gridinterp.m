function [Xq, Yq, Vq] = gridinterp(x, y, v)
% [Xq, Yq, Vq] = gridinterp(x, y, v) interp for 3-column data

xq = unique(x, 'sorted');
yq = unique(y, 'sorted');

[xq, yq] = meshgrid(xq, yq);

[Xq, Yq, Vq] = griddata(x, y, v, xq, yq, 'natural');

end
