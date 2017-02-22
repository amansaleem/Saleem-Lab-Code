function errorarea_as(X, Y, E, colour, alpha_val)
% Usage: errorarea_as(X, Y, E, colour, alpha_val)
if nargin<4
    colour = 'k';
end
if nargin<5
    alpha_val = 0.3;
end
tmp = ~isnan(Y) & ~isnan(E);
[~,PATCHHANDLE] = errorarea(X(tmp), Y(tmp), E(tmp));
hold on;
tmp = ~isnan(Y);
[LINEHANDLE] = plot(X(tmp), Y(tmp));
set(LINEHANDLE,'color',colour,'linewidth',1.5);
set(PATCHHANDLE, 'FaceColor',colour,'FaceAlpha', alpha_val);