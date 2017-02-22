function hax_new = captureSubplot(subplot_array, fignum)

if nargin>0
    figure;
    subplot(subplot_array(1), subplot_array(2), subplot_array(3));
    pos = get(gca, 'Position');
    specPos=1;
    close
else
    specPos=0;
end

    

hax_old = gca;
if nargin<2
    newfig = figure;
else
    newfig = figure(fignum);
end

hax_new = copyobj(hax_old, newfig);

if specPos
    set(hax_new, 'Position', pos);
else
    set(hax_new, 'Position', get(0, 'DefaultAxesPosition'));
end