% ContinuousSpikeSorter('M110201_BALL', 110207, 5)
% ContinuousSpikeSorter('M110201_BALL', 110207, 4)
% ContinuousSpikeSorter('M110201_BALL', 110207, 1, 1)
% LoadContinuousVRData('M110201_BALL', 110207, 4, 0, 0) - Ignore, useless
% LoadContinuousVRData('M110201_BALL', 110207, 5, 0, 0)
% LoadContinuousVRData('M110201_BALL', 110207, 1, 1, 0)
% % % % % Photodiode error on original % % % % % % %


ContinuousSpikeSorter('M110201_BALL', 110208, 4)
ContinuousSpikeSorter('M110201_BALL', 110208, 1, 1)
ContinuousSpikeSorter('M110201_BALL', 110208, 2, 1)
%
LoadContinuousVRData('M110201_BALL', 110208, 4, 0, 0)
LoadContinuousVRData('M110201_BALL', 110208, 1, 1, 0)
LoadContinuousVRData('M110201_BALL', 110208, 2, 1, 0)


% ContinuousSpikeSorter('M110201_BALL', 110209, 2)
% ContinuousSpikeSorter('M110201_BALL', 110209, 1, 1)
% ContinuousSpikeSorter('M110201_BALL', 110209, 2, 1)
%


% ContinuousSpikeSorter('M110201_BALL', 110210, 1)
% ContinuousSpikeSorter('M110201_BALL', 110210, 2)
% ContinuousSpikeSorter('M110201_BALL', 110210, 1, 1)
% ContinuousSpikeSorter('M110201_BALL', 110210, 2, 1)
%
% LoadContinuousVRData('M110201_BALL', 110210, 1, 0, 0)
% LoadContinuousVRData('M110201_BALL', 110210, 2, 0, 0)
% LoadContinuousVRData('M110201_BALL', 110210, 1, 1, 0)
% LoadContinuousVRData('M110201_BALL', 110210, 2, 1, 0)


% ContinuousSpikeSorter('M110201_BALL', 110212, 1)
% ContinuousSpikeSorter('M110201_BALL', 110212, 1, 1)
%
LoadContinuousVRData('M110201_BALL', 110209, 2, 0, 0)
LoadContinuousVRData('M110201_BALL', 110209, 1, 1, 0)
LoadContinuousVRData('M110201_BALL', 110209, 2, 1, 0)

LoadContinuousVRData('M110201_BALL', 110212, 1, 0, 0)
LoadContinuousVRData('M110201_BALL', 110212, 1, 1, 0)


ContinuousSpikeSorter('M110209_BALL', 110217, 1)
ContinuousSpikeSorter('M110209_BALL', 110217, 2, 1)

LoadContinuousVRData('M110209_BALL', 110217, 1, 0, 0)
LoadContinuousVRData('M110209_BALL', 110217, 2, 1, 0)



[M110201_BALL_110211_1, M110201_BALL_110211_1_I, M110201_BALL_110211_1_R] = getVRPlaceMaps('M110201_BALL', 110211, 1, 40);
[M110201_BALL_110211_R1, M110201_BALL_110211_R1_I, M110201_BALL_110211_R1_R] = getVRPlaceMaps('M110201_BALL', 110211, 1, 40, 1);
[M110201_BALL_110211_R2, M110201_BALL_110211_R2_I, M110201_BALL_110211_R2_R] = getVRPlaceMaps('M110201_BALL', 110211, 2, 40, 1);
close all
[M110201_BALL_110210_1, M110201_BALL_110210_1_I, M110201_BALL_110210_1_R] = getVRPlaceMaps('M110201_BALL', 110210, 1, 40);
[M110201_BALL_110210_2, M110201_BALL_110210_2_I, M110201_BALL_110210_2_R] = getVRPlaceMaps('M110201_BALL', 110210, 2, 40);
[M110201_BALL_110210_R1, M110201_BALL_110210_R1_I, M110201_BALL_110210_R1_R] = getVRPlaceMaps('M110201_BALL', 110210, 1, 40, 1);
[M110201_BALL_110210_R2, M110201_BALL_110210_R2_I, M110201_BALL_110210_R2_R] = getVRPlaceMaps('M110201_BALL', 110210, 2, 40, 1);
close all
[M110201_BALL_110209_1, M110201_BALL_110209_1_I, M110201_BALL_110209_1_R] = getVRPlaceMaps('M110201_BALL', 110209, 2, 40);
[M110201_BALL_110209_R1, M110201_BALL_110209_R1_I, M110201_BALL_110209_R1_R] = getVRPlaceMaps('M110201_BALL', 110209, 1, 40, 1);
[M110201_BALL_110209_R2, M110201_BALL_110209_R2_I, M110201_BALL_110209_R2_R] = getVRPlaceMaps('M110201_BALL', 110209, 2, 40, 1);
close all
[M110201_BALL_110212_1, M110201_BALL_110212_1_I, M110201_BALL_110212_1_R] = getVRPlaceMaps('M110201_BALL', 110212, 1, 40);
[M110201_BALL_110212_R1, M110201_BALL_110212_R1_I, M110201_BALL_110212_R1_R] = getVRPlaceMaps('M110201_BALL', 110212, 1, 40, 1);
close all
[M110209_BALL_110217_1, M110209_BALL_110217_1_I, M110209_BALL_110217_1_R] = getVRPlaceMaps('M110209_BALL', 110217, 1, 40);
[M110209_BALL_110217_R2, M110209_BALL_110217_R2_I, M110209_BALL_110217_R2_R] = getVRPlaceMaps('M110209_BALL', 110217, 2, 40, 1);
% close all

clear DEMO DIRS ans xlims ylims PICK

save VRdataparam

figure;
hold on;
plot(M110201_BALL_110210_1_I, M110201_BALL_110210_R1_I, '.k', 'markersize', 30);
% plot(M110201_BALL_110210_1_I, M110201_BALL_110210_R2_I, '.b', 'markersize', 30);

plot(M110201_BALL_110211_1_I, M110201_BALL_110211_R1_I, '.g', 'markersize', 30);
% plot(M110201_BALL_110211_1_I, M110201_BALL_110211_R2_I, '.r', 'markersize', 30);

plot(M110209_BALL_110217_1_I, M110209_BALL_110217_R1_I, '.m', 'markersize', 30);

xlims = get(gca,'XLim'); ylims = get(gca, 'YLim');
axis([min(xlims(1),ylims(1)), max(xlims(2),ylims(2)), min(xlims(1),ylims(1)), max(xlims(2),ylims(2))]); axis square
line([min(xlims(1),ylims(1)), max(xlims(2),ylims(2))], [min(xlims(1),ylims(1)), max(xlims(2),ylims(2))])
xlabel('I_{Active}', 'fontsize', 24)
ylabel('I_{Passive}', 'fontsize', 24)

figure;
hold on;
plot(1./M110201_BALL_110210_1_R, 1./M110201_BALL_110210_R1_R, '.k', 'markersize', 30);
% plot(1./M110201_BALL_110210_1_R, 1./M110201_BALL_110210_R2_R, '.b', 'markersize', 30);

plot(1./M110201_BALL_110211_1_R, 1./M110201_BALL_110211_R1_R, '.g', 'markersize', 30);
% plot(1./M110201_BALL_110211_1_R, 1./M110201_BALL_110211_R2_R, '.r', 'markersize', 30);
plot(1./M110209_BALL_110217_1_R, 1./M110209_BALL_110217_R1_R, '.m', 'markersize', 30);

xlims = get(gca,'XLim'); ylims = get(gca, 'YLim');
axis([min(xlims(1),ylims(1)), max(xlims(2),ylims(2)), min(xlims(1),ylims(1)), max(xlims(2),ylims(2))]); axis square
line([min(xlims(1),ylims(1)), max(xlims(2),ylims(2))], [min(xlims(1),ylims(1)), max(xlims(2),ylims(2))])
xlabel('R_{Active}', 'fontsize', 24)
ylabel('R_{Passive}', 'fontsize', 24)