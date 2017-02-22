function es = loadESwithTheta(animal, iseries, exp_list, shank_list, area)

if nargin<5
    area = 'CA1';
end
switch animal
    case 'M130918_BALL'
        iseries = 1030;
        switch area
            case 'V1'
                cell_list = 190:321;
            case 'CA1'
                cell_list = 1:189;
        end
    case 'M130920_BALL'
        iseries = 1025;
        switch area
            case 'V1'
                cell_list = 153:249;
            case 'CA1'
                cell_list = 1:152;
        end
    otherwise
        cell_list = [];
end

%% Load the spiking and theta data
SetDirs;
es = VRLoadMultipleExpts(animal, iseries, exp_list,'SPIKES_THETA',[18 22],shank_list);
if ~isempty(cell_list)
    es.spikeTrain = es.spikeTrain(:,cell_list);
end