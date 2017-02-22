function ContinuousSpikeSorter(animal, iseries, iexp, replay_in)

global PICK
global REPLAY

if nargin<4
    REPLAY = 0;
else
    REPLAY = replay_in;
end

PICK.expt = [];
PICK.chans = [];
PICK.picker = 1;
PICK.spikesorter = 3;


PICK.animal = animal;
PICK.iseries = iseries;
PICK.iexp = iexp;


PICK.exptinfos = [];
PICK.protocol = [];

%% Initialize the spike sorter
PICK.spikesorter = FigSpikeSorter;
set(PICK.spikesorter,'HandleVisibility','callback');
set(PICK.spikesorter,'Visible','Off','WindowStyle','Normal','RendererMode','auto','resize','off'); % 'CloseRequestFcn',@DenyInteractiveCloseRequest);

%% set the values of animal, series, and experiment on figSpikeSorter
set( findobj(PICK.spikesorter,'Tag','txtAnimal'), 		'String', PICK.animal );
set( findobj(PICK.spikesorter,'Tag','txtSeries'), 		'String', num2str(PICK.iseries) );
set( findobj(PICK.spikesorter,'Tag','txtExperiment'),	'String', num2str(PICK.iexp) );

set(PICK.spikesorter,'Visible','On');

%%
% set(PICK.spikesorter,'HandleVisibility','On');
% FigSpikeSorter_callbacks LoadCerebusTraces
FigSpikeSorter_callbacks LoadContinuousData
