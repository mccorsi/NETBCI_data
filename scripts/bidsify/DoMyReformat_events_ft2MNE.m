function DoMyReformat_events_ft2MNE(filename2load, list_events_id)
%% Short function to reformat event files to be readable by MNE afterwards
% Author: Marie-Constance Corsi
% Last modification: July, 22nd 2025

load(strcat(filename2load, '.mat'))

list_events = {event.value}.';
idx = find(ismember(list_events, list_events_id));

tmp_event = event(idx);
event_id = tmp_event(1).value;
nb_trials = length(tmp_event(1).sample);
mat_events = zeros(nb_trials,2);
mat_events(:,1) = tmp_event(1).sample;
save(strcat(filename2load, '_', event_id, '.mat'), 'mat_events', 'event_id');
clearvars mat_event nb_trials event_id


event_id = tmp_event(2).value;
nb_trials = length(tmp_event(2).sample);
mat_events = zeros(nb_trials,2);
mat_events(:,1) = tmp_event(2).sample;
save(strcat(filename2load, '_', event_id, '.mat'), 'mat_events', 'event_id');
clearvars mat_event nb_trials event_id






end
