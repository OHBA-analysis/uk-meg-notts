%
% Computes state lifetimes, interval times, fractional occupancies and switching rate
%
load([dirs.results '/hmm'], 'hmm');
load([dirs.vars '/prepFiles'], 'T');

lifetimes = getStateLifeTimes(hmm.vpath, T, hmm.train, 5);
intervals = getStateIntervalTimes(hmm.vpath, T, hmm.train, 5);
fractionalOccupancies = getFractionalOccupancy (hmm.Gamma, T, hmm.train);
switchingRate = getSwitchingRate(hmm.Gamma, T, hmm.train);

save([dirs.results '/state_statistics'], ...
     'lifetimes', 'intervals', 'fractionalOccupancies', 'switchingRate');

clearvars -except dirs freqRange hmm_options nEmbeddings nStates nSubjectsToDo session;
