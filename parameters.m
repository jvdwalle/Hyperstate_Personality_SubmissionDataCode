% Initial parameters
% For adults: 
% These models consider only adults from 2008 with personality data
% Def: WA = Wandering albatross, 5 = model version 5, F = Female, s = survival, b = breeding probability, bs = breeding success, a = adult, mean = only mean parameters, 2008 = only data from
% 2008, SCALE = personality data standardized (mean = 0, sd = 1)
load("WA_5Fs_a_mean_2008_SCALE.mat"); 
load("WA_5Fb_a_mean_2008_SCALE.mat");
load("WA_5Fbs_a_mean_2008_SCALE.mat");
load("WA_5Ms_a_mean_2008_SCALE.mat");
load("WA_5Mb_a_mean_2008_SCALE.mat");
load("WA_5Mbs_a_mean_2008_SCALE.mat");


% Juveniles parameters
% These models consider only juveniles (pre-recruited individuals) from age 1 to 16 and removing the last cohorts to avoid bias due to the low detection probability of juveniles
F_juv_mean = load("F_juv_mean_cohortremoved.mat");
F_juv_mean = F_juv_mean.F_juv_mean;

M_juv_mean = load("M_juv_mean_cohortremoved.mat");
M_juv_mean = M_juv_mean.M_juv_mean;
