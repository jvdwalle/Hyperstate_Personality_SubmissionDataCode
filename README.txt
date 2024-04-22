README

# Description of the datasets and code
These are datasets and Matlab code accompanying the manuscript entitled: Hyperstate matrix model reveals the influence of personality on demography, submitted to Ecological Monographs


# Description of the files
# Data:
WA_5Fs_a_mean_2008_SCALE.mat ; Dataset for boldness-specific estimates of survival for adult females
WA_5Fb_a_mean_2008_SCALE.mat ; Dataset for boldness-specific estimates of breeding probability for adult females
WA_5Fbs_a_mean_2008_SCALE.mat ; Dataset for boldness-specific estimates of breeding success probability for adult females
WA_5Ms_a_mean_2008_SCALE.mat ; Dataset for boldness-specific estimates of survival for adult males
WA_5Mb_a_mean_2008_SCALE.mat ; Dataset for boldness-specific estimates of breeding probability for adult males
WA_5Mbs_a_mean_2008_SCALE.mat ; Dataset for boldness-specific estimates of breeding success probability for adult males
F_juv_mean_cohortremoved.mat ; Dataset for female juvenile vital rates
M_juv_mean_cohortremoved.mat ; Dataset for female juvenile vital rates
male_score.txt; boldness scores of males
female_score.txt; boldness scores of females

# Code:
Main_Hyper.mlx ; This is the main file to follow the different steps

BD_proj_mat.m; Function to build the block matrices. From Roth and Caswell (2016)
Emat.m; Function to build the E matrix. From Roth and Caswell (2016)
get_init_dist.m; Function to build an initial vector of population
hyper_state_matrix.m; Function to build the hyperstate matrix. From Roth and Caswell (2016)
Imat.m; Function tu build the I matrix. From Roth and Caswell (2016)
Info_init_6perso.m; Code for generating initial parameters
invlogit.m; invert logit function
logit.m; logit funtion
parameter_personalityF.m; Function to generate the personality-specific parameters to go in the population model for females.
parameter_personalityM.m; Function to generate the personality-specific parameters to go in the population model for males.
parameters.m; Code to load the data
phenotype_class.m; Function to get the mid_point of the phenotypic classes
popmat_hyper.m; Function to generate the dimension-specific matrices
Qmat.m; Function to build the Q matrix. From Roth and Caswell (2016)

sensitivity_lambdaSigmaF_hyper_6perso.m; Code to run the sensitivity analyses on survival for females
sensitivity_lambdaBeta_hyper_F_6perso.m; Code to run the sensitivity analyses on breeding probability for females
sensitivity_lambdaGamma_hyper_F_6perso.m; Code to run the sensitivity analyses on breeding success for females

sensitivity_lambdaSigmaM_hyper_6perso.m; Code to run the sensitivity analyses on survival for males
sensitivity_lambdaBeta_hyper_M_6perso.m; Code to run the sensitivity analyses on breeding probability for males
sensitivity_lambdaGamma_hyper_M_6perso.m; Code to run the sensitivity analyses on breeding success for males

sensitivity_lambdaSigmaF_hyper_6perso_Uniform.m; Code to run the sensitivity analyses on survival for females when boldness distribution is uniform
sensitivity_lambdaBeta_hyper_F_6perso_Uniform.m; Code to run the sensitivity analyses on breeding probability for females when boldness distribution is uniform
sensitivity_lambdaGamma_hyper_F_6perso_Uniform.m; Code to run the sensitivity analyses on breeding success for females when boldness distribution is uniform

sensitivity_lambdaSigmaF_hyper_6perso_Extreme.m; Code to run the sensitivity analyses on survival for females under extreme cases of selection and heritability
sensitivity_lambdaBeta_hyper_F_6perso_Extreme.m; Code to run the sensitivity analyses on breeding probability for females under extreme cases of selection and heritability
sensitivity_lambdaGamma_hyper_F_6perso_Extreme.m; Code to run the sensitivity analyses on breeding success for females under extreme cases of selection and heritability

Uniform_SSD_Creation.R; Code to generate a uniform distribution of boldness
vecperm_hyp.m; Function to build the vec permutation matrix. From Roth and Caswell (2016)



