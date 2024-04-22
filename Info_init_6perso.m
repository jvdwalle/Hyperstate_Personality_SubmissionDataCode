% Number of dimensions
m=3;

% Age classes
pb_w=16; % number of juvenile ages 
w=31; % number of age classes

% Breeding state classes
b=6; % number of breeding states

% Personality classes
g = 6; % number of personality classes

% Other useful information:
first_class = 6; % first PB age when transition to adult reproductive states are possible
rho=0.5; % this is the offspring sex ratio
tau= 0.196; % this is the heritability
intercept = 0;
slope = tau/2; % This is the slope of the mother daughter regression

% Personality classes
min_pheno = -3;
max_pheno = 3;
phenotype = phenotype_class(min_pheno, max_pheno, g); 
PERSONALITY = phenotype;



[Ye Ee]    = discretize(min_pheno:max_pheno,g); % discretize the possible values within the limits and for a given number of bins g
edges_e    = Ee; % these are the edges of the bins for the phenotypes
binwidth_e = edges_e(2) - edges_e(1); % Bin width
midpoint_e = edges_e(1:end-1) +  (binwidth_e/2);




%% Initial parameters useful for sensitivity analyses
s1=w; % size of dimension 1
s2=b; % size of dimension 2
s3=g; % size of dimension 3

siz = [w,b,g];
s = prod(siz);

