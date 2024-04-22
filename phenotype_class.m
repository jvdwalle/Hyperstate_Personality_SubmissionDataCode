function phenotype = phenotype_class(min, max, nclass)
% Function to get the mid_point of the phenotypic classes
[Y,E] = discretize(min:max,nclass); % discretize the possible values within the limits and for a given number of bins g
edges_p=E; % these are the edges of the bins
binwidth_p = edges_p(2) - edges_p(1); % Bin width
phenotype = edges_p(1:end-1) +  (binwidth_p/2); % Mid phenotype in each bin
end

