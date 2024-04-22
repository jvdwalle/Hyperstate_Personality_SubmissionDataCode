function SSD = get_init_dist(g, data)

% This function allows calculating the proportion of the population within each class

% Create classes
T = array2table(data); % transform into a dataframe

[T.class, E] = discretize(data, g);

space = zeros(g,1);
mid_point = zeros(g,1);

for i=1:g
    space(i) = E(i+1) - E(i);
    mid_point(i) = E(i) + space(i)/2;
end

SSD = zeros(g,3); % col1 = class number (from 1 to 50), col2 = mid_point, col3 = number of individuals within each class
SSD(:,1) = 1:1:g;
SSD(:,2) = mid_point;

for i=1:g
    Try = find(T.class == i);
    SSD(i,3) = length(Try);
end

% Put it into proportions
n_id = length(data);
SSD(:,3) = SSD(:,3)/n_id;
sum(SSD(:,3)); % Needs to be 1
end

