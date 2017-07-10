function [coords labels] = data_generator(numcenters,numobservations)
% this code generates three dimensional clouds of datapoints with a set
% number of centers and 

% INPUT
%   numcenters - integer value that represents the number of centroids to
%               be generated
%   numobservations - integer value that represents the number of
%               datapoints within each centroid
% OUTPUT
%   coords - generated (x,y,z) data
%   labels - centroid number that each (x,y,z) tuple comes from

% TODO
% -make numobservations accept a vector of values so different centroids
%    can have different numbers of points...

if nargin < 2
    numobservations = 100;
end
x = zeros(numcenters,1);
y = zeros(numcenters,1);

while any(abs(diff(x)) < 3) || any(abs(diff(y)) < 3)
x = rand(numcenters,1)*15;
y = rand(numcenters,1)*15;
z = rand(numcenters,1)*15;
end

coords = [];
labels = [];

for i = 1:numcenters

    
mu = [x(i),y(i),z(i)];


 Sigma = rand(3);
 [U,ignore] = eig((Sigma+Sigma')/4); % (Don't really have to divide by 4)
 Sigma = U*diag(abs(randn(3,1)))*U'./4;
 

R = chol(Sigma);
coords =[coords; repmat(mu,numobservations,1) + randn(numobservations,3)*R;];
labels = [labels; repmat(i,numobservations,1)];
end
r = randperm(size(coords,1));
coords = coords(r,:);
labels = labels(r);
end
