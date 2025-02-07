function [closestPoints, minDistances] = closestPointOnTriangle_vectorized(p, vertices)

N = size(p, 2); % P: 2xN matrix of points
M = size(vertices, 1); % vertices: Mx3x2 matrix of triangle vertices

% Extract vertices of the current triangle
A = vertices(:, 1, :); % Mx1x2 matrix
B = vertices(:, 2, :);
C = vertices(:, 3, :);

% Reshape P to 1xNx2 so it can be broadcasted over M triangles
P = permute(p, [3, 2, 1]); % 1xNx2 matrix

% Precompute edge vectors
AB = B - A; % Mx1x2
BC = C - B; % Mx1x2
CA = A - C; % Mx1x2

% Precompute dot products
dotAB = sum(AB .* AB, 3); % Mx1 matrix
dotBC = sum(BC .* BC, 3);
dotCA = sum(CA .* CA, 3);

% Precompute vectors from P to each vertex (now MxNx2 matrices)
AP = P - A; % MxNx2 matrix
BP = P - B;
CP = P - C;

% Project AP onto AB, BP onto BC, and CP onto CA
tAB = max(0, min(1, sum(AP .* AB, 3) ./ dotAB)); % MxN matrix
tBC = max(0, min(1, sum(BP .* BC, 3) ./ dotBC));
tCA = max(0, min(1, sum(CP .* CA, 3) ./ dotCA));

% Calculate closest points on the edges
closestPointAB = A + tAB .* AB; % MxNx2 matrix
closestPointBC = B + tBC .* BC;
closestPointCA = C + tCA .* CA;
closestPoint_tot = [closestPointAB; closestPointBC; closestPointCA];

% % Compute squared distances to avoid sqrt until necessary
% dist2AB = sum((P - closestPointAB).^2, 3); % MxN matrix
% dist2BC = sum((P - closestPointBC).^2, 3);
% dist2CA = sum((P - closestPointCA).^2, 3);
% % Stack the distances for easier comparison
% d=[dist2AB; dist2BC; dist2CA]; %3MXN

d = [sum((P - closestPointAB).^2, 3); sum((P - closestPointBC).^2, 3); sum((P - closestPointCA).^2, 3)];

[minDistances, min_Id] = min(d);      %1XN 
minDistances = sqrt(minDistances);

% The closest points
closestPoints = zeros(2,N);

for i=1:N
    closestPoints(:,i) = closestPoint_tot(min_Id(i),i,:);
end


end