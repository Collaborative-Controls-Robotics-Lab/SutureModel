function [internal_vertices, M] = triangles_vertices_delaunay(vertices)

vertices = [vertices ; vertices(1,:)]; %connecting last point to first 

% Create a delaunay triangulation object
DT = delaunayTriangulation(vertices);

% Create a polyshape object to represent the non-convex polygon
polygon = polyshape(vertices(:,1), vertices(:,2));

% Get the centroid of each triangle in the triangulation
tri_centroids = incenter(DT);  % Returns centroids of all triangles

% Check if the centroids are inside the polygon using 'isinterior'
TF = isinterior(polygon, tri_centroids(:,1), tri_centroids(:,2));

% Extract the triangles that are inside the polygon
internal_triangles = DT.ConnectivityList(TF, :);

% Now extract the vertex coordinates of the internal triangles (Mx3x2 matrix)
M = size(internal_triangles, 1); % Number of internal triangles
internal_vertices = zeros(M, 3, 2); % Initialize Mx3x2 matrix

for i = 1:M
    % Get the indices of the vertices for each triangle
    tri_indices = internal_triangles(i, :);
    
    % Store the x and y coordinates of the vertices for this triangle
    internal_vertices(i, :, :) = DT.Points(tri_indices, :);
end
end