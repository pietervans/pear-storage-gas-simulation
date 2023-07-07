% Generation of 2D meshes using MESH2D matlab library
% To install, open Add-on menu in Matlab GUI and search for mesh2d

uniform_edge_length = 2/1000; % unit [m]
write_output = false;

pts = readmatrix('pear_contour.txt'); % List of (r,z) on outer contour of pear

% Scale points to model realistic pear dimensions: 0.03m radius, 0.12m height
range_r = max(pts(:,1)) - min(pts(:,1));
range_z = max(pts(:,2)) - min(pts(:,2));
min_z = min(pts(:,2));
pts = [pts(:,1)*0.03/range_r  (pts(:,2)-min_z)*0.12/range_z];

% Uniform mesh
[vert,~,tria] = refine2(pts,[],[],[],uniform_edge_length);

patch('faces',tria(:,1:3),'vertices',vert, ...
    'facecolor','w', ...
    'edgecolor',[.2,.2,.2]) ;

if (write_output)
    writematrix(vert, 'data/points.txt')
    writematrix(tria, 'data/triangles.txt')
end
