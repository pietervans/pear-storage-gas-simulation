function TR = get_triangulation(type)
% Return triangulation of pear. For info about returned data structure:
% https://nl.mathworks.com/help/matlab/ref/triangulation.html

% type == 'adaptive_[fine|rough]' or type == 'uniform_[5|3|1|0p5|0p25]mm']
pts =           readmatrix(strcat('pear_meshes/data/points_',   type,'.txt'));
triangle_conn = readmatrix(strcat('pear_meshes/data/triangles_',type,'.txt'));

TR0 = triangulation(triangle_conn, pts);
for i = 1:size(triangle_conn,1)
    n = faceNormal(TR0, i); % equal to [0 0 +/-1] --> cannot have -1
    if n(3)==-1
        triangle_conn(i,:) = flip(triangle_conn(i,:));
    end
end


TR = triangulation(triangle_conn, pts);
end
