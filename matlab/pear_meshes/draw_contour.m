figure
imshow('pear_cropped.jpg')
[x,y] = ginput;
r = x;
z = -y;
% plot(r,z)

rz = [r z];
writematrix(rz, 'pear_contour_new.txt')
