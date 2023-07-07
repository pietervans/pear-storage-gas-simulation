function [points, w] = get_quadrature_points_weights(order)
    % From "The Finite Element Method for Engineers" (2001, Wiley), p. 203
    
    pts_standard = [0 0; 1 0; 0 1];
    switch order
        case 1
            points = [1/3 1/3 1/3]*pts_standard;
            w = 1/2*1;
        case 2
            points = [
                1/2 1/2 0;
                0   1/2 1/2;
                1/2 0   1/2
            ]*pts_standard;
            w = 1/2*[1/3; 1/3; 1/3];
        case 3
            points = [
                1/3 1/3 1/3;
                0.6 0.2 0.2;
                0.2 0.6 0.2;
                0.2 0.2 0.6
            ]*pts_standard;
            w = 1/2*[-27/48; 25/48; 25/48; 25/48];
        otherwise
            error('Given order not supported!')
    end
end
