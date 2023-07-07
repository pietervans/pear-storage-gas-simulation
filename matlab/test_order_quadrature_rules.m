

approximations = zeros(4,1);
for k=1:4
    [points, w] = get_quadrature_points_weights(k);
    approximations(k) = quadrature(@f, points, w);
end

res_matlab = integral2(@f, 0, 1, 0, @(xi) 1-xi, 'RelTol', 1e-16);
errors = abs(approximations-res_matlab)

function res = f(xi,eta)
    res = exp(xi).*exp(eta);
end

function res = quadrature(f, points, w)
    n = size(points,1);
    eval = zeros(n,1);
    for i=1:n
        eval(i) = f(points(i,1),points(i,2));
    end
    res = w'*eval;
end
