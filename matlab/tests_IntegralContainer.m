
rz = [
    0 0;
    1 0;
    0 1;
];

test_associativity(@(i,j,rz) IntegralContainer.integral2d1(i,j,rz,1,1), rz, 3)
test_associativity(@IntegralContainer.integral2d2, rz, 3)
test_associativity(@IntegralContainer.integral1d2, rz, 2)

% Result should be invariant with respect to translation in z-direction
test_shift_invariance2d(@(i,j,rz) IntegralContainer.integral2d1(i,j,rz,1,1), rz, 3)
test_shift_invariance2d(@IntegralContainer.integral2d2, rz, 3)
test_shift_invariance2d(@IntegralContainer.integral1d2, rz, 2)
test_shift_invariance1d(@IntegralContainer.integral2d3, rz, 3)
test_shift_invariance1d(@IntegralContainer.integral1d1, rz, 2)

test_correctness(rz)

function test_associativity(integral_func, rz, max_index)
    for i=1:max_index
        for j=1:max_index
            res1 = integral_func(i,j,rz);
            res2 = integral_func(j,i,rz);
            assert(abs(res1-res2) < 1e-7)
        end
    end
end


function test_shift_invariance2d(integral_func, rz, max_index)
    rz_shifted = rz;
    rz_shifted(:,2) = rz_shifted(:,2)+0.42;
    for i=1:max_index
        for j=1:max_index
            res0 = integral_func(i,j,rz);
            res1 = integral_func(i,j,rz_shifted);
            assert(abs(res0-res1) < 1e-7);
        end
    end
end

function test_shift_invariance1d(integral_func, rz, max_index)
    rz_shifted = rz;
    rz_shifted(:,2) = rz_shifted(:,2)+0.42;
    for j=1:max_index
        res0 = integral_func(j,rz);
        res1 = integral_func(j,rz_shifted);
        assert(abs(res0-res1) < 1e-7);
    end
end

function test_correctness(rz)
    % Verify results for some parameters and reference triangle
    res_fun = IntegralContainer.integral2d1(1,1,rz,1,1);
    res_matlab = integral2(@(r,z) r*2, 0, 1, 0, @(z) 1-z);
    assert(abs(res_fun-res_matlab) < 1e-7)
    res_fun = IntegralContainer.integral2d1(1,2,rz,1,1);
    res_matlab = integral2(@(r,z) r*(-1), 0, 1, 0, @(z) 1-z);
    assert(abs(res_fun-res_matlab) < 1e-7)


    res_fun = IntegralContainer.integral2d2(1,1,rz);
    res_matlab = integral2(@(r,z) r.*(1-r-z).^2, 0, 1, 0, @(z) 1-z);
    assert(abs(res_fun-res_matlab) < 1e-7)
    res_fun = IntegralContainer.integral2d2(1,2,rz);
    res_matlab = integral2(@(r,z) r.*(1-r-z).*r, 0, 1, 0, @(z) 1-z);
    assert(abs(res_fun-res_matlab) < 1e-7)

    res_fun = IntegralContainer.integral2d3(1,rz);
    res_matlab = integral2(@(r,z) r.*(1-r-z), 0, 1, 0, @(z) 1-z);
    assert(abs(res_fun-res_matlab) < 1e-7)
    res_fun = IntegralContainer.integral2d3(2,rz);
    res_matlab = integral2(@(r,z) r.*r, 0, 1, 0, @(z) 1-z);
    assert(abs(res_fun-res_matlab) < 1e-7)

    res_fun = IntegralContainer.integral1d1(1,rz);
    res_matlab = integral(@(r) r.*(1-r), 0, 1);
    assert(abs(res_fun-res_matlab) < 1e-7)
    res_fun = IntegralContainer.integral1d1(2,rz);
    res_matlab = integral(@(r) r.*r, 0, 1);
    assert(abs(res_fun-res_matlab) < 1e-7)

    res_fun = IntegralContainer.integral1d2(1,1,rz);
    res_matlab = integral(@(r) r.*(1-r).^2, 0, 1);
    assert(abs(res_fun-res_matlab) < 1e-7)
    res_fun = IntegralContainer.integral1d2(1,2,rz);
    res_matlab = integral(@(r) r.*(1-r).*r, 0, 1);
    assert(abs(res_fun-res_matlab) < 1e-7)
end
