% This file contains functions to compute the 5 integrals
% described in theory_formulas.md

classdef IntegralContainer
    methods(Static)
        % (r,z) coordinates of triangle --> order matters for 2d integrals:
        % "flipped" triangle can result in negation of result!
        % Order of points in rz should be counter-clockwise!

        function res = integral2d1(i,j,rz,sigma_r,sigma_z)
            assert(i==1 || i==2 || i==3)
            assert(j==1 || j==2 || j==3)
            r = rz(:,1);
            det_J = det_jacobian(rz);
            res = sum(r)/6*(sigma_r*partial_r_phi(i,rz)*partial_r_phi(j,rz) ...
                + sigma_z*partial_z_phi(i,rz)*partial_z_phi(j,rz))*det_J;
        end

        function res = integral2d2(i,j,rz)
            assert(i==1 || i==2 || i==3)
            assert(j==1 || j==2 || j==3)
            r = rz(:,1);
            det_J = det_jacobian(rz);
            if (i==j)
                res = 1/60*(sum(r)+2*r(i))*det_J;
            else
                res = 1/120*(sum(r)+r(i)+r(j))*det_J;
            end
        end

        function res = integral2d3(j,rz)
            assert(j==1 || j==2 || j==3)
            r = rz(:,1);
            det_J = det_jacobian(rz);
            res = (sum(r)+r(j))*det_J/24;
        end

        function res = integral1d1(j,rz)
            % Edge outer points must be first 2 points in rz!
            assert(j==1 || j==2)
            r = rz(:,1);
            norm_E = norm(rz(1,:)-rz(2,:));
            res = (r(1)+r(2)+r(j))/6*norm_E;
        end

        function res = integral1d2(i,j,rz)
            % Edge outer points must be first 2 points in rz!
            assert(i==1 || i==2)
            assert(j==1 || j==2)
            r = rz(:,1);
            norm_E = norm(rz(1,:)-rz(2,:));
            if (i==j)
                res = 1/12*(r(1)+r(2)+2*r(i))*norm_E;
            else
                res = 1/12*(r(1)+r(2))*norm_E;
            end
        end

        function res = compute_alpha(i,m,tr,triangle_ind,c,CONSTS)
            % Inputs:
            %  - i,m: integers with 1 <= i,m <= 3
            %  - tr: triangulation object
            %  - triangle_ind: integer
            %  - c: array with 2M coefficients
            %  - CONSTS: struct
            [M, inds, r, det_J] = get_values_integrals_greek(tr,triangle_ind);
            Ni = str2func(strcat('N',int2str(i)));
            Nm = str2func(strcat('N',int2str(m)));
            integrand = @(xi,eta) ...
                (r(1)*N1(xi,eta)+r(2)*N2(xi,eta)+r(3)*N3(xi,eta))*Ni(xi,eta)*Nm(xi,eta) ...
                / (1+1/CONSTS.K_mv*(c(M+inds(1))*N1(xi,eta)+c(M+inds(2))*N2(xi,eta)+c(M+inds(3))*N3(xi,eta))) ...
                / (CONSTS.K_mu+c(inds(1))*N1(xi,eta)+c(inds(2))*N2(xi,eta)+c(inds(3))*N3(xi,eta))^2;
            res = CONSTS.V_mu*CONSTS.K_mu*det_J*quadrature_standard_triangle(integrand);
        end

        function res = compute_beta(i,m,tr,triangle_ind,c,CONSTS)
           [M, inds, r, det_J] = get_values_integrals_greek(tr,triangle_ind);
            Ni = str2func(strcat('N',int2str(i)));
            Nm = str2func(strcat('N',int2str(m)));
            integrand = @(xi,eta) ...
                (r(1)*N1(xi,eta)+r(2)*N2(xi,eta)+r(3)*N3(xi,eta))*(c(inds(1))*N1(xi,eta)+c(inds(2))*N2(xi,eta)+c(inds(3))*N3(xi,eta))*Ni(xi,eta)*Nm(xi,eta) ...
                / (1+1/CONSTS.K_mv*(c(M+inds(1))*N1(xi,eta)+c(M+inds(2))*N2(xi,eta)+c(M+inds(3))*N3(xi,eta)))^2 ...
                / (CONSTS.K_mu+c(inds(1))*N1(xi,eta)+c(inds(2))*N2(xi,eta)+c(inds(3))*N3(xi,eta));
            res = -CONSTS.V_mu/CONSTS.K_mv*det_J*quadrature_standard_triangle(integrand);
        end

        function res = compute_gamma(i,m,tr,triangle_ind,c,CONSTS)
            [~, inds, r, det_J] = get_values_integrals_greek(tr,triangle_ind);
            Ni = str2func(strcat('N',int2str(i)));
            Nm = str2func(strcat('N',int2str(m)));
            integrand = @(xi,eta) ...
                (r(1)*N1(xi,eta)+r(2)*N2(xi,eta)+r(3)*N3(xi,eta))*Ni(xi,eta)*Nm(xi,eta) ...
                / (1+1/CONSTS.K_mfu*(c(inds(1))*N1(xi,eta)+c(inds(2))*N2(xi,eta)+c(inds(3))*N3(xi,eta)))^2 ;
            res = CONSTS.V_mfv/CONSTS.K_mfu*det_J*quadrature_standard_triangle(integrand);
        end

        function res = compute_delta(m,tr,triangle_ind,c,CONSTS)
            % Inputs:
            %  - m: integers with 1 <= m <= 3
            %  - tr: triangulation object
            %  - triangle_ind: integer
            %  - c: array with 2M coefficients
            %  - CONSTS: struct
            [M, inds, r, det_J] = get_values_integrals_greek(tr,triangle_ind);
            Nm = str2func(strcat('N',int2str(m)));
            integrand = @(xi,eta) ...
                (r(1)*N1(xi,eta)+r(2)*N2(xi,eta)+r(3)*N3(xi,eta)) ...
                *(c(inds(1))*N1(xi,eta)+c(inds(2))*N2(xi,eta)+c(inds(3))*N3(xi,eta))*Nm(xi,eta) ...
                / (1+1/CONSTS.K_mv*(c(M+inds(1))*N1(xi,eta)+c(M+inds(2))*N2(xi,eta)+c(M+inds(3))*N3(xi,eta))) ...
                / (CONSTS.K_mu+c(inds(1))*N1(xi,eta)+c(inds(2))*N2(xi,eta)+c(inds(3))*N3(xi,eta));
            res = CONSTS.V_mu*det_J*quadrature_standard_triangle(integrand);
        end

        function res = compute_epsilon(m,tr,triangle_ind,c,CONSTS)
            [~, inds, r, det_J] = get_values_integrals_greek(tr,triangle_ind);
            Nm = str2func(strcat('N',int2str(m)));
            integrand = @(xi,eta) ...
                (r(1)*N1(xi,eta)+r(2)*N2(xi,eta)+r(3)*N3(xi,eta))*Nm(xi,eta) ...
                / (1+1/CONSTS.K_mfu*(c(inds(1))*N1(xi,eta)+c(inds(2))*N2(xi,eta)+c(inds(3))*N3(xi,eta)));
            res = CONSTS.V_mfv*det_J*quadrature_standard_triangle(integrand);
        end
    end
end


function res = partial_r_phi(i, rz)
det_J = det_jacobian(rz);
z = rz(:,2);
switch i
    case 1
        res = (z(2)-z(3))/det_J;
    case 2
        res = (z(3)-z(1))/det_J;
    case 3
        res = -(z(2)-z(1))/det_J;
end
end

function res = partial_z_phi(i, rz)
det_J = det_jacobian(rz);
r = rz(:,1);
switch i
    case 1
        res = (r(3)-r(2))/det_J;
    case 2
        res = -(r(3)-r(1))/det_J;
    case 3
        res = (r(2)-r(1))/det_J;
end
end

function J_det = det_jacobian(rz)
r = rz(:,1); z = rz(:,2);
J = [r(2)-r(1) r(3)-r(1); z(2)-z(1) z(3)-z(1)];
J_det = abs(det(J));
end

function res = N1(xi,eta)
    res = 1-xi-eta;
end
function res = N2(xi,~)
    res = xi;
end
function res = N3(~,eta)
    res = eta;
end

function [M, inds, r, det_J] = get_values_integrals_greek(tr,triangle_ind)
    M = size(tr.Points,1);
    inds = tr.ConnectivityList(triangle_ind,:); % Indices of triangle points
    rz = tr.Points(inds,:);                     % Triangle points
    r = rz(:,1);
    det_J = det_jacobian(rz);
end

function res = quadrature_standard_triangle(g)
    % Approximate the integral of the function handle g(xi,eta) on the
    % standard triangle in the (xi,eta)-plane using Gaussian quadrature.
    % See http://people.ucalgary.ca/~aknigh/fea/fea/triangles/num_ex.html
    res = -27/96*g(1/3,1/3) + 25/96*(g(1/5,1/5) + g(1/5,3/5) + g(3/5,1/5));
end
