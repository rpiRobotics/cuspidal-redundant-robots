classdef sew_abb
    % SEW angle used by ABB for YuMi

    properties
        e_r
    end

    methods
        function obj = sew_abb(e_r)
            obj.e_r = e_r;
        end

        function [psi, e_r, p_0W, e_x, h_4_0] = fwd_kin(obj, kin, q)
            h_4_0 = rot(kin.H(:,1), q(1)) * rot(kin.H(:,2), q(2)) * rot(kin.H(:,3), q(3)) * kin.H(:,4);
            [~, ~, p_0W] = fwdkin_inter(kin, q, 7);

            % e_r = kin.H(:,1); % Unit vector along axis 1
            e_r = obj.e_r;
            w = p_0W; % Position vector of WCP: Intersection of axis 7 and common normal of axes 6 and 7
            % p_0S = kin.P(:,1);
            % p_0S = kin.P(:,1)+ rot(kin.H(:,1), q(1))*kin.P(:,2) + rot(kin.H(:,1), q(1)) * rot(kin.H(:,2), q(2)) * kin.P(:,3);
            p_0S = kin.P(:,1);
            e_3 = h_4_0; % Unit vector along axis 4
            
            e_x = cross(e_r, w-p_0S) / norm(cross(e_r, w-p_0S));
            
            psi = atan2( ...
                sign(dot(e_3, e_r))*norm(cross(e_x, e_3)),...
                dot(e_x, e_3));

        end

        function [J_psi, sign_term] = jacobian(obj, kin, q)
            [~, e_r, p_0W, e_x, h_4_0] = fwd_kin(obj, kin, q); % TODO move e_r

            % Jacobian of p_0W wrt q
            kin_W = kin;
            kin_W.joint_type = zeros([1 6]);
            J_W_full = robotjacobian(kin_W, q);
            J_W = [J_W_full(4:6, :) zeros(3, 1)];
            
            % Jacobian of R_03*h_4 wrt q
            kin_h4.joint_type = [0 0 0];
            kin_h4.H = kin.H(:, 1:3);
            kin_h4.P = [zeros(3, 3) kin.H(:,4)];
            J_E_full = robotjacobian(kin_h4, q);
            J_h4 = [J_E_full(4:6, :) zeros(3, 4)];

            % Jacobian of e_x wrt p_0W
            J_x_W = obj.jacobian_x_W(e_r, p_0W, e_x);


            alpha = -sign(dot(h_4_0, e_r)) / norm(cross(e_x, h_4_0));

            % Jacobian of arm angle wrt R_03*h_4
            J_psi_h4 = alpha * e_x';

            % Jacobian of arm angle wrt p_0W
            J_psi_W = alpha * h_4_0' * J_x_W;

            % Jacobian of arm angle wrt q
            J_psi = J_psi_W * J_W + J_psi_h4 * J_h4;
            
            % d/dt sign(t) is undefined at t=0
            sign_term = dot(h_4_0, e_r);
        end

        % Jacobian of e_x wrt p_0W
        function J_x_W = jacobian_x_W(obj, e_r, p_0W, e_x)
            num = cross(e_r, e_x)* e_x';
            den = norm(cross(e_r, e_x)'* p_0W);
            J_x_W = num/den;
        end

        function [e_r, e_x] = inv_kin(obj, p_0W)
            e_r = obj.e_r;
            e_x = cross(e_r, p_0W) / norm(cross(e_r, p_0W));
        end
    end
end