classdef sew_sign
% SEW angle for ABB YuMi using definition from
% "Singularities of ABB's YuMi 7-DOF robot arm" by
% M. Asgari, I.A. Bonev, and C. Gosselin.
%
% This definition does NOT match the definition used by ABB in RobotStudio

    properties
        e_r
    end

    methods
        function obj = sew_abb(e_r)
            obj.e_r = e_r;
        end

        function [psi, p_SW, e_x, h_4_0] = fwd_kin(obj, kin, q)
            h_4_0 = rot(kin.H(:,1), q(1)) * rot(kin.H(:,2), q(2)) * rot(kin.H(:,3), q(3)) * kin.H(:,4);
            [~, ~, p_0W] = fwdkin_inter(kin, q, 7);

            p_0S = kin.P(:,1);
            p_SW = p_0W - p_0S;
            
            e_x = cross(obj.e_r, p_SW) / norm(cross(obj.e_r, p_SW));
            
            psi = atan2( ...
                sign(dot(h_4_0, obj.e_r))*norm(cross(e_x, h_4_0)),...
                dot(e_x, h_4_0));
        end

        function psi = fwd_kin_sp2(obj, kin, q)
            % R_03 * h_4 = R(e_x, alpha)*R(e_y, psi)*e_x
            % Where -pi/2 < alpha < pi/2

            h_4_0 = rot(kin.H(:,1), q(1)) * rot(kin.H(:,2), q(2)) * rot(kin.H(:,3), q(3)) * kin.H(:,4);
            [~, ~, p_0W] = fwdkin_inter(kin, q, 7);

            p_0S = kin.P(:,1);
            p_SW = p_0W - p_0S;
            
            e_x = cross(obj.e_r, p_SW) / norm(cross(obj.e_r, p_SW));
            e_y = cross(e_x, obj.e_r);

            [t_alpha, t_psi] = subproblem.sp_2(h_4_0, e_x, -e_x, e_y);
            if abs(t_alpha(1))<pi/2
                psi = t_psi(1);
            else
                psi = t_psi(end);
            end

        end

        function psi = fwd_kin_sp4(obj, kin, q)
            % e_x'* R_03 * h_4 = e_x'*R(e_y, psi)*e_x
            % where sign(psi) = sign(e_r' * R_03 * h_4)
            
            h_4_0 = rot(kin.H(:,1), q(1)) * rot(kin.H(:,2), q(2)) * rot(kin.H(:,3), q(3)) * kin.H(:,4);
            [~, ~, p_0W] = fwdkin_inter(kin, q, 7);

            p_0S = kin.P(:,1);
            p_SW = p_0W - p_0S;
            
            e_x = cross(obj.e_r, p_SW) / norm(cross(obj.e_r, p_SW));
            e_y = cross(e_x, obj.e_r);

            t_psi = subproblem.sp_4(e_x, e_x, e_y, e_x'* h_4_0);
            if sign(dot(h_4_0, obj.e_r)) * t_psi(1) > 0
                psi = t_psi(1);
            else
                psi = t_psi(end);
            end


        end

        function [J_psi, sign_term] = jacobian(obj, kin, q)
            [~, p_SW, e_x, h_4_0] = fwd_kin(obj, kin, q);

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
            J_x_W = obj.jacobian_x_W(obj.e_r, p_SW, e_x);


            alpha = -sign(dot(h_4_0, obj.e_r)) / norm(cross(e_x, h_4_0));

            % Jacobian of arm angle wrt R_03*h_4
            J_psi_h4 = alpha * e_x';

            % Jacobian of arm angle wrt p_0W
            J_psi_W = alpha * h_4_0' * J_x_W;

            % Jacobian of arm angle wrt q
            J_psi = J_psi_W * J_W + J_psi_h4 * J_h4;
            
            % d/dt sign(t) is undefined at t=0
            sign_term = dot(h_4_0, obj.e_r);
        end

        % Jacobian of e_x wrt p_0W
        function J_x_W = jacobian_x_W(obj, e_r, p_SW, e_x)
            num = cross(e_r, e_x)* e_x';
            den = norm(cross(e_r, e_x)'* p_SW);
            J_x_W = num/den;
        end

        function [e_r, e_x] = inv_kin(obj, p_0W)
            e_r = obj.e_r;
            e_x = cross(e_r, p_0W) / norm(cross(e_r, p_0W));
        end
    end
end