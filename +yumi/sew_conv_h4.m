classdef sew_conv_h4 < sew_conv
    methods
        function psi = fwd_kin_q(obj, q, kin)
            [~, ~, P_SEW] = fwdkin_inter(kin, q, [1 4 7]);
            
            h_4_0 = rot(kin.H(:,1), q(1)) * rot(kin.H(:,2), q(2)) * rot(kin.H(:,3), q(3)) * kin.H(:,4);
            
            psi = obj.fwd_kin(P_SEW(:,1), P_SEW(:,1) + h_4_0, P_SEW(:,3));
        end

        function J_psi = full_J_psi(obj, q, kin)
            % TODO This is not very efficient

            % kinematic chain for h_4
            kin_E = kin;
            kin_E.joint_type = [0 0 0];
            zv = [0;0;0];
            kin_E.P = [zv zv zv kin.H(:,4)];

            kin_W = kin;
            kin_W.joint_type = [0 0 0 0 0 0];
            
            [~, ~, P_SEW] = fwdkin_inter(kin,q, [1 4 7]);
            h_4_0 = rot(kin.H(:,1), q(1)) * rot(kin.H(:,2), q(2)) * rot(kin.H(:,3), q(3)) * kin.H(:,4);
            [J_psi_E, J_psi_W] = obj.jacobian(P_SEW(:,1), P_SEW(:,1) + h_4_0, P_SEW(:,3));
            
            J_E = robotjacobian(kin_E, q);
            J_E = [J_E(4:6, :) zeros(3, 4)];
            
            J_W = robotjacobian(kin_W, q);
            J_W = [J_W(4:6, :) zeros(3,1)];
            
            J_psi = J_psi_E*J_E + J_psi_W*J_W;
        end

        function J = J_aug(obj, q, kin)
            J_6x7 = robotjacobian(kin, q);
            J_psi = obj.full_J_psi(q, kin);
            J = [J_6x7; J_psi];
        end
    end
end