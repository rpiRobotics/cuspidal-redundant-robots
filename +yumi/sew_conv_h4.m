classdef sew_conv_h4 < sew_conv
    methods
        function psi = fwd_kin_q(obj, q, kin)
            [~, ~, P_SEW] = fwdkin_inter(kin, q, [1 4 7]);
            
            h_4_0 = rot(kin.H(:,1), q(1)) * rot(kin.H(:,2), q(2)) * rot(kin.H(:,3), q(3)) * kin.H(:,4);
            
            psi = obj.fwd_kin(P_SEW(:,1), P_SEW(:,1) + h_4_0, P_SEW(:,3));
        end
    end
end