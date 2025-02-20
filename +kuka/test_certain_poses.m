kin = robot_kin.kuka;

q = rand_angle([7 1]);
q(3) = 0;
% [R, p] = fwdkin(kin, q);


J = robotjacobian(kin, q);
J_psi = [1 0 0 0 0 0 0];
J_aug = [J; J_psi]; 
det(J_aug)