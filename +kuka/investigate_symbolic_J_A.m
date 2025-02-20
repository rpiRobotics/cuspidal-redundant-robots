kin = robot_kin.kuka;
kin.P(:, [1 end]) = 0;

syms q [7 1] real


J = simplify(cuspidal_3R.robotjacobian_sym(kin, q))
%%

% Factors: C_3, S_6, S^2_4
J_psi = [1 0 0 0 0 0 0];
J_aug = [J; J_psi];

D = simplify(det(J_aug))

%%
% Factors: S_2 S_3 S^2_4 S_6
J_psi = [0 1 0 0 0 0 0];
J_aug = [J; J_psi];

D = simplify(det(J_aug))
%%
% Complicated factors, but zero when S_6=0 or S_4=0
J_psi = [0 0 1 0 0 0 0];
J_aug = [J; J_psi];

D = simplify(det(J_aug))

%% det(J_A) = 0 as expected
J_psi = [0 0 0 1 0 0 0];
J_aug = [J; J_psi];

D = simplify(det(J_aug))