kin_7 = define_yumi;

kin = fwdkin_partial(kin_7, rand_angle, 1);

q = rand_angle([6 1]);
% q = zeros([6 1]); q(5) = pi/2;
[R, p] = fwdkin(kin, q)

Q = IK.IK_4_6_intersecting(R, p, kin)

[R_t, p_t] = fwdkin(kin, Q(:,1));
[R p] - [R_t p_t]
