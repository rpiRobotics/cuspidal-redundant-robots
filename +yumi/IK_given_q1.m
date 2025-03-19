function Q = IK_given_q1(R_0T, p_0T, kin_7, q1)

[kin, R_7T] = fwdkin_partial(kin_7, q1, 1);

R_07 = R_0T * R_7T';

Q_6 = IK.IK_4_6_intersecting_mex(R_07, p_0T, kin);

Q = [q1*ones(1, width(Q_6)); Q_6];

end