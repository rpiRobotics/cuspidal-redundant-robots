kin = define_yumi;
kin.P(:,end) = [36;0;0];



for i = 1:1e3
SEW = yumi.sew_abb(rand_normal_vec);
q = rand_angle([7 1]);

psi_A = SEW.fwd_kin(kin, q);
psi_B = SEW.fwd_kin_sp2(kin, q);
psi_C = SEW.fwd_kin_sp2(kin, q);

assert( (norm(wrapToPi(psi_A - psi_B))) < 1e-12 )
assert( (norm(wrapToPi(psi_A - psi_C))) < 1e-12 )
end

%%
SEW = yumi.sew_abb([0;0;1]);
q = zeros([7 1]);

psi_A = SEW.fwd_kin(kin, q)
psi_B = SEW.fwd_kin_sp2(kin, q)
psi_C = SEW.fwd_kin_sp4(kin, q)