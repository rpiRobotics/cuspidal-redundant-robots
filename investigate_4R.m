zv = [0;0;0];
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];

kin.H = [ez ey ez ey];
kin.P = [zv ex 2*ex+0.2*ey 1.5*ex 0.2*ex];
kin.joint_type = [0 0 0 0];

x_list = linspace(1, 4, 100);

q4_list = linspace(-pi, pi, 25);
q1_list = [];
q2_list = [];
q3_list = [];
x_list_display = [];

for i_x = 1:numel(x_list)
    p = ex*x_list(i_x);
    for i = 1:numel(q4_list)
        kin_i = fwdkin_partial(kin, q4_list(i), 4);
        Q = cuspidal_3R.IK(p, kin_i);
        q1_list = [q1_list Q(1,:)];
        q2_list = [q2_list Q(2,:)];
        q3_list = [q3_list Q(3,:)];
        x_list_display = [x_list_display ones(size(Q(3,:)))*x_list(i_x)];
    
    end
end

%%
scatter3(x_list_display, wrapToPi(q2_list+1)-1, q3_list, [], q3_list, '.')
xlabel("x")
ylabel("q_2")
zlabel("q_3")
view(107, 13)
axis tight

%%
scatter3(x_list_display,  wrapToPi(q2_list-1)+1, wrapTo2Pi(q1_list-1)+1, [], wrapTo2Pi(q1_list-1)+1, '.')
xlabel("x")
ylabel("q_2")
zlabel("q_1")
view(107, 13)
% axis equal

%%
scatter(x_list_display, wrapTo2Pi(q1_list-1)+1)