% Test cases from RobotStudio when e_r = e_y
% Controller > Configuration > Motion > Robot > Arm-Angle Reference Direction

q_1 = zeros([7 1]);
p_0T_1 = [341.5; 0; 598];
RPY_1 = deg2rad([0; 90; 0]);
psi_1_z = deg2rad(0);
psi_1_y = deg2rad(90);

q_2 = deg2rad(20)*ones([7 1]);
p_0T_2 = [292.46; 234.73; 348.75];
RPY_2 = deg2rad([148.67; 16.13; 175.45]);
psi_2_z = deg2rad(6.86);
psi_2_y = deg2rad(105.54);

q_3 = deg2rad(30)*ones([7 1]);
p_0T_3 = [165.49; 290.86; 255.37];
RPY_3 = deg2rad([157.34; -19.80; -175.99]);
psi_3_z = deg2rad(14.50);
psi_3_y = deg2rad(99.26);

q_4 = deg2rad(-40)*ones([7 1]);
p_0T_4 = [-219.17; 12.62; 796.12];
RPY_4 = deg2rad([13.67; -3.39; -152.29]);
psi_4_z = deg2rad(106.67);
psi_4_y = deg2rad(15.68);

q_5 = deg2rad(-80)*ones([7 1]);
p_0T_5 = [-126.60; 499.23; 441.88]; 
RPY_5 = deg2rad([5.30; 17.08; 34.27]);
psi_5_z = deg2rad(99.01);
psi_5_y = deg2rad(-26);


% New test cases

q_6 = deg2rad([104   -31     1   -18   -55   117    69])';
p_0T_6 = [32.02; 91.43; 763.94]; 
RPY_6 = deg2rad([157.59; 43.43; 170.70]);
psi_6_z = deg2rad(19.58);
psi_6_y = deg2rad(-165.56);

q_7 = deg2rad([162    1  -23  -88  241   98  -44])';
p_0T_7 = [45.84; 32.22; 846.74]; 
RPY_7 = deg2rad([-98.58; 43.30; -76.21]);
psi_7_z = deg2rad(80.07);
psi_7_y = deg2rad(-131.03);


q_8 = deg2rad([-20  -139  -144    77  -263   128  -104])';
p_0T_8 = [-47.12; -136.89; 270.77]; 
RPY_8 = deg2rad([-106.88; -60.36; 122.56]);
psi_8_z = deg2rad(139.68);
psi_8_y = deg2rad(-15.92);


q_9 = deg2rad([67  -142   121    78   -38   113  -194])';
p_0T_9 = [-100.23; 3.25; 338.01]; 
RPY_9 = deg2rad([43.05; -59.81; 44.29]);
psi_9_z = deg2rad(-146.69);
psi_9_y = deg2rad(119.74);


q_10 = deg2rad([-58 -138  -35   70 -283  133  207])';
p_0T_10 = [-123.66; 141.01; 425.04]; 
RPY_10 = deg2rad([-118.25; -37.89; -10.09]);
psi_10_z = deg2rad(165.86);
psi_10_y = deg2rad(50.38);


psi_vec_robotstudio_z = [psi_1_z
psi_2_z
psi_3_z
psi_4_z
psi_5_z
psi_6_z
psi_7_z
psi_8_z
psi_9_z
psi_10_z];


psi_vec_robotstudio_y = [psi_1_y
psi_2_y
psi_3_y
psi_4_y
psi_5_y
psi_6_y
psi_7_y
psi_8_y
psi_9_y
psi_10_y];

Q = [q_1 q_2 q_3 q_4 q_5 q_6 q_7 q_8 q_9 q_10];

%% Compare p_0T and RPY
kin = define_yumi;
% kin.P(:,end) = [36;0;0];

[R_t_1, p_t_1] = fwdkin(kin, q_1);
[R_t_2, p_t_2] = fwdkin(kin, q_2);
[R_t_3, p_t_3] = fwdkin(kin, q_3);
[R_t_4, p_t_4] = fwdkin(kin, q_4);
[R_t_5, p_t_5] = fwdkin(kin, q_5);
[R_t_6, p_t_6] = fwdkin(kin, q_6);
[R_t_7, p_t_7] = fwdkin(kin, q_7);
[R_t_8, p_t_8] = fwdkin(kin, q_8);
[R_t_9, p_t_9] = fwdkin(kin, q_9);
[R_t_10, p_t_10] = fwdkin(kin, q_10);

P_t = [p_t_1 p_t_2 p_t_3 p_t_4 p_t_5 p_t_6 p_t_7 p_t_8 p_t_9 p_t_10];
P_robotstudio = [p_0T_1 p_0T_2 p_0T_3 p_0T_4 p_0T_5 p_0T_6 p_0T_7 p_0T_8 p_0T_9 p_0T_10];
round(P_t - P_robotstudio, 2)

%%
R_7T = round(rot([0;1;0], pi/2));
[RPY_1A, RPY_1B] = R_to_RPY(R_t_1*R_7T);
[RPY_2A, RPY_2B] = R_to_RPY(R_t_2*R_7T);
[RPY_3A, RPY_3B] = R_to_RPY(R_t_3*R_7T);
[RPY_4A, RPY_4B] = R_to_RPY(R_t_4*R_7T);
[RPY_5A, RPY_5B] = R_to_RPY(R_t_5*R_7T);
[RPY_6A, RPY_6B] = R_to_RPY(R_t_6*R_7T);
[RPY_7A, RPY_7B] = R_to_RPY(R_t_7*R_7T);
[RPY_8A, RPY_8B] = R_to_RPY(R_t_8*R_7T);
[RPY_9A, RPY_9B] = R_to_RPY(R_t_9*R_7T);
[RPY_10A, RPY_10B] = R_to_RPY(R_t_10*R_7T);

RPY_t = [RPY_1A RPY_2A RPY_3A RPY_4A RPY_5A RPY_6A RPY_7A RPY_8A RPY_9A RPY_10A]
RPY_robotstudio = [RPY_1 RPY_2 RPY_3 RPY_4 RPY_5 RPY_6 RPY_7 RPY_8 RPY_9 RPY_10]
round(rad2deg(RPY_t-RPY_robotstudio), 2)

%% Test RPY function
% R = rot(e_3, Y)*rot(e_2, P)*rot(e_1, Y)
e_1 = [1;0;0];
e_2 = [0;1;0];
e_3 = [0;0;1];

R = rot(rand_normal_vec, rand_angle);
[RPY_A, RPY_B] = R_to_RPY(R);
rot(e_3, RPY_A(3))*rot(e_2, RPY_A(2))*rot(e_1, RPY_A(1))-R
rot(e_3, RPY_B(3))*rot(e_2, RPY_B(2))*rot(e_1, RPY_B(1))-R
%% Compare SEW angles using e_r = e_z

kin = define_yumi;
SEW = yumi.sew_abb([0;0;1]);

psi_sign_z = [SEW.fwd_kin(kin, q_1)
SEW.fwd_kin(kin, q_2)
SEW.fwd_kin(kin, q_3)
SEW.fwd_kin(kin, q_4)
SEW.fwd_kin(kin, q_5)
SEW.fwd_kin(kin, q_6)
SEW.fwd_kin(kin, q_7)
SEW.fwd_kin(kin, q_8)
SEW.fwd_kin(kin, q_9)
SEW.fwd_kin(kin, q_10)];

round(rad2deg(psi_sign_z - psi_vec_robotstudio_z), 2)

%% Compare SEW angles using e_r = e_y

kin = define_yumi;
SEW = yumi.sew_abb([0;1;0]);

psi_sign_y = [SEW.fwd_kin(kin, q_1)
SEW.fwd_kin(kin, q_2)
SEW.fwd_kin(kin, q_3)
SEW.fwd_kin(kin, q_4)
SEW.fwd_kin(kin, q_5)
SEW.fwd_kin(kin, q_6)
SEW.fwd_kin(kin, q_7)
SEW.fwd_kin(kin, q_8)
SEW.fwd_kin(kin, q_9)
SEW.fwd_kin(kin, q_10)];

round(rad2deg(psi_sign_y - psi_vec_robotstudio_y), 2)

%% Compare CONVENTIONAL SEW angles using e_r = e_z

kin = define_yumi;
e_r = [0;0;1];

psi_vec = [SEW_conv_h4_fwd_kin(q_1, kin, e_r)
SEW_conv_h4_fwd_kin(q_2, kin, e_r)
SEW_conv_h4_fwd_kin(q_3, kin, e_r)
SEW_conv_h4_fwd_kin(q_4, kin, e_r)
SEW_conv_h4_fwd_kin(q_5, kin, e_r)
SEW_conv_h4_fwd_kin(q_6, kin, e_r)
SEW_conv_h4_fwd_kin(q_7, kin, e_r)
SEW_conv_h4_fwd_kin(q_8, kin, e_r)
SEW_conv_h4_fwd_kin(q_9, kin, e_r)
SEW_conv_h4_fwd_kin(q_10, kin, e_r)];
psi_conv_z = wrapToPi(psi_vec+pi/2);

rad2deg([psi_conv_z psi_vec_robotstudio_z])
round(rad2deg(psi_conv_z - psi_vec_robotstudio_z), 2)

%% Compare CONVENTIONAL SEW angles using e_r = e_y

kin = define_yumi;
e_r = [0;1;0];

psi_vec = [SEW_conv_h4_fwd_kin(q_1, kin, e_r)
SEW_conv_h4_fwd_kin(q_2, kin, e_r)
SEW_conv_h4_fwd_kin(q_3, kin, e_r)
SEW_conv_h4_fwd_kin(q_4, kin, e_r)
SEW_conv_h4_fwd_kin(q_5, kin, e_r)
SEW_conv_h4_fwd_kin(q_6, kin, e_r)
SEW_conv_h4_fwd_kin(q_7, kin, e_r)
SEW_conv_h4_fwd_kin(q_8, kin, e_r)
SEW_conv_h4_fwd_kin(q_9, kin, e_r)
SEW_conv_h4_fwd_kin(q_10, kin, e_r)];
psi_conv_y = wrapToPi(psi_vec+pi/2);

round(rad2deg(psi_conv_y - psi_vec_robotstudio_y), 2)

%% temp
SEW = yumi.sew_abb( [0;1;0]);

SEW_conv = rad2deg(SEW_conv_h4_fwd_kin(q_9, kin, [0;1;0]))+90
SEW = rad2deg(SEW.fwd_kin(kin, q_9))

%%
q_disp = q_10;

q_deg = rad2deg(yumi.q2ABB(q_disp));
q_deg(1:6)
q_deg(7)

[~, p ] = fwdkin(kin, q_disp)


%%
% find random angles in joint limits with different SEW angles
N = 1e5;
psi_A_list = NaN([N 1]);
psi_B_list = NaN([N 1]); 
Q = NaN([7 N]);

[kin, q_min, q_max] = define_yumi;

e_r = [0;0;1];
% e_r = [0;1;0];
SEW = yumi.sew_abb(e_r);

for i = 1:N
    % q = rand_angle([7 1]);
    q = q_min + (q_max-q_min) .* rand([7 1]);
    q_deg = rad2deg(q);
    q_deg = round(q_deg);
    q = deg2rad(q_deg);
    
    psi_A_list(i) = SEW_conv_h4_fwd_kin(q, kin, e_r)+pi/2;
    psi_B_list(i) = SEW.fwd_kin(kin, q);
    Q(:,i) = q;
end

p = plot(wrapToPi(psi_A_list), wrapToPi(psi_B_list), '.');
xlim([-pi, pi])
ylim([-pi, pi])
xlabel("Conventional SEW angle")
ylabel("sign() SEW angle")
title("N = "+N+", with joint limits")
p.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow("i",1:N);

%%

errs_deg = rad2deg(sort(abs(wrapToPi(psi_A_list - psi_B_list))));

semilogy(round(errs_deg, 2));
ylabel("Error (DEG)")
xlim([0 N])

% sum (round(errs_deg, 2) <= 0.1) / N

%% Make a table with results

T = table();

T.ID = (1:10)';

T.q1_deg = rad2deg(Q(1,:)');
T.q2_deg = rad2deg(Q(2,:)');
T.q3_deg = rad2deg(Q(3,:)');
T.q4_deg = rad2deg(Q(4,:)');
T.q5_deg = rad2deg(Q(5,:)');
T.q6_deg = rad2deg(Q(6,:)');
T.q7_deg = rad2deg(Q(7,:)');

T.psi_z_RS_deg = round(rad2deg(psi_vec_robotstudio_z), 2);
T.psi_z_conv_deg = round(rad2deg(psi_conv_z), 2);
T.psi_z_sign_deg = round(rad2deg(psi_sign_z), 2);



T.psi_y_RS_deg = round(rad2deg(psi_vec_robotstudio_y), 2);
T.psi_y_conv_deg = round(rad2deg(psi_conv_y), 2);
T.psi_y_sign_deg = round(rad2deg(psi_sign_y), 2);
T


%% Format the table for LaTeX

% Open a file for writing (change the filename if desired)
fid = fopen('latex_table.txt','w');

% Write the LaTeX table header
fprintf(fid, '\\begin{tabular}{@{}rrrrrrrrrrrrrr@{}}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, '   &    &    &    &    &    &    &    & \\multicolumn{3}{c}{ez} & \\multicolumn{3}{c}{ey} \\\\ \\cmidrule(r){9-11} \\cmidrule(l){12-14}\n');
fprintf(fid, '\\# & \\(q_1\\) & \\(q_2\\) & \\(q_3\\) & \\(q_4\\) & \\(q_5\\) & \\(q_6\\) & \\(q_7\\) & RS   & conv   & sign   & RS   & conv   & sign   \\\\ \\midrule\n');

% Loop over each row of T and write a formatted row.
for i = 1:height(T)
    % Format the ID column as integer (no decimals)
    id_str = sprintf('%d', T.ID(i));

    % Format numeric values using %8.2f (right-aligned, 2 decimal places)
    q1_str         = sprintf('%8.2f', T.q1_deg(i));
    q2_str         = sprintf('%8.2f', T.q2_deg(i));
    q3_str         = sprintf('%8.2f', T.q3_deg(i));
    q4_str         = sprintf('%8.2f', T.q4_deg(i));
    q5_str         = sprintf('%8.2f', T.q5_deg(i));
    q6_str         = sprintf('%8.2f', T.q6_deg(i));
    q7_str         = sprintf('%8.2f', T.q7_deg(i));
    psi_z_RS_str   = sprintf('%8.2f', T.psi_z_RS_deg(i));
    psi_z_conv_str = sprintf('%8.2f', T.psi_z_conv_deg(i));
    psi_y_RS_str   = sprintf('%8.2f', T.psi_y_RS_deg(i));
    psi_y_conv_str = sprintf('%8.2f', T.psi_y_conv_deg(i));

    % Bold formatting for psi_z_sign_deg and psi_y_sign_deg in rows 6–10
    if i > 5
        psi_z_sign_str = sprintf('\\textbf{%8.2f}', T.psi_z_sign_deg(i));
        psi_y_sign_str = sprintf('\\textbf{%8.2f}', T.psi_y_sign_deg(i));
    else
        psi_z_sign_str = sprintf('%8.2f', T.psi_z_sign_deg(i));
        psi_y_sign_str = sprintf('%8.2f', T.psi_y_sign_deg(i));
    end

    % Write the row to file (separating columns with ' & ' and ending with '\\')
    fprintf(fid, '%s  & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\\n', ...
        id_str, q1_str, q2_str, q3_str, q4_str, q5_str, q6_str, q7_str, ...
        psi_z_RS_str, psi_z_conv_str, psi_z_sign_str, psi_y_RS_str, psi_y_conv_str, psi_y_sign_str);
end

% Write the table footer
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');

% Close the file
fclose(fid);





%%

function [RPY_1, RPY_2] = R_to_RPY(R)
% R = rot(e_3, Y)*rot(e_2, P)*rot(e_1, R)
e_1 = [1;0;0];
e_2 = [0;1;0];
e_3 = [0;0;1];

% R*e_1 = rot(e_3, Y)*rot(e_2, P)*e_1

[Y_t, P_t] = subproblem.sp_2(R*e_1, e_1, -e_3, e_2);

Y_1 = Y_t(1); P_1 = P_t(1);
Y_2 = Y_t(end); P_2 = P_t(end);

% R*e_2 = rot(Y, e_3)*rot(e_2, P)*rot(e_1, R)*e_2
%  rot(e_1, R)*e_2 = (rot(e_3, Y)*rot(e_2, P))' * R*e_2
R_1 = subproblem.sp_1(e_2, (rot(e_3, Y_1)*rot(e_2, P_1))' * R*e_2, e_1);
R_2 = subproblem.sp_1(e_2, (rot(e_3, Y_2)*rot(e_2, P_2))' * R*e_2, e_1);

RPY_1 = [R_1; P_1; Y_1];
RPY_2 = [R_2; P_2; Y_2];
end

function psi = SEW_conv_h4_fwd_kin(q, kin, e_r)
    SEW_conv = sew_conv(e_r);
    [~, ~, P_SEW] = fwdkin_inter(kin, q, [1 4 7]);

    % Calculate R_03 h_4
    h_4_0 = rot(kin.H(:,1), q(1)) * rot(kin.H(:,2), q(2)) * rot(kin.H(:,3), q(3)) * kin.H(:,4);

    psi = SEW_conv.fwd_kin(P_SEW(:,1), P_SEW(:,1)+h_4_0, P_SEW(:,3));
end