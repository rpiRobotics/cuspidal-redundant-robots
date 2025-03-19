ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];



diagrams.setup(); hold on

diagrams.dot(zv);
diagrams.sphere(zv, 1);
diagrams.arrow(-ex, ex);
diagrams.line(-ey, ey);
diagrams.arrow(-ez, ez);
diagrams.sphere_arc(zv, ez, 1, zv, 1);
diagrams.sphere_arc(zv, ez, 1, zv, 1, true, lineStyle=":");

diagrams.sphere_arc(zv, ey, 1, zv, 1);
diagrams.sphere_arc(zv, ey, 1, zv, 1, true, lineStyle=":");

% diagrams.arc(rot(ez, pi/4)*ex, ex, zv, pi, color=diagrams.colors.dark_green);


diagrams.arrow(zv, rot(ex, 0.7*pi)*rot(ez, pi/4)*ex);

diagrams.arrow(zv, rot(ex, pi/6)*ez);
diagrams.line(rot(ex, pi/6)*ez, ey*ey'*rot(ex, pi/6)*ez, LineStyle="--", LineWidth=1);

diagrams.text(rot(ex, 0.7*pi)*rot(ez, pi/4)*ex, "$R_{03}h_4$");

diagrams.text(ez, "$e_r$");
diagrams.text(ex, "$e_x$", verticalAlignment="bottom", horizontalAlignment="right");
diagrams.text(rot(ex, pi/6)*ez, "$e_{SW}$", verticalAlignment="bottom", horizontalAlignment="center", margin=0);

diagrams.angle_arc(ex, -ey, zv, pi/4);
diagrams.angle_arc(ex,  ey, zv, pi/4);
diagrams.angle_arc(rot(-ey, pi/4)*ex, ex, zv, 0.2*pi);
diagrams.angle_arc(rot(ey, pi/4)*ex, -ex, zv, 0.8*pi);
diagrams.text(rot(-ey, pi/4)*ex, "$\psi^{ABB}_+$", margin=20);
diagrams.text(rot(-ey, pi/4)*ex, "$\alpha_+$", verticalAlignment="bottom", horizontalAlignment="right", margin=0);
diagrams.text(rot(ey, pi/4)*ex, "$\psi^{ABB}_-$");
diagrams.text(rot(ey, pi/4)*ex, "$\alpha_-$", verticalAlignment="bottom", horizontalAlignment="right", margin=0);

view(30, 30)
diagrams.redraw(); hold off
%%
diagrams.save()