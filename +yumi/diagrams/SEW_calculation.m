ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];

psi = deg2rad(90+45);

er = rot(-ey, deg2rad(50))*ex;

h4 = rot(ez, psi)*rot(-ey, deg2rad(30))*ex;

diagrams.setup(); hold on


diagrams.arrow(zv, ez, color=diagrams.colors.blue);
diagrams.text(ez, "$\vec e_{SW} = \vec e_{z,C}$", align="<v");

diagrams.arrow(ez, 1.5*ez);
diagrams.text(1.5*ez, "$\vec p_{SW}$");

diagrams.arrow(zv, ex, color=diagrams.colors.red);
diagrams.text(ex, "$\vec e_{x,C}$")

diagrams.arrow(zv, ey, color=diagrams.colors.green);
diagrams.text(ey, "$\vec e_{y,C}$", align="v");

diagrams.line(zv, -ey, LineStyle=":", color=diagrams.colors.green);

diagrams.arrow(zv, er);
diagrams.text(er, "$\vec e_r$", align="<")
diagrams.line(er, ex'*er*ex, LineStyle=":");

diagrams.arrow(zv, h4);
diagrams.line(zv, h4-h4'*ez*ez, LineStyle=":")
diagrams.line(h4, h4-h4'*ez*ez, LineStyle=":")
diagrams.text(h4, "$\vec h_4$", align="v")

diagrams.angle_arc(0.5*ex, ez, zv, psi);

diagrams.angle_arc(0.75*(-ey), ez, zv, psi+pi/2);

diagrams.text(0.6*(h4-h4'*ez*ez),  "$\psi^{conv}$", align="v", margin=6)
diagrams.text(0.84*(h4-h4'*ez*ez), "$\psi^{ABB}$", align="v", margin=15)

diagrams.dot(zv)

view(-30, 30);

diagrams.redraw(); hold off

%%
diagrams.save();