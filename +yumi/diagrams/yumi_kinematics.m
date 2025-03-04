kin = define_yumi;

diagrams.setup(); hold on

diagrams.robot_plot(kin, zeros([7 1]), auto_scale=true);
diagrams.redraw(); hold off

%%
diagrams.save()