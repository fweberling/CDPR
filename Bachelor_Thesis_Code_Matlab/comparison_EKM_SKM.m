%Definition of Input Parameters for Model
%format long;

robot=Robot;
robot.a_1=[-0.6;-0.4;0];
robot.a_2=[0.6;-0.4;0];
robot.a_3=[-0.6;0.4;0];
robot.a_4=[0.6;0.4;0];
robot.b_1=[-sqrt(2)/10;sqrt(2)/10;0];
robot.b_2=[sqrt(2)/10;sqrt(2)/10;0];
robot.b_3=[-sqrt(2)/10;-sqrt(2)/10;0];
robot.b_4=[sqrt(2)/10;-sqrt(2)/10;0];
robot.radius=0.02;
robot.phi_x=0;
robot.phi_y=0;
robot.phi_z=pi/9;
robot.position=[0.1;-0.15;0];
robot.beta_x_1=0;
robot.beta_x_2=0;
robot.beta_x_3=0;
robot.beta_x_4=0;
robot.beta_y_1=0;
robot.beta_y_2=0;
robot.beta_y_3=0;
robot.beta_y_4=0;

robot.boundary_1=[-0.6;-0.4;0];
robot.boundary_2=[0.6;-0.4;0];
robot.boundary_3=[-0.6;0.4;0];
robot.boundary_4=[0.6;0.4;0];

robot.rot_direc_1=-1;
robot.rot_direc_2=1;
robot.rot_direc_3=1;
robot.rot_direc_4=-1;

robot.radius=0.02;
q_1=inverse_kinematics_EKM(robot)-inverse_kinematics_SKM(robot);

robot.radius=0.018;
q_2=inverse_kinematics_EKM(robot)-inverse_kinematics_SKM(robot);
robot.radius=0.016;
q_3=inverse_kinematics_EKM(robot)-inverse_kinematics_SKM(robot);
robot.radius=0.014;
q_4=inverse_kinematics_EKM(robot)-inverse_kinematics_SKM(robot);
robot.radius=0.012;
q_5=inverse_kinematics_EKM(robot)-inverse_kinematics_SKM(robot);
robot.radius=0.010;
q_6=inverse_kinematics_EKM(robot)-inverse_kinematics_SKM(robot);
robot.radius=0.008;
q_7=inverse_kinematics_EKM(robot)-inverse_kinematics_SKM(robot);
robot.radius=0.006;
q_8=inverse_kinematics_EKM(robot)-inverse_kinematics_SKM(robot);
robot.radius=0.004;
q_9=inverse_kinematics_EKM(robot)-inverse_kinematics_SKM(robot);
robot.radius=0.002;
q_10=inverse_kinematics_EKM(robot)-inverse_kinematics_SKM(robot);
robot.radius=0.000;
q_11=inverse_kinematics_EKM(robot)-inverse_kinematics_SKM(robot);

%delta ||q||_2 in cm
delta_q=[norm(q_1)*10^2;norm(q_2)*10^2;norm(q_3)*10^2;norm(q_4)*10^2;norm(q_5)*10^2;norm(q_6)*10^2;...
    norm(q_7)*10^2;norm(q_8)*10^2;norm(q_9)*10^2;norm(q_10)*10^2;norm(q_11)*10^2];

%radius of drum in cm
radius=[0.02*10^2;0.018*10^2;0.016*10^2;0.014*10^2;0.012*10^2;0.010*10^2;0.008*10^2;...
    0.006*10^2;0.004*10^2;0.002*10^2;0.000*10^2];

figure

h=plot(radius,delta_q,'.');

grid on;
%h=findobj(gcf,'type','line');
%set(h,'linewidth',30);
x_axis_label=xlabel('r_d [cm]');
y_axis_label=ylabel('||\Delta q||_2 [cm]');
heading=title('Comparison SKM and EKM with Decreasing Radius ');

% font weight
set(heading,'FontWeight','bold');
set(x_axis_label,'FontWeight','bold');
set(y_axis_label,'FontWeight','bold');

% font size
set(heading,'FontSize',12);
set(x_axis_label,'FontSize',12);
set(y_axis_label,'FontSize',12);


set(gca,'FontWeight','bold');
set(gca,'FontSize',10);

%converting fig to latex
matlab2tikz('comparison_SKM_EKM_radius.tex');
