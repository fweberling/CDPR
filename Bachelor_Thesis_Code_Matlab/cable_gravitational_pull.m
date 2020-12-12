%maximales Durchhängen des Seils
format long;
% cable length SKM
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
robot.phi_z=-pi/9;
robot.position=[-0.4;-0.4;0];

q=inverse_kinematics_SKM(robot);

L=q(4);  %pose[-0.4;-0.4;-pi/9];
% pretension
S=5;

mu=0.0013;
g=9.81;

% homogenous load due to cable mass
q_0=mu*g;

% maximum slack at L/2
w_max=q_0*L^2/(8*S);

u = -0.5*10^(-3):0.00001:2*10^(-3); 

%linear functon

lin=(2*w_max)/L*u+1;

syms u_krit;

coshypersolver=cosh(u_krit);
linsolver=(2*w_max)/L*u_krit+1;

figure


plot(u,cosh(u),'r'); 
hold on;
plot(u,lin,'b');
grid on;
x_axis_label=xlabel('u [-]');
y_axis_label=ylabel('f(u) [-]');
h=findobj(gcf,'type','line');
set(h,'linewidth',2);
legend({'f_{lin}(u)','cosh(u)'}, 'Location','northwest')
heading=title('Determination of u_{krit}');

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

solve(linsolver==coshypersolver,u_krit)
%%
%determine graphically

u_krit=0.00173;

% Horizontal pull

H=q_0*L/(2*u_krit);


% real cable length determined through catenary line

L_real=2*H/q_0*sinh((q_0/H)*L/2);

% approximate real cable length

L_approx=L+8/3*w_max^2/L;

%difference in cable length

L_Delta=L_real-L;
%%
%converting figure to latex
matlab2tikz('Determination_of_u_krit.tex');
