% SKM Levenberg-Marquardt with circle method (SKM LM circle method)

%Definition of Input Parameters for Model
%format long;

% generating object of class Robot
robot=Robot;

% actuators' position K_0
robot.a_1=[-0.6;-0.4;0];
robot.a_2=[0.6;-0.4;0];
robot.a_3=[-0.6;0.4;0];
robot.a_4=[0.6;0.4;0];

% cable attachment position in K_P
robot.b_1=[-sqrt(2)/10;sqrt(2)/10;0];
robot.b_2=[sqrt(2)/10;sqrt(2)/10;0];
robot.b_3=[-sqrt(2)/10;-sqrt(2)/10;0];
robot.b_4=[sqrt(2)/10;-sqrt(2)/10;0];

% drums' radius
robot.radius=0.02;

% rotation direction of actuators
robot.rot_direc_1=-1;
robot.rot_direc_2=1;
robot.rot_direc_3=1;
robot.rot_direc_4=-1;

% pose of CDPR
robot.phi_x=0;%pelvic drop
robot.phi_y=0;%anterior tilt
robot.phi_z=pi/9;%rotation
robot.position=[0.1;-0.15;0];%position

% number of increment steps n-1  (20-->19)
robot.steps=15;

n=robot.steps-1;

%definition of workspace's size
matrix_x_components=0.6;
matrix_y_components=0.6;

increment_step_x=matrix_x_components/n;
increment_step_y=matrix_y_components/n;

% array of x- and y-components
x=-0.3:increment_step_x:0.3;
y=-0.3:increment_step_y:0.3;

% number of components in y
m=numel(y);

%% Kalibrierung
% Levenberg-Marquardt (SKM) with x_0=[0;0;0] method
robot.L_SKM=inverse_kinematics_SKM(robot);
robot.x_0=0;
robot.y_0=0;
robot.phi_0=0;
robot.k_max=100; %maximum number of iterations
robot.tau=10^-3; % proposing that x_0 is a good guess (10^-3)
robot.epsilon_1=10^-6; %proposed 10^-6
robot.epsilon_2=10^-6; %proposed 10^-6
[X_LM,~]=Levenberg_Marquardt_SKM(robot);

% Levenberg-Marquardt (EKM) with x_0=[0;0;0] method
robot.L_EKM=inverse_kinematics_EKM(robot);
robot.x_0=0;
robot.y_0=0;
robot.phi_0=0;
robot.k_max=100; %maximum number of iterations
robot.tau=10^-3; % proposing that x_0 is a good guess (10^-3)
robot.epsilon_1=10^-6; %proposed 10^-6
robot.epsilon_2=10^-6; %proposed 10^-6
[X_LM,~]=Levenberg_Marquardt_EKM(robot);

% Dog Leg (SKM) with x_0=[0;0;0] method
robot.L_SKM=inverse_kinematics_SKM(robot);
robot.x_0=0; %starting point x-coordinate
robot.y_0=0; %starting point y-coordinate
robot.phi_0=0; %starting pint phi-coordinate
robot.delta=0.7; %radius of trust region proposed 0.7
robot.k_max=100; %maximum number of iterations
robot.epsilon_1=10^-5; %proposed 10^-5
robot.epsilon_2=10^-5; %proposed 10^-5
robot.epsilon_3=10^-5; %proposed 10^-5
[X_DL,~]=Dog_Leg_SKM(robot);

% Dog Leg (SKM) with x_0=[0;0;0] method
robot.L_EKM=inverse_kinematics_EKM(robot);
robot.x_0=0; %starting point x-coordinate
robot.y_0=0; %starting point y-coordinate
robot.phi_0=0; %starting pint phi-coordinate
robot.delta=0.7; %radius of trust region proposed 0.7
robot.k_max=100; %maximum number of iterations
robot.epsilon_1=10^-5; %proposed 10^-5
robot.epsilon_2=10^-5; %proposed 10^-5
robot.epsilon_3=10^-5; %proposed 10^-5
[X_DL,~]=Dog_Leg_EKM(robot);

box_method(robot);
circle_method(robot);


%% 

% Ideal Position error matrix
position_accuracy_matrix_ideal=zeros(robot.steps);


% Position deviation Matrix
for iterations_1=1:m
    
    
    y_vector=y(1,iterations_1);
    for iterations_2=1:robot.steps
        robot.position=[x(1,iterations_2);y_vector;0];
        
        % Levenberg-Marquardt (SKM) with circle method for initial estimate
        robot.L_SKM=inverse_kinematics_EKM(robot);
        estimated_position=circle_method(robot);
        robot.x_0=estimated_position(1);
        robot.y_0=estimated_position(2);
        robot.phi_0=0;
        robot.k_max=100; %maximum number of iterations
        robot.tau=10^-3; % proposing that x_0 is a good guess (10^-3)
        robot.epsilon_1=10^-5; %proposed 10^-5
        robot.epsilon_2=10^-5; %proposed 10^-5
        [X_LM,~]=Levenberg_Marquardt_SKM(robot);
        
        position_accuracy_array(iterations_2)=...
            sqrt((X_LM(1)-x(1,iterations_2))^2 ...
            +(X_LM(2)-y(1,iterations_1))^2);
        
    end
    position_accuracy_matrix_real(iterations_1,:)=position_accuracy_array;
end

% position deviation matrix in millimeter 
position_deviation_matrix=10^3*abs(position_accuracy_matrix_ideal...
    -position_accuracy_matrix_real);

% average error in position
number_of_components=numel(position_deviation_matrix);
sum=0;
[row, column]=size(position_deviation_matrix);
for itera_1=1:column
    for itera_2=1:row
        
        sum=sum+abs(position_deviation_matrix(itera_2,itera_1));
    end
end

position_average=sum/number_of_components;

% Ideal orientation matrix
phi_matrix_ideal=zeros(robot.steps);


for iter_1=1:m
    
    for iter_2=1:m
        phi_matrix_ideal(iter_2,iter_1)=phi_matrix_ideal(iter_2,iter_1)+robot.phi_z;
    end
end

% orientation deviation matrix
for iterations_1=1:m
    
    
    y_vector=y(1,iterations_1);
    for iterations_2=1:robot.steps
        robot.position=[x(1,iterations_2);y_vector;0];
        
        % Levenberg-Marquardt (SKM) with circle method for initial estimate
        robot.L_SKM=inverse_kinematics_EKM(robot);
        estimated_position=circle_method(robot);
        robot.x_0=estimated_position(1);
        robot.y_0=estimated_position(2);
        robot.phi_0=0;
        robot.k_max=100; %maximum number of iterations
        robot.tau=10^-3; % proposing that x_0 is a good guess (10^-3)
        robot.epsilon_1=10^-5; %proposed 10^-5
        robot.epsilon_2=10^-5; %proposed 10^-5
        [X_LM,~]=Levenberg_Marquardt_SKM(robot);
        
        phi_accuracy_array(iterations_2)=X_LM(3);
        
    end
    phi_accuracy_matrix_real(iterations_1,:)=phi_accuracy_array;
end

%orientation deviation matrix in degree
orientation_deviation_matrix=180/pi*abs(phi_matrix_ideal-phi_accuracy_matrix_real);

% average of error in phi_z
number_of_components=numel(orientation_deviation_matrix);
sum=0;
[row, column]=size(orientation_deviation_matrix);
for itera_1=1:column
    for itera_2=1:row
        
        sum=sum+abs(orientation_deviation_matrix(itera_2,itera_1));
    end
end

phi_average=sum/number_of_components;

% Iterations efficiency matrix
for iterations_1=1:m
    
    
    y_vector=y(1,iterations_1);
    for iterations_2=1:robot.steps
        robot.position=[x(1,iterations_2);y_vector;0];
        
        % Levenberg-Marquardt (SKM) with circle method for initial estimate
        robot.L_SKM=inverse_kinematics_EKM(robot);
        estimated_position=circle_method(robot);
        robot.x_0=estimated_position(1);
        robot.y_0=estimated_position(2);
        robot.phi_0=0;
        robot.k_max=100; %maximum number of iterations
        robot.tau=10^-3; % proposing that x_0 is a good guess (10^-3)
        robot.epsilon_1=10^-5; %proposed 10^-5
        robot.epsilon_2=10^-5; %proposed 10^-5
        [~,iterations_LM]=Levenberg_Marquardt_SKM(robot);
        
        iterations_array(iterations_2)=iterations_LM;
        
    end
    iterations_matrix(iterations_1,:)=iterations_array;
end

% average of executed iterations
number_of_components=numel(iterations_matrix);
sum=0;
[row, column]=size(iterations_matrix);
for itera_1=1:column
    for itera_2=1:row
        
        sum=sum+abs(iterations_matrix(itera_2,itera_1));
    end
end

iterations_average=sum/number_of_components;

%duration matrix
for iterations_1=1:m
    
    
    y_vector=y(1,iterations_1);
    for iterations_2=1:robot.steps
        robot.position=[x(1,iterations_2);y_vector;0];
        
        % Levenberg-Marquardt (SKM) with circle method for initial estimate
        robot.L_SKM=inverse_kinematics_EKM(robot);
        
        robot.phi_0=0;
        robot.k_max=100; %maximum number of iterations
        robot.tau=10^-3; % proposing that x_0 is a good guess (10^-3)
        robot.epsilon_1=10^-5; %proposed 10^-5
        robot.epsilon_2=10^-5; %proposed 10^-5
        
        tic
        estimated_position=circle_method(robot);
        robot.x_0=estimated_position(1);
        robot.y_0=estimated_position(2);
        [X_LM,iterations_LM]=Levenberg_Marquardt_SKM(robot);
        
        duration_array(iterations_2)=toc;
        
    end
    duration_matrix(iterations_1,:)=duration_array;
end

%duration matrix in milliseconds
duration_matrix=10^3*duration_matrix;

% duration average
number_of_components=numel(duration_matrix);
sum=0;
[row, column]=size(duration_matrix);
for itera_1=1:column
    for itera_2=1:row
        
        sum=sum+abs(duration_matrix(itera_2,itera_1));
    end
end

duration_average=sum/number_of_components;



% contourplots of position, orientation, iteration and duration devian
% matrix
[X,Y]=meshgrid(x,y);

figure(2)

[ha, pos]=tight_subplot(2,2,[.02 .02],[.01 .01],[0.2 0.2]);
for ii = 1:4;
    axes(ha(ii));
    if ii==1
        contourf(X,Y,position_deviation_matrix,10);
        c=colorbar;
        
        hold on;
        
        caxis([0 25]);
        
        % area where human moves
        PolyX = 0.1*[cosd(0:1:360) 1];
        PolyY = 0.1*[sind(0:1:360) 0];
        
        fill(PolyX,PolyY,'w','Edgecolor', 'w','FaceAlpha',0.1);
        hold on;
        
        % title and label
        heading=title('Position Accuracy');
        x_axis_label=xlabel('x [m]');
        y_axis_label=ylabel('y [m]');
        colorbar_label=ylabel(c, 'x');
        
        colorbar_label=ylabel(c, '\Delta r [mm]');
        
        % font weight
        set(heading,'FontWeight','bold');
        set(x_axis_label,'FontWeight','bold');
        set(y_axis_label,'FontWeight','bold');
        set(colorbar_label,'FontWeight','bold');
        
        % font size
        set(heading,'FontSize',12);
        set(x_axis_label,'FontSize',12);
        set(y_axis_label,'FontSize',12);
        set(colorbar_label,'FontSize',12);
        
        set(gca,'FontWeight','bold');
        set(gca,'FontSize',10);
        axis equal;
        axis([-0.30 0.3 -0.3 0.3]);
        
    end
    if ii==2
        % orientation accuracy
        %subplot(2,2,2);
        contourf(X,Y,orientation_deviation_matrix,10);
        c=colorbar;
        
        hold on;
        
        caxis([0 2.5]);
        
        % area where human moves
        PolyX = 0.1*[cosd(0:1:360) 1];
        PolyY = 0.1*[sind(0:1:360) 0];
        
        fill(PolyX,PolyY,'w','Edgecolor', 'w','FaceAlpha',0.1);
        hold on;
        
        % title and label
        heading=title('Orientation Accuracy');
        x_axis_label=xlabel('x [m]');
        y_axis_label=ylabel('y [m]');
        
        
        colorbar_label=ylabel(c, '\Delta \phi_z [^\circ]');
        
        % font weight
        set(heading,'FontWeight','bold');
        set(x_axis_label,'FontWeight','bold');
        set(y_axis_label,'FontWeight','bold');
        set(colorbar_label,'FontWeight','bold');
        
        % font size
        set(heading,'FontSize',12);
        set(x_axis_label,'FontSize',12);
        set(y_axis_label,'FontSize',12);
        set(colorbar_label,'FontSize',12);
        
        set(gca,'FontWeight','bold');
        set(gca,'FontSize',10);
        axis equal;
        axis([-0.30 0.3 -0.3 0.3]);
        
    end
    
    if ii==3
        %iteration efficiency
        %subplot(2,2,3);
        contourf(X,Y,iterations_matrix,15);
        
        
        c=colorbar;%Create Colorbar
        caxis([2 7]);
        set(c, 'YTick', linspace(2, 7, 6));
        
        hold on;
        
        
        % area where human moves
        PolyX = 0.1*[cosd(0:1:360) 1];
        PolyY = 0.1*[sind(0:1:360) 0];
        
        fill(PolyX,PolyY,'w','Edgecolor', 'w','FaceAlpha',0.1);
        hold on;
        
        % title and label
        heading=title('Iteration''s Efficiency');
        x_axis_label=xlabel('x [m]');
        y_axis_label=ylabel('y [m]');
        
        colorbar_label=ylabel(c, '\Delta k [-]');
        
        % font weight
        set(heading,'FontWeight','bold');
        set(x_axis_label,'FontWeight','bold');
        set(y_axis_label,'FontWeight','bold');
        set(colorbar_label,'FontWeight','bold');
        
        % font size
        set(heading,'FontSize',12);
        set(x_axis_label,'FontSize',12);
        set(y_axis_label,'FontSize',12);
        set(colorbar_label,'FontSize',12);
        
        set(gca,'FontWeight','bold');
        set(gca,'FontSize',10);
        axis equal;
        axis([-0.30 0.3 -0.3 0.3]);
        
    end
    
    if ii==4
        %duration matrix
        %subplot(2,2,4);
        contourf(X,Y,duration_matrix,10);
        c=colorbar;
        
        hold on;
        
        caxis([0 6]);
        
        % area where human moves
        PolyX = 0.1*[cosd(0:1:360) 1];
        PolyY = 0.1*[sind(0:1:360) 0];
        
        fill(PolyX,PolyY,'w','Edgecolor', 'w','FaceAlpha',0.1);
        hold on;
        
        % title and label
        heading=title('Duration');
        x_axis_label=xlabel('x [m]');
        y_axis_label=ylabel('y [m]');
        colorbar_label=ylabel(c, 'x');
        
        colorbar_label=ylabel(c, '\Delta t [ms]');
        
        % font weight
        set(heading,'FontWeight','bold');
        set(x_axis_label,'FontWeight','bold');
        set(y_axis_label,'FontWeight','bold');
        set(colorbar_label,'FontWeight','bold');
        
        % font size
        set(heading,'FontSize',12);
        set(x_axis_label,'FontSize',12);
        set(y_axis_label,'FontSize',12);
        set(colorbar_label,'FontSize',12);
        
        set(gca,'FontWeight','bold');
        set(gca,'FontSize',10);
        axis equal;
        axis([-0.30 0.3 -0.3 0.3]);
        
    end
end

disp('position average');
disp(position_average);
disp('orientation average');
disp(phi_average);
disp('iteration average');
disp(iterations_average);
disp('duration average');
disp(duration_average);

% Output Maxima in Konsole

disp('maximum position error');
disp(max(max(position_deviation_matrix)));
disp('maximum orientation error');
disp(max(max(orientation_deviation_matrix)));
disp('maximum iterations');
disp(max(max(iterations_matrix)));
disp('maximum duration');
disp(max(max(duration_matrix)));

%%
% convert figure to latex
matlab2tikz('SKM_LM_circle_method_based_EKM.tex');