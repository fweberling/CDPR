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
robot.phi_z=0;
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



robot.phi_upperboundary=(1/9)*pi;
robot.phi_lowerboundary=(-1/9)*pi;


% number of increment steps n-1  (20-->19)
robot.steps=30;

n=robot.steps-1;

%definition of workspace's size
matrix_x_components=0.6;
matrix_y_components=0.6;

increment_step_x=matrix_x_components/n;
increment_step_y=matrix_y_components/n;

x=-0.3:increment_step_x:0.3;
y=-0.3:increment_step_y:0.3;

phi_components=robot.phi_upperboundary-robot.phi_lowerboundary;
increment_step_phi=phi_components/n;

phi_z_vector=robot.phi_lowerboundary:increment_step_phi:robot.phi_upperboundary;

m=numel(y);

%rotation
robot.phi_z=0;

position_accuracy_matrix_ideal=zeros(robot.steps);

%% circle method

for iterations_1=1:m
    
    
    y_vector=y(1,iterations_1);
    for iterations_2=1:robot.steps
        robot.position=[x(1,iterations_2);y_vector;0];
        robot.L_SKM=inverse_kinematics_SKM(robot);
        
        %estimated position circle method
        position_c=circle_method(robot);
        
        position_accuracy_array(iterations_2)=...
            sqrt((position_c(1)-x(1,iterations_2))^2 ...
            +(position_c(2)-y(1,iterations_1))^2);
        
    end
    position_accuracy_matrix_real(iterations_1,:)=position_accuracy_array;
end

deviation_matrix=abs(position_accuracy_matrix_ideal-position_accuracy_matrix_real);

[X,Y]=meshgrid(x,y);

figure

contourf(X,Y,deviation_matrix,15);
colorbar();

caxis([0 0.12]);
hold on;


% position accuracy average deviation
number_of_components=numel(deviation_matrix);
sum=0;
[row, column]=size(deviation_matrix);
for itera_1=1:column
    for itera_2=1:row
        
        sum=sum+abs(deviation_matrix(itera_2,itera_1));
    end
end

position_average_circle_method=sum/number_of_components;

%% box method

for iterations_1=1:m
    
    
    y_vector=y(1,iterations_1);
    for iterations_2=1:robot.steps
        robot.position=[x(1,iterations_2);y_vector;0];
        robot.L_SKM=inverse_kinematics_SKM(robot);
        
        %estimated position circle method
        position_b=box_method(robot);
        
        position_accuracy_array(iterations_2)=...
            sqrt((position_b(1)-x(1,iterations_2))^2 ...
            +(position_b(2)-y(1,iterations_1))^2);
        
    end
    position_accuracy_matrix_real(iterations_1,:)=position_accuracy_array;
end

deviation_matrix=abs(position_accuracy_matrix_ideal-position_accuracy_matrix_real);

[X,Y]=meshgrid(x,y);

figure

contourf(X,Y,deviation_matrix,15);
colorbar();

caxis([0 0.12]);
hold on;


% position accuracy average deviation
number_of_components=numel(deviation_matrix);
sum=0;
[row, column]=size(deviation_matrix);
for itera_1=1:column
    for itera_2=1:row
        
        sum=sum+abs(deviation_matrix(itera_2,itera_1));
    end
end

position_average_box_method=sum/number_of_components;


