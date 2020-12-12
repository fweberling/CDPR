function [ cable_vector ] = l_standardmodel( a,r,b,phi_x,phi_y,phi_z)
% cable vector Standard Kinematic Model
% a= vector to anchor point on the robot base(centre of winch); b=vector to
% anchor point on the platform, r= vector representing position of the
% platform fixed frame K_P w.r.t the world coordinate frame K_0;
% phi_x=rotation angle x-axis; phi_y=rotation angle y_axis; phi_z=rotation
% angle z-axis

cable_vector=a-r-R_Kardan(phi_x,phi_y,phi_z)*b;


end

