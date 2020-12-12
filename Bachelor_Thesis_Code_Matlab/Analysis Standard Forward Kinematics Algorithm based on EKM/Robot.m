classdef Robot
    % Properties' definition:
    %   general parameters for all calculations:
    %       vector a_i: defines position of a winch's centerpoint, where i is
    %           (1;2;3;4), 1 is the lower left winch, 2 is the lower right winch, 3
    %           is the upper left winch, 4 is the upper right winch
    %       vector b_i: defines the cable's attachment point on the platform, where i is
    %           (1;2;3;4), 1 is the lower left attachment, 2 is the lower right attachment, 3
    %           is the upper left attachment, 4 is the upper right attachment
    %       radius: defines the radius of all winches
    %       angle phi_x: defines orientation of platform (=0), equals 0 because
    %           of planar robot
    %       angle phi_y: defines orientation of platform (=0), equals 0 because
    %           of planar robot
    %       angle phi_z: defines orientation of platform, modifiable parameter
    %       vector position: defines position of platform's mid point in
    %           [x,y,z]-coordinates, where z equals 0 due to the planar robot
    %       angle beta_x_i: defines the orientation of the winch, equals 0 due
    %           to the planar robot
    %       angle beta_y_i: defines the orientation of the winch, equals 0 due
    %           to the planar robot
    %   special parameters for contour plots
    %       vector boundary_i: defines the workspace's boundary for the
    %           contour plots of the structure matrix
    %       steps: defines the number of increment for the contour plots of the
    %           structure matrix
    %       angle phi_upperboundary: defines the platform's maximum rotation
    %           (counter-clockwise) for contour plots of the structure matrix
    %       angle phi_lowerboundary: defines the platform's minimum rotation
    %           (clockwise) for contour plots
    %   special parameters for forward kinematics
    %       vector l: vector containing length of each cable derived from
    %           inverse kinematics
    %       coordinate x: x-coordinate of position is calculated by
    %           algorithms
    %       coordinate y: y-coordinate of position is calculated by
    %           algorithms
    %       coordinate phi: phi-coordinate, referencing the rotation, is
    %           calculated by algorithms
    %       coordinate x_0: x-coordinate of starting point for algorithms
    %       coordinate y_0: y-coordinate of starting point for algorithms
    %       coordinate phi_0: phi-coordinate of starting point for
    %           algorithms
    %       coordinate x_new: x-coordinate of position calculated during
    %           each algorithm's iteration
    %       coordinate y_new: y-coordinate of position calculated during
    %           each algorithm's iteration
    %       coordinate phi_new: phi-coordinate, referencing the rotation,
    %           calculated during each algorithm's iteration
    %       k-max: maximum number of iterations for each algorithm
    %       tau: factor of initial damping parameter for LM-algorithm,
    %           proposing x_0 is a good guess -> tau~10^-3
    %       epsilon_1: applied by LM- and DL-algorithm, proposed 10^-8
    %       epsilon_2: applied by LM- and DL-algorithm, proposed 10^-8
    %       epsilon_3: only applied by DL-algorithm, proposed 10^-8
    %       delta: radius of trust region for DL-algorithm
    %       vector f_min: contains the minimal necessary cable forces
    %       vector f_max: contains the maximal possible cable forces
    %       vector w: describes the wrench (all external forces and
    %           torques)
    % Methods's description:
    
    
    
    
    properties
        % general variables for all calculations
        a_1
        a_2
        a_3
        a_4
        b_1
        b_2
        b_3
        b_4
        radius
        phi_x
        phi_y
        phi_z
        position
        beta_x_1
        beta_x_2
        beta_x_3
        beta_x_4
        beta_y_1
        beta_y_2
        beta_y_3
        beta_y_4
        
        % special variables for contour plots
        boundary_1
        boundary_2
        boundary_3
        boundary_4
        steps
        phi_upperboundary
        phi_lowerboundary
        
        % Parameters for Levenberg-Marquardt algorithm
        L_SKM                                    % also liable for winch model
        x
        y
        phi
        
        r_x
        r_y
        r_phi
        
        x_0
        y_0
        phi_0
        x_new
        y_new
        phi_new
        k_max
        tau
        epsilon_1
        epsilon_2
        epsilon_3
        delta
        
        % Parameters for CF-method
        f_min
        f_max
        w
        
        % Parameters for EKM model
        rot_direc_1              % motor rotation direction
        rot_direc_2
        rot_direc_3
        rot_direc_4
        
        L_EKM
        
    end
    
    methods
        
        
        
        % Standard Kinematic Model (SKM)
        
        function cable_vector= inverse_kinematics_SKM( obj)
            % cable vector Standard Kinematic Model
            % a= vector to anchor point on the robot base(centre of winch); b=vector to
            % anchor point on the platform, r= vector representing position of the
            % platform fixed frame K_P w.r.t the world coordinate frame K_0;
            % phi_x=rotation angle x-axis; phi_y=rotation angle y_axis; phi_z=rotation
            % angle z-axis
            
            vector_L_1=obj.a_1-obj.position-R_Kardan(obj.phi_x,obj.phi_y,obj.phi_z)*obj.b_1;
            vector_L_2=obj.a_2-obj.position-R_Kardan(obj.phi_x,obj.phi_y,obj.phi_z)*obj.b_2;
            vector_L_3=obj.a_3-obj.position-R_Kardan(obj.phi_x,obj.phi_y,obj.phi_z)*obj.b_3;
            vector_L_4=obj.a_4-obj.position-R_Kardan(obj.phi_x,obj.phi_y,obj.phi_z)*obj.b_4;
            
            cable_vector=[norm(vector_L_1);norm(vector_L_2);norm(vector_L_3);norm(vector_L_4)];
            
        end
        function forwardkinematics_vector_function=vector_function_SKM(obj)
            
            f_1=sqrt((obj.a_1(1)-obj.r_x-cos(obj.r_phi)*obj.b_1(1)+sin(obj.r_phi)*obj.b_1(2))^2 ...
                +(obj.a_1(2)-obj.r_y-sin(obj.r_phi)*obj.b_1(1)-cos(obj.r_phi)*obj.b_1(2))^2)-obj.L_SKM(1);
            
            f_2=sqrt((obj.a_2(1)-obj.r_x-cos(obj.r_phi)*obj.b_2(1)+sin(obj.r_phi)*obj.b_2(2))^2 ...
                +(obj.a_2(2)-obj.r_y-sin(obj.r_phi)*obj.b_2(1)-cos(obj.r_phi)*obj.b_2(2))^2)-obj.L_SKM(2);
            
            f_3=sqrt((obj.a_3(1)-obj.r_x-cos(obj.r_phi)*obj.b_3(1)+sin(obj.r_phi)*obj.b_3(2))^2 ...
                +(obj.a_3(2)-obj.r_y-sin(obj.r_phi)*obj.b_3(1)-cos(obj.r_phi)*obj.b_3(2))^2)-obj.L_SKM(3);
            
            f_4=sqrt((obj.a_4(1)-obj.r_x-cos(obj.r_phi)*obj.b_4(1)+sin(obj.r_phi)*obj.b_4(2))^2 ...
                +(obj.a_4(2)-obj.r_y-sin(obj.r_phi)*obj.b_4(1)-cos(obj.r_phi)*obj.b_4(2))^2)-obj.L_SKM(4);
            
            forwardkinematics_vector_function=[f_1;f_2;f_3;f_4];
        end
        function forwardkinematics_F=objective_function_SKM(obj)
            
            f=vector_function_SKM(obj);
            
            forwardkinematics_F=1/2*(f(1)^2+f(2)^2+f(3)^2+f(4)^2);
            
        end
        function depiction=contourplot_objective_function_SKM(obj)
            
            %cable length's of robot calculated by SKM inverse kinematics
            obj.L_SKM=inverse_kinematics_SKM(obj);
            
            n=obj.steps-1;
            
            %definition of the workspace's size
            matrix_x_components=obj.boundary_2(1,1)-obj.boundary_1(1,1);
            matrix_y_components=obj.boundary_3(2,1)-obj.boundary_1(2,1);
            
            increment_step_x=matrix_x_components/n;
            increment_step_y=matrix_y_components/n;
            
            x=obj.boundary_1(1):increment_step_x:obj.boundary_2(1);
            y=obj.boundary_1(2):increment_step_y:obj.boundary_3(2);
            
            m=numel(y);
            
            
            for iterations_1=1:m
                r_y=y(1,iterations_1);
                obj.r_y=r_y;
                obj.r_phi=obj.phi_z;
                for iterations_2=1:obj.steps
                    r_x=x(1,iterations_2);
                    obj.r_x=r_x;
                    F(iterations_2)=objective_function_SKM(obj);
                end
                
                
                F_matrix(iterations_1,:)=F;
                
                
            end
            
            [X,Y]=meshgrid(x,y);
            
            
            
            pic=figure;
            contourf(X,Y,F_matrix,30);
            c=colorbar();
            heading=title('Minimum of Objective Function');
            x_axis_label=xlabel('x');
            y_axis_label=ylabel('y');
            colorbar_label=ylabel(c, 'Objective Function F');
            
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
            
            
            
            %pic=contour3(X,Y,F_matrix,obj.steps);
            
            
            
            depiction=pic;
            
        end
        function forwardkinematics_J=jacobian_matrix_SKM(obj)
            
            %standard cable length's x- and y- components
            
            %standard cable length's x- components
            L_SKM_1_x=obj.a_1(1)-obj.x-cos(obj.phi)*obj.b_1(1)+sin(obj.phi)*obj.b_1(2);
            L_SKM_2_x=obj.a_2(1)-obj.x-cos(obj.phi)*obj.b_2(1)+sin(obj.phi)*obj.b_2(2);
            L_SKM_3_x=obj.a_3(1)-obj.x-cos(obj.phi)*obj.b_3(1)+sin(obj.phi)*obj.b_3(2);
            L_SKM_4_x=obj.a_4(1)-obj.x-cos(obj.phi)*obj.b_4(1)+sin(obj.phi)*obj.b_4(2);
            
            %standard cable length's y- components
            L_SKM_1_y=obj.a_1(2)-obj.y-sin(obj.phi)*obj.b_1(1)-cos(obj.phi)*obj.b_1(2);
            L_SKM_2_y=obj.a_2(2)-obj.y-sin(obj.phi)*obj.b_2(1)-cos(obj.phi)*obj.b_2(2);
            L_SKM_3_y=obj.a_3(2)-obj.y-sin(obj.phi)*obj.b_3(1)-cos(obj.phi)*obj.b_3(2);
            L_SKM_4_y=obj.a_4(2)-obj.y-sin(obj.phi)*obj.b_4(1)-cos(obj.phi)*obj.b_4(2);
            
            %standard cable length's x- and y- components derived by phi
            del_L_SKM_1_x_del_phi=sin(obj.phi)*obj.b_1(1)+cos(obj.phi)*obj.b_1(2);
            del_L_SKM_2_x_del_phi=sin(obj.phi)*obj.b_2(1)+cos(obj.phi)*obj.b_2(2);
            del_L_SKM_3_x_del_phi=sin(obj.phi)*obj.b_3(1)+cos(obj.phi)*obj.b_3(2);
            del_L_SKM_4_x_del_phi=sin(obj.phi)*obj.b_4(1)+cos(obj.phi)*obj.b_4(2);
            
            del_L_SKM_1_y_del_phi=-cos(obj.phi)*obj.b_1(1)+sin(obj.phi)*obj.b_1(2);
            del_L_SKM_2_y_del_phi=-cos(obj.phi)*obj.b_2(1)+sin(obj.phi)*obj.b_2(2);
            del_L_SKM_3_y_del_phi=-cos(obj.phi)*obj.b_3(1)+sin(obj.phi)*obj.b_3(2);
            del_L_SKM_4_y_del_phi=-cos(obj.phi)*obj.b_4(1)+sin(obj.phi)*obj.b_4(2);
            
            %standard cable length
            L_SKM_1=sqrt(L_SKM_1_x^2+L_SKM_1_y^2);
            L_SKM_2=sqrt(L_SKM_2_x^2+L_SKM_2_y^2);
            L_SKM_3=sqrt(L_SKM_3_x^2+L_SKM_3_y^2);
            L_SKM_4=sqrt(L_SKM_4_x^2+L_SKM_4_y^2);
            
            %partial derivations for jacobian matrix
            
            del_f_1_del_x=-L_SKM_1_x/L_SKM_1;
            del_f_2_del_x=-L_SKM_2_x/L_SKM_2;
            del_f_3_del_x=-L_SKM_3_x/L_SKM_3;
            del_f_4_del_x=-L_SKM_4_x/L_SKM_4;
            
            del_f_1_del_y=-L_SKM_1_y/L_SKM_1;
            del_f_2_del_y=-L_SKM_2_y/L_SKM_2;
            del_f_3_del_y=-L_SKM_3_y/L_SKM_3;
            del_f_4_del_y=-L_SKM_4_y/L_SKM_4;
            
            del_f_1_del_phi=(L_SKM_1_x*del_L_SKM_1_x_del_phi+L_SKM_1_y*del_L_SKM_1_y_del_phi)/L_SKM_1;
            del_f_2_del_phi=(L_SKM_2_x*del_L_SKM_2_x_del_phi+L_SKM_2_y*del_L_SKM_2_y_del_phi)/L_SKM_2;
            del_f_3_del_phi=(L_SKM_3_x*del_L_SKM_3_x_del_phi+L_SKM_3_y*del_L_SKM_3_y_del_phi)/L_SKM_3;
            del_f_4_del_phi=(L_SKM_4_x*del_L_SKM_4_x_del_phi+L_SKM_4_y*del_L_SKM_4_y_del_phi)/L_SKM_4;
            
            %jacobian matrix
            forwardkinematics_J=...
                [
                del_f_1_del_x del_f_1_del_y del_f_1_del_phi; ...
                del_f_2_del_x del_f_2_del_y del_f_2_del_phi; ...
                del_f_3_del_x del_f_3_del_y del_f_3_del_phi; ...
                del_f_4_del_x del_f_4_del_y del_f_4_del_phi; ...
                ];
        end
        function [forwardkinematics_LM, iterations]=Levenberg_Marquardt_SKM(obj)
            %fixed parameters
            k=0;
            
            %starting point
            
            obj.x=obj.x_0;
            obj.y=obj.y_0;
            obj.phi=obj.phi_0;
            
            obj.r_x=obj.x;
            obj.r_y=obj.y;
            obj.r_phi=obj.phi;
            
            J=jacobian_matrix_SKM(obj);
            f=vector_function_SKM(obj);
            
            A=J'*J;
            g=J'*f;
            
            %             if (norm(obj.x_new-obj.x)<=obj.epsilon_2*(norm(obj.x)+obj.epsilon_2))
            %                 found=true;
            %             else
            %                 found=false;
            %             end
            
            
            if (norm(g,Inf)<=obj.epsilon_1)
                found=true;
            else
                found=false;
            end
            
            %mue=damping_parameter,
            a_ii=diag(A);
            mue=obj.tau*max(a_ii); %jacobian_matrix w.r.t x_0
            v=2;
            while ( ~found && k<obj.k_max)
                k=k+1;
                I=diag(diag(A));
                
                M=A+mue*I;
                
                M_det=M(1,1)*M(2,2)*M(3,3)+M(1,2)*M(2,3)*M(3,1)+M(1,3)*M(2,1)*M(3,2)...
                    -M(3,1)*M(2,2)*M(1,3)-M(3,2)*M(2,3)*M(1,1)-M(3,3)*M(2,1)*M(1,2);
                
                M_inv=1/M_det*[M(2,2)*M(3,3)-M(2,3)*M(3,2) M(1,3)*M(3,2)-M(1,2)*M(3,3) M(1,2)*M(2,3)-M(1,3)*M(2,2);...
                    M(2,3)*M(3,1)-M(2,1)*M(3,3) M(1,1)*M(3,3)-M(1,3)*M(3,1) M(1,3)*M(2,1)-M(1,1)*M(2,3);...
                    M(2,1)*M(3,2)-M(2,2)*M(3,1) M(1,2)*M(3,1)-M(1,1)*M(3,2) M(1,1)*M(2,2)-M(1,2)*M(2,1);];
                
                h_lm=-M_inv*g; %h_lm=-(A+mue*I)\g; %eye(3)
                
                if (norm(h_lm)<=obj.epsilon_2*(norm([obj.x; obj.y; obj.phi])...
                        +obj.epsilon_2))
                    found=true;
                    
                else
                    
                    F_old=1/2*(f'*f);
                    %F_old=objective_function_SKM(obj);
                    
                    obj.x_new=obj.x+h_lm(1,1);
                    obj.y_new=obj.y+h_lm(2,1);
                    obj.phi_new=obj.phi+h_lm(3,1);
                    
                    obj.r_x=obj.x_new;
                    obj.r_y=obj.y_new;
                    obj.r_phi=obj.phi_new;
                    
                    f_new=vector_function_SKM(obj);
                    
                    F_new=1/2*(f_new'*f_new);
                    %F_new=objective_function_SKM(obj);
                    
                    obj.r_x=obj.x;
                    obj.r_y=obj.y;
                    obj.r_phi=obj.phi;
                    
                    rho=(F_old-F_new)/(1/2*h_lm'*(mue*h_lm-g));
                    
                    if (rho>0)
                        
                        obj.x=obj.x_new;
                        obj.y=obj.y_new;
                        obj.phi=obj.phi_new;
                        obj.r_x=obj.x;
                        obj.r_y=obj.y;
                        obj.r_phi=obj.phi;
                        
                        
                        J_new=jacobian_matrix_SKM(obj);
                        %f_new=vector_function_SKM(obj);
                        
                        A=J_new'*J_new;
                        
                        g=J_new'*f_new;
                        
                        if (norm(g,Inf)<=obj.epsilon_1)
                            found=true;
                            
                        else
                            found=false;
                        end
                        mue=mue*max([1/3 (1-(2*rho-1)^3)]);% possibly mue_factor missing
                        v=2;
                        
                    else
                        
                        mue=mue*v;
                        v=2*v;
                    end
                end
            end
            
            X=[obj.x;obj.y;obj.phi];
            
            forwardkinematics_LM=X;
            
            iterations=k;
        end
        function [forwardkinematics_DL,iterations]=Dog_Leg_SKM(obj)
            k=0;
            
            %starting point
            obj.x=obj.x_0;
            obj.y=obj.y_0;
            obj.phi=obj.phi_0;
            
            obj.r_x=obj.x;
            obj.r_y=obj.y;
            obj.r_phi=obj.phi;
            
            J=jacobian_matrix_SKM(obj);
            f=vector_function_SKM(obj);
            g=J'*f;
            
            if((norm(f,Inf)<=obj.epsilon_3))|| ...
                    (norm(g,Inf)<=obj.epsilon_1)
                found=true;
            else
                found=false;
            end
            
            while ( (~found) && (k<obj.k_max))
                
                k=k+1;
                
                %J=jacobian_matrix_SKM(obj);
                %f=vector_function_SKM(obj);
                
                alpha=(norm(g)^2)/(norm(J*g)^2);
                
                h_sd=-g;
                
                M=J'*J;
                
                M_det=M(1,1)*M(2,2)*M(3,3)+M(1,2)*M(2,3)*M(3,1)+M(1,3)*M(2,1)*M(3,2)...
                    -M(3,1)*M(2,2)*M(1,3)-M(3,2)*M(2,3)*M(1,1)-M(3,3)*M(2,1)*M(1,2);
                
                M_inv=1/M_det*[M(2,2)*M(3,3)-M(2,3)*M(3,2) M(1,3)*M(3,2)-M(1,2)*M(3,3) M(1,2)*M(2,3)-M(1,3)*M(2,2);...
                    M(2,3)*M(3,1)-M(2,1)*M(3,3) M(1,1)*M(3,3)-M(1,3)*M(3,1) M(1,3)*M(2,1)-M(1,1)*M(2,3);...
                    M(2,1)*M(3,2)-M(2,2)*M(3,1) M(1,2)*M(3,1)-M(1,1)*M(3,2) M(1,1)*M(2,2)-M(1,2)*M(2,1);];
                
                h_gn=-M_inv*g;%h_gn=(J'*J)\(-J'*f);
                
                %trust region and Dog Leg step
                a=alpha*h_sd;
                b=h_gn;
                
                %F=objective_function_SKM(obj);
                F=1/2*(f'*f);
                
                %step
                if(norm(h_gn)<=obj.delta)
                    h_dl=h_gn;
                    LinModel=F;
                    
                elseif(norm(a)>=obj.delta)
                    h_dl=(obj.delta/norm(h_sd))*h_sd;
                    LinModel=obj.delta*(2*norm(alpha*g)-obj.delta)/(2*alpha);
                    
                else
                    c=a'*(b-a);
                    d2_a2=obj.delta^2-norm(a)^2;
                    b_a_2=norm(b-a)^2;
                    
                    if(c<=0)
                        beta=-c+sqrt(c^2+b_a_2*d2_a2)/b_a_2;
                    else
                        beta=d2_a2/(c+sqrt(c^2+b_a_2*d2_a2));
                    end
                    h_dl=a+beta*(b-a);
                    %beta chosen so that norm(h_dl)=delta
                    
                    LinModel=1/2*alpha*(1-beta)^2*norm(g)^2+beta*(2-beta)*F;
                end
                
                if(norm(h_dl)<=obj.epsilon_2*(norm([obj.x;obj.y;obj.phi])+obj.epsilon_2))
                    found=true;
                else
                    obj.x_new=obj.x+h_dl(1,1);
                    obj.y_new=obj.y+h_dl(2,1);
                    obj.phi_new=obj.phi+h_dl(3,1);
                    
                    obj.r_x=obj.x_new;
                    obj.r_y=obj.y_new;
                    obj.r_phi=obj.phi_new;
                    
                    f_new=vector_function_SKM(obj);
                    %F_new=objective_function_SKM(obj);
                    F_new=1/2*(f_new'*f_new);
                    
                    obj.r_x=obj.x;
                    obj.r_y=obj.y;
                    obj.r_phi=obj.phi;
                    
                    %gain ratio rho
                    rho=(F-F_new)/LinModel;
                    
                    if(rho>0)
                        obj.x=obj.x_new;
                        obj.y=obj.y_new;
                        obj.phi=obj.phi_new;
                        
                        obj.r_x=obj.x;
                        obj.r_y=obj.y;
                        obj.r_phi=obj.phi;
                        
                        J=jacobian_matrix_SKM(obj);
                        f=f_new;
                        
                        g=J'*f;
                        if((norm(f,Inf)<=obj.epsilon_3))|| ...
                                (norm(g,Inf)<=obj.epsilon_1)
                            found=true;
                        else
                            found=false;
                        end
                    end
                    
                    if(rho>0.75)
                        obj.delta=max(obj.delta,3*norm(h_dl));
                    elseif(rho<0.25)
                        obj.delta=obj.delta/2;
                        if(obj.delta<=obj.epsilon_2*...
                                (norm([obj.x;obj.y;obj.phi])+obj.epsilon_2))
                            found=true;
                        else
                            found=false;
                        end
                    end
                    
                end
                
            end
            forwardkinematics_DL=[obj.x;obj.y;obj.phi];
            iterations=k;
            
        end
        
        % Extended Kinematic Model (EKM)
        
        function [q_new,q_f_new,q_delta]=inverse_kinematics_EKM(obj)
            % starting position
            
            % 0 stands for start position
            
            % calculation of cable vectors (standard model)
            vector_L_1_0=l_standardmodel(obj.a_1,[0;0;0],obj.b_1,...
                0,0,0);
            vector_L_2_0=l_standardmodel(obj.a_2,[0;0;0],obj.b_2,...
                0,0,0);
            vector_L_3_0=l_standardmodel(obj.a_3,[0;0;0],obj.b_3,...
                0,0,0);
            vector_L_4_0=l_standardmodel(obj.a_4,[0;0;0],obj.b_4,...
                0,0,0);
            
            %calculation angle epsilon_0
            beta_0_L_1=atan(vector_L_1_0(2)/vector_L_1_0(1));
            beta_0_L_2=atan(vector_L_2_0(2)/vector_L_2_0(1));
            beta_0_L_3=atan(vector_L_3_0(2)/vector_L_3_0(1));
            beta_0_L_4=atan(vector_L_4_0(2)/vector_L_4_0(1));
            
            
            % new position
            
            % calculation of cable vectors (standard model)
            vector_L_1_new=l_standardmodel(obj.a_1,obj.position,obj.b_1,...
                obj.phi_x,obj.phi_y,obj.phi_z);
            vector_L_2_new=l_standardmodel(obj.a_2,obj.position,obj.b_2,...
                obj.phi_x,obj.phi_y,obj.phi_z);
            vector_L_3_new=l_standardmodel(obj.a_3,obj.position,obj.b_3,...
                obj.phi_x,obj.phi_y,obj.phi_z);
            vector_L_4_new=l_standardmodel(obj.a_4,obj.position,obj.b_4,...
                obj.phi_x,obj.phi_y,obj.phi_z);
            
            % cable lengths (euclidean norm)
            L_1_SKM_new=norm(vector_L_1_new);
            L_2_SKM_new=norm(vector_L_2_new);
            L_3_SKM_new=norm(vector_L_3_new);
            L_4_SKM_new=norm(vector_L_4_new);
            
            % free cable lengths
            L_f_1_new=sqrt(L_1_SKM_new^2-obj.radius^2);
            L_f_2_new=sqrt(L_2_SKM_new^2-obj.radius^2);
            L_f_3_new=sqrt(L_3_SKM_new^2-obj.radius^2);
            L_f_4_new=sqrt(L_4_SKM_new^2-obj.radius^2);
            
            % vector q_f_new
            q_f_new=[L_f_1_new;L_f_2_new;L_f_3_new;L_f_4_new;];
            
            %calculation angle delta
            delta_1_new=acos(obj.radius/L_1_SKM_new);
            delta_2_new=acos(obj.radius/L_2_SKM_new);
            delta_3_new=acos(obj.radius/L_3_SKM_new);
            delta_4_new=acos(obj.radius/L_4_SKM_new);
            
            % q_delta_new
            q_delta=[delta_1_new;delta_2_new;delta_3_new;...
                delta_4_new;];
            
            % cable vector l_i_new in new coordinate frame
            L_SKM_beta_1=R_Kardan(0,0,-beta_0_L_1)*vector_L_1_new;
            L_SKM_beta_2=R_Kardan(0,0,-beta_0_L_2)*vector_L_2_new;
            L_SKM_beta_3=R_Kardan(0,0,-beta_0_L_3)*vector_L_3_new;
            L_SKM_beta_4=R_Kardan(0,0,-beta_0_L_4)*vector_L_4_new;
            
            % angle rho between cable_vector_l_i_start and cable_vector_l_i_new
            theta_1=atan(L_SKM_beta_1(2)/L_SKM_beta_1(1));
            theta_2=atan(L_SKM_beta_2(2)/L_SKM_beta_2(1));
            theta_3=atan(L_SKM_beta_3(2)/L_SKM_beta_3(1));
            theta_4=atan(L_SKM_beta_4(2)/L_SKM_beta_4(1));
            
            % cable length wrapped around winch
            L_delta_1=(pi-delta_1_new)*obj.radius;
            L_delta_2=(pi-delta_2_new)*obj.radius;
            L_delta_3=(pi-delta_3_new)*obj.radius;
            L_delta_4=(pi-delta_4_new)*obj.radius;
            
            % cable length correcting factor
            L_corr_1=-obj.rot_direc_1*obj.radius*theta_1;
            L_corr_2=-obj.rot_direc_2*obj.radius*theta_2;
            L_corr_3=-obj.rot_direc_3*obj.radius*theta_3;
            L_corr_4=-obj.rot_direc_4*obj.radius*theta_4;
            
            % calculation of cable lengths in new position
            L_1_new=L_f_1_new+L_delta_1+L_corr_1;
            L_2_new=L_f_2_new+L_delta_2+L_corr_2;
            L_3_new=L_f_3_new+L_delta_3+L_corr_3;
            L_4_new=L_f_4_new+L_delta_4+L_corr_4;
            % l_corrected is not identical with Jean's model
            
            % inverse kinematic starting position
            q_new=[L_1_new; L_2_new; L_3_new; L_4_new;];
            
            
            
        end
        function forwardkinematics_vector_function=vector_function_EKM(obj)
            
            
            
            
            %start configuration
            
            % calculation of cable vectors (standard model)
            vector_L_1_0=l_standardmodel(obj.a_1,[0;0;0],obj.b_1,0,0,0);
            vector_L_2_0=l_standardmodel(obj.a_2,[0;0;0],obj.b_2,0,0,0);
            vector_L_3_0=l_standardmodel(obj.a_3,[0;0;0],obj.b_3,0,0,0);
            vector_L_4_0=l_standardmodel(obj.a_4,[0;0;0],obj.b_4,0,0,0);
            
            % calculation epsilon start
            beta_0_L_1=atan(vector_L_1_0(2)/(vector_L_1_0(1)));
            beta_0_L_2=atan(vector_L_2_0(2)/(vector_L_2_0(1)));
            beta_0_L_3=atan(vector_L_3_0(2)/(vector_L_3_0(1)));
            beta_0_L_4=atan(vector_L_4_0(2)/(vector_L_4_0(1)));
            
            
            %cable new configuration
            
            % calculation of cable vectors (standard model)
            vector_L_1_new=l_standardmodel(obj.a_1,[obj.r_x;obj.r_y;0],obj.b_1,...
                obj.phi_x,obj.phi_y,obj.r_phi);
            vector_L_2_new=l_standardmodel(obj.a_2,[obj.r_x;obj.r_y;0],obj.b_2,...
                obj.phi_x,obj.phi_y,obj.r_phi);
            vector_L_3_new=l_standardmodel(obj.a_3,[obj.r_x;obj.r_y;0],obj.b_3,...
                obj.phi_x,obj.phi_y,obj.r_phi);
            vector_L_4_new=l_standardmodel(obj.a_4,[obj.r_x;obj.r_y;0],obj.b_4,...
                obj.phi_x,obj.phi_y,obj.r_phi);
            
            % cable vector l_i_new in new coordinate frame
            L_SKM_beta_1=R_Kardan(0,0,-beta_0_L_1)*vector_L_1_new;
            L_SKM_beta_2=R_Kardan(0,0,-beta_0_L_2)*vector_L_2_new;
            L_SKM_beta_3=R_Kardan(0,0,-beta_0_L_3)*vector_L_3_new;
            L_SKM_beta_4=R_Kardan(0,0,-beta_0_L_4)*vector_L_4_new;
            
            % angle rho between cable_vector_l_i_start and cable_vector_l_i_new
            theta_1=atan(L_SKM_beta_1(2)/L_SKM_beta_1(1));
            theta_2=atan(L_SKM_beta_2(2)/L_SKM_beta_2(1));
            theta_3=atan(L_SKM_beta_3(2)/L_SKM_beta_3(1));
            theta_4=atan(L_SKM_beta_4(2)/L_SKM_beta_4(1));
            
            
            L_SKM_1_x=obj.a_1(1)-obj.r_x-cos(obj.r_phi)*obj.b_1(1)+sin(obj.r_phi)*obj.b_1(2);
            L_SKM_2_x=obj.a_2(1)-obj.r_x-cos(obj.r_phi)*obj.b_2(1)+sin(obj.r_phi)*obj.b_2(2);
            L_SKM_3_x=obj.a_3(1)-obj.r_x-cos(obj.r_phi)*obj.b_3(1)+sin(obj.r_phi)*obj.b_3(2);
            L_SKM_4_x=obj.a_4(1)-obj.r_x-cos(obj.r_phi)*obj.b_4(1)+sin(obj.r_phi)*obj.b_4(2);
            
            L_SKM_1_y=obj.a_1(2)-obj.r_y-sin(obj.r_phi)*obj.b_1(1)-cos(obj.r_phi)*obj.b_1(2);
            L_SKM_2_y=obj.a_2(2)-obj.r_y-sin(obj.r_phi)*obj.b_2(1)-cos(obj.r_phi)*obj.b_2(2);
            L_SKM_3_y=obj.a_3(2)-obj.r_y-sin(obj.r_phi)*obj.b_3(1)-cos(obj.r_phi)*obj.b_3(2);
            L_SKM_4_y=obj.a_4(2)-obj.r_y-sin(obj.r_phi)*obj.b_4(1)-cos(obj.r_phi)*obj.b_4(2);
            
            %free cable lengths
            L_f_1=sqrt((L_SKM_1_x)^2+(L_SKM_1_y)^2-obj.radius^2);
            L_f_2=sqrt((L_SKM_2_x)^2+(L_SKM_2_y)^2-obj.radius^2);
            L_f_3=sqrt((L_SKM_3_x)^2+(L_SKM_3_y)^2-obj.radius^2);
            L_f_4=sqrt((L_SKM_4_x)^2+(L_SKM_4_y)^2-obj.radius^2);
            
            %cable length wrapped around winch
            L_delta_1=obj.radius*(pi-acos(obj.radius/sqrt((L_SKM_1_x)^2 ...
                +(L_SKM_1_y)^2)));
            L_delta_2=obj.radius*(pi-acos(obj.radius/sqrt((L_SKM_2_x)^2 ...
                +(L_SKM_2_y)^2)));
            L_delta_3=obj.radius*(pi-acos(obj.radius/sqrt((L_SKM_3_x)^2 ...
                +(L_SKM_3_y)^2)));
            L_delta_4=obj.radius*(pi-acos(obj.radius/sqrt((L_SKM_4_x)^2 ...
                +(L_SKM_4_y)^2)));
            %cable length correcting term
            L_corr_1=-obj.rot_direc_1*obj.radius*theta_1;
            L_corr_2=-obj.rot_direc_2*obj.radius*theta_2;
            L_corr_3=-obj.rot_direc_3*obj.radius*theta_3;
            L_corr_4=-obj.rot_direc_4*obj.radius*theta_4;
            
            %vector function components
            f_1=L_f_1+L_delta_1+L_corr_1-obj.L_EKM(1);
            f_2=L_f_2+L_delta_2+L_corr_2-obj.L_EKM(2);
            f_3=L_f_3+L_delta_3+L_corr_3-obj.L_EKM(3);
            f_4=L_f_4+L_delta_4+L_corr_4-obj.L_EKM(4);
            %vector function
            forwardkinematics_vector_function=[f_1;f_2;f_3;f_4];
        end
        function forwardkinematics_F=objective_function_EKM(obj)
            
            f=vector_function_EKM(obj);
            
            forwardkinematics_F=1/2*(f(1)^2+f(2)^2+f(3)^2+f(4)^2);
            
        end
        function depiction=contourplot_objective_function_EKM(obj)
            
            %cable length's of robot calculated by EKM inverse kinematics
            obj.L_EKM=inverse_kinematics_EKM(obj);
            
            n=obj.steps-1;
            
            %definition of the workspace's size
            matrix_x_components=obj.boundary_2(1,1)-obj.boundary_1(1,1);
            matrix_y_components=obj.boundary_3(2,1)-obj.boundary_1(2,1);
            
            increment_step_x=matrix_x_components/n;
            increment_step_y=matrix_y_components/n;
            
            x=obj.boundary_1(1):increment_step_x:obj.boundary_2(1);
            y=obj.boundary_1(2):increment_step_y:obj.boundary_3(2);
            
            m=numel(y);
            
            
            for iterations_1=1:m
                r_y=y(1,iterations_1);
                obj.r_y=r_y;
                obj.r_phi=obj.phi_z;
                for iterations_2=1:obj.steps
                    r_x=x(1,iterations_2);
                    obj.r_x=r_x;
                    F(iterations_2)=objective_function_EKM(obj);
                end
                
                
                F_matrix(iterations_1,:)=F;
                
                
            end
            
            [X,Y]=meshgrid(x,y);
            
            pic=contourf(X,Y,F_matrix,obj.steps);
            
            %pic=contour3(X,Y,F_matrix,obj.steps);
            
            colorbar();
            
            depiction=pic;
            
        end
        function forwardkinematics_J=jacobian_matrix_EKM(obj)
            
            %start configuration determination angle epsilon_0
            
            % calculation of cable vectors (standard model)
            vector_L_1_0=l_standardmodel(obj.a_1,[0;0;0],obj.b_1,0,0,0);
            vector_L_2_0=l_standardmodel(obj.a_2,[0;0;0],obj.b_2,0,0,0);
            vector_L_3_0=l_standardmodel(obj.a_3,[0;0;0],obj.b_3,0,0,0);
            vector_L_4_0=l_standardmodel(obj.a_4,[0;0;0],obj.b_4,0,0,0);
            
            % calculation epsilon start
            beta_0_L_1=atan(vector_L_1_0(2)/(vector_L_1_0(1)));
            beta_0_L_2=atan(vector_L_2_0(2)/(vector_L_2_0(1)));
            beta_0_L_3=atan(vector_L_3_0(2)/(vector_L_3_0(1)));
            beta_0_L_4=atan(vector_L_4_0(2)/(vector_L_4_0(1)));
            
            
            % standard cable length's x and y components
            L_SKM_1_x=(obj.a_1(1)-obj.x-cos(obj.phi)*obj.b_1(1)+sin(obj.phi)*obj.b_1(2));
            L_SKM_2_x=(obj.a_2(1)-obj.x-cos(obj.phi)*obj.b_2(1)+sin(obj.phi)*obj.b_2(2));
            L_SKM_3_x=(obj.a_3(1)-obj.x-cos(obj.phi)*obj.b_3(1)+sin(obj.phi)*obj.b_3(2));
            L_SKM_4_x=(obj.a_4(1)-obj.x-cos(obj.phi)*obj.b_4(1)+sin(obj.phi)*obj.b_4(2));
            L_SKM_1_y=(obj.a_1(2)-obj.y-sin(obj.phi)*obj.b_1(1)-cos(obj.phi)*obj.b_1(2));
            L_SKM_2_y=(obj.a_2(2)-obj.y-sin(obj.phi)*obj.b_2(1)-cos(obj.phi)*obj.b_2(2));
            L_SKM_3_y=(obj.a_3(2)-obj.y-sin(obj.phi)*obj.b_3(1)-cos(obj.phi)*obj.b_3(2));
            L_SKM_4_y=(obj.a_4(2)-obj.y-sin(obj.phi)*obj.b_4(1)-cos(obj.phi)*obj.b_4(2));
            
            %standard cable length
            L_SKM_1=sqrt(L_SKM_1_x^2+L_SKM_1_y^2);
            L_SKM_2=sqrt(L_SKM_2_x^2+L_SKM_2_y^2);
            L_SKM_3=sqrt(L_SKM_3_x^2+L_SKM_3_y^2);
            L_SKM_4=sqrt(L_SKM_4_x^2+L_SKM_4_y^2);
            
            %standard cable length's x and y components derived by phi
            del_L_SKM_1_x_del_phi=(sin(obj.phi)*obj.b_1(1)+cos(obj.phi)*obj.b_1(2));
            del_L_SKM_2_x_del_phi=(sin(obj.phi)*obj.b_2(1)+cos(obj.phi)*obj.b_2(2));
            del_L_SKM_3_x_del_phi=(sin(obj.phi)*obj.b_3(1)+cos(obj.phi)*obj.b_3(2));
            del_L_SKM_4_x_del_phi=(sin(obj.phi)*obj.b_4(1)+cos(obj.phi)*obj.b_4(2));
            del_L_SKM_1_y_del_phi=(-cos(obj.phi)*obj.b_1(1)+sin(obj.phi)*obj.b_1(2));
            del_L_SKM_2_y_del_phi=(-cos(obj.phi)*obj.b_2(1)+sin(obj.phi)*obj.b_2(2));
            del_L_SKM_3_y_del_phi=(-cos(obj.phi)*obj.b_3(1)+sin(obj.phi)*obj.b_3(2));
            del_L_SKM_4_y_del_phi=(-cos(obj.phi)*obj.b_4(1)+sin(obj.phi)*obj.b_4(2));
            
            %free cable length for Extended Kinematic Model
            L_f_1=sqrt((L_SKM_1_x)^2+(L_SKM_1_y)^2-obj.radius^2);
            L_f_2=sqrt((L_SKM_2_x)^2+(L_SKM_2_y)^2-obj.radius^2);
            L_f_3=sqrt((L_SKM_3_x)^2+(L_SKM_3_y)^2-obj.radius^2);
            L_f_4=sqrt((L_SKM_4_x)^2+(L_SKM_4_y)^2-obj.radius^2);
            
            %derivations of the free cable length
            
            %free cable length derived by x
            del_L_f_1_del_x=-L_SKM_1_x/L_f_1;
            del_L_f_2_del_x=-L_SKM_2_x/L_f_2;
            del_L_f_3_del_x=-L_SKM_3_x/L_f_3;
            del_L_f_4_del_x=-L_SKM_4_x/L_f_4;
            
            %free cable length derived by y
            del_L_f_1_del_y=-L_SKM_1_y/L_f_1;
            del_L_f_2_del_y=-L_SKM_2_y/L_f_2;
            del_L_f_3_del_y=-L_SKM_3_y/L_f_3;
            del_L_f_4_del_y=-L_SKM_4_y/L_f_4;
            
            %free cable length derived by phi
            del_L_f_1_del_phi=(L_SKM_1_x*del_L_SKM_1_x_del_phi+L_SKM_1_y*del_L_SKM_1_y_del_phi)...
                /L_f_1;
            del_L_f_2_del_phi=(L_SKM_2_x*del_L_SKM_2_x_del_phi+L_SKM_2_y*del_L_SKM_2_y_del_phi)...
                /L_f_2;
            del_L_f_3_del_phi=(L_SKM_3_x*del_L_SKM_3_x_del_phi+L_SKM_3_y*del_L_SKM_3_y_del_phi)...
                /L_f_3;
            del_L_f_4_del_phi=(L_SKM_4_x*del_L_SKM_4_x_del_phi+L_SKM_4_y*del_L_SKM_4_y_del_phi)...
                /L_f_4;
            
            %derivations of cable length wrapped around winch
            
            %cable length wrapped around winch derived by x
            del_L_delta_1_del_x=obj.radius^2*L_SKM_1_x/...
                sqrt((1-(obj.radius/L_SKM_1)^2)*(L_SKM_1_x^2+L_SKM_1_y^2)^3);
            del_L_delta_2_del_x=obj.radius^2*L_SKM_2_x/...
                sqrt((1-(obj.radius/L_SKM_2)^2)*(L_SKM_2_x^2+L_SKM_2_y^2)^3);
            del_L_delta_3_del_x=obj.radius^2*L_SKM_3_x/...
                sqrt((1-(obj.radius/L_SKM_3)^2)*(L_SKM_3_x^2+L_SKM_3_y^2)^3);
            del_L_delta_4_del_x=obj.radius^2*L_SKM_4_x/...
                sqrt((1-(obj.radius/L_SKM_4)^2)*(L_SKM_4_x^2+L_SKM_4_y^2)^3);
            
            %cable length wrapped around winch derived by y
            del_L_delta_1_del_y=obj.radius^2*L_SKM_1_y/...
                sqrt((1-(obj.radius/L_SKM_1)^2)*(L_SKM_1_x^2+L_SKM_1_y^2)^3);
            del_L_delta_2_del_y=obj.radius^2*L_SKM_2_y/...
                sqrt((1-(obj.radius/L_SKM_2)^2)*(L_SKM_2_x^2+L_SKM_2_y^2)^3);
            del_L_delta_3_del_y=obj.radius^2*L_SKM_3_y/...
                sqrt((1-(obj.radius/L_SKM_3)^2)*(L_SKM_3_x^2+L_SKM_3_y^2)^3);
            del_L_delta_4_del_y=obj.radius^2*L_SKM_4_y/...
                sqrt((1-(obj.radius/L_SKM_4)^2)*(L_SKM_4_x^2+L_SKM_4_y^2)^3);
            
            %cable length wrapped around winch derived by phi
            del_L_delta_1_del_phi=-obj.radius^2*(L_SKM_1_x*del_L_SKM_1_x_del_phi...
                +L_SKM_1_y*del_L_SKM_1_y_del_phi)/...
                sqrt((1-(obj.radius/L_SKM_1)^2)*(L_SKM_1_x^2+L_SKM_1_y^2)^3);
            del_L_delta_2_del_phi=-obj.radius^2*(L_SKM_2_x*del_L_SKM_2_x_del_phi...
                +L_SKM_2_y*del_L_SKM_2_y_del_phi)/...
                sqrt((1-(obj.radius/L_SKM_2)^2)*(L_SKM_2_x^2+L_SKM_2_y^2)^3);
            del_L_delta_3_del_phi=-obj.radius^2*(L_SKM_3_x*del_L_SKM_3_x_del_phi...
                +L_SKM_3_y*del_L_SKM_3_y_del_phi)/...
                sqrt((1-(obj.radius/L_SKM_3)^2)*(L_SKM_3_x^2+L_SKM_3_y^2)^3);
            del_L_delta_4_del_phi=-obj.radius^2*(L_SKM_4_x*del_L_SKM_4_x_del_phi...
                +L_SKM_4_y*del_L_SKM_4_y_del_phi)/...
                sqrt((1-(obj.radius/L_SKM_4)^2)*(L_SKM_4_x^2+L_SKM_4_y^2)^3);
            
            %derivations of cable length correcting term
            
            %rotatated standard cable length's x- and y-components
            L_SKM_beta_1_x=cos(-beta_0_L_1)*L_SKM_1_x-sin(-beta_0_L_1)*L_SKM_1_y;
            L_SKM_beta_2_x=cos(-beta_0_L_2)*L_SKM_2_x-sin(-beta_0_L_2)*L_SKM_2_y;
            L_SKM_beta_3_x=cos(-beta_0_L_3)*L_SKM_3_x-sin(-beta_0_L_3)*L_SKM_3_y;
            L_SKM_beta_4_x=cos(-beta_0_L_4)*L_SKM_4_x-sin(-beta_0_L_4)*L_SKM_4_y;
            L_SKM_beta_1_y=sin(-beta_0_L_1)*L_SKM_1_x+cos(-beta_0_L_1)*L_SKM_1_y;
            L_SKM_beta_2_y=sin(-beta_0_L_2)*L_SKM_2_x+cos(-beta_0_L_2)*L_SKM_2_y;
            L_SKM_beta_3_y=sin(-beta_0_L_3)*L_SKM_3_x+cos(-beta_0_L_3)*L_SKM_3_y;
            L_SKM_beta_4_y=sin(-beta_0_L_4)*L_SKM_4_x+cos(-beta_0_L_4)*L_SKM_4_y;
            
            %rotated standard cable length's x-components derived by x
            del_L_SKM_beta_1_x_del_x=-cos(-beta_0_L_1);
            del_L_SKM_beta_2_x_del_x=-cos(-beta_0_L_2);
            del_L_SKM_beta_3_x_del_x=-cos(-beta_0_L_3);
            del_L_SKM_beta_4_x_del_x=-cos(-beta_0_L_4);
            
            %rotated standard cable length's y-components derived by x
            del_L_SKM_beta_1_y_del_x=-sin(-beta_0_L_1);
            del_L_SKM_beta_2_y_del_x=-sin(-beta_0_L_2);
            del_L_SKM_beta_3_y_del_x=-sin(-beta_0_L_3);
            del_L_SKM_beta_4_y_del_x=-sin(-beta_0_L_4);
            
            %rotated standard cable length's x-components derived by y
            del_L_SKM_beta_1_x_del_y=sin(-beta_0_L_1);
            del_L_SKM_beta_2_x_del_y=sin(-beta_0_L_2);
            del_L_SKM_beta_3_x_del_y=sin(-beta_0_L_3);
            del_L_SKM_beta_4_x_del_y=sin(-beta_0_L_4);
            
            %rotated standard cable length's y-components derived by y
            del_L_SKM_beta_1_y_del_y=-cos(-beta_0_L_1);
            del_L_SKM_beta_2_y_del_y=-cos(-beta_0_L_2);
            del_L_SKM_beta_3_y_del_y=-cos(-beta_0_L_3);
            del_L_SKM_beta_4_y_del_y=-cos(-beta_0_L_4);
            
            %rotated standard cable length's x-components derived by phi
            del_L_SKM_beta_1_x_del_phi=cos(-beta_0_L_1)*del_L_SKM_1_x_del_phi...
                -sin(-beta_0_L_1)*del_L_SKM_1_y_del_phi;
            del_L_SKM_beta_2_x_del_phi=cos(-beta_0_L_2)*del_L_SKM_2_x_del_phi...
                -sin(-beta_0_L_2)*del_L_SKM_2_y_del_phi;
            del_L_SKM_beta_3_x_del_phi=cos(-beta_0_L_3)*del_L_SKM_3_x_del_phi...
                -sin(-beta_0_L_3)*del_L_SKM_3_y_del_phi;
            del_L_SKM_beta_4_x_del_phi=cos(-beta_0_L_4)*del_L_SKM_4_x_del_phi...
                -sin(-beta_0_L_4)*del_L_SKM_4_y_del_phi;
            
            %rotated standard cable length's y-components derived by phi
            del_L_SKM_beta_1_y_del_phi=sin(-beta_0_L_1)*del_L_SKM_1_x_del_phi...
                +cos(-beta_0_L_1)*del_L_SKM_1_y_del_phi;
            del_L_SKM_beta_2_y_del_phi=sin(-beta_0_L_2)*del_L_SKM_2_x_del_phi...
                +cos(-beta_0_L_2)*del_L_SKM_2_y_del_phi;
            del_L_SKM_beta_3_y_del_phi=sin(-beta_0_L_3)*del_L_SKM_3_x_del_phi...
                +cos(-beta_0_L_3)*del_L_SKM_3_y_del_phi;
            del_L_SKM_beta_4_y_del_phi=sin(-beta_0_L_4)*del_L_SKM_4_x_del_phi...
                +cos(-beta_0_L_4)*del_L_SKM_4_y_del_phi;
            
            %cable length correcting term derived by x
            del_L_corr_1_del_x=-obj.rot_direc_1*obj.radius*...
                (L_SKM_beta_1_x*del_L_SKM_beta_1_y_del_x-L_SKM_beta_1_y*del_L_SKM_beta_1_x_del_x)/...
                ((1+(L_SKM_beta_1_y/L_SKM_beta_1_x))*L_SKM_beta_1_x)^2;
            del_L_corr_2_del_x=-obj.rot_direc_2*obj.radius*...
                (L_SKM_beta_2_x*del_L_SKM_beta_2_y_del_x-L_SKM_beta_2_y*del_L_SKM_beta_2_x_del_x)/...
                ((1+(L_SKM_beta_2_y/L_SKM_beta_2_x))*L_SKM_beta_2_x)^2;
            del_L_corr_3_del_x=-obj.rot_direc_3*obj.radius*...
                (L_SKM_beta_3_x*del_L_SKM_beta_3_y_del_x-L_SKM_beta_3_y*del_L_SKM_beta_3_x_del_x)/...
                ((1+(L_SKM_beta_3_y/L_SKM_beta_3_x))*L_SKM_beta_3_x)^2;
            del_L_corr_4_del_x=-obj.rot_direc_4*obj.radius*...
                (L_SKM_beta_4_x*del_L_SKM_beta_4_y_del_x-L_SKM_beta_4_y*del_L_SKM_beta_4_x_del_x)/...
                ((1+(L_SKM_beta_4_y/L_SKM_beta_4_x))*L_SKM_beta_4_x)^2;
            
            %cable length correcting term derived by y
            del_L_corr_1_del_y=-obj.rot_direc_1*obj.radius*...
                (L_SKM_beta_1_x*del_L_SKM_beta_1_y_del_y-L_SKM_beta_1_y*del_L_SKM_beta_1_x_del_y)/...
                ((1+(L_SKM_beta_1_y/L_SKM_beta_1_x))*L_SKM_beta_1_x)^2;
            del_L_corr_2_del_y=-obj.rot_direc_2*obj.radius*...
                (L_SKM_beta_2_x*del_L_SKM_beta_2_y_del_y-L_SKM_beta_2_y*del_L_SKM_beta_2_x_del_y)/...
                ((1+(L_SKM_beta_2_y/L_SKM_beta_2_x))*L_SKM_beta_2_x)^2;
            del_L_corr_3_del_y=-obj.rot_direc_3*obj.radius*...
                (L_SKM_beta_3_x*del_L_SKM_beta_3_y_del_y-L_SKM_beta_3_y*del_L_SKM_beta_3_x_del_y)/...
                ((1+(L_SKM_beta_3_y/L_SKM_beta_3_x))*L_SKM_beta_3_x)^2;
            del_L_corr_4_del_y=-obj.rot_direc_4*obj.radius*...
                (L_SKM_beta_4_x*del_L_SKM_beta_4_y_del_y-L_SKM_beta_4_y*del_L_SKM_beta_4_x_del_y)/...
                ((1+(L_SKM_beta_4_y/L_SKM_beta_4_x))*L_SKM_beta_4_x)^2;
            
            %cable length correcting term derived by phi
            del_L_corr_1_del_phi=-obj.rot_direc_1*obj.radius*...
                (L_SKM_beta_1_x*del_L_SKM_beta_1_y_del_phi-L_SKM_beta_1_y*del_L_SKM_beta_1_x_del_phi)/...
                ((1+(L_SKM_beta_1_y/L_SKM_beta_1_x))*L_SKM_beta_1_x)^2;
            del_L_corr_2_del_phi=-obj.rot_direc_2*obj.radius*...
                (L_SKM_beta_2_x*del_L_SKM_beta_2_y_del_phi-L_SKM_beta_2_y*del_L_SKM_beta_2_x_del_phi)/...
                ((1+(L_SKM_beta_2_y/L_SKM_beta_2_x))*L_SKM_beta_2_x)^2;
            del_L_corr_3_del_phi=-obj.rot_direc_3*obj.radius*...
                (L_SKM_beta_3_x*del_L_SKM_beta_3_y_del_phi-L_SKM_beta_3_y*del_L_SKM_beta_3_x_del_phi)/...
                ((1+(L_SKM_beta_3_y/L_SKM_beta_3_x))*L_SKM_beta_3_x)^2;
            del_L_corr_4_del_phi=-obj.rot_direc_4*obj.radius*...
                (L_SKM_beta_4_x*del_L_SKM_beta_4_y_del_phi-L_SKM_beta_4_y*del_L_SKM_beta_4_x_del_phi)/...
                ((1+(L_SKM_beta_4_y/L_SKM_beta_4_x))*L_SKM_beta_4_x)^2;
            
            
            %partial derivations for jacobian matrix
            
            del_f_1_del_x=del_L_f_1_del_x+del_L_delta_1_del_x+del_L_corr_1_del_x;
            del_f_2_del_x=del_L_f_2_del_x+del_L_delta_2_del_x+del_L_corr_2_del_x;
            del_f_3_del_x=del_L_f_3_del_x+del_L_delta_3_del_x+del_L_corr_3_del_x;
            del_f_4_del_x=del_L_f_4_del_x+del_L_delta_4_del_x+del_L_corr_4_del_x;
            
            del_f_1_del_y=del_L_f_1_del_y+del_L_delta_1_del_y+del_L_corr_1_del_y;
            del_f_2_del_y=del_L_f_2_del_y+del_L_delta_2_del_y+del_L_corr_2_del_y;
            del_f_3_del_y=del_L_f_3_del_y+del_L_delta_3_del_y+del_L_corr_3_del_y;
            del_f_4_del_y=del_L_f_4_del_y+del_L_delta_4_del_y+del_L_corr_4_del_y;
            
            del_f_1_del_phi=del_L_f_1_del_phi+del_L_delta_1_del_phi+del_L_corr_1_del_phi;
            del_f_2_del_phi=del_L_f_2_del_phi+del_L_delta_2_del_phi+del_L_corr_2_del_phi;
            del_f_3_del_phi=del_L_f_3_del_phi+del_L_delta_3_del_phi+del_L_corr_3_del_phi;
            del_f_4_del_phi=del_L_f_4_del_phi+del_L_delta_4_del_phi+del_L_corr_4_del_phi;
            
            %jacobian matrix
            
            forwardkinematics_J=...
                [
                del_f_1_del_x del_f_1_del_y del_f_1_del_phi; ...
                del_f_2_del_x del_f_2_del_y del_f_2_del_phi; ...
                del_f_3_del_x del_f_3_del_y del_f_3_del_phi; ...
                del_f_4_del_x del_f_4_del_y del_f_4_del_phi; ...
                ];
        end
        function [forwardkinematics_LM, iterations]=Levenberg_Marquardt_EKM(obj)
            %fixed parameters
            k=0;
            
            %starting point
            
            obj.x=obj.x_0;
            obj.y=obj.y_0;
            obj.phi=obj.phi_0;
            
            obj.r_x=obj.x;
            obj.r_y=obj.y;
            obj.r_phi=obj.phi;
            
            J=jacobian_matrix_EKM(obj);
            f=vector_function_EKM(obj);
            
            A=J'*J;
            g=J'*f;
            
            %             if (norm(obj.x_new-obj.x)<=obj.epsilon_2*(norm(obj.x)+obj.epsilon_2))
            %                 found=true;
            %             else
            %                 found=false;
            %             end
            
            
            if (norm(g,Inf)<=obj.epsilon_1)
                found=true;
            else
                found=false;
            end
            
            %mue=damping_parameter,
            a_ii=diag(A);
            mue=obj.tau*max(a_ii); %jacobian_matrix w.r.t x_0
            v=2;
            while ( ~found && k<obj.k_max)
                k=k+1;
                I=diag(diag(A));
                
                M=A+mue*I;
                
                M_det=M(1,1)*M(2,2)*M(3,3)+M(1,2)*M(2,3)*M(3,1)+M(1,3)*M(2,1)*M(3,2)...
                    -M(3,1)*M(2,2)*M(1,3)-M(3,2)*M(2,3)*M(1,1)-M(3,3)*M(2,1)*M(1,2);
                
                M_inv=1/M_det*[M(2,2)*M(3,3)-M(2,3)*M(3,2) M(1,3)*M(3,2)-M(1,2)*M(3,3) M(1,2)*M(2,3)-M(1,3)*M(2,2);...
                    M(2,3)*M(3,1)-M(2,1)*M(3,3) M(1,1)*M(3,3)-M(1,3)*M(3,1) M(1,3)*M(2,1)-M(1,1)*M(2,3);...
                    M(2,1)*M(3,2)-M(2,2)*M(3,1) M(1,2)*M(3,1)-M(1,1)*M(3,2) M(1,1)*M(2,2)-M(1,2)*M(2,1);];
                
                h_lm=-M_inv*g; %h_lm=-(A+mue*I)\g; %eye(3)
                
                if (norm(h_lm)<=obj.epsilon_2*(norm([obj.x; obj.y; obj.phi])...
                        +obj.epsilon_2))
                    found=true;
                    
                    
                else
                    
                    %F_old=objective_function_EKM(obj);
                    F_old=1/2*(f'*f);
                    
                    
                    obj.x_new=obj.x+h_lm(1,1);
                    obj.y_new=obj.y+h_lm(2,1);
                    obj.phi_new=obj.phi+h_lm(3,1);
                    
                    obj.r_x=obj.x_new;
                    obj.r_y=obj.y_new;
                    obj.r_phi=obj.phi_new;
                    
                    f_new=vector_function_EKM(obj);
                    %F_new=objective_function_EKM(obj);
                    F_new=1/2*(f_new'*f_new);
                    
                    obj.r_x=obj.x;
                    obj.r_y=obj.y;
                    obj.r_phi=obj.phi;
                    
                    rho=(F_old-F_new)/(1/2*h_lm'*(mue*h_lm-g));
                    
                    if (rho>0)
                        
                        obj.x=obj.x_new;
                        obj.y=obj.y_new;
                        obj.phi=obj.phi_new;
                        obj.r_x=obj.x;
                        obj.r_y=obj.y;
                        obj.r_phi=obj.phi;
                        
                        
                        J_new=jacobian_matrix_EKM(obj);
                        
                        A=J_new'*J_new;
                        
                        g=J_new'*f_new;
                        
                        if (norm(g,Inf)<=obj.epsilon_1)
                            found=true;
                            
                        else
                            found=false;
                        end
                        mue=mue*max([1/3 (1-(2*rho-1)^3)]);
                        v=2;
                        
                    else
                        
                        mue=mue*v;
                        v=2*v;
                    end
                end
            end
            
            X=[obj.x;obj.y;obj.phi];
            
            forwardkinematics_LM=X;
            
            iterations=k;
            
            
            
            
        end
        function [forwardkinematics_DL,iterations]=Dog_Leg_EKM(obj)
            k=0;
            
            %starting point
            obj.x=obj.x_0;
            obj.y=obj.y_0;
            obj.phi=obj.phi_0;
            
            obj.r_x=obj.x;
            obj.r_y=obj.y;
            obj.r_phi=obj.phi;
            
            J=jacobian_matrix_EKM(obj);
            f=vector_function_EKM(obj);
            g=J'*f;
            
            if((norm(f,Inf)<=obj.epsilon_3))|| ...
                    (norm(g,Inf)<=obj.epsilon_1)
                found=true;
            else
                found=false;
            end
            
            while ( (~found) && (k<obj.k_max))
                
                k=k+1;
                
                %J=jacobian_matrix_SKM(obj);
                %f=vector_function_SKM(obj);
                
                alpha=(norm(g)^2)/(norm(J*g)^2);
                
                h_sd=-g;
                
                M=J'*J;
                
                M_det=M(1,1)*M(2,2)*M(3,3)+M(1,2)*M(2,3)*M(3,1)+M(1,3)*M(2,1)*M(3,2)...
                    -M(3,1)*M(2,2)*M(1,3)-M(3,2)*M(2,3)*M(1,1)-M(3,3)*M(2,1)*M(1,2);
                
                M_inv=1/M_det*[M(2,2)*M(3,3)-M(2,3)*M(3,2) M(1,3)*M(3,2)-M(1,2)*M(3,3) M(1,2)*M(2,3)-M(1,3)*M(2,2);...
                    M(2,3)*M(3,1)-M(2,1)*M(3,3) M(1,1)*M(3,3)-M(1,3)*M(3,1) M(1,3)*M(2,1)-M(1,1)*M(2,3);...
                    M(2,1)*M(3,2)-M(2,2)*M(3,1) M(1,2)*M(3,1)-M(1,1)*M(3,2) M(1,1)*M(2,2)-M(1,2)*M(2,1);];
                
                h_gn=-M_inv*g;%h_gn=(J'*J)\(-J'*f);
                
                %trust region and Dog Leg step
                a=alpha*h_sd;
                b=h_gn;
                
                %F=objective_function_EKM(obj);
                F=1/2*(f'*f);
                
                %step
                if(norm(h_gn)<=obj.delta)
                    h_dl=h_gn;
                    LinModel=F;
                    
                elseif(norm(a)>=obj.delta)
                    h_dl=(obj.delta/norm(h_sd))*h_sd;
                    LinModel=obj.delta*(2*norm(alpha*g)-obj.delta)/(2*alpha);
                    
                else
                    c=a'*(b-a);
                    d2_a2=obj.delta^2-norm(a)^2;
                    b_a_2=norm(b-a)^2;
                    
                    if(c<=0)
                        beta=-c+sqrt(c^2+b_a_2*d2_a2)/b_a_2;
                    else
                        beta=d2_a2/(c+sqrt(c^2+b_a_2*d2_a2));
                    end
                    h_dl=a+beta*(b-a);
                    %beta chosen so that norm(h_dl)=delta
                    
                    LinModel=1/2*alpha*(1-beta)^2*norm(g)^2+beta*(2-beta)*F;
                end
                
                if(norm(h_dl)<=obj.epsilon_2*(norm([obj.x;obj.y;obj.phi])+obj.epsilon_2))
                    found=true;
                else
                    obj.x_new=obj.x+h_dl(1,1);
                    obj.y_new=obj.y+h_dl(2,1);
                    obj.phi_new=obj.phi+h_dl(3,1);
                    
                    obj.r_x=obj.x_new;
                    obj.r_y=obj.y_new;
                    obj.r_phi=obj.phi_new;
                    
                    f_new=vector_function_EKM(obj);
                    %F_new=objective_function_EKM(obj);
                    F_new=1/2*(f_new'*f_new);
                    
                    obj.r_x=obj.x;
                    obj.r_y=obj.y;
                    obj.r_phi=obj.phi;
                    
                    %gain ratio rho
                    rho=(F-F_new)/LinModel;
                    
                    if(rho>0)
                        obj.x=obj.x_new;
                        obj.y=obj.y_new;
                        obj.phi=obj.phi_new;
                        
                        obj.r_x=obj.x;
                        obj.r_y=obj.y;
                        obj.r_phi=obj.phi;
                        
                        J=jacobian_matrix_EKM(obj);
                        f=f_new;
                        
                        g=J'*f;
                        if((norm(f,Inf)<=obj.epsilon_3))|| ...
                                (norm(g,Inf)<=obj.epsilon_1)
                            found=true;
                        else
                            found=false;
                        end
                    end
                    
                    if(rho>0.75)
                        obj.delta=max(obj.delta,3*norm(h_dl));
                    elseif(rho<0.25)
                        obj.delta=obj.delta/2;
                        if(obj.delta<=obj.epsilon_2*...
                                (norm([obj.x;obj.y;obj.phi])+obj.epsilon_2))
                            found=true;
                        else
                            found=false;
                        end
                    end
                    
                end
                
            end
            forwardkinematics_DL=[obj.x;obj.y;obj.phi];
            iterations=k;
            
        end
        
        % Initial Estimation Algorithms for Forward Kinematics
        
        function forwardkinematics_estimatedpostition=box_method(obj)
            
            % euclidean norm of vector b_i in 2 dimensions
            b_1_norm=norm([obj.b_1(1);obj.b_1(2)]);
            b_2_norm=norm([obj.b_2(1);obj.b_2(2)]);
            b_3_norm=norm([obj.b_3(1);obj.b_3(2)]);
            b_4_norm=norm([obj.b_4(1);obj.b_4(2)]);
            
            % vector a_i in 2 dimensions
            a_1_2D=[obj.a_1(1);obj.a_1(2)];
            a_2_2D=[obj.a_2(1);obj.a_2(2)];
            a_3_2D=[obj.a_3(1);obj.a_3(2)];
            a_4_2D=[obj.a_4(1);obj.a_4(2)];
            
            % sphere's radius
            r_1_sphere=obj.L_SKM(1)+b_1_norm;
            r_2_sphere=obj.L_SKM(2)+b_2_norm;
            r_3_sphere=obj.L_SKM(3)+b_3_norm;
            r_4_sphere=obj.L_SKM(4)+b_4_norm;
            
            % determination r_i_low
            r_1_low=a_1_2D-r_1_sphere*[1;1];
            r_2_low=a_2_2D-r_2_sphere*[1;1];
            r_3_low=a_3_2D-r_3_sphere*[1;1];
            r_4_low=a_4_2D-r_4_sphere*[1;1];
            
            r_low_x_vector=[r_1_low(1);r_2_low(1);r_3_low(1);r_4_low(1)];
            r_low_y_vector=[r_1_low(2);r_2_low(2);r_3_low(2);r_4_low(2)];
            
            % determination r_i_high
            r_1_high=a_1_2D+r_1_sphere*[1;1];
            r_2_high=a_2_2D+r_2_sphere*[1;1];
            r_3_high=a_3_2D+r_3_sphere*[1;1];
            r_4_high=a_4_2D+r_4_sphere*[1;1];
            
            r_high_x_vector=[r_1_high(1);r_2_high(1);r_3_high(1);r_4_high(1)];
            r_high_y_vector=[r_1_high(2);r_2_high(2);r_3_high(2);r_4_high(2)];
            
            % greatest component of r_low_x_vector
            r_low_x=max(r_low_x_vector);
            
            % greatest component of r_low_y_vector
            r_low_y=max(r_low_y_vector);
            
            % smallest component of r_low_x_vector
            r_high_x=min(r_high_x_vector);
            
            % smallest component of r_low_y_vector
            r_high_y=min(r_high_y_vector);
            
            % definition r_low and r_high
            r_low=[r_low_x;r_low_y];
            r_high=[r_high_x;r_high_y];
            
            %center point of bounding box around estimated pose
            forwardkinematics_estimatedpostition=1/2*(r_low+r_high);
            
        end
        function midpoint_est=circle_method(obj)
            % circle method
            
            % euclidean norm (length) of vectors b_i
            b_1_length=norm(obj.b_1);
            b_2_length=norm(obj.b_2);
            b_3_length=norm(obj.b_3);
            b_4_length=norm(obj.b_4);
            
            % radius of estimation circles
            radius_1=obj.L_SKM(1)+b_1_length;
            radius_2=obj.L_SKM(2)+b_2_length;
            radius_3=obj.L_SKM(3)+b_3_length;
            radius_4=obj.L_SKM(4)+b_4_length;
            
            %intersection points winch 1 and winch 2
            
            % coordinate frame at A_1
            c12_x_s=((obj.a_2(1)-obj.a_1(1))^2+(radius_1)^2 ...
                -(radius_2)^2)/(2*(obj.a_2(1)-obj.a_1(1)));
            c12_y_1_s=sqrt((radius_1)^2-(c12_x_s^2));
            
            % transition to world coordinate frame
            c12_x=obj.a_1(1)+c12_x_s;
            c12_y_1=c12_y_1_s+obj.a_1(2);
            
            intersec_12=[c12_x;c12_y_1];
            
            %intersection points winch 1 and 3
            
            % coordinate frame at A_1
            c13_y_s=((obj.a_3(2)-obj.a_1(2))^2+(radius_1)^2 ...
                -(radius_3)^2)/(2*(obj.a_3(2)-obj.a_1(2)));
            c13_x_1_s=sqrt((radius_1)^2-(c13_y_s^2));
            
            % transition to world coordinate frame
            c13_y=obj.a_1(2)+c13_y_s;
            c13_x_1=obj.a_1(1)+c13_x_1_s;
            
            intersec_13=[c13_x_1;c13_y];
            
            %intersection points winch 3 and 4
            
            % coordinate frame at A_3
            c34_x_s=((obj.a_4(1)-obj.a_3(1))^2+(radius_3)^2 ...
                -(radius_4)^2)/(2*(obj.a_4(1)-obj.a_3(1)));
            c34_y_2_s=-sqrt((radius_3)^2-(c34_x_s^2));
            
            % transition to world coordinate frame
            c34_x=obj.a_3(1)+c34_x_s;
            c34_y_2=c34_y_2_s+obj.a_3(2);
            
            intersec_34=[c34_x;c34_y_2];
            
            %intersection points winch 2 and 4
            
            % coordinate frame at A_2
            c24_y_s=((obj.a_4(2)-obj.a_2(2))^2+(radius_2)^2 ...
                -(radius_4)^2)/(2*(obj.a_4(2)-obj.a_2(2)));
            c24_x_2_s=-sqrt((radius_2)^2-(c24_y_s^2));
            
            % transition to world coordinate frame
            c24_y=obj.a_2(2)+c24_y_s;
            c24_x_2=obj.a_2(1)+c24_x_2_s;
            
            intersec_24=[c24_x_2;c24_y];
            
            midpoint_x=1/2*(intersec_13(1)+intersec_24(1));
            midpoint_y=1/2*(intersec_12(2)+intersec_34(2));
            
            midpoint_est=[midpoint_x; midpoint_y];
            
        end
        
        
        
        
    end
    
    
end



