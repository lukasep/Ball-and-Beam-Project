% EE222: Nonlinear Systems
% Lab Project Phase I: Simulations
% Soomi Lee, Arvind Kruthiventy, Emily Lukas

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Our Controller Interface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef studentControllerInterface < matlab.System
    properties (Access = private)
        %% You can add values that you want to store and updae while running your controller.
        % For more information of the supported data type, see
        % https://www.mathworks.com/help/simulink/ug/data-types-supported-by-simulink.html
        t_prev = -1;
        theta_d = 0;
        V_servo = 0; 
        % our feedback linearization control function
        control_func;
        
        %% CHOOSE WHICH CONTROLLER TO USE:
        fl_lin = false; % true to use Controller 1, false to use Controller 2

        %% 
        % deriv gain
        L = [22.9338 , 1.0388 ; 
             130.7798 , 12.977; 
             0.9570 , 23.0662; 
             10.9830 , 132.2187]; 
        % L = ones(4, 2) * 100; 
     
        ol_est = [0; 0.00; 0; 0]; 
        ol_est_dot = [0 ; 0.00; 0; 0]; 
        
        % LQR params
        Q = zeros(4,4);
        R = 0.05;
        lie1_func;
        lie2_func;
        lie3_func;
        lie4_func;
        y_func;
        K;
        A = [0, 1, 0, 0; 0, 0, 1 , 0;0, 0, 0,1; 0, 0, 0, 0 ]; 
        B = [0; 0; 0; 1]; 
        C = [1, 0, 0, 0; 0, 0, 1, 0]; 
        Ad = zeros(4,4); 
        Bd = zeros(4,1);
        Cd = zeros(2, 4); 
        new_est = [-0.19, 0, 0, 0]
        v = 0 
        init_state = [0.19, 0, 0, 0]'; 
        
    end
    methods(Access = protected)

        function V_servo = stepImpl(obj, t, ball_pos, ball_vel, beam_ang, beam_ang_vel)
            % This is the main function called every iteration. You have to implement
            % the controller in this function, bu you are not allowed to
            % change the signature of this function.
            % Input arguments:
            %   t: current time
            %   ball_pos: position of the ball provided by the ball position sensor (m)
            %   ball_vel: velocity of the ball (m/s)
            %   beam_angle: servo motor angle provided by the encoder of the motor (rad)
            %   beam_angular_vel: angular velocity of the beam (rad/s)
            % Output:
            %   V_servo: voltage to the servo input
            % Safety Params
            beam_ang_min = -pi/4;
            beam_ang_max = pi/4;
            
            % Extract reference trajectory at the current timestep.
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);
            dt = t - obj.t_prev;

            % Observer 
            dt = t - obj.t_prev;

            % Choose which controller to use
            if obj.fl_lin == true
                error_x1 = obj.y_func(ball_pos, ball_vel, beam_ang, beam_ang_vel) - p_ball_ref;
                error_x2 = obj.lie1_func(ball_pos, ball_vel, beam_ang, beam_ang_vel) - v_ball_ref;
                error_x3 = obj.lie2_func(ball_pos, ball_vel, beam_ang, beam_ang_vel) - a_ball_ref;
                error_x4 = obj.lie3_func(ball_pos, ball_vel, beam_ang, beam_ang_vel) - 0;
                error = [error_x1; error_x2; error_x3; error_x4];
             
                % Apply the LQR gain matrix to calculate the virtual control input
                % The negative sign ensures negative feedback (reducing error)
                % obj.K is the optimal gain matrix calculated by LQR
                % v is then the virtua input for the linearized system
                K = obj.Feedback_LQR(dt); 
                v = K * error;
                V_servo = obj.control_func(ball_pos, ball_vel, beam_ang, beam_ang_vel, v);
                
            % disp(state_vec)
            else
                % LQR for system linearized around trajectory
                
                A = zeros(4,4);
                A(1:3,2:4) = eye(3);
                B = [0; 0; 0; 1];
                
                state = [ball_pos; ball_vel; beam_ang; beam_ang_vel];
                ref = [p_ball_ref; v_ball_ref; 0; 0]; 
                error_x1 = ball_pos - p_ball_ref; 
                error_x2 = ball_vel - v_ball_ref; 
                error_x3 = beam_ang ; 
                error_x4 = beam_ang_vel; 
                error = [error_x1; error_x2; error_x3; error_x4];
                K = obj.lqr_linear(obj.init_state, ref,dt);
                v = -K * error;
            
                V_servo =v;
            end

            % Observer called here, runs after controller because we need
            % feedback linearizations virtual control for estimation
            obj.luenberger_obs(ball_pos, beam_ang, v, dt);  
            obj.init_state = obj.ol_est;
            % ball_pos = obj.ol_est(1); 
            % ball_vel = obj.ol_est(2); 
            % beam_ang = obj.ol_est(3); 
            % beam_ang_vel = obj.ol_est(4);
            
            % V_servo = v; 
            % Restict to the safe region as in the original script (I
            % tightened this a bit)
            if beam_ang > beam_ang_max
                V_servo = min(V_servo, 10 * (beam_ang_max - beam_ang));
            elseif beam_ang < beam_ang_min
                V_servo = max(V_servo, 10 * (beam_ang_min - beam_ang));
            end

            obj.t_prev = t;
        end
    end

    methods(Access = public)
        function obj = studentControllerInterface
            % Initialize controller by computing the feedback
            % disp(eig(obj.Ad - obj.L * obj.Cd))
            % disp(obj.Bd)
            obj = setupFeedbackLinearization(obj);
            %obj.fl_lin = false; 
            if obj.fl_lin == false
                fprintf('Use Controller 2\n')
            else
                fprintf('Use Controller 1\n')
            end

        end
        
        function [V_servo, theta_d] = stepController(obj, t, ball_pos, ball_vel, beam_angle, beam_angular_vel)
            V_servo = stepImpl(obj, t, ball_pos, ball_vel, beam_angle, beam_angular_vel);
            theta_d = obj.theta_d;
            obj.V_servo = V_servo; 
        end
    end


    methods(Access = private)
        
        function obj = setupFeedbackLinearization(obj)
            %% Feedback Linearization:
            % Note that we had to change this from symbolic variables to
            % the way they are inputted due to simulink's inability to run
            % with symbolic vars. 
            %% Given parameters from the PDF
            ball_rad = 0.0254;   % Ball radius [m]
            beam_len = 0.4255;   % Beam length [m]
            g_val = 9.81;         % Gravitational acceleration [m/s^2]
            servo_gain = 1.5;       % Servo motor gain [rad/(V*s)]
            tau_val = 0.025;  % Motor time constant [s]
            
            % Constants for simplicity
            const_1 = 5*g_val*ball_rad/(7*beam_len);      % const_1 = 5*g*rg/(7*L)
            const_2 = (5/7)*(ball_rad/beam_len)^2;        % const_2 = 5/7*(rg/L)^2
            
            %% Computed Lie functions
            obj.lie1_func = @(x1,x2,x3,x4) x2;
            obj.lie2_func = @(x1,x2,x3,x4) const_1*sin(x3) - const_2*((beam_len/2) - x1)*x4.^2.*cos(x3).^2;
            obj.lie3_func = @(x1,x2,x3,x4) const_2*x2*x4.^2.*cos(x3).^2 + const_1*x4.*cos(x3) + ...
                2*const_2*((beam_len/2)-x1)*x4.^3.*cos(x3).*sin(x3) + (2*const_2/tau_val)*((beam_len/2)-x1)*x4.^2.*cos(x3).^2;
            
            % lie4 = phi(x) + psi(x) * u
            % phi(x)
            phi_func = @(x1,x2,x3,x4) ( (-2*const_2*x4.^3.*cos(x3).*sin(x3) - (2*const_2/tau_val)*x4.^2.*cos(x3).^2).*x2 ...
                  + (const_2*x4.^2.*cos(x3).^2).*(const_1*sin(x3) - const_2*((beam_len/2)-x1)*x4.^2.*cos(x3).^2) ...
                  + (-2*const_2*x2.*x4.^2.*cos(x3).*sin(x3) - const_1*x4.*sin(x3) ...
                     + 2*const_2*((beam_len/2)-x1)*x4.^3.*(cos(x3).^2 - sin(x3).^2) ...
                     - (4*const_2/tau_val)*((beam_len/2)-x1)*x4.^2.*cos(x3).*sin(x3)).*x4 ...
                  + (2*const_2*x2.*x4.*cos(x3).^2 + const_1*cos(x3) ...
                     + 6*const_2*((beam_len/2)-x1)*x4.^2.*cos(x3).*sin(x3) ...
                     + (4*const_2/tau_val)*((beam_len/2)-x1)*x4.*cos(x3).^2) .* (-x4/tau_val) );
            
            % psi(x)
            psi_func = @(x1,x2,x3,x4) ( 2*const_2*x2.*x4.*cos(x3).^2 + const_1*cos(x3) ...
                     + 6*const_2*((beam_len/2)-x1)*x4.^2.*cos(x3).*sin(x3) ...
                     + (4*const_2/tau_val)*((beam_len/2)-x1)*x4.*cos(x3).^2 )*(servo_gain/tau_val);

            obj.lie4_func = @(x1,x2,x3,x4,u) phi_func(x1,x2,x3,x4) + psi_func(x1,x2,x3,x4)*u;
            
            %% Feedback linearization control law u = (v - phi(x)) / psi(x).
            obj.control_func = @(x1,x2,x3,x4,v) (v - phi_func(x1,x2,x3,x4)) / psi_func(x1,x2,x3,x4);
            
            %% Output function: y = x1.
            obj.y_func = @(x1,x2,x3,x4) x1;
        end
        
        function [K] = Feedback_LQR(obj,dt)
            % LQR setup
            
            % Define system matrices
            A = zeros(4,4);
            A(1:3,2:4) = eye(3);
            B = [0; 0; 0; 1]; 
            
            obj.Q(1,1) = 550; % made this 550 for method 1
            obj.Q(2,2) = 50;  % made this 50 for method 1
            obj.R = 1;        % made this 1 for method 1
            
            [Ad, Bd] = obj.cont2discrete(A,B, dt); 
            % Compute LQR gain
            % obj.K = lqr(A, B, obj.Q, obj.R);
            K = obj.ARESolve(Ad, Bd, obj.Q, obj.R); % compatible with simulink
        end

        function K = lqr_linear(obj, state, ref, dt)
            % State equations
            %% 2. Define system dynamics: dx/dt = f(x) + g(x)*u
            % Nonlinear dynamics of the ball and beam system
            rg = 0.0254;
            L = 0.4255;
            g = 9.81;
            tau = 0.025;
            k = 1.5; 
            
            x1 = state(1); 
            x2 = state(2); 
            x3 = state(3); 
            x4 = state(4); 
            
            jacb = [0, 1, 0 , 0;
                    ((5/7)*x1)*((rg/L)^2)*(x4^2)*cos(x3)^2, 0, (5*g*rg/(7*L))*cos(x3)+(5/7) * ((L/2) - x1) * ((rg/L)^2) * x4^2 * sin(2*x3), (10/7) * ((L/2) - x1) * ((rg/L)^2)*x4*cos(x3)^2 ; 
                    0, 0, 0, 1; 
                    0, 0, 0, -1/tau]; 
            B = [0; 0; 0; -k/tau];
            

            % K = obj.ARESolve(jacb, B, Q, R); 
            obj.Q(1,1) = 500;  % Weight on position error
            obj.Q(2, 2) = 100;
            obj.R = 1; 
            obj.Q(3,3) = 0; 
            obj.Q(4,4) = 0; 
            % K = lqr(jacb, B, obj.Q, obj.R);
            [Ad, Bd] = obj.cont2discrete(jacb,B, dt); 
            K = obj.ARESolve(Ad, Bd, obj.Q, obj.R);
              
        end

        function obj = luenberger_obs(obj, ball_pos, beam_ang,v, dt)
            % The setup for the Luenberger observer
            y_t = [ball_pos; beam_ang]; 
            obj.ol_est_dot = obj.A* obj.ol_est + obj.B * v + obj.L *(y_t - (obj.C * obj.ol_est));
            new_est = obj.ol_est_dot; 
            obj.ol_est = obj.ol_est + dt * new_est; 
        end  

        function K = ARESolve(obj, A, B, Q, R)
            Pk = Q; 
            for i = 1:400
                Pk = A' * Pk * A + Q - A' * Pk * B * inv(R + B' * Pk * B) * B' * Pk * A; 
                % disp(Pk)
            end
           
            K = -inv(R + B' * Pk * B)*B'*Pk*A; 
        end 
        
        function [Ad, Bd] = cont2discrete(obj, A, B, dt)
            Ad = expm(A * dt); 
            B_func = @(t) expm(A.*t) * B; 
            Bd = integral(B_func, 0, dt, ArrayValued=true);
        end 
        
    end
end