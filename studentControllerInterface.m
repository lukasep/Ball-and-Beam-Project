% EE222: Nonlinear Systems
% Lab Project Phase I: Simulations
% Soomi Lee, Arvind Kruthiventy, Emily Lukas

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Our Controller Interface:
% TRY 1: Just feedback linearization. Didn't work super well alone. I tuned
% the values of kp and kd a lot but they just weren't great.
% TRY 2: Added LQR. Does work with score of 0.91!
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
        
        % tune gains (this was from try 1)
        % kp = 0.1;  % prop gain
        % kd = 
        % 
        % 0.6;  
        
        
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
            %observer 
 
            

            %LQR for system linearized around trajectory
            dt = t - obj.t_prev;
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
            %Observer called here, runs after controller because we need
            %feedback linearizations virtual control for estimation
            obj.luoberger_obs(ball_pos, beam_ang, v, dt);  
            obj.init_state = obj.ol_est;
            ball_pos = obj.ol_est(1); 
            ball_vel = obj.ol_est(2); 
            beam_ang = obj.ol_est(3); 
            beam_ang_vel = obj.ol_est(4);
            
            error_x1 = obj.y_func(ball_pos, ball_vel, beam_ang, beam_ang_vel) - p_ball_ref;
            error_x2 = obj.lie1_func(ball_pos, ball_vel, beam_ang, beam_ang_vel) - v_ball_ref;
            error_x3 = obj.lie2_func(ball_pos, ball_vel, beam_ang, beam_ang_vel) - a_ball_ref;
            error_x4 = obj.lie3_func(ball_pos, ball_vel, beam_ang, beam_ang_vel) - 0;
            error = [error_x1; error_x2; error_x3; error_x4];
            
            
            % Apply the LQR gain matrix to calculate the virtual control input
            % The negative sign ensures negative feedback (reducing error)
            % obj.K is the optimal gain matrix calculated by LQR
            % v is then the virtua input for the linearized system
            
            % v = - obj.K * error;
            % disp(state_vec)
            % Apply the feedback linearization func
            
            % V_servo = obj.control_func(ball_pos, ball_vel, beam_ang, beam_ang_vel, v);
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
            % Feedback Linearization: try 1
            % disp(eig(obj.Ad - obj.L * obj.Cd))
            % disp(obj.Bd)
            obj = setupFeedbackLinearization(obj);
            
            % Set up LQR: try 2
            obj = setupLQR(obj);
            
        end
        
        function [V_servo, theta_d] = stepController(obj, t, ball_pos, ball_vel, beam_angle, beam_angular_vel)
            V_servo = stepImpl(obj, t, ball_pos, ball_vel, beam_angle, beam_angular_vel);
            theta_d = obj.theta_d;
            obj.V_servo = V_servo; 
        end
    end

    methods(Access = private)
        % function error = get_fl_error(state, ref)
        % 
        % end 
        function obj = setupFeedbackLinearization(obj)
            %% Feedback Linearization: TRY #1
            % 1. Define symbolic variables
            % State variables based on the ball_and_beam.pdf:
            %   ball_pos (x1) - Position of ball on beam [m]
            %   ball_vel (x2) - Velocity of ball [m/s]
            %   beam_angle (x3) - Angle of beam [rad]
            %   beam_angular_vel (x4) - Angular velocity of beam [rad/s]

            % Input and parameters:
            %   u - Control input (voltage to servo) [V]
            %   g - Gravitational acceleration [m/s^2]
            %   rg - Radius of ball [m]
            %   L - Length of beam [m]
            %   K - Servo motor gain [rad/(V*s)]
            %   tau - Servo motor time constant [s]
            %   v - Virtual input after linearization

            syms x1 x2 x3 x4 u g rg L K tau v real

            %% 2. Define system dynamics: dx/dt = f(x) + g(x)*u
            % Nonlinear dynamics of the ball and beam system
            f = [x2;  % dx1/dt = x2 (ball velocity)
                (5*g*rg/(7*L)) * sin(x3) - (5/7) * ((L/2) - x1) * ((rg/L)^2) * x4^2 * cos(x3)^2;  % dx2/dt (ball acceleration)
                x4;  % dx3/dt = x4 (beam angular velocity)
                -x4/tau];  % dx4/dt (beam angular acceleration (without control for now))

            % Control input vector (only affects x4 aka beam angular acceleration)
            u_vec = [0; 0; 0; K/tau];

            % Combined dynamics: dx/dt = f(x) + g(x)*u
            state_vector = [x1; x2; x3; x4];

            %% 3. Feedback Linearization Procedure
            % Define output (we want to control the ball position on the baem)
            y = x1;

            % Step 1: Compute Lie derivs
            lie1 = simplify(jacobian(y, state_vector) * (f + u_vec*u));
            lie2 = simplify(jacobian(lie1, state_vector) * (f + u_vec*u));
            lie3 = simplify(jacobian(lie2, state_vector) * (f + u_vec*u));

            % Step 2: Get rid of u off of the 3rd and 4th derivatives
            lie3 = simplify(expand(lie3) - (-(5*K*rg^2*u*x4*cos(x3)^2)/(7*L*tau) + (10*K*rg^2*u*x1*x4*cos(x3)^2)/(7*L^2*tau)));
            % I had to manually su

            % get lie4 without u
            lie4 = jacobian(lie3, state_vector) * (f + u_vec*u);
            lie4_exp = expand(lie4);

            % Step 3: get control input u to achieve v = lie4
            control_law = simplify(solve(lie4 == v, u));

            %% 4. Given parameters from the PDF
            ball_rad = 0.0254;   % Ball radius [m]
            beam_len = 0.4255;   % Beam length [m]
            g_val = 9.81;         % Gravitational acceleration [m/s^2]
            servo_gain = 1.5;       % Servo motor gain [rad/(V*s)]
            tau_val = 0.025;  % Motor time constant [s]

            % Subs params in
            u_real = subs(control_law, {g, rg, L, K, tau}, {g_val, ball_rad, beam_len, servo_gain, tau_val});
            f = subs(f, {g, rg, L, K, tau}, {g_val, ball_rad, beam_len, servo_gain, tau_val});
            lie1 = subs(lie1, {g, rg, L, K, tau}, {g_val, ball_rad, beam_len, servo_gain, tau_val});
            lie2 = subs(lie2, {g, rg, L, K, tau}, {g_val, ball_rad, beam_len, servo_gain, tau_val});
            lie3 = subs(lie3, {g, rg, L, K, tau}, {g_val, ball_rad, beam_len, servo_gain, tau_val});

            obj.control_func = matlabFunction(u_real, 'Vars', [state_vector' v]);
            obj.y_func = matlabFunction(y, 'Vars', state_vector');
            obj.lie1_func = matlabFunction(lie1, 'Vars', state_vector');
            obj.lie2_func = matlabFunction(lie2, 'Vars', state_vector');
            obj.lie3_func = matlabFunction(lie3, 'Vars', state_vector');
        end

        function obj = setupLQR(obj)
            % LQR setup
            
            % Define system matrices
            A = zeros(4,4);
            A(1:3,2:4) = eye(3);
            B = [0; 0; 0; 1];

            % Compute LQR gain
            % obj.K = lqr(A, B, obj.Q, obj.R);
            % obj.K = obj.ARESolve(A, B, obj.Q, obj.R); 
        end
        function K = lqr_linear(obj, state, ref, dt)
            %state equations

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
            Q = zeros(4,4);
            Q(1, 1) = 5;
            R = 0.001; 
            sys = ss(jacb,B, zeros(2, 4), zeros(2, 1)); 
            sysd = c2d(sys, 0.01); 
            Ad = sysd.A; 
            Bd = sysd.B; 

            % K = obj.ARESolve(jacb, B, Q, R); 
            obj.Q(1,1) = 300;  % Weight on position error
            obj.Q(2, 2) = 100;
            obj.R = 10; 
            obj.Q(3,3) = 0; 
            obj.Q(4,4) = 0; 
            % K = lqr(jacb, B, obj.Q, obj.R);
            K = obj.ARESolve(Ad, Bd, obj.Q, obj.R);
              
        end
        function obj = luoberger_obs(obj, ball_pos, beam_ang,v, dt)
 
            y_t = [ball_pos; beam_ang]; 
            obj.ol_est_dot = obj.A* obj.ol_est + obj.B * v + obj.L *(y_t - (obj.C * obj.ol_est));
            new_est = obj.ol_est_dot; 
            obj.ol_est = obj.ol_est + dt * new_est; 
        end  
        function K = ARESolve(obj, A, B, Q, R)
            Pk = Q; 
            for i = 1:100
                Pk = A' * Pk * A + Q - A' * Pk * B * inv(R + B' * Pk * B) * B' * Pk * A; 
                % disp(Pk)
            end
           
            K = -inv(R + B' * Pk * B)*B'*Pk*A; 
        end 
        function obj = ekf(obj, ball)
        end
    end
end
