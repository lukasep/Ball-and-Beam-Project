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
        % our feedback linearization control function
        control_func;

        % tune gains (this was from try 1)
        % kp = 0.1;  % prop gain
        % kd = 0.6;    % deriv gain

        % LQR params
        Q = zeros(4,4);
        R = 0.001;
        lie1_func;
        lie2_func;
        lie3_func;
        lie4_func;
        y_func;
        K;

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
            %   V_servo: voltage to the servo input.

            %% Feedback Linearization Controller with LQR (try 1 and 2)
            % Safety Params
            beam_ang_min = -pi/4;
            beam_ang_max = pi/4;

            % Extract reference trajectory at the current timestep.
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);

            % Calc the errors
            error_x1 = obj.y_func(ball_pos, ball_vel, beam_ang, beam_ang_vel) - p_ball_ref;
            error_x2 = obj.lie1_func(ball_pos, ball_vel, beam_ang, beam_ang_vel) - v_ball_ref;
            error_x3 = obj.lie2_func(ball_pos, ball_vel, beam_ang, beam_ang_vel) - a_ball_ref;
            error_x4 = obj.lie3_func(ball_pos, ball_vel, beam_ang, beam_ang_vel) - 0;
            error = [error_x1; error_x2; error_x3; error_x4];

            % Apply the LQR gain matrix to calculate the virtual control input
            % The negative sign ensures negative feedback (reducing error)
            % obj.K is the optimal gain matrix calculated by LQR
            % v is then the virtua input for the linearized system
            v = - obj.K * error;

            % Apply the feedback linearization func
            V_servo = obj.control_func(ball_pos, ball_vel, beam_ang, beam_ang_vel, v);

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
            obj = setupFeedbackLinearization(obj);

            % Set up LQR: try 2
            obj = setupLQR(obj);
        end

        function [V_servo, theta_d] = stepController(obj, t, ball_pos, ball_vel, beam_angle, beam_angular_vel)
            V_servo = stepImpl(obj, t, ball_pos, ball_vel, beam_angle, beam_angular_vel);
            theta_d = obj.theta_d;
        end
    end

    methods(Access = private)
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
            obj.Q(1,1) = 1;  % Weight on position error

            % Define system matrices
            A = zeros(4,4);
            A(1:3,2:4) = eye(3);
            B = [0; 0; 0; 1];

            % Compute LQR gain
            obj.K = lqr(A, B, obj.Q, obj.R);
        end
    end
end
