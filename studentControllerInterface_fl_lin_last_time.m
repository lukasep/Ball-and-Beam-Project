% EE222: Nonlinear Systems
% Lab Project Phase I: Simulations
% Soomi Lee, Arvind Kruthiventy, Emily Lukas
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Our Controller Interface:
% % Notes: 4/8 only the first method works well (FL_LIN = TRUE)
% % We could probably make the second one better but I was too lazy to find
% % the proper gains. 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
classdef studentControllerInterface < matlab.System
    properties (Access = private)
        % Timing variables
        t_prev = -1;
        % Control signals
        theta_d = 0;
        V_servo = 0;
        % Controller selection flag
        fl_lin = true;  % Using Feedback Lin
        % fl_lin = false; % Using Trajectory LQR
        % System parameters
        ball_rad = 0.0254;   % Ball radius [m]
        beam_len = 0.4255;   % Beam length [m]
        g_val = 9.81;        % Gravitational acceleration [m/s^2]
        servo_gain = 1.5;    % Servo motor gain [rad/(V*s)]
        tau_val = 0.025;     % Motor time constant [s]
        const_1 = 0;
        const_2 = 0;
        % State estimate
        state_estimate = [0; 0; 0; 0];
        % LQR parameters
        Q = zeros(4,4);
        R = 1;
        K = [];
        % System matrices
        A = [0, 1, 0, 0;
            0, 0, 1, 0;
            0, 0, 0, 1;
            0, 0, 0, 0];
        B = [0; 0; 0; 1];
        C = [1, 0, 0, 0; 0, 0, 1, 0];
        % Observer gain matrix
%         observer_gain =[22.9338 , 1.0388 ;
%             130.7798 , 12.977;
%             0.9570 , 23.0662;
%             10.9830 , 132.2187];
%         observer_gain = [    2.5476,    0.2454;
%     1.6221  ,  1.3068;
%     0.2454  ,  2.4524;
%     0.3068   , 1.5029]; 
        observer_gain = [9.8792    2.9513;
           16.1168    8.9874;
            2.9152    8.6208;
            7.8843   12.6327;]
    end
    methods(Access = protected)
        function [V_servo, estimated_ball_pos, estimated_ball_vel, estimated_beam_ang, estimated_beam_ang_vel] = stepImpl(obj, t, ball_pos, beam_ang)
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
            coder.extrinsic("lqr")
            beam_ang_min = -pi/4;
            beam_ang_max = pi/4;
            % Calculate time step
            dt = t - obj.t_prev;
            if obj.t_prev < 0
                dt = 0.001;
            end
            % Get reference trajectory
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);
            % Measurements available from sensors
            measurements = [ball_pos; beam_ang];
            % Get the previous control input
            control_input = obj.V_servo;
            if isempty(control_input)
                control_input = 0;
            end
            % Update state estimate using the Luenberger observer
            obj.state_estimate = obj.runObserver(dt, obj.state_estimate, measurements, control_input);
            % Extract state estimates
            estimated_ball_pos = obj.state_estimate(1);
            estimated_ball_vel = obj.state_estimate(2);
            estimated_beam_ang = obj.state_estimate(3);
            estimated_beam_ang_vel = obj.state_estimate(4);
            % Use estimated states for control
            if obj.fl_lin == true
                % Feedback linearization controller
                % Compute the error state
                error_x1 = estimated_ball_pos - p_ball_ref;
                error_x2 = estimated_ball_vel - v_ball_ref;
                error_x3 = obj.lie2(estimated_ball_pos, estimated_ball_vel, estimated_beam_ang, estimated_beam_ang_vel) - a_ball_ref;
                error_x4 = obj.lie3(estimated_ball_pos, estimated_ball_vel, estimated_beam_ang, estimated_beam_ang_vel) - 0;
                error = [error_x1; error_x2; error_x3; error_x4];
                % This K is using the hard coded matrices from the function
                % set below
       
                disp([error_x3, error_x4])
                v = -obj.K * error;
                % Apply feedback linearization to get the actual control
%                 V_servo = obj.computeControl(estimated_ball_pos, estimated_ball_vel, estimated_beam_ang, estimated_beam_ang_vel, v);
                V_servo = obj.computeControl(estimated_ball_pos, estimated_ball_vel, estimated_beam_ang, estimated_beam_ang_vel, v);
            else % this doesn't work rn
                % LQR for system linearized around trajectory
                % Calculate error state
                error_x1 = estimated_ball_pos - p_ball_ref;
                error_x2 = estimated_ball_vel - v_ball_ref;
                error_x3 = estimated_beam_ang;
                error_x4 = estimated_beam_ang_vel;
                error = [error_x1; error_x2; error_x3; error_x4];
                % Calculate Jacobian at current state
                x1 = estimated_ball_pos;
                x2 = estimated_ball_vel;
                x3 = estimated_beam_ang;
                x4 = estimated_beam_ang_vel;
                L = obj.beam_len;
                rg = obj.ball_rad;
                jacb = [0, 1, 0, 0;
                    ((5/7)*x1)*((rg/L)^2)*(x4^2)*cos(x3)^2, 0, (5*obj.g_val*rg/(7*L))*cos(x3)+(5/7)*((L/2)-x1)*((rg/L)^2)*x4^2*sin(2*x3), (10/7)*((L/2)-x1)*((rg/L)^2)*x4*cos(x3)^2;
                    0, 0, 0, 1;
                    0, 0, 0, -1/obj.tau_val];
                B_control = [0; 0; 0; -obj.servo_gain/obj.tau_val];
                % Set LQR weights
                obj.Q(1,1) = 100;  % Position error weight
                obj.Q(2,2) = 10;  % Velocity error weight
                obj.Q(3,3) = 1;    % Angle weight
                obj.Q(4,4) = 1;    % Angular velocity weight
                obj.R = 2;         % Control effort weight
                % Calculate LQR gain
                % K = lqr(jacb, B_control, obj.Q, obj.R)
                K = [20, 25.1525, 13.0233, 2.6315];
                % Apply LQR to get control
                V_servo = -K * error;
            end
            % Apply safety limits to servo voltage
            if beam_ang > beam_ang_max
                V_servo = min(V_servo, 10 * (beam_ang_max - beam_ang));
            elseif beam_ang < beam_ang_min
                V_servo = max(V_servo, 10 * (beam_ang_min - beam_ang));
            end
            obj.V_servo = V_servo;
            obj.t_prev = t;
        end
    end
    methods(Access = public)
        function obj = studentControllerInterface
            % Constructor to initialize the controller
            obj.const_1 = 5*obj.g_val*obj.ball_rad/(7*obj.beam_len);
            obj.const_2 = (5/7)*(obj.ball_rad/obj.beam_len)^2;
            % % Set LQR weights and compute gains
            % obj.Q(1,1) = 10;  % Position error weight
            % obj.Q(2,2) = 70;   % Velocity error weight
            % obj.R = 0.1;         % Control effort weight
            % Calculate LQR gain for the Feedback Linearization method
            %obj.K = lqr(obj.A, obj.B, obj.Q, obj.R)
            obj.K = [23.4521   29.3025   17.2402    5.8720]; % MATLAB SINE WAVE-- Score: 0.84
%             obj.K = [100.0000   89.3115   36.3827    8.5303]; % MATLAB SQUARE WAVE-- Score: 3.33
%             obj.K = [120   89.3115   20    5]; % Better
%             obj.K = [120   80   20   2]; % SHITTY
%             obj.K = [120   90   0   0]; % Shittier
%             obj. K = [3.0002   6.0101    8.0129   4.0423]; % SIMULINK SINE WAVE-- Score: 1.5
            %obj.K = [3.0102   6.1143    7.9920   4.0120]; % SIMULINK SQUARE WAVE-- 4.37
            
            
            % Display which controller is active
            if obj.fl_lin
                disp('Using Feedback Linearization Controller');
            else
                disp('Using LQR Controller');
            end
        end
        function [V_servo, theta_d] = stepController(obj, t, ball_pos, beam_angle)
            % Interface method for the controller - matches the expected call signature
            % from run_matlab_ball_and_beam
            V_servo = stepImpl(obj, t, ball_pos, beam_angle);
            theta_d = obj.theta_d;
        end
    end
    methods(Access = private)
        function result = lie1(obj, x1, x2, x3, x4)
            % First Lie derivative
            result = x2;
        end
        function result = lie2(obj, x1, x2, x3, x4)
            % Second Lie derivative
            result = obj.const_1 * sin(x3) - obj.const_2 * ((obj.beam_len/2) - x1) * x4^2 * cos(x3)^2;
        end
        function result = lie3(obj, x1, x2, x3, x4)
            % Third Lie derivative
            result = obj.const_2 * x2 * x4^2 * cos(x3)^2 + ...
                obj.const_1 * x4 * cos(x3) + ...
                2 * obj.const_2 * ((obj.beam_len/2) - x1) * x4^3 * cos(x3) * sin(x3) + ...
                (2 * obj.const_2 / obj.tau_val) * ((obj.beam_len/2) - x1) * x4^2 * cos(x3)^2;
        end
        function result = phi(obj, x1, x2, x3, x4)
            % phi(x) part of the control law
            result = (-2 * obj.const_2 * x4^3 * cos(x3) * sin(x3) - ...
                (2 * obj.const_2 / obj.tau_val) * x4^2 * cos(x3)^2) * x2 + ...
                (obj.const_2 * x4^2 * cos(x3)^2) * (obj.const_1 * sin(x3) - ...
                obj.const_2 * ((obj.beam_len/2) - x1) * x4^2 * cos(x3)^2) + ...
                (-2 * obj.const_2 * x2 * x4^2 * cos(x3) * sin(x3) - obj.const_1 * x4 * sin(x3) + ...
                2 * obj.const_2 * ((obj.beam_len/2) - x1) * x4^3 * (cos(x3)^2 - sin(x3)^2) - ...
                (4 * obj.const_2 / obj.tau_val) * ((obj.beam_len/2) - x1) * x4^2 * cos(x3) * sin(x3)) * x4 + ...
                (2 * obj.const_2 * x2 * x4 * cos(x3)^2 + obj.const_1 * cos(x3) + ...
                6 * obj.const_2 * ((obj.beam_len/2) - x1) * x4^2 * cos(x3) * sin(x3) + ...
                (4 * obj.const_2 / obj.tau_val) * ((obj.beam_len/2) - x1) * x4 * cos(x3)^2) * (-x4 / obj.tau_val);
        end
        function result = psi(obj, x1, x2, x3, x4)
            % psi(x) part of the control law
            result = (2 * obj.const_2 * x2 * x4 * cos(x3)^2 + obj.const_1 * cos(x3) + ...
                6 * obj.const_2 * ((obj.beam_len/2) - x1) * x4^2 * cos(x3) * sin(x3) + ...
                (4 * obj.const_2 / obj.tau_val) * ((obj.beam_len/2) - x1) * x4 * cos(x3)^2) * ...
                (obj.servo_gain / obj.tau_val);
        end
        function result = computeControl(obj, x1, x2, x3, x4, v)
            % Feedback linearization control law: u = (v - phi(x))/psi(x)
            result = (v - obj.phi(x1, x2, x3, x4)) / obj.psi(x1, x2, x3, x4);
        end
        function state_estimate = runObserver(obj, time_step, previous_estimate, measurements, control_input)
            % Luenberger Observer implementation
            % Extract state variables for clarity
%             A = [0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 0, 0]; 
%             C = [1, 0, 0, 0; 0, 0, 1, 0]; 
%             desired_poles = [-10, -5, -2.25, -1.25]
%             L = place(A', C', desired_poles)' 
            ball_position = previous_estimate(1);
            ball_velocity = previous_estimate(2);
            beam_angle = previous_estimate(3);
            beam_angular_velocity = previous_estimate(4);
            % Calculate the nonlinear system dynamics
            state_derivative = zeros(4, 1);
            state_derivative(1) = ball_velocity;
            state_derivative(2) = obj.const_1 * sin(beam_angle) - ...
                (obj.const_2 * ((obj.beam_len/2) - ball_position) * beam_angular_velocity^2 * cos(beam_angle)^2);
            state_derivative(3) = beam_angular_velocity;
            state_derivative(4) = (-beam_angular_velocity + obj.servo_gain * control_input) / obj.tau_val;
            % Add correction term based on measurement error
            measurement_error = measurements - obj.C * previous_estimate;
            state_derivative = state_derivative + obj.observer_gain * measurement_error;
            % Update the state estimate using forward Euler integration
            state_estimate = previous_estimate + state_derivative * time_step;
        end
    end
end
% classdef studentControllerInterface < matlab.System
%     properties (Access = private)
%         % EKF state and covariance
%         state_estimate = zeros(4,1);   % [p; v; theta; theta_dot]
%         P              = eye(4);       % Initial covariance
%         
%         % EKF process and measurement noise
%         Q = diag([100, 10, 1, 1]);
%         R = diag([2, 2]);
%         
%         % Physical parameters
%         ball_rad   = 0.0254;
%         beam_len   = 0.4255;
%         g_val      = 9.81;
%         servo_gain = 1.5;
%         tau_val    = 0.025;
%         
%         % Output map (we measure ball_pos and beam_ang)
%         C = [1, 0, 0, 0;
%              0, 0, 1, 0];
%         
%         % Controller gains
%         fl_lin = true;                  % Toggle FL vs. LQR
%         K_fl   = [120   89.3115   20    5];
% %         K_fl = [3.0002   6.0101    8.0129   4.0423];
%         K_lqr  = [20, 25.1525, 13.0233, 2.6315]; 
%         
%         % Bookkeeping variables
%         t_prev  = -1;
%         V_servo = 0;
%         theta_d = 0;
%         
%         % Lie constants for feedback linearization
%         const_1 = 0;
%         const_2 = 0;
%     end
%     
%     methods
%         function obj = studentControllerInterface
%             % Initialize constants for Lie derivatives
%             obj.const_1 = 5 * obj.g_val * obj.ball_rad / (7 * obj.beam_len);
%             obj.const_2 = (5/7) * (obj.ball_rad / obj.beam_len)^2;
%             
%             if obj.fl_lin
%                 disp('Using Feedback Linearization + LQR-v');
%             else
%                 disp('Using pure LQR around origin');
%             end
%         end
%     end
%     
%     methods (Access = protected)
%         function [V_servo, estimated_ball_pos, estimated_ball_vel, estimated_beam_ang, estimated_beam_ang_vel] = stepImpl(obj, t, ball_pos, beam_ang)
%             % Time step calculation
%             dt = t - obj.t_prev;
%             if obj.t_prev < 0, dt = 1e-3; end
%             
%             % EKF prediction and update
%             z = [ball_pos; beam_ang];
%             u = obj.V_servo;
%             [obj.state_estimate, obj.P] = ...
%                 obj.EKF(dt, obj.state_estimate, obj.P, z, u);
%             
%             % Control computation
%             [p_ref, v_ref, a_ref] = get_ref_traj(t);
%             xhat = obj.state_estimate;
%             err  = [ xhat(1)-p_ref;
%                      xhat(2)-v_ref;
%                      xhat(3)-a_ref;
%                      xhat(4) ];
% 
%  estimated_ball_pos=xhat(1);
%                 estimated_ball_vel=xhat(2);
%                 estimated_beam_ang=xhat(3);
%                 estimated_beam_ang_vel=xhat(4);
%                  
%             if obj.fl_lin
%                 % Feedback linearization control law
%                 v = -obj.K_fl * err;
%                 V_servo = obj.computeControl(xhat(1), xhat(2), xhat(3), xhat(4), v);
%             else
%                 V_servo = -obj.K_lqr * err;
%             end
%             
%             % Beam angle saturation
%             th_max = pi/4;  th_min = -th_max;
%             if ball_pos > th_max
%                 V_servo = min(V_servo, 10*(th_max - ball_pos));
%             elseif ball_pos < th_min
%                 V_servo = max(V_servo, 10*(th_min - ball_pos));
%             end
%             
%             % Update internal state
%             obj.V_servo = V_servo;
%             obj.t_prev  = t;
%         end
%     end
%     
%     methods (Access = public)
%         function [V_servo, theta_d] = stepController(obj, t, ball_pos, beam_ang)
%             V_servo = stepImpl(obj, t, ball_pos, beam_ang);
%             theta_d = obj.theta_d;  % Not used in EKF + LQR mode
%         end
%     end
%     
%     methods (Access = private)
%         function [x_upd, P_upd] = EKF(obj, dt, x_prev, P_prev, z, u)
%             % EKF prediction step
%             p   = x_prev(1);  v = x_prev(2);
%             th  = x_prev(3);  dth = x_prev(4);
%             
%             p_dot   = v;
%             v_dot   = obj.const_1 * sin(th);
%             th_dot  = dth;
%             dth_dot = (-dth + obj.servo_gain * u) / obj.tau_val;
%             x_pred  = x_prev + dt * [p_dot; v_dot; th_dot; dth_dot];
%             
%             % Covariance prediction F = I + A*dt, A = df/dx
%             F = eye(4);
%             F(1,2) = dt;
%             F(2,3) = dt * obj.const_1 * cos(th);
%             F(3,4) = dt;
%             F(4,4) = 1 - dt / obj.tau_val;
%             P_pred = F * P_prev * F.' + obj.Q;
%             
%             % EKF update step
%             H      = obj.C;
%             z_pred = H * x_pred;
%             y      = z - z_pred;
%             S      = H * P_pred * H.' + obj.R;
%             Kk     = P_pred * H.' / S;
%             
%             x_upd = x_pred + Kk * y;
%             P_upd = (eye(4) - Kk * H) * P_pred;
%         end
%         
%         function y = lie2(obj, xhat)
%             x1 = xhat(1); x3 = xhat(3); x4 = xhat(4);
%             term1 = obj.const_1 * sin(x3);
%             term2 = obj.const_2 * ((obj.beam_len/2) - x1) * x4^2 * cos(x3)^2;
%             y = term1 - term2;
%         end
%         
%         function result = phi(obj, x1, x2, x3, x4)
%             % Computes phi(x) for nonlinear feedback law
%             result = (-2 * obj.const_2 * x4^3 * cos(x3) * sin(x3) - ...
%                 (2 * obj.const_2 / obj.tau_val) * x4^2 * cos(x3)^2) * x2 + ...
%                 (obj.const_2 * x4^2 * cos(x3)^2) * (obj.const_1 * sin(x3) - ...
%                 obj.const_2 * ((obj.beam_len/2) - x1) * x4^2 * cos(x3)^2) + ...
%                 (-2 * obj.const_2 * x2 * x4^2 * cos(x3) * sin(x3) - obj.const_1 * x4 * sin(x3) + ...
%                 2 * obj.const_2 * ((obj.beam_len/2) - x1) * x4^3 * (cos(x3)^2 - sin(x3)^2) - ...
%                 (4 * obj.const_2 / obj.tau_val) * ((obj.beam_len/2) - x1) * x4^2 * cos(x3) * sin(x3)) * x4 + ...
%                 (2 * obj.const_2 * x2 * x4 * cos(x3)^2 + obj.const_1 * cos(x3) + ...
%                 6 * obj.const_2 * ((obj.beam_len/2) - x1) * x4^2 * cos(x3) * sin(x3) + ...
%                 (4 * obj.const_2 / obj.tau_val) * ((obj.beam_len/2) - x1) * x4 * cos(x3)^2) * (-x4 / obj.tau_val);
%         end
%         function result = psi(obj, x1, x2, x3, x4)
%             % Computes psi(x) for nonlinear feedback law
%             result = (2 * obj.const_2 * x2 * x4 * cos(x3)^2 + obj.const_1 * cos(x3) + ...
%                 6 * obj.const_2 * ((obj.beam_len/2) - x1) * x4^2 * cos(x3) * sin(x3) + ...
%                 (4 * obj.const_2 / obj.tau_val) * ((obj.beam_len/2) - x1) * x4 * cos(x3)^2) * ...
%                 (obj.servo_gain / obj.tau_val);
%         end
%         function result = computeControl(obj, x1, x2, x3, x4, v)
%             % Feedback linearization control law
%             result = (v - obj.phi(x1, x2, x3, x4)) / obj.psi(x1, x2, x3, x4);
%         end
%     end
% end