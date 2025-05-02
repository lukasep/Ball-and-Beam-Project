classdef studentControllerInterface < matlab.System
    properties (Access = private)
        % EKF state and covariance
        state_estimate = zeros(4,1);   % [p; v; theta; theta_dot]
        P              = eye(4);       % Initial covariance

        % EKF process and measurement noise
%         Q = diag([10, 100, 1, 1]); % best
%         R = diag([0.01, 0.01]); % best
        Q = 10 * diag([1, 1, 0.1, 0.1]);
        R = diag([0.01, 0.01]);

        
%         Q = diag([100, 100, 10, 1]);
%         R = diag([0.01, 0.01]);


        % Physical parameters
        ball_rad   = 0.0254;
        beam_len   = 0.4255;
        g_val      = 9.81;
        servo_gain = 1.5;
        tau_val    = 0.025;

        % Output map (ball_pos and beam_ang)
        C = [1, 0, 0, 0;
            0, 0, 1, 0];

        % Controller gains
        fl_lin = true;                  % Toggle FL vs. LQR
%         K_fl   = [120   150   20    5];
        K_fl = [40, 20, 13.0233, 2.6315]; % best 2
%         K_fl   = [120   121   20    5]; % best
%         K_fl = [23.4521   29.3025   17.2402    5.8720];
        %         K_fl = [3.0002   6.0101    8.0129   4.0423];
        K_lqr  = [20, 25.1525, 13.0233, 2.6315];

        % Bookkeeping variables
        t_prev  = -1;
        V_servo = 0;
        theta_d = 0;

        % Lie constants for feedback linearization
        const_1 = 0;
        const_2 = 0;
    end

    methods
        function obj = studentControllerInterface
            % Initialize constants for Lie derivatives
            obj.const_1 = 5 * obj.g_val * obj.ball_rad / (7 * obj.beam_len);
            obj.const_2 = (5/7) * (obj.ball_rad / obj.beam_len)^2;

            if obj.fl_lin
                disp('Using Feedback Linearization + LQR-v');
            else
                disp('Using pure LQR around origin');
            end
        end
    end

    methods (Access = protected)
        function [V_servo, estimated_ball_pos, estimated_ball_vel, estimated_beam_ang, estimated_beam_ang_vel] = stepImpl(obj, t, ball_pos, beam_ang)
            % Time step calculation
            dt = t - obj.t_prev;
            if obj.t_prev < 0, dt = 1e-3; end

            % EKF prediction and update
            z = [ball_pos; beam_ang];
            u = obj.V_servo;
            [obj.state_estimate, obj.P] = ...
                obj.EKF(dt, obj.state_estimate, obj.P, z, u);

            % Control computation
            [p_ref, v_ref, a_ref] = get_ref_traj(t);
            xhat = obj.state_estimate;
            err  = [xhat(1)-p_ref;
                    xhat(2)-v_ref;
                    xhat(3)-a_ref;
                    xhat(4) ];

            estimated_ball_pos      =   xhat(1);
            estimated_ball_vel      =   xhat(2);
            estimated_beam_ang      =   xhat(3);
            estimated_beam_ang_vel  =   xhat(4);

            if obj.fl_lin
                % Feedback linearization control law
                v = -obj.K_fl * err;
                V_servo = obj.computeControl(xhat(1), xhat(2), xhat(3), xhat(4), v);
            else
                V_servo = -obj.K_lqr * err;
            end

            % Beam angle saturation
            th_max = pi/4;  th_min = -th_max;
            if ball_pos > th_max
                V_servo = min(V_servo, 10*(th_max - ball_pos));
            elseif ball_pos < th_min
                V_servo = max(V_servo, 10*(th_min - ball_pos));
            end

            % Update internal state
            obj.V_servo = V_servo;
            obj.t_prev  = t;
        end
    end

    methods (Access = public)
        function [V_servo, theta_d] = stepController(obj, t, ball_pos, beam_ang)
            V_servo = stepImpl(obj, t, ball_pos, beam_ang);
            theta_d = obj.theta_d;  % Not used in EKF + LQR mode
        end
    end

    methods (Access = private)
        function [x_upd, P_upd] = EKF(obj, dt, x_prev, P_prev, z, u)
            % EKF prediction step
            p   = x_prev(1);  v = x_prev(2);
            th  = x_prev(3);  dth = x_prev(4);

            p_dot   = v;
            v_dot   = obj.const_1 * sin(th);
%             v_dot   = obj.const_1 * sin(th) - obj.const_2 * (obj.beam_len - p) * dth^2 * cos(th)^2;
            th_dot  = dth;
            dth_dot = (-dth + obj.servo_gain * u) / obj.tau_val;
            x_pred  = x_prev + dt * [p_dot; v_dot; th_dot; dth_dot];

            % Covariance prediction F = I + A*dt, A = df/dx
            F = eye(4);
            F(1,2) = dt;
            F(2,3) = dt * obj.const_1 * cos(th);
            F(3,4) = dt;
            F(4,4) = 1 - dt / obj.tau_val;
            P_pred = F * P_prev * F.' + obj.Q;

            % EKF update step
            H      = obj.C;
            z_pred = H * x_pred;
            y      = z - z_pred;
            S      = H * P_pred * H.' + obj.R;
            Kk     = P_pred * H.' / S;

            x_upd = x_pred + Kk * y;
            P_upd = (eye(4) - Kk * H) * P_pred;
        end

        function y = lie2(obj, xhat)
            x1 = xhat(1); x3 = xhat(3); x4 = xhat(4);
            term1 = obj.const_1 * sin(x3);
            term2 = obj.const_2 * ((obj.beam_len/2) - x1) * x4^2 * cos(x3)^2;
            y = term1 - term2;
        end

        function result = phi(obj, x1, x2, x3, x4)
            % Computes phi(x) for nonlinear feedback law
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
            % Computes psi(x) for nonlinear feedback law
            result = (2 * obj.const_2 * x2 * x4 * cos(x3)^2 + obj.const_1 * cos(x3) + ...
                6 * obj.const_2 * ((obj.beam_len/2) - x1) * x4^2 * cos(x3) * sin(x3) + ...
                (4 * obj.const_2 / obj.tau_val) * ((obj.beam_len/2) - x1) * x4 * cos(x3)^2) * ...
                (obj.servo_gain / obj.tau_val);
        end
        function result = computeControl(obj, x1, x2, x3, x4, v)
            % Feedback linearization control law
            result = (v - obj.phi(x1, x2, x3, x4)) / obj.psi(x1, x2, x3, x4);
        end
    end
end

