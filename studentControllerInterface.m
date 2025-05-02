classdef studentControllerInterface < matlab.System
    properties (Constant, Access = private)
        C = [1 0 0 0; 0 0 1 0];
        input = eye(2);
        A = [0 1 0 0;0 0 1 0;0 0 0 1;0 0 0 0];
        B = [0 0 0 1]';
        Q = diag([1000, 400, 0.1, 0.1]);
        R = diag(0.055);
        Ki = 2;
        K_fl = [23.4521, 29.3025, 17.2402, 5.8720]; % this was generated for MATLAB SINE WAVE but ended up working for all 
    end
    properties (Access = private)
        t_prev = 0;
        dt_prev = 0.001;
        control_input = 0;
        state_estimate = [0; 0; 0; 0];
        cov_est = diag([1e-2 1e-2 1e-2 1e-2]);
        cumulative_error = 0;
        const_1 = 0;
        const_2 = 0;
        % Flag to use feedback linearization (true) or DLQR (false)
        use_FL = true;
    end
    methods(Access = protected)
        function [V_servo, est_pos, est_vel, est_ang, est_ang_vel, cumulative_error_out] = stepImpl(obj, t, p_ball, theta)
            r_b = 0.0254;   
            L    = 0.4255; 
            g           = 9.81;   
            K    = 1.5;   
            tau      = 0.025; 
            beam_ang_min = -pi/4;
            beam_ang_max = pi/4;
            if obj.const_1 == 0
                obj.const_1 = 5*g*r_b/(7*L);
                obj.const_2 = (5/7)*(r_b/L)^2;
            end
            if t == obj.t_prev
                dt = obj.dt_prev;
            else
                dt = t - obj.t_prev;
            end

            % EKF
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);
            process_noise_cov      = diag([0.02, 0.055, 0.01, 0.1]);
            measurement_noise_cov  = diag([1, 0.09]);
            prev_state   = obj.state_estimate;
            prev_control = obj.control_input;
            c3 = cos(prev_state(3));
            s3 = sin(prev_state(3));
            state_deriv = [prev_state(2);
                               (5*g/7)*(r_b/L)*s3 - ...
                               (5/7)*(L/2 - prev_state(1))*(r_b/L)^2 * prev_state(4)^2 * c3^2;
                               prev_state(4);
                               -prev_state(4)/tau + (K/tau)*prev_control];
            predicted_state = prev_state + dt * state_deriv;
            A_pred = [1, dt, 0, 0;
                      (5*dt*r_b^2*prev_state(4)^2*c3^2)/(7*L^2), 1, ...
                      (5*dt*r_b*(2*L*g*c3 + ...
                       L*r_b*prev_state(4)^2*sin(2*prev_state(3)) - ...
                       2*r_b*prev_state(1)*prev_state(4)^2*sin(2*prev_state(3))))/(14*L^2), ...
                      -(5*dt*r_b^2*prev_state(4)*c3^2*(L - 2*prev_state(1)))/(7*L^2);
                      0, 0, 1, dt;
                      0, 0, 0, 1 - dt/tau];
            P_plus = A_pred * obj.cov_est * A_pred' + eye(4) * process_noise_cov * eye(4)';
            z     = [p_ball; theta];     
            y_Pred = obj.C * predicted_state;                    
            S     = obj.C * P_plus * obj.C' + obj.input * measurement_noise_cov * obj.input';
            K     = P_plus * obj.C' / S;                        
            updated_state = predicted_state + K * (z - y_Pred);  
            % Covariance correction (Joseph form for numerical stability)
            updated_cov = (eye(4) - K * obj.C) * P_plus * (eye(4) - K * obj.C)' + ...
                         K * obj.input * measurement_noise_cov * obj.input' * K';
            obj.state_estimate     = updated_state;
            obj.cov_est = updated_cov;

            % Control calculation - choose between feedback linearization and LQR
            if obj.use_FL
                % FL
                error_x1 = updated_state(1) - p_ball_ref;
                error_x2 = updated_state(2) - v_ball_ref;
                error_x3 = obj.lie2(updated_state(1), updated_state(2), updated_state(3), updated_state(4)) - a_ball_ref;
                error_x4 = obj.lie3(updated_state(1), updated_state(2), updated_state(3), updated_state(4)) - 0;
                error = [error_x1; error_x2; error_x3; error_x4];
                v = -obj.K_fl * error;
                u_final = obj.computeControl(updated_state(1), updated_state(2), updated_state(3), updated_state(4), v);
            else
                % DLQR
                nonlinear_gain = 5/7 * g * r_b / L;
                f_x = -nonlinear_gain/tau * updated_state(4) * cos(updated_state(3)) - ...
                      nonlinear_gain * updated_state(4)^2 * sin(updated_state(3));
                g_x = nonlinear_gain * K/tau * cos(updated_state(3));
                Ad = expm(obj.A * dt);
                n_states = size(obj.A,1);
                if rank(obj.A) == n_states
                    Bd = obj.A \ (Ad - eye(n_states)) * obj.B;
                else
                    N_steps = 100;
                    tau_vec = linspace(0, dt, N_steps);
                    d_tau   = dt/(N_steps - 1);
                    Bd     = zeros(n_states,1);
                    for i = 1:N_steps
                        Bd = Bd + expm(obj.A * tau_vec(i)) * obj.B * d_tau;
                    end
                end

                max_iter = 100000;
                P_mat   = eye(4);
                K_lqr   = zeros(1,4);
                for iter = 1:max_iter
                    P_new = Ad'*P_mat*Ad - (Ad'*P_mat*Bd)*(obj.R + Bd'*P_mat*Bd)^-1*(Bd'*P_mat*Ad) + obj.Q;
                    K_new = (obj.R + Bd'*P_mat*Bd)\(Bd'*P_mat*Ad);
                    if isequal(K_lqr, K_new)
                        break;
%                     elseif iter == maxIter
%                         warning('LQR did not converge');
                    end
                    P_mat = P_new;
                    K_lqr = K_new;
                end
                e = [updated_state(1) - p_ball_ref;
                     updated_state(2) - v_ball_ref;
                     nonlinear_gain * sin(updated_state(3)) - a_ball_ref;
                     nonlinear_gain * updated_state(4) * cos(updated_state(3))];
                u_n = 1/g_x * (-f_x - K_lqr * e); % torque control correction
                obj.cumulative_error = obj.cumulative_error + (updated_state(1) - p_ball_ref) * dt;
                obj.cumulative_error = max(min(obj.cumulative_error, 1), -1);
                u_int    = obj.Ki * obj.cumulative_error;  % integral correction
                u_final  = u_n - u_int;
            end
            obj.control_input = u_final;
            V_servo               = u_final;
            if theta > beam_ang_max
                V_servo = min(V_servo, 10 * (beam_ang_max - theta));
            elseif theta < beam_ang_min
                V_servo = max(V_servo, 10 * (beam_ang_min - theta));
            end
            
            est_pos     = updated_state(1);    
            est_vel     = updated_state(2);     
            est_ang        = updated_state(3);       
            est_ang_vel = updated_state(4);      
            cumulative_error_out    = obj.cumulative_error;    
            obj.t_prev      = t;
            obj.dt_prev = dt;
        end
    end
    methods(Access = public)
        function [V_servo, theta_d] = stepController(obj, t, p_ball, theta)
            [V_servo, ~, ~, ~, ~, ~] = stepImpl(obj, t, p_ball, theta);
            theta_d = 0; 
        end
    end
    
    methods(Access = private)
        function result = lie1(obj, x1, x2, x3, x4)
            % First Lie derivative
            result = x2;
        end
        
        function result = lie2(obj, x1, x2, x3, x4)
            % Second Lie derivative
            L = 0.4255;
            result = obj.const_1 * sin(x3) - obj.const_2 * ((L/2) - x1) * x4^2 * cos(x3)^2;
        end
        
        function result = lie3(obj, x1, x2, x3, x4)
            % Third Lie derivative
            L = 0.4255;
            tau = 0.025;
            result = obj.const_2 * x2 * x4^2 * cos(x3)^2 + ...
                obj.const_1 * x4 * cos(x3) + ...
                2 * obj.const_2 * ((L/2) - x1) * x4^3 * cos(x3) * sin(x3) + ...
                (2 * obj.const_2 / tau) * ((L/2) - x1) * x4^2 * cos(x3)^2;
        end
        
        function result = phi(obj, x1, x2, x3, x4)
            % phi(x) part of the control law
            L = 0.4255;
            tau = 0.025;
            result = (-2 * obj.const_2 * x4^3 * cos(x3) * sin(x3) - ...
                (2 * obj.const_2 / tau) * x4^2 * cos(x3)^2) * x2 + ...
                (obj.const_2 * x4^2 * cos(x3)^2) * (obj.const_1 * sin(x3) - ...
                obj.const_2 * ((L/2) - x1) * x4^2 * cos(x3)^2) + ...
                (-2 * obj.const_2 * x2 * x4^2 * cos(x3) * sin(x3) - obj.const_1 * x4 * sin(x3) + ...
                2 * obj.const_2 * ((L/2) - x1) * x4^3 * (cos(x3)^2 - sin(x3)^2) - ...
                (4 * obj.const_2 / tau) * ((L/2) - x1) * x4^2 * cos(x3) * sin(x3)) * x4 + ...
                (2 * obj.const_2 * x2 * x4 * cos(x3)^2 + obj.const_1 * cos(x3) + ...
                6 * obj.const_2 * ((L/2) - x1) * x4^2 * cos(x3) * sin(x3) + ...
                (4 * obj.const_2 / tau) * ((L/2) - x1) * x4 * cos(x3)^2) * (-x4 / tau);
        end
        
        function result = psi(obj, x1, x2, x3, x4)
            % psi(x) part of the control law
            L = 0.4255;
            tau = 0.025;
            K = 1.5;
            result = (2 * obj.const_2 * x2 * x4 * cos(x3)^2 + obj.const_1 * cos(x3) + ...
                6 * obj.const_2 * ((L/2) - x1) * x4^2 * cos(x3) * sin(x3) + ...
                (4 * obj.const_2 / tau) * ((L/2) - x1) * x4 * cos(x3)^2) * ...
                (K / tau);
        end
        
        function result = computeControl(obj, x1, x2, x3, x4, v)
            % Feedback linearization control law: u = (v - phi(x))/psi(x)
            result = (v - obj.phi(x1, x2, x3, x4)) / obj.psi(x1, x2, x3, x4);
        end
    end
end
