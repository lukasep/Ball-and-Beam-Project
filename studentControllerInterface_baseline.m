classdef studentControllerInterface < matlab.System
    properties (Access = private)
        t_prev = -1;
        V_servo = 0;
        ball_rad = 0.0254;   
        beam_len = 0.4255;   
        g_val = 9.81;        
        servo_gain = 1.5;    
        tau_val = 0.025;     
        state_estimate = [0; 0; 0; 0];
        Q = zeros(4,4);
        R = 1;
        K = [];
        A_lin;
        B_lin;
        C = [1, 0, 0, 0;
             0, 0, 1, 0];
        observer_gain =[22.9338 , 1.0388;
                        130.7798, 12.977;
                        0.9570  , 23.0662;
                        10.9830 , 132.2187];
    end
    methods(Access = protected)
        function V_servo = stepImpl(obj, t, ball_pos, beam_ang)
            beam_ang_min = -pi/4;
            beam_ang_max = pi/4;
            dt = t - obj.t_prev;
            if obj.t_prev < 0, dt = 0.001; end
            [p_ball_ref, v_ball_ref, ~] = get_ref_traj(t);
            measurements = [ball_pos; beam_ang];
            control_input = obj.V_servo;
            obj.state_estimate = obj.runObserver(dt, obj.state_estimate, measurements, control_input);
            error = [obj.state_estimate(1) - p_ball_ref;
                     obj.state_estimate(2) - v_ball_ref;
                     obj.state_estimate(3);
                     obj.state_estimate(4)];
            V_servo = -obj.K * error;
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
            const_1 = (5 * obj.g_val * obj.ball_rad) / (7 * obj.beam_len);
            obj.A_lin = [0,         1,       0,           0;
                         0,         0,  const_1,           0;
                         0,         0,       0,           1;
                         0,         0,       0,    -1/obj.tau_val];
            obj.B_lin = [0; 0; 0; obj.servo_gain/obj.tau_val];
            obj.Q(1,1) = 100;
            obj.Q(2,2) = 10;
            obj.Q(3,3) = 1;
            obj.Q(4,4) = 1;
            obj.R = 2;
            
            obj.K = [20, 25.1525, 13.0233, 2.6315];
        end
        function [V_servo, dummy] = stepController(obj, t, ball_pos, beam_ang)
            V_servo = stepImpl(obj, t, ball_pos, beam_ang);
            dummy = 0;
        end
    end
    methods(Access = private)
        function state_estimate = runObserver(obj, time_step, previous_estimate, measurements, control_input)
            ball_position    = previous_estimate(1);
            ball_velocity    = previous_estimate(2);
            beam_angle       = previous_estimate(3);
            beam_ang_vel     = previous_estimate(4);
            state_derivative = zeros(4, 1);
            state_derivative(1) = ball_velocity;
            state_derivative(2) = (5 * obj.g_val * obj.ball_rad/(7 * obj.beam_len)) * sin(beam_angle);
            state_derivative(3) = beam_ang_vel;
            state_derivative(4) = (-beam_ang_vel + obj.servo_gain * control_input) / obj.tau_val;
            measurement_error = measurements - obj.C * previous_estimate;
            state_derivative = state_derivative + obj.observer_gain * measurement_error;
            state_estimate = previous_estimate + state_derivative * time_step;
        end
    end
end