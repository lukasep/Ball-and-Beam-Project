classdef studentControllerInterface < matlab.System
    properties (Constant, Access = private)
        % Measurement matrix: maps full state [position; velocity; angle; angular vel]
        C = [1 0 0 0; 0 0 1 0];
        % Measurement noise input mapping (identity for direct measurement)
        input = eye(2);
        % Continuous-time state matrix for linearized pendulum dynamics
        A = [0 1 0 0;0 0 1 0;0 0 0 1;0 0 0 0];
        % Continuous-time input matrix for control torque on angular acceleration
        B = [0 0 0 1]';
        % LQR state cost matrix (penalize position, velocity, angle, angular vel)
        Q = diag([1000, 100, 0.1, 1e-10]);
        % LQR input cost (penalize large control effort)
        R = diag(0.05);
        % Integrator gain 
        % for eliminating steady-state error (Ki term)
        Ki = 2;
         
    end
    properties (Access = private)
        % Last time stamp (s)
        previousTime = 0;
        % Last time step duration (s)
        previousDeltaTime = 0.001;
        % Last computed control input (servo command)
        controlInput = 0;
        % Estimated state vector [pos; vel; angle; angVel]
        stateEstimate = [0; 0; 0; 0];
        % Estimated covariance matrix of state uncertainty
        estimateCovariance = diag([0.01, 0.01, 0.01, 0.01]);
        % Accumulated position error for integrator anti-windup
        cumulativeError = -0.1;
        p_old = 0;
    end
    methods(Access = protected)
        function [V_servo, estimatedPosition, estimatedVelocity, estimatedAngle, estimatedAngularVelocity, cumulativeErrorOut] = stepImpl(obj, currentTime, measuredBallPosition, measuredAngle)
            % stepImpl Run one control step: estimate state, compute control
            %  Inputs:
            %    currentTime          - current simulation time (s)
            %    measuredBallPosition - measured ball position along pendulum (m)
            %    measuredAngle        - measured pendulum angle (rad)
            %  Outputs:
            %    V_servo              - computed servo voltage command
            %    estimatedPosition    - Kalman-estimated position (m)
            %    estimatedVelocity    - Kalman-estimated velocity (m/s)
            %    estimatedAngle       - Kalman-estimated pendulum angle (rad)
            %    estimatedAngularVelocity - estimated angular velocity (rad/s)
            %    cumulativeErrorOut   - current integral of position error
            % Physical constants for pendulum-ball system
            ballRadiusMeters = 0.0254;    % ball radius (m)
            pendulumLength    = 0.4255;   % pendulum length (m)
            gravity           = 9.81;     % gravitational acceleration (m/s^2)
            controllerGain    = 1.5;      % gain mapping input voltage to torque
            timeConstant      = 0.025;    % motor time constant (s)
            % Determine time step dt, handling first-call or zero increment
            if currentTime == obj.previousTime
                deltaTime = obj.previousDeltaTime;
            else
                deltaTime = currentTime - obj.previousTime;
            end
            % Get desired reference trajectory at current time
            [refPos, refVel, refAccel] = get_ref_traj(currentTime);
            % Define process and measurement noise covariance matrices
            processNoiseCov      = diag([0.0101, 0.0501, 0.021, 0.2]);
            measurementNoiseCov  = diag([10, 0.08]);
            % Retrieve previous state estimate and control input
            prevState   = obj.stateEstimate;
            prevControl = obj.controlInput;
            % Precompute sin and cos of estimated angle for speed
            c3 = cos(prevState(3));
            s3 = sin(prevState(3));
%             p_curr = measuredBallPosition(1);
%             measuredBallPosition(1) = 0.7 * p_curr + 0.3 * obj.p_old; 
%             obj.p_old = measuredBallPosition(1); 
            % Continuous-time state derivative under nonlinear dynamics
            stateDerivative = [prevState(2);
                               (5*gravity/7)*(ballRadiusMeters/pendulumLength)*s3 - ...
                               (5/7)*(pendulumLength/2 - prevState(1))*(ballRadiusMeters/pendulumLength)^2 * prevState(4)^2 * c3^2;
                               prevState(4);
                               -prevState(4)/timeConstant + (controllerGain/timeConstant)*prevControl];
            % Euler integration to predict next state
            predictedState = prevState + deltaTime * stateDerivative;
            % Linearized discretized A matrix for covariance prediction
            A_pred = [1, deltaTime, 0, 0;
                      (5*deltaTime*ballRadiusMeters^2*prevState(4)^2*c3^2)/(7*pendulumLength^2), 1, ...
                      (5*deltaTime*ballRadiusMeters*(2*pendulumLength*gravity*c3 + ...
                       pendulumLength*ballRadiusMeters*prevState(4)^2*sin(2*prevState(3)) - ...
                       2*ballRadiusMeters*prevState(1)*prevState(4)^2*sin(2*prevState(3))))/(14*pendulumLength^2), ...
                      -(5*deltaTime*ballRadiusMeters^2*prevState(4)*c3^2*(pendulumLength - 2*prevState(1)))/(7*pendulumLength^2);
                      0, 0, 1, deltaTime;
                      0, 0, 0, 1 - deltaTime/timeConstant];
            % A-priori covariance update
            P_plus = A_pred * obj.estimateCovariance * A_pred' + eye(4) * processNoiseCov * eye(4)';
            % Kalman filter update step
            z     = [measuredBallPosition; measuredAngle];           % measurement vector
            yPred = obj.C * predictedState;                         % predicted measurement
            S     = obj.C * P_plus * obj.C' + obj.input * measurementNoiseCov * obj.input';
            K     = P_plus * obj.C' / S;                                % Kalman gain
            updatedState = predictedState + K * (z - yPred);        % state correction
            % Covariance correction (Joseph form for numerical stability)
            updatedCov = (eye(4) - K * obj.C) * P_plus * (eye(4) - K * obj.C)' + ...
                         K * obj.input * measurementNoiseCov * obj.input' * K';
            % Save updated estimates
            obj.stateEstimate     = updatedState;
            obj.estimateCovariance = updatedCov;
            % Nonlinear dynamics feedforward terms for control law
            nonlinearGain = 5/7 * gravity * ballRadiusMeters / pendulumLength;
            f_x = -nonlinearGain/timeConstant * updatedState(4) * cos(updatedState(3)) - ...
                  nonlinearGain * updatedState(4)^2 * sin(updatedState(3));
            g_x = nonlinearGain * controllerGain/timeConstant * cos(updatedState(3));
            % Compute discrete-time system matrices for LQR
            Ad = expm(obj.A * deltaTime);
            nStates = size(obj.A,1);
            if rank(obj.A) == nStates
                % Analytical input matrix for full-rank continuous A
                Bd = obj.A \ (Ad - eye(nStates)) * obj.B;
            else
                % Numerical integration fallback for rank-deficient A
                Nsteps = 100;
                tauVec = linspace(0, deltaTime, Nsteps);
                dtau   = deltaTime/(Nsteps - 1);
                Bd     = zeros(nStates,1);
                for i = 1:Nsteps
                    Bd = Bd + expm(obj.A * tauVec(i)) * obj.B * dtau;
                end
            end
            % Discrete-time LQR via iterative Riccati loop
            maxIter = 100000;
            P_mat   = eye(4);
            K_lqr   = zeros(1,4);
            for iter = 1:maxIter
                P_new = Ad'*P_mat*Ad - (Ad'*P_mat*Bd)*(obj.R + Bd'*P_mat*Bd)^-1*(Bd'*P_mat*Ad) + obj.Q;
                K_new = (obj.R + Bd'*P_mat*Bd)\(Bd'*P_mat*Ad);
                if isequal(K_lqr, K_new)
                    break;
                elseif iter == maxIter
%                     warning('LQR did not converge');
                end
                P_mat = P_new;
                K_lqr = K_new;
            end
            % Form LQR error vector and nominal control
            e = [updatedState(1) - refPos;
                 updatedState(2) - refVel;
                 nonlinearGain * sin(updatedState(3)) - refAccel;
                 nonlinearGain * updatedState(4) * cos(updatedState(3))];
            u_nominal = 1/g_x * (-f_x - K_lqr * e);
            % Integrator anti-windup: accumulate and clamp error
            obj.cumulativeError = obj.cumulativeError + (updatedState(1) - refPos) * deltaTime;
            obj.cumulativeError = max(min(obj.cumulativeError, 1), -1);
            u_int    = obj.Ki * obj.cumulativeError;  % integral correction
            u_final  = u_nominal - u_int;
            % Store final control input
            obj.controlInput = u_final;
            % Outputs to Simulink
            V_servo               = u_final;                 % servo voltage command
            estimatedPosition     = updatedState(1);         % estimated ball position
            estimatedVelocity     = updatedState(2);         % estimated ball velocity
            estimatedAngle        = updatedState(3);         % estimated pendulum angle
            estimatedAngularVelocity = updatedState(4);      % estimated pendulum angular velocity
            cumulativeErrorOut    = obj.cumulativeError;     % integrator state
            % Update time history for next step
            obj.previousTime      = currentTime;
            obj.previousDeltaTime = deltaTime;
        end
    end
    methods(Access = public)
        function [V_servo, desiredAngle] = stepController(obj, currentTime, measuredBallPosition, measuredAngle)
            % stepController Public wrapper to call stepImpl
            V_servo    = stepImpl(obj, currentTime, measuredBallPosition, measuredAngle);
            desiredAngle = obj.desiredAngle;  % placeholder for future desired angle output
        end
    end
end