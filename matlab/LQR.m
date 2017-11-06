% Class for calculating linear feedback law
%
% Discrete time linear quadratic regulator
% Calculates the optimal linear feedback law in various diff. situations
%
% Inputs:
% Q,R - Weight matrices
% Qf - Penalty of final state
% A,B - Cts. Matrices of the system differential equation
% C - Observation matrix
% 
% Functions - Infinite horizon or finite horizon, as well as time
% varying/invariant
%
% Outputs:
% K - Feedback gain
%

classdef LQR < handle

    properties   
        % Flag for linear time invariant system
        LTI
        % State penalty matrix
        Q
        % Input penalty matrix
        R
        % Time horizon
        N
        % discrete A and B matrices
        Ad, Bd
        % observation matrix
        C
        % final state penalty matrix
        Qf
    end

    methods
        
        %% constructor for LQR
        function obj = LQR(Q,R,Qf,A,B,C,N,h,discretized)
            
            obj.N = N;
            obj.Q = C'*Q*C;
            obj.R = R;
            obj.Qf = C'*Qf*C;
            obj.LTI = length(size(A)) == 2;
            
            % get the model matrices
            if ~discretized && obj.LTI
                warning('LQR only discretizes LTI systems!');
                % get Ad, Bd matrices
                obj.discretizeMatrices(A,B,h);
            else
                % discrete matrices provided
                obj.Ad = A;
                obj.Bd = B;
            end
            obj.C = C;

        end
        
        %% Get discrete matrices for dLQR
        function discretizeMatrices(obj,A,B,h)
            
            dimu = size(B,2);
            dimx = size(B,1);
            % trick to get discrete time versions
            Mat = [A, B; zeros(dimu, dimx + dimu)];
            MD = expm(h * Mat);
            obj.Ad = MD(1:dimx,1:dimx);
            obj.Bd = MD(1:dimx,dimx+1:end);
        end
        
        %% Check for output controllability in LTI systems
        % make sure the system is controllable/reachable
        % otherwise give an error
        function assertOutputControllability(obj)
            
            if (obj.LTI)
                dimu = size(obj.Bd,2);
                dimx = size(obj.Bd,1);
                dimy = size(obj.C,1);
                % construct controllability Kalman matrix
                K = zeros(dimx,dimx*dimu);
                for i = 0:dimx-1
                    K(:,(dimu*i+1):dimu*(i+1)) = obj.C*(obj.Ad^i)*obj.Bd;
                end
                assert(rank(K) == dimy, 'System is not controllable!');
            end
        end
        
        %% Infinite Horizon Regulator, assumed to be LTI
        function K = computeInfHorizon(obj)
            
            % inf horizon            
            Pinf = obj.Q;
            % iterate till convergence
            delta = 1e-3;
            diff = 1;
            i = 0;
            while diff > delta
                Ppre = Pinf;
                Pinf = dynamic_riccati(Pinf,obj.Q,obj.R,obj.Ad,obj.Bd);
                matDiff = Ppre - Pinf;
                diff = norm(matDiff,'fro');
                i = i+1;
            end
            fprintf('Iterated dynamic riccati equation %d times\n', i);
            K = -(obj.R + obj.Bd'*Pinf*obj.Bd)\(obj.Bd'*Pinf*obj.Ad);
        end
        
        %% LTI, simplest discrete-time case, i.e. no tracking
        function K = computeFinHorizonLTI(obj)
            
            % finite horizon
            dimu = size(obj.Bd,2);
            dimx = size(obj.Bd,1);
            Q = obj.Q;
            R = obj.R;
            A = obj.Ad;
            B = obj.Bd;
            Qf = obj.Qf;
            N = obj.N;

            P = zeros(dimx,dimx,N+1);
            P(:,:,end) = Qf;

            % calculate optimal K
            K = zeros(dimu,dimx,N);            
            for i = N:-1:1
                P(:,:,i) = dynamic_riccati(P(:,:,i+1),Q,R,A,B);
                K(:,:,i) = -(R + B'*P(:,:,i+1)*B)\(B'*P(:,:,i+1)*A);
            end
        end
        
        %% LTV for finite horizon
        function [K,P] = computeFinHorizonLTV(obj)
            
            % Arrays A and B hold time-varying matrices
            dimu = size(obj.Bd,2);
            dimx = size(obj.Bd,1);
            Q = obj.Q;
            R = obj.R;
            A = obj.Ad;
            B = obj.Bd;
            Qf = obj.Qf;
            N = obj.N;
            
            P = zeros(dimx,dimx,N+1);
            P(:,:,end) = Qf;

            % calculate optimal K
            K = zeros(dimu,dimx,N);
            for i = N:-1:1
                P(:,:,i) = dynamic_riccati(P(:,:,i+1),Q,R,A(:,:,i),B(:,:,i));
                K(:,:,i) = -(R + B(:,:,i)'*P(:,:,i+1)*B(:,:,i))\...
                            (B(:,:,i)'*P(:,:,i+1)*A(:,:,i));
            end
        end
        
        %% Discrete-time finite-horizon trajectory tracking
        
        % discrete-time finite-horizon trajectory tracking
        % 
        % Outputs:
        % Kbar is K and uff combined in a n+1 x n+1 feedback/feedforward
        % matrix
        function Kbar = computeFinHorizonTracking(obj,s)
                        
            dimu = size(obj.Bd,2);
            dimx = size(obj.Bd,1);
            A = obj.Ad;
            B = obj.Bd;
            N = obj.N;
            C = obj.C;
            s = C'*((C*C')\s);
            s = [s, zeros(dimx,1)];
            
            % form the time varying matrices Abar and Bbar
            Abar = zeros(dimx+1,dimx+1,N);
            Bbar = zeros(dimx+1,dimu,N);
            for i = 1:N
                Abar(:,:,i) = [A(:,:,i), A(:,:,i)*s(:,i) - s(:,i+1); ...
                               zeros(1,dimx), 1];
                Bbar(:,:,i) = [B(:,:,i); zeros(1,dimu)];
            end
            
            obj.Ad = Abar;
            obj.Bd = Bbar;
            % modify Q and Qf
            obj.Q = diag([diag(obj.Q);0]);
            obj.Qf = diag([diag(obj.Qf);0]);
            
            Kbar = obj.computeFinHorizonLTV();
            
        end
        
        % The error/state-form of LQR tracking is used
        % 
        % 
        % Outputs:
        % K and uff calculated using the error form
        % Hack in error-form:
        % Ks is actually added to uff instead! I don't see why not!
        function [K,uff] = computeFinHorizonTracking2(obj,s,error_form)
            
            dimu = size(obj.Bd,2);
            dimx = size(obj.Bd,1);
            Q = obj.Q; 
            R = obj.R;
            Qf = obj.Qf;
            A = obj.Ad;
            B = obj.Bd;
            N = obj.N;
            C = obj.C;
            s = C'*((C*C')\s);
            
            % value matrix entries
            P = zeros(dimx,dimx,N+1);
            P(:,:,end) = Qf;
            nu = zeros(dimx,N+1);
            nu(:,end) = -Qf*s(:,end);

            % calculate optimal K and uff
            K = zeros(dimu,dimx,N);
            uff = zeros(dimu,N);
            
            for i = N:-1:1
                P(:,:,i) = dynamic_riccati(P(:,:,i+1),Q,R,A(:,:,i),B(:,:,i));
                K(:,:,i) = -(R + B(:,:,i)'*P(:,:,i+1)*B(:,:,i))\(B(:,:,i)'*P(:,:,i+1)*A(:,:,i));
                nu(:,i) = (A(:,:,i) + B(:,:,i)*K(:,:,i))'*nu(:,i+1) - Q*s(:,i);
                %uff(:,i) = (R + B(:,:,i)'*P(:,:,i+1)*B(:,:,i))\B(:,:,i)'*...
                %            (P(:,:,i+1)*s(:,i));
                uff(:,i) = -(R + B(:,:,i)'*P(:,:,i+1)*B(:,:,i))\(B(:,:,i)'*nu(:,i));
                if (error_form)
                    uff(:,i) = uff(:,i) + K(:,:,i)*s(:,i);
                end
                %}
            end
            
        end
        
        % The error/state-form of LQR tracking is used
        % 
        % 
        % Outputs:
        % K and uff calculated using the error form
        % Hack in error-form:
        % Ks is actually added to uff instead! I don't see why not!
        function [K,uff] = computeTrackingWithDist(obj,s,d)
            
            dimu = size(obj.Bd,2);
            dimx = size(obj.Bd,1);
            Q = obj.Q; 
            R = obj.R;
            Qf = obj.Qf;
            A = obj.Ad;
            B = obj.Bd;
            N = obj.N;
            C = obj.C;
            s = C'*((C*C')\s);
            d = C'*((C*C')\d);
            
            % value matrix entries
            P = zeros(dimx,dimx,N);
            P(:,:,end) = Qf;
            nu = zeros(dimx,N);
            b = zeros(dimx,N);
            nu(:,end) = -Qf*s(:,end);
            b(:,end) = -Qf*d(:,end);
            
            % calculate optimal K and uff
            K = zeros(dimu,dimx,N-1);
            uff = zeros(dimu,N-1);
            
            for i = N-1:-1:1
                P(:,:,i) = dynamic_riccati(P(:,:,i+1),Q,R,A(:,:,i),B(:,:,i));
                K(:,:,i) = -(R + B(:,:,i)'*P(:,:,i+1)*B(:,:,i))\(B(:,:,i)'*P(:,:,i+1)*A(:,:,i));
                nu(:,i) = (A(:,:,i) + B(:,:,i)*K(:,:,i))'*nu(:,i+1) - Q*s(:,i);
                b(:,i) = (A(:,:,i) + B(:,:,i)*K(:,:,i))'*b(:,i+1) - P(:,:,i)*d(:,i);
                uff(:,i) = -(R + B(:,:,i)'*P(:,:,i+1)*B(:,:,i))\B(:,:,i)'*(nu(:,i+1) - b(:,i+1));
                uff(:,i) = uff(:,i) + K(:,:,i)*s(:,i);
            end
        end
        
    end
end