% Ran Ma
% 08/12/2021
%
% Anand, L. and Kothari, M., 1996. A computational procedure for
% rate-independent crystal plasticity. Journal of the Mechanics and
% Physics of Solids, 44(4), pp.525-558.
%
classdef CP_Anand_rate_independent < handle
    properties
        % Input parameters
        K;
        nu;
        hardening;
        Euler_phi1;
        Euler_phi;
        Euler_phi2;
        
        % NaCl Slip system
        % ALL slip systems with different orientation
        % Do NOT include the same slip system
        % with opposite normal or shear direction
        %     0: rock salt
        %     1: FCC
        %     2: beta-HMX
        slip_type;
        num_slip_system;
        slip_normal;
        slip_direct;
        
        % Yield stress
        tauy_n1;
        tauy_n;
        
        % history
        Fp_n1;
        Fp_n;
        
        new_active_set; % std::vector<unsigned>
        old_active_set; % std::vector<unsigned>
        
        % slip system
        P_Schmid; % 0.5 * ( s x m + m x s ) in Voigt notation
        P_Schmid_full; % s x m
        
        % constants
        tol = 1.0e-9;
    end
    
    properties (Dependent)
        youngs;
        mu;
        lambda;
    end
    
    methods
        % constructor without default input value
        function obj = CP_Anand_rate_independent(s_type)
            obj.slip_type = s_type;
            if s_type == 1 % FCC
                obj.num_slip_system = 12;
                obj.slip_normal = zeros(obj.num_slip_system, 3);
                obj.slip_direct = zeros(obj.num_slip_system, 3);
                f = 1/sqrt(2.0D0);

                obj.slip_direct(1,1)=0;
                obj.slip_direct(1,2)=-f;
                obj.slip_direct(1,3)=f;
                obj.slip_direct(2,1)=f;
                obj.slip_direct(2,2)=0;
                obj.slip_direct(2,3)=-f;
                obj.slip_direct(3,1)=-f;
                obj.slip_direct(3,2)=f;
                obj.slip_direct(3,3)=0;
                obj.slip_direct(4,1)=0;
                obj.slip_direct(4,2)=f;
                obj.slip_direct(4,3)=f;
                obj.slip_direct(5,1)=f;
                obj.slip_direct(5,2)=0;
                obj.slip_direct(5,3)=f;
                obj.slip_direct(6,1)=f;
                obj.slip_direct(6,2)=-f;
                obj.slip_direct(6,3)=0;
                obj.slip_direct(7,1)=0;
                obj.slip_direct(7,2)=-f;
                obj.slip_direct(7,3)=f;
                obj.slip_direct(8,1)=-f;
                obj.slip_direct(8,2)=0;
                obj.slip_direct(8,3)=-f;
                obj.slip_direct(9,1)=f;
                obj.slip_direct(9,2)=f;
                obj.slip_direct(9,3)=0;
                obj.slip_direct(10,1)=0;
                obj.slip_direct(10,2)=f;
                obj.slip_direct(10,3)=f;
                obj.slip_direct(11,1)=f;
                obj.slip_direct(11,2)=0;
                obj.slip_direct(11,3)=-f;
                obj.slip_direct(12,1)=-f;
                obj.slip_direct(12,2)=-f;
                obj.slip_direct(12,3)=0;

                f = 1/sqrt(3.0D0);
                obj.slip_normal(1,1)=f;
                obj.slip_normal(1,2)=f;
                obj.slip_normal(1,3)=f;
                obj.slip_normal(2,1)=f;
                obj.slip_normal(2,2)=f;
                obj.slip_normal(2,3)=f;
                obj.slip_normal(3,1)=f;
                obj.slip_normal(3,2)=f;
                obj.slip_normal(3,3)=f;
                obj.slip_normal(4,1)=-f;
                obj.slip_normal(4,2)=-f;
                obj.slip_normal(4,3)=f;
                obj.slip_normal(5,1)=-f;
                obj.slip_normal(5,2)=-f;
                obj.slip_normal(5,3)=f;
                obj.slip_normal(6,1)=-f;
                obj.slip_normal(6,2)=-f;
                obj.slip_normal(6,3)=f;
                obj.slip_normal(7,1)=-f;
                obj.slip_normal(7,2)=f;
                obj.slip_normal(7,3)=f;
                obj.slip_normal(8,1)=-f;
                obj.slip_normal(8,2)=f;
                obj.slip_normal(8,3)=f;
                obj.slip_normal(9,1)=-f;
                obj.slip_normal(9,2)=f;
                obj.slip_normal(9,3)=f;
                obj.slip_normal(10,1)=f;
                obj.slip_normal(10,2)=-f;
                obj.slip_normal(10,3)=f;
                obj.slip_normal(11,1)=f;
                obj.slip_normal(11,2)=-f;
                obj.slip_normal(11,3)=f;
                obj.slip_normal(12,1)=f;
                obj.slip_normal(12,2)=-f;
                obj.slip_normal(12,3)=f;
            else
                error('>>> Error: Unknown slip type!');
            end
        end
        
        function initialize(obj)
            % Parse Parameters
            obj.K = 36.7e3; % bulk modulus
            obj.nu = 0.276; % poisson ratio
            tau_ini = 100; % critical resolved shear stress
            hab = 100; % hardening parameter
            
            % initial Euler angle
            obj.Euler_phi1 = 0.01;
            obj.Euler_phi  = 0;
            obj.Euler_phi2 = 0;
            
            % active set
            obj.old_active_set = zeros(2 * obj.num_slip_system + 1, 1);
            obj.new_active_set = zeros(2 * obj.num_slip_system + 1, 1);
            
            % hardening
            nslip = obj.num_slip_system;
            obj.tauy_n  = repmat(tau_ini, nslip, 1);
            obj.tauy_n1 = repmat(tau_ini, nslip, 1);
            obj.hardening = repmat(hab, nslip, nslip);
            
            % Plastic Internal Variable / Strain
            obj.Fp_n1   = eye(3);
            obj.Fp_n    = eye(3);
            
            obj.update_slip_plane();
        end
        
        function update_slip_plane(obj)
            rot_matrix = obj.update_R_matrix();
            
            % initialize slip plane from crystal to spatial
            nslip = obj.num_slip_system;
            n20   = 2 * nslip;
            P_local = zeros(6,n20);
            P_full  = zeros(3,3,n20);
            for ii = 1:nslip
                m_normal = transpose(obj.slip_normal(ii,1:3));
                s_direct = transpose(obj.slip_direct(ii,1:3));
                % Normalize
                m_normal = m_normal / norm(m_normal);
                s_direct = s_direct / norm(s_direct);
                % Assemble Schmid tensor alpha
                P_Schmid_tmp = 0.5*rot_matrix*(m_normal*s_direct' + s_direct*m_normal')*transpose(rot_matrix);
                P_full(:,:,ii) = rot_matrix*(s_direct*m_normal')*transpose(rot_matrix);
                P_local(1,ii) = P_Schmid_tmp(1,1);
                P_local(2,ii) = P_Schmid_tmp(2,2);
                P_local(3,ii) = P_Schmid_tmp(3,3);
                P_local(4,ii) = P_Schmid_tmp(1,2)+P_Schmid_tmp(2,1);
                P_local(5,ii) = P_Schmid_tmp(1,3)+P_Schmid_tmp(3,1);
                P_local(6,ii) = P_Schmid_tmp(2,3)+P_Schmid_tmp(3,2);
                
                P_local(1:6,ii+nslip) = -P_local(1:6,ii);
                P_full(:,:,ii+nslip) = -P_full(:,:,ii);
            end
            obj.P_Schmid = P_local;
            obj.P_Schmid_full = P_full;
        end
        
        function save_state(obj)
            obj.Fp_n   = obj.Fp_n1;    % plastic deformation gradient
            obj.tauy_n = obj.tauy_n1;  % hardening variable
            obj.old_active_set = obj.new_active_set;
        end
        
        function [S, C] = elasticity(obj, E)
            I2 = [1; 1; 1; 0; 0; 0];
            I2I2 = I2 * transpose(I2);
            I4 = diag([1,1,1,0.5,0.5,0.5]);
            C = obj.lambda*I2I2 + 2*obj.mu*I4;
            E_v = obj.e2v( E );
            S_v = C * E_v;
            S = obj.v2s( S_v );
        end
        
        function T_n1 = update(obj, Fn1)
            obj.tauy_n1 = obj.tauy_n;
            obj.new_active_set = obj.old_active_set;

            % Assign crystal slip system
            nslip = obj.num_slip_system;
            n20   = 2 * nslip;
            n21   = n20 + 1;
            P_Schmid_loc = obj.P_Schmid;
            
            % step 1: trial state
            Fe_tr = Fn1 / obj.Fp_n;
            Ce_tr = transpose(Fe_tr) * Fe_tr;
            Ee_tr = 0.5 * ( Ce_tr - eye(3) );
            C_S_v = zeros(6, n20);
            for ii = 1 : n20
                Ce_tr_P = Ce_tr * obj.P_Schmid_full(:,:,ii);
                Ce_tr_P = 0.5 * ( Ce_tr_P + transpose( Ce_tr_P ) );
                C_S_v(:,ii) = obj.s2v( obj.elasticity( Ce_tr_P ) );
            end
            
            % step 2: calculate the trial stress
            obj.Fp_n1 = obj.Fp_n;
            Te_tr = obj.elasticity( Ee_tr );
            Te_tr_v = obj.s2v( Te_tr );
            Te_n1 = obj.v2s( Te_tr_v );
            T_n1 = Fe_tr * Te_n1 * transpose(Fe_tr) / det(Fn1);
            
            
%% *************************************************************************
            % step 3: calculate the trial resolved shear stress
            % step 4: Identify potentially active systems (F > 0)
            tau_100 = repmat(obj.tauy_n, 2, 1);
            trial_active_set = zeros(n21, 1);
            for ii=1:n20
                trial_F = transpose( P_Schmid_loc(:,ii) ) * Te_tr_v - tau_100(ii);
                if ( trial_F > obj.tol)
                    trial_active_set(n21) = trial_active_set(n21) + 1;
                    temp = trial_active_set(n21);
                    trial_active_set(temp) = ii;
                end
            end
            
            if (trial_active_set(n21) < 1)
                return;
            end
            
            % Else Plastic
            
            update_II  = 1;
            residual_all = transpose( P_Schmid_loc ) * Te_tr_v - tau_100;
            
            % Active set iteration
            active_set = obj.old_active_set;
            for S_iter = 0:14
                
                % Check the size of J_active set
                update_I = active_set(n21) >= 1;
                
                % Start while loop
                while update_I
                    
                    % residual, replace with NR iteration
                    n_active = active_set(n21);
                    residual = residual_all(active_set(1 : n_active));
                    
                    % Jacobian matrix for local return mapping
                    jacobian = zeros(n_active, n_active);
                    for ii = 1 : n_active
                        for jj = 1 : n_active
                            % find the No. of slip system in active set
                            alpha = active_set(ii);
                            beta = active_set(jj);
                            currslip_1 = mod(alpha-1, nslip) + 1;
                            currslip_2 = mod(beta-1, nslip) + 1;
                            jacobian(ii,jj) = transpose(P_Schmid_loc(:,alpha)) * C_S_v(:,beta)...
                                            + obj.hardening(currslip_1, currslip_2);

                        end
                    end

                    % Inverse of Jacobian
                    [U,S,V] = svd(jacobian);
                    jacobian_inv = obj.update_jacobian_inverse(S,U,V);
                    Gamma = jacobian_inv * residual;
                    
                    % step 6: compute plastic deformation gradient
                    % step 7: normalize Fp
                    dFp = eye(3);
                    for ii = 1 : n_active
                        alpha = active_set(ii);
                        dFp = dFp + Gamma(ii) * obj.P_Schmid_full(:,:,alpha);
                    end
                    obj.Fp_n1 = dFp * obj.Fp_n;
                    obj.Fp_n1 = obj.Fp_n1 / nthroot(det(obj.Fp_n1), 3);

                    % step 8: compute Fe and stress
                    Te = Te_tr_v;
                    for ii = 1 : n_active
                        alpha = active_set(ii);
                        Te = Te - Gamma(ii) * C_S_v(:,alpha);
                    end
                    Te_n1 = obj.v2s( Te );
                    Fe_n1 = Fn1 / obj.Fp_n1;
                    T_n1 = Fe_n1 * Te_n1 * transpose(Fe_n1) / det(Fn1);
                        
                    
                    % Update I: if Gamma_min < 0, drop it and reconstruct the active set
                    n_active = active_set(n21);
                    [Gamma_min, drop_slip_system] = min(Gamma(1:n_active));
                    update_I = ( Gamma_min < 0 );
                    update_II = ~update_I;
                    if ( update_I )
                        active_set(drop_slip_system:n_active-1) = ...
                            active_set(drop_slip_system+1:n_active);
                        active_set(n_active) = 0;
                        active_set(n21) = active_set(n21) - 1;
                    end
                    
                end % end while loop for Update I
                
                % update II
                if update_II
                    % Find the slip systems not in the active set
                    check_F = transpose(P_Schmid_loc) * obj.s2v( Te_n1 ) - tau_100;
                    n_active = active_set(n21);
                    temp = active_set(1:n_active);
                    check_F(temp) = 0;
                    
                    % Find the maximum F
                    [F_max, add_slip_system] = max(check_F);
                    
                    % Check the maximum F
                    update_I = F_max > obj.tol;
                    if update_I
                       active_set(n21) = active_set(n21) + 1;
                       n_active = active_set(n21);
                       active_set(n_active) = add_slip_system;
                    else
                        break;
                    end
                end % End Update II
            end % End active set iteration

%% *************************************************************************
            
            % step 9: Update the variables
            obj.new_active_set = active_set;
            Gamma_all = zeros(obj.num_slip_system, 2);
            Gamma_all(active_set(1:n_active)) = Gamma;
            obj.tauy_n1 = obj.tauy_n + obj.hardening * sum(Gamma_all, 2);
            
        end
        
        function euler_matrix = update_R_matrix(obj)
            % Note: R matrix here is different from the one from MTEX
            % In fact, R = transpose(R_mtex)
            % Reason: https://mtex-toolbox.github.io/MTEXvsBungeConvention.html
            cos_phi1 = cos(obj.Euler_phi1);
            sin_phi1 = sin(obj.Euler_phi1);
            cos_phi  = cos(obj.Euler_phi);
            sin_phi  = sin(obj.Euler_phi);
            cos_phi2 = cos(obj.Euler_phi2);
            sin_phi2 = sin(obj.Euler_phi2);
            
            phi1_matrix = zeros(3);
            phi_matrix  = zeros(3);
            phi2_matrix = zeros(3);
            
            % Assign each element (Miehe and Schroder, 2001)
            % rotation w.r.t. z axis
            phi1_matrix(1,1) = cos_phi1;
            phi1_matrix(1,2) = sin_phi1;
            phi1_matrix(1,3) = 0.0;
            phi1_matrix(2,1) = -sin_phi1;
            phi1_matrix(2,2) = cos_phi1;
            phi1_matrix(2,3) = 0.0;
            phi1_matrix(3,1) = 0.0;
            phi1_matrix(3,2) = 0.0;
            phi1_matrix(3,3) = 1.0;
            
            % rotation w.r.t. x axis
            phi_matrix(1,1) = 1.0;
            phi_matrix(1,2) = 0.0;
            phi_matrix(1,3) = 0.0;
            phi_matrix(2,1) = 0.0;
            phi_matrix(2,2) = cos_phi;
            phi_matrix(2,3) = sin_phi;
            phi_matrix(3,1) = 0.0;
            phi_matrix(3,2) = -sin_phi;
            phi_matrix(3,3) = cos_phi;
            
            % rotation w.r.t. z axis
            phi2_matrix(1,1) = cos_phi2;
            phi2_matrix(1,2) = sin_phi2;
            phi2_matrix(1,3) = 0.0;
            phi2_matrix(2,1) = -sin_phi2;
            phi2_matrix(2,2) = cos_phi2;
            phi2_matrix(2,3) = 0.0;
            phi2_matrix(3,1) = 0.0;
            phi2_matrix(3,2) = 0.0;
            phi2_matrix(3,3) = 1.0;
            
            euler_matrix = phi2_matrix * phi_matrix * phi1_matrix;
        end
        
        function value = get.youngs(obj)
            value = 3 * obj.K * (1 - 2*obj.nu);
        end
        
        function value = get.mu(obj)
            value = 0.5 * obj.youngs / (1 + obj.nu);
        end
        
        function value = get.lambda(obj)
            value = obj.K - 2*obj.mu/3;
        end
        
    end
    
    methods (Static)
        
        function jacobian_inv = update_jacobian_inverse(S, U, V)
            s_size = length(S);
            S_inv = zeros(s_size, s_size);
            for ii = 1:s_size
                if S(ii,ii)/S(1,1) < 1.0e-6
                    S_inv(ii,ii) = 0;
                else
                    S_inv(ii,ii) = 1 / S(ii,ii);
                end
            end
            jacobian_inv = V * S_inv * transpose(U);
        end
        
        function v = e2v(t)
            v = [t(1,1); t(2,2); t(3,3); ( t(1,2) + t(2,1) );
                ( t(1,3) + t(3,1) ); ( t(2,3) + t(3,2) )];
        end
        
        function t = v2s(v)
            t = [v(1), v(4), v(5);
                 v(4), v(2), v(6);
                 v(5), v(6), v(3)];
        end
        
        function v = s2v(t)
            v = [t(1,1); t(2,2); t(3,3); 0.5 * ( t(1,2) + t(2,1) );
                0.5 * ( t(1,3) + t(3,1) ); 0.5 * ( t(2,3) + t(3,2) )];
        end
        
        function RV = mm10_RT2RVE(RT)
%             implicit none
%             double precision, dimension(3,3), intent(in) :: RT
%             double precision, dimension(6,6), intent(out) :: RV
%
%           voigt notation: 11, 22, 33, 12, 13, 23
%           for strain type tensor only
            RV(1,1)=RT(1,1)^2;
            RV(1,2)=RT(1,2)^2;
            RV(1,3)=RT(1,3)^2;
            RV(1,4)=RT(1,1)*RT(1,2);
            RV(1,5)=RT(1,1)*RT(1,3);
            RV(1,6)=RT(1,3)*RT(1,2);
            RV(2,1)=RT(2,1)^2;
            RV(2,2)=RT(2,2)^2;
            RV(2,3)=RT(2,3)^2;
            RV(2,4)=RT(2,1)*RT(2,2);
            RV(2,5)=RT(2,1)*RT(2,3);
            RV(2,6)=RT(2,3)*RT(2,2);
            RV(3,1)=RT(3,1)^2;
            RV(3,2)=RT(3,2)^2;
            RV(3,3)=RT(3,3)^2;
            RV(3,4)=RT(3,1)*RT(3,2);
            RV(3,5)=RT(3,1)*RT(3,3);
            RV(3,6)=RT(3,3)*RT(3,2);
            RV(4,1)=2*RT(1,1)*RT(2,1);
            RV(4,2)=2*RT(1,2)*RT(2,2);
            RV(4,3)=2*RT(1,3)*RT(2,3);
            RV(4,4)=RT(1,1)*RT(2,2)+RT(2,1)*RT(1,2);
            RV(4,5)=RT(1,1)*RT(2,3)+RT(1,3)*RT(2,1);
            RV(4,6)=RT(1,2)*RT(2,3)+RT(1,3)*RT(2,2);
            RV(5,1)=2*RT(1,1)*RT(3,1);
            RV(5,2)=2*RT(1,2)*RT(3,2);
            RV(5,3)=2*RT(1,3)*RT(3,3);
            RV(5,4)=RT(1,1)*RT(3,2)+RT(1,2)*RT(3,1);
            RV(5,5)=RT(1,1)*RT(3,3)+RT(3,1)*RT(1,3);
            RV(5,6)=RT(1,2)*RT(3,3)+RT(1,3)*RT(3,2);
            RV(6,1)=2*RT(2,1)*RT(3,1);
            RV(6,2)=2*RT(3,2)*RT(2,2);
            RV(6,3)=2*RT(2,3)*RT(3,3);
            RV(6,4)=RT(2,1)*RT(3,2)+RT(2,2)*RT(3,1);
            RV(6,5)=RT(2,1)*RT(3,3)+RT(2,3)*RT(3,1);
            RV(6,6)=RT(2,2)*RT(3,3)+RT(3,2)*RT(2,3);
            %
            %             return
        end
    end
end