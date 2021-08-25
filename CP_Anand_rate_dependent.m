% Ran Ma
% 03/20/2020
% Anand, L., 2004. Single-crystal elasto-viscoplasticity:
% application to texture evolution in polycrystalline metals at large
% strains. Computer methods in applied mechanics and engineering,
% 193(48-51), pp.5359-5383.
%
classdef CP_Anand_rate_dependent < handle
    properties
        % Input parameters
        K;
        nu;
        tau_ini;
        hardening;
        prm_Q;    % Activation energy
        prm_C;    % Fitting parameter
        exponent; % Power-law parameter
        Euler_phi1;
        Euler_phi;
        Euler_phi2;
        
        % Computed Parameters
        %   Elastic moduli
        lambda;
        mu;
        youngs;
        viscosity;
        
        %   Temperature, time increment
        temperature;
        dt;
        
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
        F_n1;   F_n;     % deformation gradient
        Fp_n1;  Fp_n;    % plastic deformation gradient
        T_n1;   T_n;     % Cauchy stress
        Te_n1;  Te_n;    % stress in intermediate configuration
        
        % algorithmic tangent stiffness
        C_stiff;
        
        % slip system
        P_Schmid; % 0.5 * ( s x m + m x s ) in Voigt notation
        P_Schmid_full; % s x m
        
        % constants
        tol = 1.0e-10;
    end
    
    methods
        % constructor without default input value
        function obj = CP_Anand_rate_dependent(s_type)
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
            obj.tau_ini = 3; % critical resolved shear stress
%             obj.hardening = 100; % hardening parameter
            obj.hardening = 0; % hardening parameter
            obj.prm_Q = 14; % activation energy
            obj.prm_C = 1.0; % shape factor
            obj.exponent = 10.0;
            
            % initial Euler angle
            obj.Euler_phi1 = pi/3;
            obj.Euler_phi  = 0;
            obj.Euler_phi2 = 0;
            
            % Elastic Tangent Operator
            obj.youngs = 3 * obj.K * (1 - 2*obj.nu);
            obj.mu = 0.5 * obj.youngs / (1 + obj.nu);
            obj.lambda = obj.K - 2*obj.mu/3;
            
            % Plastic Internal Variable / Strain
            obj.F_n1    = eye(3);
            obj.F_n     = eye(3);
            obj.Fp_n1   = eye(3);
            obj.Fp_n    = eye(3);
            obj.T_n1    = zeros(3);
            obj.T_n     = zeros(3);
            obj.Te_n1   = zeros(3);
            obj.Te_n    = zeros(3);
            obj.tauy_n  = obj.tau_ini;
            obj.tauy_n1 = obj.tau_ini;
            
            obj.update_slip_plane();
            obj.save_state();
        end
        
        function update(obj, F_in)
            if length(F_in) == 2
                strain_n1 = eye(3);
                strain_n1(1:2, 1:2) = F_in;
                obj.run_3D_update(strain_n1);
            else
                strain_n1 = F_in;
                obj.run_3D_update(strain_n1);
            end
        end
        
        function update_parameters(obj, pass_temp, pass_time)
            obj.temperature = pass_temp + 273.0;
            obj.dt = pass_time;
        end
        
        function update_slip_plane(obj)
            rot_matrix = obj.update_R_matrix();
            
            % initialize slip plane from crystal to spatial
            nslip = obj.num_slip_system;
            P_local = zeros(6,nslip);
            P_full  = zeros(3,3,nslip);
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
            end
            obj.P_Schmid = P_local;
            obj.P_Schmid_full = P_full;
        end
        
        function save_state(obj)
            obj.F_n    = obj.F_n1;     % deformation gradient
            obj.Fp_n   = obj.Fp_n1;    % plastic deformation gradient
            obj.T_n    = obj.T_n1;     % Cauchy stress
            obj.Te_n   = obj.Te_n1;    % stress in intermediate configuration
            obj.tauy_n = obj.tauy_n1;  % hardening variable
        end
        
        function sigma = stress(obj)
            sigma = obj.T_n1;
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
        
        function [dGamma, dGamma_dtau] = crystal_plasticity(obj, tau)
            tau_abs = abs( tau / obj.tauy_n1 );
            dGamma = obj.dt * obj.viscosity * ( tau_abs ) .^ (obj.exponent) .* sign(tau);
            a = obj.dt * obj.viscosity * obj.exponent / obj.tauy_n1;
            dGamma_dtau = a * ( tau_abs ) .^ (obj.exponent - 1);
        end
        
        function [temp, Q66] = run_3D_update(obj, strain_n1)
            obj.F_n1 = strain_n1;
            obj.tauy_n1 = obj.tauy_n;
            I4 = eye(6);
%             I4(4,4) = 0.5; I4(5,5) = 0.5; I4(6,6) = 0.5;
            prm_R    = 1.986e-3; % Gas constant (kcal/(mol.K))
            % this viscosity is different from the original one in deal.ii
            obj.viscosity = (obj.prm_C*exp((-1.0)*obj.prm_Q/(prm_R*obj.temperature)));
            
            % trial state
            Fe_tr = strain_n1 / obj.Fp_n;
            Ce_tr = transpose(Fe_tr) * Fe_tr;
            Ee_tr = 0.5 * ( Ce_tr - eye(3) );
            Te_tr = obj.elasticity( Ee_tr );
            C_S_v = zeros(6,obj.num_slip_system);
            for ii = 1 : obj.num_slip_system
                Ce_tr_P = Ce_tr * obj.P_Schmid_full(:,:,ii);
                Ce_tr_P = 0.5 * ( Ce_tr_P + transpose( Ce_tr_P ) );
                C_S_v(:,ii) = obj.s2v( obj.elasticity( Ce_tr_P ) );
            end
            
            % transfer to Voigt notation
%             Ce_tr_v = obj.e2v( Ce_tr );
%             Ee_tr_v = obj.e2v( Ee_tr );
%             Te_n_v  = obj.s2v( obj.Te_n );
            Te_tr_v = obj.s2v( Te_tr );
            
            % prepare for Newton iteration
            Te_n1_v = obj.s2v( obj.Te_n );
            Te_tmp = Inf(6,1);
            Te_n1_norm = norm( Te_n1_v );
            if Te_n1_norm < obj.tol
                Te_n1_norm = 1;
            end
            Te_res = norm( Te_n1_v - Te_tmp ) / Te_n1_norm;
            
            while Te_res > obj.tol
                % store old stress from previous iteration
                Te_tmp = Te_n1_v;
                
                % first level of iterations
                % fix hardening variable, update stress
                while true
                    
                    % update slip increment
                    tau = transpose(obj.P_Schmid) * Te_n1_v;
                    [dGamma, dGamma_dtau] = obj.crystal_plasticity(tau);
                    dGamma_dTe = repmat(dGamma_dtau',6,1) .* obj.P_Schmid;
                    
                    % compute residual and Jacobian
                    % ATTENTION: J_n is not the actuall stiffness
                    %            J_n(:,4:6) is multiplied by two
                    %            This is because it will be multiplied
                    %            by stress instead of strain!
                    G_n = Te_n1_v - Te_tr_v + C_S_v * dGamma;
                    J_n = I4 + C_S_v * transpose(dGamma_dTe);
                    
                    % finite difference Jacobian
%                     J_fd = zeros(6,6);
%                     tol_fd = 1.0e-6;
%                     for ii = 1:6
%                         Te_fd = Te_n1_v;
%                         Te_fd(ii) = Te_fd(ii) + tol_fd;
%                         tau_fd = transpose(obj.P_Schmid) * Te_fd;
%                         [dGamma_fd, ~] = obj.crystal_plasticity(tau_fd);
%                         G_fd = Te_fd - Te_tr_v + C_S_v * dGamma_fd;
%                         J_fd(:,ii) = (G_fd - G_n) / tol_fd;
%                     end
                    
                    if ( norm(G_n) / Te_n1_norm < obj.tol )
                        break;
                    end
                    
                    % update stress
                    Te_n1_v = Te_n1_v - J_n \ G_n;
                    Te_n1_norm = norm( Te_n1_v );
                    if Te_n1_norm < obj.tol
                        Te_n1_norm = 1;
                    end
                end
                
                % second level of iterations
                % fix stress, update hardening variable
                tau = transpose(obj.P_Schmid) * Te_n1_v;
                obj.tauy_n1 = obj.tauy_n + obj.hardening * sum(abs(dGamma));
                [dGamma, dGamma_dtau] = obj.crystal_plasticity(tau);
                dGamma_dtau_m = dGamma_dtau .* abs(tau) / obj.tauy_n1;
                res_tauy = obj.tauy_n1 - obj.tauy_n - obj.hardening * sum(abs(dGamma));
                J_res = 1 + obj.hardening * sum(abs(dGamma_dtau_m));
                while abs(res_tauy)/obj.tauy_n1 > obj.tol
                    obj.tauy_n1 = obj.tauy_n1 - res_tauy / J_res;
                    [dGamma, dGamma_dtau] = obj.crystal_plasticity(tau);
                    dGamma_dtau_m = dGamma_dtau .* abs(tau) / obj.tauy_n1;
                    res_tauy = obj.tauy_n1 - obj.tauy_n - obj.hardening * sum(abs(dGamma));
                    J_res = 1 + obj.hardening * sum(abs(dGamma_dtau_m));
                end
                
                % update global residual
                Te_res = norm( Te_n1_v - Te_tmp ) / Te_n1_norm;
            end
            
            % compute plastic deformation gradient
            dFp = eye(3);
            for ii = 1 : obj.num_slip_system
                dFp = dFp + dGamma(ii) * obj.P_Schmid_full(:,:,ii);
            end
            obj.Fp_n1 = dFp * obj.Fp_n;
            obj.Fp_n1 = obj.Fp_n1 / nthroot(det(obj.Fp_n1), 3);
            
            % compute Cauchy stress and stiffness
            obj.Te_n1 = obj.v2s( Te_n1_v );
            Fe_n1 = strain_n1 / obj.Fp_n1;
            obj.T_n1 = Fe_n1 * obj.Te_n1 * transpose(Fe_n1) / det(strain_n1);
            
            temp = obj.Te_n1;
            
            % compute algorithmic stiffness
            % =============================================================
            % attention: Ran add these lines to update dGamma_dtau
            % =============================================================
            nslip = obj.num_slip_system;
            ccc = zeros(nslip);
            bbb = (obj.tauy_n1) ^ (obj.exponent-1) * ( obj.exponent * obj.hardening * sign( dGamma ) );
            tauym = (obj.tauy_n1) ^ (obj.exponent);
            for ii = 1 : nslip
                ccc(ii,:) = dGamma(ii) * bbb;
                ccc(ii,ii) = ccc(ii,ii) + tauym;
            end
            aaa = obj.dt * obj.viscosity * obj.exponent * abs(tau) .^ (obj.exponent-1);
            dGamma_dtau = ccc \ aaa;
            dGamma_dTe = repmat(dGamma_dtau',6,1) .* obj.P_Schmid;
            J_n = I4 + C_S_v * transpose(dGamma_dTe);
            % =============================================================
            % attention: Ran add these lines to update dGamma_dtau
            % =============================================================
            Q66 = obj.stiffness(dGamma, dGamma_dtau, J_n);
            
        end
        
        function Q66 = stiffness(obj, dGamma, dGamma_dtau, K66)
            % input : dGamma_dTau, dGamma, K99
            % output: W
            % class property: dF, Fe, C, S0, Te
            
            % useful matrix
            idx = [1,2,3,1,1,2; 1,2,3,2,3,3];
            m6_m6 = diag([1,1,1,2,2,2]);
            
            % prepare intermediate variables
            nslip = obj.num_slip_system;
            dF = obj.F_n1 / obj.F_n;
            Fe_n1 = obj.F_n1 / obj.Fp_n1;
            Fe_n = obj.F_n / obj.Fp_n;
            [~, Ce66] = obj.elasticity(zeros(3));
            [dR, dU] = poldec(dF);
            
            % step (1): calculate L(9*9)
            L66 = zeros(6);
            UFe = dU * Fe_n;
            for aa = 1:6
                ii = idx(1,aa);
                jj = idx(2,aa);
                for bb = 1:6
                    kk = idx(1,bb);
                    ll = idx(2,bb);
%                     L99(aa,bb) = Fe_n1(kk,ii)*UFe(ll,jj) + UFe(kk,ii)*Fe_n1(ll,jj);
                    L66(aa,bb) = 0.5 * ( Fe_n(ll,ii)*UFe(kk,jj) + ...
                        UFe(kk,ii)*Fe_n(ll,jj) + Fe_n(kk,ii)*UFe(ll,jj) ...
                        + UFe(ll,ii)*Fe_n(kk,jj) );
                end
            end
            
            % step (2): map Ce from lattice coordinate to global coordinate
            % Ce(6*6) already in intermediate configuration
            
            % step (3): calculate D(9*9)
            D66 = 0.5 * Ce66 * m6_m6 * L66;
            
            % step (4): calculate G(9*9*nslip) and J(9*9*nslip)
            % step (5): calculate B(9*nslip)
            B6 = zeros(6, nslip);
            for alpha = 1 : nslip
                S0 = obj.P_Schmid_full(:,:,alpha);
                G66 = zeros(6);
                B33 = 0.5 * dGamma_dtau(alpha) * ( S0 + transpose(S0) );
                for bb = 1:6
                    L33 = obj.v2s( L66(:,bb) );
                    G33 = L33 * S0 + transpose(S0) * L33;
                    G66(:,bb) = obj.s2v( G33 );
                end
                B6(:,alpha) = obj.s2v(B33);
                J66 = 0.5 * Ce66 * m6_m6 * G66;
                % compute D = D - dGamma * J
                D66 = D66 - dGamma(alpha) * J66;
            end
            
            % step (6): calculate K(9*9) and Q(9*9)
            % ATTENTION: K66 is not the actuall stiffness
            %            K66(:,4:6) is multiplied by two
            %            This is because it will be multiplied
            %            by stress instead of strain!
            Q66 = K66 \ D66;
            
            % step (7): calculate S(9*9)
            R6 = transpose(Q66) * m6_m6 * B6;
            
% Q66 = [50306.0334501981,27103.6948298331,32746.1383619720,-1237.33334476128,-0.0261060506545618,-0.0388503451631550;
% 27178.9353504914,54073.9901879306,28782.3959823186,1415.66663600656,0.0285957924006652,0.0209844586152030;
% 32758.9150278973,28741.5739881114,48608.6233723482,-302.402436846982,-0.00310380610102357,0.0172178049595573;
% -1097.49705501372,1335.59366766178,-280.896531990038,16252.8242242388,0.00859827764543297,0.000613553652328847;
% -6.69547236300250e-08,-6.43439322898849e-08,-6.52834564704117e-08,-6.74503701400462e-08,16367.3386978482,467.030699351862;
% -3.98020276997997e-08,-4.34359746270834e-08,-4.63186702440687e-08,-3.82046671738894e-08,525.650903725398,15742.6693294381];
% R6 = [-0.295184151686450,0.0792060808569662,0.215255447411674,-0.0658309657396790,0.110913585375805,-0.192309977996654;
% 0.0417774636316903,0.0211411936845978,-0.0626739678591290,-0.0843553361700827,-0.0547479140212119,-0.0313829378764230;
% 2.40271757294796e-09,-2.79048323371876e-09,3.84476213971461e-10,6.75058159759173e-09,1.68766038845501e-09,5.45313308885858e-09;
% -0.295184151686450,0.0792060808569662,0.215255447411674,-0.0658309657396790,-0.110913974036534,0.192307931911999;
% -0.0417774636316903,-0.0211411936845978,0.0626739678591290,0.0843553361700827,-0.0547492024575969,-0.0313831022912778;
% 2.40271757294796e-09,-2.79048323371876e-09,3.84476213971461e-10,6.75058159759173e-09,-1.68753708873317e-09,-5.45179267795128e-09;
% -0.00131078008182002,-0.0270179036978232,0.0284657468968638,-0.0342534523080269,-0.0114933834137790,0.0199282862262763;
% -0.0558260046710731,0.188412747638536,-0.133090975938758,-0.0244615120727784,-0.0929632995190217,-0.0532882683617307;
% 2.40271757294796e-09,-2.79048323371876e-09,3.84476213971461e-10,6.75058159759173e-09,5.57641743209501e-09,-1.28895048237190e-09;
% -0.00131078008182002,-0.0270179036978232,0.0284657468968638,-0.0342534523080269,0.0114934774060606,-0.0199277099741769;
% -0.0558260046710731,0.188412747638536,-0.133090975938758,-0.0244615120727784,0.0929615268132331,0.0532877355942699;
% 2.40271757294796e-09,-2.79048323371876e-09,3.84476213971461e-10,6.75058159759173e-09,-5.57500683048501e-09,1.28901146368267e-09];
% R6 = transpose(R6);
            
            % step (7): calculate S(9*9)
            Gamma_S0 = zeros(3);
            for alpha = 1 : nslip
                Gamma_S0 = Gamma_S0 + dGamma(alpha) * obj.P_Schmid_full(:,:,alpha);
            end
            FeGammaS0 = Fe_n * ( eye(3) - Gamma_S0 );
            S96 = zeros(9,6);
            for aa = 1:9
                ii = floor( (aa-1) / 3 ) + 1;
                jj = mod( aa-1, 3 )+1;
                for bb = 1:6
                    kk = idx( 1, bb );
                    ll = idx( 2, bb );
                    % S99(aa,bb) = dR(ii,kk) * FeGammaS0(ll,jj);
                    S96(aa,bb) = 0.5 * ( dR(ii,kk) * FeGammaS0(ll,jj) +...
                        dR(ii,ll) * FeGammaS0(kk,jj) );
                end
            end
            RUF = dF * Fe_n;
            RUFS99 = zeros(9,nslip);
            for alpha = 1 : nslip
                RUFS33 = RUF * obj.P_Schmid_full(:,:,alpha);
                RUFS99(:,alpha) = reshape(transpose(RUFS33),9,1);
            end
            S96 = S96 - RUFS99 * transpose(R6);
            
            % step (8): calculate W(9*9)
            detF = det(Fe_n1);
            Te = obj.Te_n1;
            FTF33 = Fe_n1 * Te * transpose(Fe_n1);
            Fe_inv_9 = reshape( inv(Fe_n1), 1, 9 ); % transpose of inverse
            W66 = zeros(6);
            for bb = 1:6
                S33 = transpose( reshape( S96(:,bb), 3, 3 ) );
                Q33 = obj.v2s( Q66(:,bb) );
                SF  = Fe_inv_9 * S96(:,bb);
                W33 = S33 * Te * transpose(Fe_n1) + Fe_n1 * Q33 * transpose(Fe_n1) ...
                    + Fe_n1 * Te * transpose(S33) - SF * FTF33;
                W66(:,bb) = obj.s2v( W33 );
            end
            
            % reduce Wijkl to WIJ
            obj.C_stiff = W66 / detF;
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
        
    end
    
    methods (Static)
        
        function t = v2e(v)
            v(4:6) = 0.5 * v(4:6);
            t = [v(1), v(4), v(5);
                 v(4), v(2), v(6);
                 v(5), v(6), v(3)];
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