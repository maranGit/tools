% Ran Ma
% 11/29/2019
% multi-phase-field crystal plasticity
% Voigt notation: 11, 22, 33, 12, 13, 23
%
% Na, S., & Sun, W. (2018). Computational thermomechanics of crystalline rock, 
% Part I: A combined multi-phase-field/crystal plasticity approach for single
% crystal simulations. Computer Methods in Applied Mechanics and Engineering,
% 338, 657-691.
%
% update 12/6/2019: reduce to 1 phase field for fracture
% update 12/9/2019: update to my phase field fracture formulation
% update 12/9/2019: twinning phase field
%
classdef PhaseFieldMultiCrystal < handle
    properties
        % Input parameters
        K;
        nu;
        tau_ini;
        hardening;
        prm_Q;    % Activation energy
        prm_C;    % Fitting parameter
        exponent; % Power-law parameter
%         Euler_Theta;
%         Euler_Phi;
        Euler_phi1;
        Euler_phi;
        Euler_phi2;
        verbose = 0; % bool
        
        % Computed Parameters
        %   Elastic moduli
        lambda;
        mu;
        youngs;
        viscosity;
        
        %   Temperature, time increment
        temperature;
        dt;
        
        % Stresses, Strains, Moduli
        new_stress; % SymmetricTensor<2,3>
        old_stress; % SymmetricTensor<2,3>
        
        new_strain; % SymmetricTensor<2,3>
        old_strain; % SymmetricTensor<2,3>
        
        new_elastic_strain; % SymmetricTensor<2,3>
        old_elastic_strain; % SymmetricTensor<2,3>
        old_old_elastic_strain; % SymmetricTensor<2,3>
        
        new_plastic_strain; % SymmetricTensor<2,3>
        old_plastic_strain; % SymmetricTensor<2,3>
        
        new_equiv_plastic_strain;
        old_equiv_plastic_strain;
        
        new_cto; % SymmetricTensor<4,3>
        elastic_cto; % SymmetricTensor<4,3>
        
        % Internal Variables
        new_plastic_slip;
        old_plastic_slip;
        old_old_plastic_slip;
        
        new_plastic_slip_i;
        old_plastic_slip_i;
        old_old_plastic_slip_i;
        
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
        
        new_active_set; % std::vector<unsigned>
        old_active_set; % std::vector<unsigned>
        
        % Yield stress
        new_tau;
        old_tau;
        
        % phase field fracture
        Gc_c;
        l0_c;
        k_c;
        Hn_c;
        Hn1_c;
        g_c;
        gp_c;
        gpp_c;
        
        % phase field twinning
        alpha_t; % constant in representative function
        gamma0_t; % magnitude of twinning shear
        Gc_t; % twin boundary energy
        l0_t; % equilibrium boundary thickness
        A_t; % double-well energy parameter
        k_t; % gradient energy parameter
        phi_t; % representative function
        phip_t; % 1st order derivative of representative function
        phipp_t; % 2nd order derivative of representative function
        sm_t; % twinning strain mode
        strain_twin; % twinning strain
        tau_t; % shear stress on the twinning plane
        svec_t = [1; 0; 0]; % shear direction of twinning plane
        mvec_t = [0; 1; 0]; % normal direction of twinning plane
        Q_t; % rotation matrix of twinning deformation
        
        % constants
        tol = 1.0e-9;
    end
    
    methods
        % constructor without default input value
        function obj = PhaseFieldMultiCrystal(s_type)
            obj.slip_type = s_type;
            if s_type == 0
                obj.num_slip_system = 6;
                obj.slip_normal = [ 1.0, 1.0, 0.0;
                    1.0,-1.0, 0.0;
                    1.0, 0.0, 1.0;
                    1.0, 0.0,-1.0;
                    0.0, 1.0, 1.0;
                    0.0, 1.0,-1.0];
                obj.slip_direct = [1.0,-1.0, 0.0;
                    1.0, 1.0, 0.0;
                    1.0, 0.0,-1.0;
                    1.0, 0.0, 1.0;
                    0.0, 1.0,-1.0;
                    0.0, 1.0, 1.0];
            elseif s_type == 1
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
            elseif s_type == 2
                error('>>> Error: HMX slip system is not implemented yet!');
            else
                error('>>> Error: Unknown slip type!');
            end
        end
        
        function initialize(obj)
            % Parse Parameters
            obj.K = 25; % bulk modulus
            obj.nu = 0.25; % poisson ratio
            obj.tau_ini = 5.0e-3; % critical resolved shear stress
            obj.hardening = 0.1; % hardening parameter
            obj.prm_Q = 14; % activation energy
            obj.prm_C = 1.0; % shape factor
            obj.exponent = 10.0;
%             obj.Euler_Theta = 10.0 * pi / 180;
%             obj.Euler_Phi = 30.0 * pi / 180;
            obj.Euler_phi1 = 30.0 * pi / 180;
            obj.Euler_phi  = 50.0 * pi / 180;
            obj.Euler_phi2 = 10.0 * pi / 180;
            obj.verbose = 0;
            
            % Elastic Tangent Operator
            obj.youngs = 3 * obj.K * (1 - 2*obj.nu);
            obj.mu = 0.5 * obj.youngs / (1 + obj.nu);
            obj.lambda = obj.K - 2*obj.mu/3;
            I2 = [1; 1; 1; 0; 0; 0];
            I2I2 = I2 * transpose(I2);
            I4 = eye(6);
            I4(4,4) = 0.5; I4(5,5) = 0.5; I4(6,6) = 0.5;
            obj.elastic_cto = obj.lambda*I2I2 + 2*obj.mu*I4;
            
            % Plastic Internal Variable / Strain
            obj.new_plastic_strain       = zeros(3,1);
            obj.new_equiv_plastic_strain = 0;
            obj.new_plastic_slip         = 0;
            obj.new_plastic_slip_i       = zeros(obj.num_slip_system, 1);
            obj.old_plastic_slip_i       = zeros(obj.num_slip_system, 1);
            obj.old_old_plastic_slip_i   = zeros(obj.num_slip_system, 1);
            obj.new_tau                  = obj.tau_ini;
            obj.new_active_set = zeros(2 * obj.num_slip_system + 1, 1);
            
            % phase field fracture
            obj.Gc_c  = 1.15e-6;
            obj.l0_c  = 0.02;
            obj.k_c   = 1e-3;
            obj.Hn_c  = 0;
            obj.Hn1_c = 0;
            obj.g_c   = 0;
            obj.gp_c  = 0;
            obj.gpp_c = 0;
            
            % phase field twinning
            obj.Gc_t        = 1.17e-4; % mJ/(mm^2)
            obj.l0_t        = 1.0e-6;  % mm
            obj.alpha_t     = 3;
            obj.gamma0_t    = 0.1295;
            obj.A_t         = 12 * obj.Gc_t / obj.l0_t;
            obj.k_t         = 0.75 * obj.Gc_t * obj.l0_t;
            obj.phi_t       = 0;
            obj.phip_t      = 0;
            obj.phipp_t     = 0;
            obj.sm_t        = zeros(6,1);
            obj.strain_twin = zeros(6,1);
            
            % extra initiation for tensor variables
            obj.new_stress = zeros(6,1);
            obj.old_stress = zeros(6,1);
            obj.new_strain = zeros(6,1);
            obj.old_strain = zeros(6,1);
            obj.new_elastic_strain = zeros(6,1);
            obj.old_elastic_strain = zeros(6,1);
            obj.old_old_elastic_strain = zeros(6,1);
            obj.new_plastic_strain = zeros(6,1);
            obj.old_plastic_strain = zeros(6,1);
            obj.new_cto = obj.elastic_cto;
            
            obj.update_slip_plane();
            obj.save_state();
        end
        
        function set_stress(obj, in_situ_stress)
            obj.new_stress = in_situ_stress;
            obj.new_elastic_strain = obj.elastic_cto \ in_situ_stress;
            obj.save_state();
        end
        
        function update(obj, incr_strain)
            if length(incr_strain) == 3
                incr_strain_2d = zeros(6,1);
                incr_strain_2d(1) = incr_strain(1);
                incr_strain_2d(2) = incr_strain(2);
                incr_strain_2d(4) = incr_strain(3);
                strain_n1 = obj.old_strain + incr_strain_2d;
                obj.run_3D_update(strain_n1);
            else
                strain_n1 = obj.old_strain + incr_strain;
                obj.run_3D_update(strain_n1);
            end
        end
        
        function update_parameters(obj, pass_temp, pass_time)
            obj.temperature = pass_temp + 273.0;
            obj.dt = pass_time;
        end
        
        function update_phasefield(obj, pft, pfc)
            % fracture
            one       = 1.0e0;
            two       = 2.0e0;
            three     = 3.0e0;
            four      = 4.0e0;
            six       = 6.0e0;
            twelve    = 12.0e0;
            onemk     = 1 - obj.k_c;
            obj.g_c   = onemk * (1-pfc) * (1-pfc) + obj.k_c;
            obj.g_c   = onemk * (one - pfc) * (one - pfc) + obj.k_c;
            obj.gp_c  = two * onemk * (pfc - one);
            obj.gpp_c = two * onemk;
            
            % twinning
            a1 = obj.alpha_t;
            a2 = two - a1;
            a3 = a1 - three;
            pft2 = pft * pft;
            pft3 = pft * pft2;
            pft4 = pft * pft3;
            obj.phi_t   = a1*pft2 + two*a2*pft3 + a3*pft4;
            obj.phip_t  = two*a1*pft + six*a2*pft2 + four*a3*pft3;
            obj.phipp_t = two*a1 + twelve*a2*pft + twelve*a3*pft2;
            obj.strain_twin = obj.gamma0_t * obj.phi_t * obj.sm_t;
        end
        
        function update_slip_plane(obj)
            rotation_matrix = obj.update_R_matrix();
            % update twinning plane
            svec = obj.svec_t;
            mvec = obj.mvec_t;
            svec = rotation_matrix * svec;
            mvec = rotation_matrix * mvec;
            sm   = zeros(6,1);
            sm(1:3) = svec .* mvec;
            sm(4)   = svec(1)*mvec(2) + svec(2)*mvec(1);
            sm(5)   = svec(1)*mvec(3) + svec(3)*mvec(1);
            sm(6)   = svec(2)*mvec(3) + svec(3)*mvec(2);
            obj.sm_t = sm;
            % proper twinning rotation matrix
            % Clayton, Knap, 2011
            % rotate around mvec for 180 degree
            % -Q*v is the mirror of v w.r.t. the mvec plane
            obj.Q_t  = 2 * mvec * transpose(mvec) - eye(3);
        end
        
        function update_strain_energy(obj)
            
            I2 = [1; 1; 1; 0; 0; 0];
            I4 = diag([1, 1, 1, 0.5, 0.5, 0.5]);
            
            %----------------------------------------------------%
            % Decomposition                                      %
            %----------------------------------------------------%
            
            one_third = 1/3;
            % Construct stress from tensile stress + compressive stress
            tr_eps = obj.new_elastic_strain(1) + obj.new_elastic_strain(2) ...
                + obj.new_elastic_strain(3);
            tr_eps_plus  = max(0.0, tr_eps);
            tr_eps_minus = tr_eps - tr_eps_plus;
            
            dev_starin = obj.new_elastic_strain - one_third*tr_eps*I2;
            
            new_stress_plus = obj.K* tr_eps_plus*I2 + 2*obj.mu*I4*dev_starin;
            new_stress_minus = obj.K* tr_eps_minus*I2;
            
            obj.new_stress = obj.g_c * new_stress_plus + new_stress_minus;
            
            % shear stress on twinning plane
            obj.tau_t = transpose(obj.new_stress) * obj.sm_t;
            
            % Update strain history variable
            strain_energy_plus = 0.5*obj.K*tr_eps_plus*tr_eps_plus...
                + obj.mu*transpose(dev_starin)*I4*dev_starin;
            
            plastic_energy = 0.5*obj.hardening*(...
                transpose(obj.new_plastic_slip_i) * obj.new_plastic_slip_i);
            
            % Update H plus for multiphase
            strain_energy = strain_energy_plus+plastic_energy;
            if strain_energy > obj.Hn_c
                obj.Hn1_c = strain_energy;
            else
                obj.Hn1_c = obj.Hn_c;
            end
            
            obj.new_cto = obj.g_c * obj.new_cto;
        end
        
        function save_state(obj)
            % Stress
            obj.old_stress               = obj.new_stress;
            % Elastic strain
            obj.old_old_elastic_strain   = obj.old_elastic_strain;
            obj.old_elastic_strain       = obj.new_elastic_strain;
            % Others
            obj.old_plastic_strain       = obj.new_plastic_strain;
            obj.old_strain               = obj.new_strain;
            obj.old_equiv_plastic_strain = obj.new_equiv_plastic_strain;
            obj.old_tau                  = obj.new_tau;
            obj.old_active_set           = obj.new_active_set;
            % Plastic slip
            obj.old_old_plastic_slip     = obj.old_plastic_slip;
            obj.old_plastic_slip         = obj.new_plastic_slip;
            obj.old_old_plastic_slip_i   = obj.old_plastic_slip_i;
            obj.old_plastic_slip_i       = obj.new_plastic_slip_i;
            % Driving force
            obj.Hn_c                     = obj.Hn1_c;
        end
        
        function recover_state(obj)
            % Stress
            obj.new_stress               = obj.old_stress;
            % Elastic strain
            obj.new_elastic_strain       = obj.old_elastic_strain;
            obj.old_elastic_strain       = obj.old_old_elastic_strain;
            % Others
            obj.new_plastic_strain       = obj.old_plastic_strain;
            obj.new_strain               = obj.old_strain;
            obj.new_equiv_plastic_strain = obj.old_equiv_plastic_strain;
            obj.new_tau                  = obj.old_tau;
            obj.new_active_set           = obj.old_active_set;
            % Plastic slip
            obj.new_plastic_slip         = obj.old_plastic_slip;
            obj.old_plastic_slip         = obj.old_old_plastic_slip;
            obj.new_plastic_slip_i       = obj.old_plastic_slip_i;
            obj.old_plastic_slip_i       = obj.old_old_plastic_slip_i;
            % Driving force
            obj.Hn1_c                    = obj.Hn_c;
        end
        
        function sigma = stress(obj)
            sigma = obj.new_stress;
        end
        
        function sigma = full_stress(obj)
            sigma = obj.new_stress;
        end
        
        function eps = strain(obj)
            eps = obj.new_strain;
        end
        
        function eps = full_strain(obj)
            eps = obj.new_strain;
        end
        
        function vlms = volumetric_strain(obj)
            vlms = obj.new_plastic_strain(1) ...
                + obj.new_plastic_strain(2) ...
                + obj.new_plastic_strain(3);
        end
        
        function dvts = deviatoric_strain(obj)
            one_third = 1/3;
            I2 = [1; 1; 1; 0; 0; 0];
            vlms = obj.new_plastic_strain(1) ...
                + obj.new_plastic_strain(2) ...
                + obj.new_plastic_strain(3);
            dev_strain = obj.new_plastic_strain ...
                - one_third * vlms * I2;
            dvts = dev_strain(1)*dev_strain(1) ...
                + dev_strain(2)*dev_strain(2) ...
                + dev_strain(3)*dev_strain(3) ...
                + 0.5 * dev_strain(4)*dev_strain(4) ...
                + 0.5 * dev_strain(5)*dev_strain(5) ...
                + 0.5 * dev_strain(6)*dev_strain(6);
        end
        
        function eqpls = equiv_plastic_strain(obj)
            eqpls = obj.new_equiv_plastic_strain;
        end
        
        function stiff = cto(obj)
            stiff = obj.new_cto;
        end
        
        function bkm = bulk_mod(obj)
            bkm = obj.K;
        end
        
        function shm = shear_mod(obj)
            shm = obj.mu;
        end
        
        function psn = plastic_slip_new(obj)
            psn = obj.new_plastic_slip;
        end
        
        function plsl = plastic_slip_old(obj)
            plsl = obj.old_plastic_slip;
        end
        
        function plsl = plastic_slip_new_0(obj)
            plsl = obj.new_plastic_slip_i(1);
        end
        
        function plsl = plastic_slip_new_1(obj)
            plsl = obj.new_plastic_slip_i(2);
        end
        
        function plsl = plastic_slip_new_2(obj)
            plsl = obj.new_plastic_slip_i(3);
        end
        
        function plsl = plastic_slip_new_3(obj)
            plsl = obj.new_plastic_slip_i(4);
        end
        
        function plsl = plastic_slip_new_4(obj)
            plsl = obj.new_plastic_slip_i(5);
        end
        
        function plsl = plastic_slip_new_5(obj)
            plsl = obj.new_plastic_slip_i(6);
        end
        
        % Energy dissipation
        function temp = old_plastic_dissipation(obj)
            old_plastic_slip_rate = obj.old_plastic_slip - obj.old_old_plastic_slip;
            if (old_plastic_slip_rate  < 0)
                old_plastic_slip_rate = 0.0;
            end
            temp = obj.old_tau * old_plastic_slip_rate;
        end
        
        function temp = old_hardening_dissipation(obj)
            old_plastic_slip_rate_i = obj.old_plastic_slip_i - obj.old_old_plastic_slip_i;
            old_hds_i = obj.hardening * obj.old_plastic_slip_i * old_plastic_slip_rate_i;
            tmp_hds = sum(old_hds_i);
            temp = -tmp_hds;
        end
        
        function temp = old_elastic_strain_rate(obj)
            tmp_elastic_strain = obj.old_elastic_strain - obj.old_old_elastic_strain;
            temp = tmp_elastic_strain(1) + tmp_elastic_strain(2) + tmp_elastic_strain(3);
        end
        
        % Update
        function run_3D_update(obj,strain_n1)
            % Compute a trial state in which the increment is assumed to be fully elastic
%             obj.new_elastic_strain       = obj.old_elastic_strain + incr_strain;
            obj.new_elastic_strain       = strain_n1 - obj.old_plastic_strain - obj.strain_twin;
            obj.new_cto                  = obj.elastic_cto;
            obj.new_stress               = obj.new_cto*obj.new_elastic_strain;
            obj.new_plastic_strain       = obj.old_plastic_strain;
            obj.new_equiv_plastic_strain = obj.old_equiv_plastic_strain;
            obj.new_plastic_slip         = obj.old_plastic_slip;
            obj.new_plastic_slip_i       = obj.old_plastic_slip_i;
            obj.new_tau                  = obj.old_tau;
            obj.new_active_set           = obj.old_active_set;
            
            prm_R    = 1.986e-3; % Gas constant (kcal/(mol.K))
            obj.viscosity = 1.0/(obj.prm_C*exp((-1.0)*obj.prm_Q/(prm_R*obj.temperature)));
            
            % Assign crystal slip system
            nslip = obj.num_slip_system;
            n20   = 2 * nslip;
            n21   = n20 + 1;
            P_Schmid = zeros(6,n20);
            rot_matrix = obj.update_R_matrix();
            for ii = 1:nslip
                m_normal = transpose(obj.slip_normal(ii,1:3));
                s_direct = transpose(obj.slip_direct(ii,1:3));
                % Normalize
                m_normal = m_normal / norm(m_normal);
                s_direct = s_direct / norm(s_direct);
                
                % Assemble Schmid tensor alpha
                P_Schmid_tmp = 0.5*rot_matrix*(m_normal*s_direct' + s_direct*m_normal')*transpose(rot_matrix);
                P_Schmid_twn = obj.Q_t * P_Schmid_tmp * transpose( obj.Q_t );
                P_Schmid_tpt = (1-obj.phi_t)*P_Schmid_tmp + obj.phi_t*P_Schmid_twn;
                
                P_Schmid(1,ii) = P_Schmid_tpt(1,1);
                P_Schmid(2,ii) = P_Schmid_tpt(2,2);
                P_Schmid(3,ii) = P_Schmid_tpt(3,3);
                P_Schmid(4,ii) = P_Schmid_tpt(1,2)+P_Schmid_tpt(2,1);
                P_Schmid(5,ii) = P_Schmid_tpt(1,3)+P_Schmid_tpt(3,1);
                P_Schmid(6,ii) = P_Schmid_tpt(2,3)+P_Schmid_tpt(3,2);
                
                P_Schmid(1:6,ii+nslip) = -P_Schmid(1:6,ii);
            end
            
            % Identify potentially active systems (F > 0)
            trial_active_set = zeros(n21, 1);
            for ii=1:n20
                trial_F = transpose( P_Schmid(:,ii) ) * obj.new_stress - obj.new_tau;
                if ( trial_F > obj.tol)
                    trial_active_set(n21) = trial_active_set(n21) + 1;
                    temp = trial_active_set(n21);
                    trial_active_set(temp) = ii;
                end
            end
            
            % ---------------------------------------------------------
            % If size of trial_active_Set = 0, then Elastic response %
            %                                                        %
            if (trial_active_set(n21) < 1)                          %
%                 obj.new_strain = obj.new_elastic_strain + obj.new_plastic_strain;   %
                obj.new_strain = strain_n1;   %
                %new_active_set.resize(0);                             %
                obj.update_strain_energy();                                 %
                return;                                                 %
            end                                                         %
            %                                                        %
            % ---------------------------------------------------------
            % Else Plastic
            
            if obj.verbose
                if obj.old_active_set(n21) < 1
                    fprintf("No initial active set!\n");
                else
                    fprintf("Initial active set: \n");
                    format = repmat('%i ', 1, 2*nslip);
                    format = strcat(format, '\n');
                    fprintf(format, obj.old_active_set);
                end
            end
            
            % Define variables
            active_set = zeros(n21,1);
            % Gamma_min  = 0;
            update_II  = 1;
            
            tmp_elastic_strain = zeros(6,1);
            tmp_plastic_strain = zeros(6,1);
            
            % Active set iteration
            for S_iter = 0:14
                if S_iter == 0
                    active_set = obj.old_active_set;
                end
                % Check the size of J_active set
                update_I = active_set(n21) >= 1;
                
                % Start while loop
                while update_I
                    % Initial values for plastic slip
                    Gamma = zeros(active_set(n21), 1);
                    jacobian = zeros(active_set(n21),active_set(n21));
                    jacobian_inv = zeros(active_set(n21),active_set(n21));
                    residual = zeros(active_set(n21), 1);
                    % xdelta = zeros(active_set(n21), 1);
                    % hardening_moduli = zeros(active_set(n21),active_set(n21));
                    
                    tmp_elastic_strain = zeros(6,1);
                    
                    resid = norm(residual);
                    resid_ini = 1 + resid;
                    
                    % Newton iteration
                    for N_iter = 0:19
                        tmp_plastic_strain = zeros(6,1);
                        tmp_plastic_slip   = 0;
                        tmp_plastic_slip_i = zeros(obj.num_slip_system, 1);
                        % Temporary plastic strain
                        for ii = 1:active_set(n21)
                            alpha = active_set(ii);
                            tmp_plastic_strain = tmp_plastic_strain + Gamma(ii)*P_Schmid(:,alpha);
                            tmp_plastic_slip   = tmp_plastic_slip + Gamma(ii);
                            curr_slip = mod(alpha-1, obj.num_slip_system) + 1;
                            tmp_plastic_slip_i(curr_slip) = tmp_plastic_slip_i(curr_slip) + Gamma(ii);
                        end

                        tmp_elastic_strain = obj.new_elastic_strain - tmp_plastic_strain;

                        obj.new_plastic_slip   = obj.old_plastic_slip   + tmp_plastic_slip;
                        obj.new_plastic_slip_i = obj.old_plastic_slip_i + tmp_plastic_slip_i;

                        obj.new_stress         = obj.elastic_cto * tmp_elastic_strain;

                        % Isotropic hardening
                        obj.new_tau = obj.old_tau + obj.hardening*tmp_plastic_slip;
                        hardening_moduli = repmat(obj.hardening, active_set(n21), active_set(n21));

                        % Residual
                        for ii = 1 : active_set(n21)
                            if Gamma(ii) < 0
                                continue;
                            end
                            alpha = active_set(ii);
                            residual(ii) = transpose(obj.new_stress) * P_Schmid(:,alpha) ...
                                - obj.new_tau*( (obj.viscosity/obj.dt*Gamma(ii)+1)^(1/obj.exponent) );
                        end

                        resid = norm(residual);
                        if N_iter == 0
                            resid_ini = 1 + resid;
                        end

                        if resid/resid_ini < obj.tol || resid < obj.tol
                            obj.new_active_set = active_set;
                            break;
                        end

                        % Jacobian matrix for local return mapping
                        for ii = 1 : active_set(n21)
                            if Gamma(ii) < 0
                                continue;
                            end
                            for jj = 1 : active_set(n21)
                                % find the No. of slip system in active set
                                alpha = active_set(ii);
                                beta = active_set(jj);
                                tmp_beta = obj.elastic_cto*P_Schmid(:,beta);

                                jacobian(ii,jj) = transpose(P_Schmid(:,alpha)) * tmp_beta ...
                                    + hardening_moduli(ii,jj)...
                                    *((obj.viscosity/obj.dt*Gamma(ii)+1)^(1/obj.exponent));
                                if (ii == jj)
                                    jacobian(ii,jj) = jacobian(ii,jj) ...
                                        + obj.new_tau*(obj.viscosity/(obj.exponent*obj.dt))...
                                        *(((obj.viscosity/obj.dt)*Gamma(ii)+1)^(1/obj.exponent-1));
                                end

                            end
                        end

                        % Inverse of Jacobian
                        [U,S,V] = svd(jacobian);
                        % Inverse Method I - Pseudo inverse
                        jacobian_inv = obj.update_jacobian_inverse(S,U,V);
                        % Newton update: xnew = xcurrent - xdelta, xdelta = jacobian_inv*residual
                        xdelta = jacobian_inv * residual;
                        Gamma = Gamma + xdelta;
                        
                    end % end of Newton return mapping
                    
                    % Update Jacobian inverse
                    CTO_jacobian_inv = jacobian_inv;
                    
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
                    check_F = transpose(P_Schmid) * obj.new_stress - obj.new_tau;
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
            
            % Update CTO
            obj.new_cto = obj.elastic_cto;
            n_active = obj.new_active_set(n21);
            for ii = 1:n_active
                for jj = 1:n_active
                    alpha = obj.new_active_set(ii);
                    beta  = obj.new_active_set(jj);
                    tmp_Schmid_a = obj.elastic_cto * P_Schmid(:,alpha);
                    tmp_Schmid_b = obj.elastic_cto * P_Schmid(:,beta);
                    obj.new_cto = obj.new_cto - CTO_jacobian_inv(ii,jj) * tmp_Schmid_a * transpose(tmp_Schmid_b);
                end
            end
            
            % Update variables
            obj.new_elastic_strain = tmp_elastic_strain;
            obj.new_plastic_strain = obj.old_plastic_strain + tmp_plastic_strain;
            tmp_plastic_strain(4:6) = tmp_plastic_strain(4:6) * sqrt(0.5);
            obj.new_equiv_plastic_strain = obj.old_equiv_plastic_strain ...
                + sqrt( 2 / 3 ) * norm(tmp_plastic_strain);
%             obj.new_strain = obj.new_elastic_strain + obj.new_plastic_strain;
            obj.new_strain = strain_n1;
            obj.update_strain_energy();
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
    
        % A whole bunch of 'get' function
        
        function fren = get_Gc_c(obj)
            fren = obj.Gc_c;
        end
        
        function lgsc = get_l0_c(obj)
            lgsc = obj.l0_c;
        end
        
        function tmp = get_A_t(obj)
            tmp = obj.A_t;
        end
        
        function tmp = get_k_t(obj) % k = 3 * Gamma * l / 4
            tmp = obj.k_t;
        end
        
        function tmp = get_tau(obj) % return resolved shear stress
            tmp = obj.tau_t;
        end
        
        function tmp = get_gamma0_t(obj)
            tmp = obj.gamma0_t;
        end
        
        function tmp = get_phip_t(obj)
            tmp = obj.phip_t;
        end
        
        function tmp = get_phipp_t(obj)
            tmp = obj.phipp_t;
        end
        
        function tmp = get_g_c(obj)
            tmp = obj.g_c;
        end
        
        function tmp = get_gp_c(obj)
            tmp = obj.gp_c;
        end
        
        function tmp = get_gpp_c(obj)
            tmp = obj.gpp_c;
        end
        
        function tmp = get_dtau_dpft(obj)
            tmp = obj.tol;
        end
        
        function tmp = get_dtau_dpfc(obj)
            tmp = obj.tol;
        end
        
        function tmp = get_dH_dpft(obj)
            tmp = obj.tol;
        end
        
        function tmp = get_Hn1_c(obj)
            tmp = obj.Hn1_c;
        end
        
        function tmp = get_dH_dpfc(obj)
            tmp = obj.tol;
        end
        
        function tmp = get_dH_deps(obj)
            tmp = obj.tol;
        end
        
        function tmp = get_dtau_deps(obj)
            tmp = obj.tol;
        end
        
        function tmp = get_dsigma_dpft(obj)
            tmp = obj.tol;
        end
        
        function tmp = get_dsigma_dpfc(obj)
            tmp = obj.tol;
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
    end
end