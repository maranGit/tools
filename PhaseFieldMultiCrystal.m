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
% update 12/15/2019: remove yield surface, improve robustness
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
        elastic_cto_damage; % SymmetricTensor<4,3>
        
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
        % tangent stiffness
        dtau_dpft;
        dtau_dpfc;
        dH_dpft;
        dH_dpfc;
        dH_deps;
        dtau_deps;
        dsigma_dpft;
        dsigma_dpfc;
        
        % slip system
        P_Schmid_notw;
        P_Schmid_twin;
        
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
            rot_matrix = obj.update_R_matrix();
            
            % update twinning plane
            svec = obj.svec_t / norm(obj.svec_t);
            mvec = obj.mvec_t / norm(obj.mvec_t);
            svec = rot_matrix * svec;
            mvec = rot_matrix * mvec;
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
            
            % initialize slip plane from crystal to spatial
            nslip = obj.num_slip_system;
            P_Schmid = zeros(6,nslip);
            for ii = 1:nslip
                m_normal = transpose(obj.slip_normal(ii,1:3));
                s_direct = transpose(obj.slip_direct(ii,1:3));
                % Normalize
                m_normal = m_normal / norm(m_normal);
                s_direct = s_direct / norm(s_direct);
                % Assemble Schmid tensor alpha
                P_Schmid_tmp = 0.5*rot_matrix*(m_normal*s_direct' + s_direct*m_normal')*transpose(rot_matrix);
                
                P_Schmid(1,ii) = P_Schmid_tmp(1,1);
                P_Schmid(2,ii) = P_Schmid_tmp(2,2);
                P_Schmid(3,ii) = P_Schmid_tmp(3,3);
                P_Schmid(4,ii) = P_Schmid_tmp(1,2)+P_Schmid_tmp(2,1);
                P_Schmid(5,ii) = P_Schmid_tmp(1,3)+P_Schmid_tmp(3,1);
                P_Schmid(6,ii) = P_Schmid_tmp(2,3)+P_Schmid_tmp(3,2);
            end
            obj.P_Schmid_notw = P_Schmid;
            
            % initialize slip system in twinning region
            Q66 = obj.mm10_RT2RVE(obj.Q_t);
            obj.P_Schmid_twin = Q66 * P_Schmid;
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
            
%             obj.new_cto = obj.g_c * obj.new_cto;
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
            eps_tr                       = strain_n1(1) + strain_n1(2) + strain_n1(3);
            I2                           = [1; 1; 1; 0; 0; 0];
            I4                           = diag([1, 1, 1, 0.5, 0.5, 0.5]);
            I2I2                         = I2 * transpose(I2);
            if eps_tr >= 0
                obj.elastic_cto_damage   = obj.g_c * obj.elastic_cto;
            else
                obj.elastic_cto_damage   = obj.K * I2I2 + ...
                    obj.g_c * 2 * obj.mu * (I4 - I2I2 / 3);
            end
            
            prm_R    = 1.986e-3; % Gas constant (kcal/(mol.K))
            obj.viscosity = 1.0/(obj.prm_C*exp((-1.0)*obj.prm_Q/(prm_R*obj.temperature)));
            
            % Assign crystal slip system
            nslip = obj.num_slip_system;
            P_Schmid = (1-obj.phi_t) * obj.P_Schmid_notw + obj.phi_t * obj.P_Schmid_twin;
            
            % preprocessing for initial guess
            Gamma0  = zeros(nslip, 1);
            pre_iter = 10;
            tau_s = transpose(P_Schmid) * obj.elastic_cto * obj.new_elastic_strain ...
                - transpose(P_Schmid) * obj.elastic_cto * P_Schmid * Gamma0;
            tau_s = tau_s / pre_iter;
            dt_pre = obj.dt / pre_iter;
            for ii = 1:pre_iter
                tau_y  = obj.old_tau + obj.hardening*sum( abs( Gamma0 ) );
                dGamma0 = dt_pre / obj.viscosity * (abs(tau_s)*ii/tau_y).^(obj.exponent) .* sign(tau_s);
                Gamma0 = Gamma0 + dGamma0;
            end

            % actural Newton iteration
            [resid, jacobian] = obj.formR(Gamma0, P_Schmid);
            resid_ini = norm(resid);
            resid_a = resid_ini;
            resid_r = resid_a / resid_ini;
            Gamma = Gamma0;
            N_iter = 0;
            while resid_r > 1.0e-10 && resid_a > 1.0e-20 && N_iter < 20
                Gamma = Gamma + jacobian \ resid;
                [resid, jacobian] = obj.formR(Gamma, P_Schmid);
                resid_a = norm(resid);
                resid_r = resid_a / resid_ini;
                N_iter = N_iter + 1;
            end
            if N_iter == 20
                error('>>> Error: material model fails to converge!');
            end
            
            % update state variables
            tmp_plastic_strain = P_Schmid * Gamma;
            tmp_plastic_slip_i = abs( Gamma );
            tmp_plastic_slip = sum( tmp_plastic_slip_i );
            tmp_elastic_strain = obj.new_elastic_strain - tmp_plastic_strain;
            obj.new_plastic_slip   = obj.old_plastic_slip   + tmp_plastic_slip;
            obj.new_plastic_slip_i = obj.old_plastic_slip_i + tmp_plastic_slip_i;
            obj.new_stress         = obj.elastic_cto * tmp_elastic_strain;
            obj.new_tau = obj.old_tau + obj.hardening*tmp_plastic_slip;

            % Update Jacobian inverse - Pseudo inverse
            [U,S,V] = svd(jacobian);
            CTO_jacobian_inv = obj.update_jacobian_inverse(S,U,V);
            
            % Update CTO
            tau_s = transpose( P_Schmid ) * obj.new_stress;
            obj.new_cto = obj.elastic_cto_damage;
            for ii = 1:nslip
                for jj = 1:nslip
                    ptaup = obj.exponent * ( abs( tau_s(jj) ) ) ^ ( obj.exponent - 1 );
                    tmp_Schmid_a = obj.elastic_cto * P_Schmid(:,ii);
                    tmp_Schmid_b = obj.elastic_cto * P_Schmid(:,jj);
                    obj.new_cto = obj.new_cto - obj.g_c * ptaup * CTO_jacobian_inv(ii,jj) * tmp_Schmid_a * transpose(tmp_Schmid_b);
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
            
            % update coupled stiffness
            obj.update_tan_utc(CTO_jacobian_inv, Gamma);
        end
        
        function [R,J] = formR(obj, Gamma, P_Schmid)
            % attention: J = - d(R) / d(gamma)
            % verified against finite difference
            nslip = obj.num_slip_system;
            J = zeros(nslip, nslip);
            
            tmp_plastic_strain = P_Schmid * Gamma;
            tmp_plastic_slip_i = abs( Gamma );
            tmp_plastic_slip = sum( tmp_plastic_slip_i );
            tmp_elastic_strain = obj.new_elastic_strain - tmp_plastic_strain;
            fd_new_stress         = obj.elastic_cto * tmp_elastic_strain;
            % Isotropic hardening
            fd_new_tau = obj.old_tau + obj.hardening*tmp_plastic_slip;
            % Residual
            tau = transpose(P_Schmid) * fd_new_stress;
            R = ( abs(tau) ) .^ (obj.exponent) .* sign(tau) ...
                - (fd_new_tau^(obj.exponent)) * obj.viscosity / obj.dt * Gamma;
            
            % Jacobian matrix for local return mapping
            for ii = 1 : nslip
                for jj = 1 : nslip
                    tmp_beta = obj.elastic_cto*P_Schmid(:,jj);

                    J(ii,jj) = transpose(P_Schmid(:,ii)) * tmp_beta * obj.exponent * (abs(tau(ii))^(obj.exponent-1) )...
                        + obj.hardening ...
                        *obj.exponent*fd_new_tau^(obj.exponent-1)*obj.viscosity/obj.dt*Gamma(ii)*sign(Gamma(jj));
                    if (ii == jj)
                        J(ii,jj) = J(ii,jj) ...
                            + fd_new_tau^(obj.exponent)*(obj.viscosity/obj.dt);
                    end

                end
            end
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
        
        % update coupled tangent stiffness
        function update_tan_utc(obj, Dinv, Gamma)
            
            % (0) volumetric stress and deviatoric stress
            I2 = [1; 1; 1; 0; 0; 0];
            I4 = diag([1, 1, 1, 0.5, 0.5, 0.5]);
            one_third = 1/3;
            tr_eps = obj.new_elastic_strain(1) + obj.new_elastic_strain(2) ...
                + obj.new_elastic_strain(3);
            tr_eps_plus  = max(0.0, tr_eps);
            % tr_eps_minus = tr_eps - tr_eps_plus;
            dev_starin = obj.new_elastic_strain - one_third*tr_eps*I2;
            new_stress_plus = obj.K* tr_eps_plus*I2 + 2*obj.mu*I4*dev_starin;
            stress_hat = obj.elastic_cto * obj.new_elastic_strain;
            % new_stress_minus = obj.K* tr_eps_minus*I2;

            % (1) depends on original data
            % strain_n1 = obj.new_elastic_strain + obj.new_plastic_strain + obj.strain_twin;
            strain_tr = obj.new_elastic_strain + obj.new_plastic_strain - obj.old_plastic_strain;
            P_Schmid = (1-obj.phi_t) * obj.P_Schmid_notw + obj.phi_t * obj.P_Schmid_twin;
            dP_dpft = obj.phip_t * (obj.P_Schmid_twin - obj.P_Schmid_notw);
            depst_dpft = obj.gamma0_t * obj.phip_t * obj.sm_t;
            sslip = obj.new_plastic_slip_i;
%             tau_s = transpose( P_Schmid ) * obj.new_stress;
            tau_s = transpose( P_Schmid ) * stress_hat;
            ptaup = obj.exponent * ( abs( tau_s ) ) .^ ( obj.exponent - 1 );

            % (2) depends on (1)
            depstr_dpft = - depst_dpft;
            dgamma_deps = Dinv * ( repmat(ptaup,1,6) .* transpose(P_Schmid) * obj.elastic_cto );
            dgamma_dpft = Dinv * ( ptaup .* ( ...
                  transpose(dP_dpft)  * obj.elastic_cto * strain_tr ...
                + transpose(P_Schmid) * obj.elastic_cto * depstr_dpft ...
                - transpose(dP_dpft)  * obj.elastic_cto * P_Schmid * Gamma ...
                - transpose(P_Schmid) * obj.elastic_cto * dP_dpft  * Gamma ) );

            % (3) depends on (2)
            depsp_dpft = P_Schmid * dgamma_dpft + dP_dpft * Gamma;

            % (4) final results
            obj.dH_dpft = obj.hardening * transpose(sslip) * ( dgamma_dpft .* sign( Gamma ) ) ...
                - transpose( stress_hat ) * ( depsp_dpft + depst_dpft );
            obj.dH_dpfc = 0;
            obj.dH_deps = new_stress_plus - transpose(dgamma_deps) * transpose(P_Schmid) * new_stress_plus ...
                + obj.hardening * transpose(dgamma_deps) * ( sslip .* sign( Gamma ) );
            obj.dsigma_dpft = - obj.g_c * obj.elastic_cto * (depsp_dpft + depst_dpft);
            obj.dsigma_dpfc = obj.gp_c * new_stress_plus;
            obj.dtau_dpft = transpose(obj.dsigma_dpft) * obj.sm_t;
            obj.dtau_dpfc = transpose(obj.dsigma_dpfc) * obj.sm_t;
            obj.dtau_deps = transpose(obj.new_cto) * obj.sm_t;
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
            tmp = obj.dtau_dpft;
        end
        
        function tmp = get_dtau_dpfc(obj)
            tmp = obj.dtau_dpfc;
        end
        
        function tmp = get_dH_dpft(obj)
            tmp = obj.dH_dpft;
        end
        
        function tmp = get_Hn1_c(obj)
            tmp = obj.Hn1_c;
        end
        
        function tmp = get_dH_dpfc(obj)
            tmp = obj.dH_dpfc;
        end
        
        function tmp = get_dH_deps(obj)
            tmp = obj.dH_deps;
        end
        
        function tmp = get_dtau_deps(obj)
            tmp = obj.dtau_deps;
        end
        
        function tmp = get_dsigma_dpft(obj)
            tmp = obj.dsigma_dpft;
        end
        
        function tmp = get_dsigma_dpfc(obj)
            tmp = obj.dsigma_dpfc;
        end
        
    end
    
    methods (Static)
        function jacobian_inv = update_jacobian_inverse(S, U, V)
            s_size = length(S);
            S_inv = zeros(s_size, s_size);
            for ii = 1:s_size
%                 if S(ii,ii)/S(1,1) < 1.0e-6
%                     S_inv(ii,ii) = 0;
%                 else
%                     S_inv(ii,ii) = 1 / S(ii,ii);
%                 end
                S_inv(ii,ii) = 1 / S(ii,ii);
            end
            jacobian_inv = V * S_inv * transpose(U);
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