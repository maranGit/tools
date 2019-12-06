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
classdef PhaseFieldMultiCrystal < handle
    properties
        % Input parameters
        K;
        nu;
        G_c;
        l_c;
        S_c;
        tau_ini;
        hardening;
        prm_Q;    % Activation energy
        prm_C;    % Fitting parameter
        exponent; % Power-law parameter
        Euler_Theta;
        Euler_Phi;
        aniso_factor;
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
        
        kappa_reg;
        
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
        
        new_plastic_slip_0;
        old_plastic_slip_0;
        old_old_plastic_slip_0;
        
        new_plastic_slip_1;
        old_plastic_slip_1;
        old_old_plastic_slip_1;
        
        new_plastic_slip_2;
        old_plastic_slip_2;
        old_old_plastic_slip_2;
        
        new_plastic_slip_3;
        old_plastic_slip_3;
        old_old_plastic_slip_3;
        
        new_plastic_slip_4;
        old_plastic_slip_4;
        old_old_plastic_slip_4;
        
        new_plastic_slip_5;
        old_plastic_slip_5;
        old_old_plastic_slip_5;
        
        % NaCl Slip system
        num_slip_system = 6;
        slip_normal = [ 1.0, 1.0, 0.0;
            1.0,-1.0, 0.0;
            1.0, 0.0, 1.0;
            1.0, 0.0,-1.0;
            0.0, 1.0, 1.0;
            0.0, 1.0,-1.0];
        slip_direct = [1.0,-1.0, 0.0;
            1.0, 1.0, 0.0;
            1.0, 0.0,-1.0;
            1.0, 0.0, 1.0;
            0.0, 1.0,-1.0;
            0.0, 1.0, 1.0];
        
        new_active_set; % std::vector<unsigned>
        old_active_set; % std::vector<unsigned>
        
        % Yield stress
        new_tau;
        old_tau;
        
        % Degradation
        degrade;
        pf_eq;
        slip_crit;
        degrade_0;
        degrade_1;
        degrade_2;
        degrade_3;
        degrade_4;
        degrade_5;
        
        degrade_deriv;
        degrade_deriv_2;
        
        degrade_0_deriv;
        degrade_1_deriv;
        degrade_2_deriv;
        degrade_3_deriv;
        degrade_4_deriv;
        degrade_5_deriv;
        degrade_0_deriv_2;
        degrade_1_deriv_2;
        degrade_2_deriv_2;
        degrade_3_deriv_2;
        degrade_4_deriv_2;
        degrade_5_deriv_2;
        
        % Strain energy
        H_plus_updated; % bool
        new_H_plus;
        old_H_plus;
        
        new_H_plus_0;
        old_H_plus_0;
        new_H_plus_1;
        old_H_plus_1;
        new_H_plus_2;
        old_H_plus_2;
        new_H_plus_3;
        old_H_plus_3;
        new_H_plus_4;
        old_H_plus_4;
        new_H_plus_5;
        old_H_plus_5;
        
        H_plastic;
        
        % Structure tensors
        omega_0; % SymmetricTensor<2,3>
        omega_1; % SymmetricTensor<2,3>
        omega_2; % SymmetricTensor<2,3>
        omega_3; % SymmetricTensor<2,3>
        omega_4; % SymmetricTensor<2,3>
        omega_5; % SymmetricTensor<2,3>
        
        % constants
        tol = 1.0e-9;
    end
    methods
        function obj = PhaseFieldMultiCrystal()
        end
        
        function initialize(obj)
            % Parse Parameters
            obj.K = 25; % bulk modulus
            obj.nu = 0.25; % poisson ratio
            obj.G_c = 1.15e-6;
            obj.l_c = 0.02;
            obj.tau_ini = 5.0e-3; % critical resolved shear stress
            obj.hardening = 0.1; % hardening parameter
            obj.prm_Q = 14; % activation energy
            obj.prm_C = 1.0; % shape factor
            obj.exponent = 10.0;
            obj.Euler_Theta = 10.0 * pi / 180;
            obj.Euler_Phi = 30.0 * pi / 180;
            obj.aniso_factor = 40;
            obj.S_c = 0;
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
            obj.new_plastic_slip_0       = 0;
            obj.new_plastic_slip_1       = 0;
            obj.new_plastic_slip_2       = 0;
            obj.new_plastic_slip_3       = 0;
            obj.new_plastic_slip_4       = 0;
            obj.new_plastic_slip_5       = 0;
            obj.new_tau                  = obj.tau_ini;
            obj.new_active_set = zeros(2 * obj.num_slip_system + 1, 1);
            
            obj.new_H_plus      = 0;
            obj.degrade         = 0;
            obj.degrade_deriv   = 0;
            
            obj.new_H_plus_0    = 0;
            obj.degrade_0       = 0;
            obj.degrade_0_deriv = 0;
            
            obj.new_H_plus_1    = 0;
            obj.degrade_1       = 0;
            obj.degrade_1_deriv = 0;
            
            obj.new_H_plus_2    = 0;
            obj.degrade_2       = 0;
            obj.degrade_2_deriv = 0;
            
            obj.new_H_plus_3    = 0;
            obj.degrade_3       = 0;
            obj.degrade_3_deriv = 0;
            
            obj.new_H_plus_4    = 0;
            obj.degrade_4       = 0;
            obj.degrade_4_deriv = 0;
            
            obj.new_H_plus_5    = 0;
            obj.degrade_5       = 0;
            obj.degrade_5_deriv = 0;
            
            obj.H_plastic       = 0;
            
            obj.kappa_reg       = 1e-3;
            
            % End
            
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
            obj.omega_0 = zeros(6,1);
            obj.omega_1 = zeros(6,1);
            obj.omega_2 = zeros(6,1);
            obj.omega_3 = zeros(6,1);
            obj.omega_4 = zeros(6,1);
            obj.omega_5 = zeros(6,1);
            
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
                obj.run_3D_update(incr_strain_2d);
            else
                obj.run_3D_update(incr_strain);
            end
        end
        
        function update_parameters(obj, pass_temp, pass_time)
            obj.temperature = pass_temp + 273.0;
            obj.dt = pass_time;
        end
        
        function update_phasefield(obj, pf_0, pf_1, pf_2, pf_3, pf_4, pf_5)
            obj.degrade_0 = (1-pf_0)*(1-pf_0);
            obj.degrade_1 = (1-pf_1)*(1-pf_1);
            obj.degrade_2 = (1-pf_2)*(1-pf_2);
            obj.degrade_3 = (1-pf_3)*(1-pf_3);
            obj.degrade_4 = (1-pf_4)*(1-pf_4);
            obj.degrade_5 = (1-pf_5)*(1-pf_5);
            obj.degrade   = (1 - obj.kappa_reg) * obj.degrade_0...
                *obj.degrade_1*obj.degrade_2*obj.degrade_3*obj.degrade_4...
                *obj.degrade_5 + obj.kappa_reg;
            
            obj.degrade_0_deriv   = -2*(1-pf_0);
            obj.degrade_1_deriv   = -2*(1-pf_1);
            obj.degrade_2_deriv   = -2*(1-pf_2);
            obj.degrade_3_deriv   = -2*(1-pf_3);
            obj.degrade_4_deriv   = -2*(1-pf_4);
            obj.degrade_5_deriv   = -2*(1-pf_5);
            
            obj.degrade_0_deriv_2 = 2;
            obj.degrade_1_deriv_2 = 2;
            obj.degrade_2_deriv_2 = 2;
            obj.degrade_3_deriv_2 = 2;
            obj.degrade_4_deriv_2 = 2;
            obj.degrade_5_deriv_2 = 2;
            obj.degrade_deriv_2 = 2;
            
            % Update anisotropic plane
            obj.update_slip_plane();
        end
        
        function update_slip_plane(obj)
            % Slip plane
            aniso = obj.aniso_factor;
            A = zeros (6, obj.num_slip_system); % for elasticity
            
            rotation_matrix = obj.update_R_matrix();
            
            for i = 1:obj.num_slip_system
                structure = transpose( obj.slip_normal(i,1:3) );
                % Normalize
                structure = structure / norm(structure);
                Rs = rotation_matrix * structure;
                A_temp = Rs * transpose(Rs);
                A(1,i) = A_temp(1,1);
                A(2,i) = A_temp(2,2);
                A(3,i) = A_temp(3,3);
                A(4,i) = 0.5 * ( A_temp(1,2) + A_temp(2,1) );
                A(5,i) = 0.5 * ( A_temp(1,3) + A_temp(3,1) );
                A(6,i) = 0.5 * ( A_temp(2,3) + A_temp(3,2) );
            end
            
            I2 = [1; 1; 1; 0; 0; 0];
            obj.omega_0 = I2 + aniso*(I2 - A(:,1));
            obj.omega_1 = I2 + aniso*(I2 - A(:,2));
            obj.omega_2 = I2 + aniso*(I2 - A(:,3));
            obj.omega_3 = I2 + aniso*(I2 - A(:,4));
            obj.omega_4 = I2 + aniso*(I2 - A(:,5));
            obj.omega_5 = I2 + aniso*(I2 - A(:,6));
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
            
            %  new_stress_plus = 0.5*K* tr_eps_plus*eye + 2*mu*dev_starin;
            %  new_stress_minus = 0.5*K* tr_eps_minus*eye;
            %  Ran: bug?
            new_stress_plus = obj.K* tr_eps_plus*I2 + 2*obj.mu*I4*dev_starin;
            new_stress_minus = obj.K* tr_eps_minus*I2;
            
            obj.new_stress = obj.degrade*new_stress_plus + new_stress_minus;
            
            % Update strain history variable
            strain_energy_plus = 0.5*obj.K*tr_eps_plus*tr_eps_plus...
                + obj.mu*transpose(dev_starin)*I4*dev_starin;
            plastic_energy = 0.5*obj.hardening*(...
                obj.new_plastic_slip_0 * obj.new_plastic_slip_0 ...
                + obj.new_plastic_slip_1 * obj.new_plastic_slip_1 ...
                + obj.new_plastic_slip_2 * obj.new_plastic_slip_2 ...
                + obj.new_plastic_slip_3 * obj.new_plastic_slip_3 ...
                + obj.new_plastic_slip_4 * obj.new_plastic_slip_4 ...
                + obj.new_plastic_slip_5 * obj.new_plastic_slip_5);
            
            obj.H_plastic = plastic_energy;
            
            if strain_energy_plus > obj.old_H_plus
                obj.new_H_plus = strain_energy_plus;
            else
                obj.new_H_plus = obj.old_H_plus;
            end
            
            reg = 1 - obj.kappa_reg;
            
            % Update H plus for multiphase
            d0 = obj.degrade_0;
            d1 = obj.degrade_1;
            d2 = obj.degrade_2;
            d3 = obj.degrade_3;
            d4 = obj.degrade_4;
            d5 = obj.degrade_5;
            
            dd = d1*d2*d3*d4*d5;
            strain_energy = strain_energy_plus+plastic_energy;
            if reg*dd*strain_energy > obj.old_H_plus_0
                obj.new_H_plus_0 = reg*dd*strain_energy;
            else
                obj.new_H_plus_0 = obj.old_H_plus_0;
            end
            
            dd = d0*d2*d3*d4*d5;
            if reg*dd*strain_energy > obj.old_H_plus_1
                obj.new_H_plus_1 = reg*dd*strain_energy;
            else
                obj.new_H_plus_1 = obj.old_H_plus_1;
            end
            
            dd = d0*d1*d3*d4*d5;
            if reg*dd*strain_energy > obj.old_H_plus_2
                obj.new_H_plus_2 = reg*dd*strain_energy;
            else
                obj.new_H_plus_2 = obj.old_H_plus_2;
            end
            
            dd = d0*d1*d2*d4*d5;
            if reg*dd*strain_energy > obj.old_H_plus_3
                obj.new_H_plus_3 = reg*dd*strain_energy;
            else
                obj.new_H_plus_3 = obj.old_H_plus_3;
            end
            
            dd = d0*d1*d2*d3*d5;
            if reg*dd*strain_energy > obj.old_H_plus_4
                obj.new_H_plus_4 = reg*dd*strain_energy;
            else
                obj.new_H_plus_4 = obj.old_H_plus_4;
            end
            
            dd = d0*d1*d2*d3*d4;
            if reg*dd*strain_energy > obj.old_H_plus_5
                obj.new_H_plus_5 = reg*dd*strain_energy;
            else
                obj.new_H_plus_5 = obj.old_H_plus_5;
            end
            
            %  new_cto = degrade*(K*eye_dyad_eye + 2*mu*(big_eye - one_third*eye_dyad_eye)) + K*eye_dyad_eye;
            % Ran: bug?
            % obj.new_cto = obj.degrade * ( obj.K * I2I2 + 2 * obj.mu * (I4 - one_third*I2I2) );
            obj.new_cto = obj.degrade * obj.new_cto;
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
            obj.old_old_plastic_slip_0   = obj.old_plastic_slip_0;
            obj.old_old_plastic_slip_1   = obj.old_plastic_slip_1;
            obj.old_old_plastic_slip_2   = obj.old_plastic_slip_2;
            obj.old_old_plastic_slip_3   = obj.old_plastic_slip_3;
            obj.old_old_plastic_slip_4   = obj.old_plastic_slip_4;
            obj.old_old_plastic_slip_5   = obj.old_plastic_slip_5;
            obj.old_plastic_slip         = obj.new_plastic_slip;
            obj.old_plastic_slip_0       = obj.new_plastic_slip_0;
            obj.old_plastic_slip_1       = obj.new_plastic_slip_1;
            obj.old_plastic_slip_2       = obj.new_plastic_slip_2;
            obj.old_plastic_slip_3       = obj.new_plastic_slip_3;
            obj.old_plastic_slip_4       = obj.new_plastic_slip_4;
            obj.old_plastic_slip_5       = obj.new_plastic_slip_5;
            % Driving force
            obj.old_H_plus               = obj.new_H_plus;
            obj.old_H_plus_0             = obj.new_H_plus_0;
            obj.old_H_plus_1             = obj.new_H_plus_1;
            obj.old_H_plus_2             = obj.new_H_plus_2;
            obj.old_H_plus_3             = obj.new_H_plus_3;
            obj.old_H_plus_4             = obj.new_H_plus_4;
            obj.old_H_plus_5             = obj.new_H_plus_5;
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
            obj.new_plastic_slip_0       = obj.old_plastic_slip_0;
            obj.new_plastic_slip_1       = obj.old_plastic_slip_1;
            obj.new_plastic_slip_2       = obj.old_plastic_slip_2;
            obj.new_plastic_slip_3       = obj.old_plastic_slip_3;
            obj.new_plastic_slip_4       = obj.old_plastic_slip_4;
            obj.new_plastic_slip_5       = obj.old_plastic_slip_5;
            obj.old_plastic_slip         = obj.old_old_plastic_slip;
            obj.old_plastic_slip_0       = obj.old_old_plastic_slip_0;
            obj.old_plastic_slip_1       = obj.old_old_plastic_slip_1;
            obj.old_plastic_slip_2       = obj.old_old_plastic_slip_2;
            obj.old_plastic_slip_3       = obj.old_old_plastic_slip_3;
            obj.old_plastic_slip_4       = obj.old_old_plastic_slip_4;
            obj.old_plastic_slip_5       = obj.old_old_plastic_slip_5;
            % Driving force
            obj.new_H_plus               = obj.old_H_plus;
            obj.new_H_plus_0             = obj.old_H_plus_0;
            obj.new_H_plus_1             = obj.old_H_plus_1;
            obj.new_H_plus_2             = obj.old_H_plus_2;
            obj.new_H_plus_3             = obj.old_H_plus_3;
            obj.new_H_plus_4             = obj.old_H_plus_4;
            obj.new_H_plus_5             = obj.old_H_plus_5;
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
        
        function fren = fracture_energy(obj)
            fren = obj.G_c;
        end
        
        function lgsc = length_scale(obj)
            lgsc = obj.l_c;
        end
        
        function psn = plastic_slip_new(obj)
            psn = obj.new_plastic_slip;
        end
        
        function plsl = plastic_slip_old(obj)
            plsl = obj.old_plastic_slip;
        end
        
        function plsl = plastic_slip_new_0(obj)
            plsl = obj.new_plastic_slip_0;
        end
        
        function plsl = plastic_slip_new_1(obj)
            plsl = obj.new_plastic_slip_1;
        end
        
        function plsl = plastic_slip_new_2(obj)
            plsl = obj.new_plastic_slip_2;
        end
        
        function plsl = plastic_slip_new_3(obj)
            plsl = obj.new_plastic_slip_3;
        end
        
        function plsl = plastic_slip_new_4(obj)
            plsl = obj.new_plastic_slip_4;
        end
        
        function plsl = plastic_slip_new_5(obj)
            plsl = obj.new_plastic_slip_5;
        end
        
        function temp = anisotropy_omega_0(obj)
            temp = obj.omega_0;
        end
        
        function temp = anisotropy_omega_1(obj)
            temp = obj.omega_1;
        end
        
        function temp = anisotropy_omega_2(obj)
            temp = obj.omega_2;
        end
        
        function temp = anisotropy_omega_3(obj)
            temp = obj.omega_3;
        end
        
        function temp = anisotropy_omega_4(obj)
            temp = obj.omega_4;
        end
        
        function temp = anisotropy_omega_5(obj)
            temp = obj.omega_5;
        end
        
        function temp = crack_driving_force(obj)
            temp = obj.degrade_deriv * obj.new_H_plus;
        end
        
        function temp = crack_driving_force_deriv(obj)
            temp = obj.degrade_deriv_2*obj.new_H_plus;
        end
        
        function temp = crack_driving_force_0(obj)
            temp = obj.degrade_0_deriv*obj.new_H_plus_0;
        end
        
        function temp = crack_driving_force_deriv_0(obj)
            temp = obj.degrade_0_deriv_2*obj.new_H_plus_0;
        end
        
        function temp = crack_driving_force_1(obj)
            temp = obj.degrade_1_deriv*obj.new_H_plus_1;
        end
        
        function temp = crack_driving_force_deriv_1(obj)
            temp = obj.degrade_1_deriv_2*obj.new_H_plus_1;
        end
        
        function temp = crack_driving_force_2(obj)
            temp = obj.degrade_2_deriv*obj.new_H_plus_2;
        end
        
        function temp = crack_driving_force_deriv_2(obj)
            temp = obj.degrade_2_deriv_2*obj.new_H_plus_2;
        end
        
        function temp = crack_driving_force_3(obj)
            temp = obj.degrade_3_deriv*obj.new_H_plus_3;
        end
        
        function temp = crack_driving_force_deriv_3(obj)
            temp = obj.degrade_3_deriv_2*obj.new_H_plus_3;
        end
        
        function temp = crack_driving_force_4(obj)
            temp = obj.degrade_4_deriv*obj.new_H_plus_4;
        end
        
        function temp = crack_driving_force_deriv_4(obj)
            temp = obj.degrade_4_deriv_2*obj.new_H_plus_4;
        end
        
        function temp = crack_driving_force_5(obj)
            temp = obj.degrade_5_deriv*obj.new_H_plus_5;
        end
        
        function temp = crack_driving_force_deriv_5(obj)
            temp = obj.degrade_5_deriv_2*obj.new_H_plus_5;
        end
        
        function temp = strain_energy_plus(obj)
            temp = obj.new_H_plus;
        end
        
        function temp = plastic_strain_energy(obj)
            temp = obj.H_plastic;
        end
        
        function temp = update_H_plus_0(obj)
            temp = obj.new_H_plus_0;
        end
        
        function temp = update_H_plus_1(obj)
            temp = obj.new_H_plus_1;
        end
        
        function temp = update_H_plus_2(obj)
            temp = obj.new_H_plus_2;
        end
        
        function temp = update_H_plus_3(obj)
            temp = obj.new_H_plus_3;
        end
        
        function temp = update_H_plus_4(obj)
            temp = obj.new_H_plus_4;
        end
        
        function temp = update_H_plus_5(obj)
            temp = obj.new_H_plus_5;
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
            old_plastic_slip_rate_0 = obj.old_plastic_slip_0 - obj.old_old_plastic_slip_0;
            old_plastic_slip_rate_1 = obj.old_plastic_slip_1 - obj.old_old_plastic_slip_1;
            old_plastic_slip_rate_2 = obj.old_plastic_slip_2 - obj.old_old_plastic_slip_2;
            old_plastic_slip_rate_3 = obj.old_plastic_slip_3 - obj.old_old_plastic_slip_3;
            old_plastic_slip_rate_4 = obj.old_plastic_slip_4 - obj.old_old_plastic_slip_4;
            old_plastic_slip_rate_5 = obj.old_plastic_slip_5 - obj.old_old_plastic_slip_5;
            
            old_hds_0 = obj.hardening * obj.old_plastic_slip_0 * old_plastic_slip_rate_0;
            old_hds_1 = obj.hardening * obj.old_plastic_slip_1 * old_plastic_slip_rate_1;
            old_hds_2 = obj.hardening * obj.old_plastic_slip_2 * old_plastic_slip_rate_2;
            old_hds_3 = obj.hardening * obj.old_plastic_slip_3 * old_plastic_slip_rate_3;
            old_hds_4 = obj.hardening * obj.old_plastic_slip_4 * old_plastic_slip_rate_4;
            old_hds_5 = obj.hardening * obj.old_plastic_slip_5 * old_plastic_slip_rate_5;
            
            tmp_hds = old_hds_0 + old_hds_1 + old_hds_2 + old_hds_3 + old_hds_4 + old_hds_5;
            
            temp = -tmp_hds;
        end
        
        function temp = old_elastic_strain_rate(obj)
            tmp_elastic_strain = obj.old_elastic_strain - obj.old_old_elastic_strain;
            temp = tmp_elastic_strain(1) + tmp_elastic_strain(2) + tmp_elastic_strain(3);
        end
        
        % Update
        function run_3D_update(obj,incr_strain)
            % Compute a trial state in which the increment is assumed to be fully elastic
            obj.new_elastic_strain       = obj.old_elastic_strain + incr_strain;
            obj.new_cto                  = obj.elastic_cto;
            obj.new_stress               = obj.new_cto*obj.new_elastic_strain;
            obj.new_plastic_strain       = obj.old_plastic_strain;
            obj.new_equiv_plastic_strain = obj.old_equiv_plastic_strain;
            obj.new_plastic_slip         = obj.old_plastic_slip;
            obj.new_plastic_slip_0       = obj.old_plastic_slip_0;
            obj.new_plastic_slip_1       = obj.old_plastic_slip_1;
            obj.new_plastic_slip_2       = obj.old_plastic_slip_2;
            obj.new_plastic_slip_3       = obj.old_plastic_slip_3;
            obj.new_plastic_slip_4       = obj.old_plastic_slip_4;
            obj.new_plastic_slip_5       = obj.old_plastic_slip_5;
            obj.new_tau                  = obj.old_tau;
            obj.new_active_set           = obj.old_active_set;
            
            %double exponent = 3.26;
            %double prm_C    = 3.16e+8;  % Fitting parameter (s^-1)
            %double prm_Q    = 12.0;     % Activation energy (kcal/mol)
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
                
                P_Schmid(1,ii) = P_Schmid_tmp(1,1);
                P_Schmid(2,ii) = P_Schmid_tmp(2,2);
                P_Schmid(3,ii) = P_Schmid_tmp(3,3);
                P_Schmid(4,ii) = P_Schmid_tmp(1,2)+P_Schmid_tmp(2,1);
                P_Schmid(5,ii) = P_Schmid_tmp(1,3)+P_Schmid_tmp(3,1);
                P_Schmid(6,ii) = P_Schmid_tmp(2,3)+P_Schmid_tmp(3,2);
                
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
                obj.new_strain = obj.new_elastic_strain + obj.new_plastic_strain;   %
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
            % tmp_plastic_slip   = 0;
            % tmp_plastic_slip_0 = 0;
            % tmp_plastic_slip_1 = 0;
            % tmp_plastic_slip_2 = 0;
            % tmp_plastic_slip_3 = 0;
            % tmp_plastic_slip_4 = 0;
            % tmp_plastic_slip_5 = 0;
            
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
                        tmp_plastic_slip_0 = 0;
                        tmp_plastic_slip_1 = 0;
                        tmp_plastic_slip_2 = 0;
                        tmp_plastic_slip_3 = 0;
                        tmp_plastic_slip_4 = 0;
                        tmp_plastic_slip_5 = 0;
                    
                        % Temporary plastic strain
                        for ii = 1:active_set(n21)
                            alpha = active_set(ii);
                            tmp_plastic_strain = tmp_plastic_strain + Gamma(ii)*P_Schmid(:,alpha);
                            tmp_plastic_slip   = tmp_plastic_slip + Gamma(ii);

                            % if (alpha == 0 || alpha == 6)
                            if (alpha == 1 || alpha == 7)
                                tmp_plastic_slip_0 = tmp_plastic_slip_0 + Gamma(ii);
                            end

                            % if (alpha == 1 || alpha == 7)
                            if (alpha == 2 || alpha == 8)
                                tmp_plastic_slip_1 = tmp_plastic_slip_1 + Gamma(ii);
                            end

                            % if (alpha == 2 || alpha == 8)
                            if (alpha == 3 || alpha == 9)
                                tmp_plastic_slip_2 = tmp_plastic_slip_2 + Gamma(ii);
                            end

                            % if (alpha == 3 || alpha == 9)
                            if (alpha == 4 || alpha == 10)
                                tmp_plastic_slip_3 = tmp_plastic_slip_3 + Gamma(ii);
                            end

                            % if (alpha == 4 || alpha == 10)
                            if (alpha == 5 || alpha == 11)
                                tmp_plastic_slip_4 = tmp_plastic_slip_4 + Gamma(ii);
                            end

                            % if (alpha == 5 || alpha == 11)
                            if (alpha == 6 || alpha == 12)
                                tmp_plastic_slip_5 = tmp_plastic_slip_5 + Gamma(ii);
                            end
                        end

                        tmp_elastic_strain = obj.new_elastic_strain - tmp_plastic_strain;

                        obj.new_plastic_slip   = obj.old_plastic_slip   + tmp_plastic_slip;
                        obj.new_plastic_slip_0 = obj.old_plastic_slip_0 + tmp_plastic_slip_0;
                        obj.new_plastic_slip_1 = obj.old_plastic_slip_1 + tmp_plastic_slip_1;
                        obj.new_plastic_slip_2 = obj.old_plastic_slip_2 + tmp_plastic_slip_2;
                        obj.new_plastic_slip_3 = obj.old_plastic_slip_3 + tmp_plastic_slip_3;
                        obj.new_plastic_slip_4 = obj.old_plastic_slip_4 + tmp_plastic_slip_4;
                        obj.new_plastic_slip_5 = obj.old_plastic_slip_5 + tmp_plastic_slip_5;

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
            obj.new_strain = obj.new_elastic_strain + obj.new_plastic_strain;
            obj.update_strain_energy();
        end
        
        function euler_matrix = update_R_matrix(obj)
            % construct rotation matrix
            cos_theta = cos(obj.Euler_Theta);
            sin_theta = -sin(obj.Euler_Theta);
            cos_phi = cos(obj.Euler_Phi);
            sin_phi = -sin(obj.Euler_Phi);
            
            theta_matrix = zeros(3,3);
            phi_matrix = zeros(3,3);
            
            % Assign each element (Miehe and Schroder, 2001)
            theta_matrix(1,1) = cos_theta;
            theta_matrix(1,2) = 0.0;
            theta_matrix(1,3) = -sin_theta;
            theta_matrix(2,1) = 0.0;
            theta_matrix(2,2) = 1.0;
            theta_matrix(2,3) = 0.0;
            theta_matrix(3,1) = sin_theta;
            theta_matrix(3,2) = 0.0;
            theta_matrix(3,3) = cos_theta;
            
            phi_matrix(1,1) = cos_phi;
            phi_matrix(1,2) = sin_phi;
            phi_matrix(1,3) = 0.0;
            phi_matrix(2,1) = -sin_phi;
            phi_matrix(2,2) = cos_phi;
            phi_matrix(2,3) = 0.0;
            phi_matrix(3,1) = 0.0;
            phi_matrix(3,2) = 0.0;
            phi_matrix(3,3) = 1.0;
            
            euler_matrix = theta_matrix * phi_matrix;
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