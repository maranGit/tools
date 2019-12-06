clear
clc

% input material parameters
% @@@@@@@@@@@@@@@@@@@@@ input options @@@@@@@@@@@@@@@@@@@@@
% ul = rand(3,4);
ul = [0.123932277598070,0.873927405861733,0.564979570738201,0.205975515532243;
    0.490357293468018,0.270294332292698,0.640311825162758,0.947933121293169;
    0.852998155340816,0.208461358751314,0.417028951642886,0.0820712070977259;
    0.852998155340816,0.208461358751314,0.417028951642886,0.0820712070977259];
theta = 0*pi/180;
K = 102.0833; % MPa
nu = 0.1667; % MPa
alpha_t = 3;
gamma0_t = 0.1295;
% Gc_t = 1.17e+4; % mJ.mm^{-2} @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ use smaller value for debug
Gc_t = 1.17e-4; % mJ.mm^{-2}
l0_t = 0.015; % mm
Gc_c = 2.7e-3;
l0_c = 0.015;
k_c = 0.001;
% @@@@@@@@@@@@@@@@@@@@@ input options @@@@@@@@@@@@@@@@@@@@@

% @@@@@@@@@@@@@@@@@@@@@ derived @@@@@@@@@@@@@@@@@@@@@
k_t = 0.75 * Gc_t * l0_t; % mJ/mm
A_t = 12 * Gc_t / l0_t; % MPa
s_vec = [cos(theta); sin(theta)];
m_vec = [-sin(theta); cos(theta)];
sm_vec = zeros(3,1);
sm_vec(1) = s_vec(1)*m_vec(1);
sm_vec(2) = s_vec(2)*m_vec(2);
sm_vec(3) = s_vec(1)*m_vec(2) + s_vec(2)*m_vec(1);
lambda = 3*K*nu/(1+nu);
mu = 1.5*K*(1-2*nu)/(1+nu);
sm_vec_2 = sm_vec;
sm_vec_2(3) = sm_vec_2(3) * 0.5;
Gcdl0 = Gc_c / l0_c;
Gcml0 = Gc_c * l0_c;
% @@@@@@@@@@@@@@@@@@@@@ derived @@@@@@@@@@@@@@@@@@@@@


ndf = 4; % u1, u2, t, c
der = 0;
bf = 0;
ib = 0;
nel = 4;
nen = 4;
lint = 4;
xl = [0 1 1 0;
      0 0 1 1];
thick = 1;

ulres = reshape(ul,ndf*nen,1);
block_u = [1,2,5,6,9,10,13,14]';
block_t = [3,7,11,15]';
block_c = [4,8,12,16];

ElemK = zeros(ndf*nel);
ElemF = zeros(ndf*nel,1);
Nmat = zeros(2,2*nel);
Bmat = zeros(3,2*nel);

I4 = eye(3);
I4(3,3) = I4(3,3) * 0.5;
I2 = transpose( [1,1,0] );
I2I2 = I2 * I2';
Cmat = lambda * I2I2 + 2 * mu * I4;

for je = 1:lint
    
    [Wgt,litr,lits] = intpntq(je,lint,ib);
    [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
    [Qxy, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
    
    c1 = Wgt*Jdet*thick;
    
    for mm = 1:nel
        % shape functions
        Nmat(:,2*mm-1:2*mm) = [shl(mm,1)     0
            0        shl(mm,1)    ];
        % derivatives w.r.t. x_n+1^i
        Bmat(:,2*mm-1:2*mm) = [Qxy(mm,1) 0
            0         Qxy(mm,2)
            Qxy(mm,2) Qxy(mm,1)];
    end
    
    disp = ulres(block_u);
    pft_node = ulres(block_t);
    pfc_node = ulres(block_c);
    pft1 = transpose(shl) * pft_node;
    pft2 = pft1*pft1;
    pft3 = pft1*pft2;
    pft4 = pft2 * pft2;
    pft5 = pft2 * pft3;
    pft6 = pft3 * pft3;
    pfc1 = transpose(shl) * pfc_node;
    g_c = (1-k_c)*(1-pfc1)*(1-pfc1) + k_c;
    gp_c = 2*(1-k_c)*(pfc1-1);
    gpp_c = 2*(1-k_c);
    phi_t = alpha_t*pft2 + 2*(2-alpha_t)*pft3 + (alpha_t-3)*pft4;
    phip_t = 2*alpha_t*pft1 + 6*(2-alpha_t)*pft2 + 4*(alpha_t-3)*pft3;
    phipp_t = 2*alpha_t + 12*(2-alpha_t)*pft1 + 12*(alpha_t-3)*pft2;
    
    % displacement residual
    epse = Bmat * disp - phi_t * gamma0_t * sm_vec;
    stress_eff = Cmat * epse;
    stress = g_c * stress_eff;
    if epse(1)+epse(2)>=0 % volumetric tension
        Hn1_c = 0.5 * transpose(stress_eff) * epse;
    else % volumetric compression
        epsev = epse(1) + epse(2);
        epsed = epse;
        epsed(1) = epsed(1) - 0.5 * epsev;
        epsed(2) = epsed(2) - 0.5 * epsev;
        eded = epsed(1)*epsed(1) + epsed(2)*epsed(2) + 0.5*epsed(3)*epsed(3);
        Hn1_c = mu*eded;
    end
    ElemF(block_u) = ElemF(block_u) + c1 * transpose(Bmat) * stress;
    
    % twinning residual
    usm = transpose(sm_vec_2)*Bmat*disp;
    tmp1 = 2*A_t*(pft1-3*pft2+2*pft3);
    tmp2 = mu*gamma0_t*( ( alpha_t*pft2 + 2*(2-alpha_t)*pft3 + (alpha_t-3)*pft4 ) * gamma0_t - 2*usm )...
         * ( 2*alpha_t*pft1 + 6*(2-alpha_t)*pft2 + 4*(alpha_t-3)*pft3 );
    ElemF(block_t) = ElemF(block_t) + c1 * g_c * shl*(tmp1+tmp2) + c1 * g_c * 2*k_t*(Qxy*Qxy')*pft_node;
    
    % crack residual
    pft_grad = Qxy'*pft_node;
    tmp1 = Gcdl0 * pfc1 + gp_c * (Hn1_c + A_t*(pft4 - 2*pft3 + pft2) + k_t*transpose(pft_grad)*pft_grad);
    tmp2 = Gcml0 * (Qxy*Qxy')*pfc_node;
    ElemF(block_c) = ElemF(block_c) + c1 * shl*tmp1 + c1*tmp2;
    
    % displacement stiffness
    tmp1 = c1 * Bmat' * g_c * Cmat * Bmat;
    ElemK(block_u,block_u) = ElemK(block_u,block_u) + tmp1;
    
    % twinning stiffness
    tmp1 = 2*c1*A_t*(1-6*pft1+6*pft2)*(shl*shl');
    tmp2 = 2*c1*k_t*(Qxy*Qxy');
    tmp3 = 6*alpha_t*alpha_t*pft2*gamma0_t + 40*alpha_t*(2-alpha_t)*pft3*gamma0_t + 30*(alpha_t-1)*(3*alpha_t-8)*pft4*gamma0_t ...
        +84*(2-alpha_t)*(alpha_t-3)*pft5*gamma0_t + 28*(alpha_t-3)*(alpha_t-3)*pft6*gamma0_t ...
        - (4*alpha_t+24*(2-alpha_t)*pft1+24*(alpha_t-3)*pft2)*usm;
    tmp4 = c1 * mu * gamma0_t * tmp3 * (shl*shl');
    ElemK(block_t,block_t) = ElemK(block_t,block_t) + g_c * tmp1 + g_c * tmp2 + g_c * tmp4;
    
    % @@@@@@@@@@@@@ stop here @@@@@@@@@@@@@@@@@
    % crack stiffness
    tmp1 = (Gcdl0+gpp_c*(Hn1_c + A_t*(pft4 - 2*pft3 + pft2) + k_t*transpose(pft_grad)*pft_grad))*(shl*shl');
    tmp2 = Gcml0 * (Qxy*Qxy');
    ElemK(block_c,block_c) = ElemK(block_c,block_c) + c1*tmp1 + c1*tmp2;
    
    % displacement + twinning
    tmp1 = -2*mu*gamma0_t*(2*alpha_t*pft1 + 6*(2-alpha_t)*pft2 + 4*(alpha_t-3)*pft3);
    tmp2 = Bmat' * sm_vec_2 * shl';
    ElemK(block_u,block_t) = ElemK(block_u,block_t) + c1*g_c * tmp1*tmp2;
    ElemK(block_t,block_u) = ElemK(block_t,block_u) + c1*g_c * tmp1*tmp2';
    
    % displacement + crack
    dsigma_dpfc = gp_c * stress_eff;
    tmp1 = Bmat' * dsigma_dpfc * shl';
    ElemK(block_u,block_c) = ElemK(block_u,block_c) + c1*tmp1;
    ElemK(block_c,block_u) = ElemK(block_c,block_u) + c1*tmp1';
    
    % twinning + crack
    dtau_dpfc = gp_c * transpose(stress_eff) * sm_vec;
    tmp1 = ( 2*A_t*(2*pft3-3*pft2+pft1)*gp_c - gamma0_t*dtau_dpfc*phip_t ) * (shl * shl');
    tmp2 = 2*k_t*gp_c*Qxy*pft_grad*shl';
    ElemK(block_c,block_t) = ElemK(block_c,block_t) + c1*tmp1 + c1*transpose(tmp2);
    ElemK(block_t,block_c) = ElemK(block_t,block_c) + c1*tmp1 + c1*tmp2;
    
end %je


% get stiffness matrix based on finite difference
tol = 1.0e-6;
ElemK_fd = zeros(ndf*nel);

for ii = 1:ndf*nel
    
    ulres_fd = ulres;
    ulres_fd(ii) = ulres_fd(ii) + tol;
    ElemF_fd = zeros(ndf*nel,1);

for je = 1:lint
    
    [Wgt,litr,lits] = intpntq(je,lint,ib);
    [shl,shld,shls,be] = shlq(litr,lits,nel,nel,der,bf);
    [Qxy, shgs, Jdet, be, xs] = shgq(xl,nel,shld,shls,nen,bf,der,be);
    
    c1 = Wgt*Jdet*thick;
    
    for mm = 1:nel
        % shape functions
        Nmat(:,2*mm-1:2*mm) = [shl(mm,1)     0
            0        shl(mm,1)    ];
        % derivatives w.r.t. x_n+1^i
        Bmat(:,2*mm-1:2*mm) = [Qxy(mm,1) 0
            0         Qxy(mm,2)
            Qxy(mm,2) Qxy(mm,1)];
    end
    
%     disp = ulres_fd(block_u);
%     pft_node = ulres_fd(block_t);
%     pft1 = transpose(shl) * pft_node;
%     pft2 = pft1*pft1;
%     pft3 = pft1*pft2;
%     pft4 = pft2 * pft2;
%     pft5 = pft2 * pft3;
%     pft6 = pft3 * pft3;
    disp = ulres_fd(block_u);
    pft_node = ulres_fd(block_t);
    pfc_node = ulres_fd(block_c);
    pft1 = transpose(shl) * pft_node;
    pft2 = pft1*pft1;
    pft3 = pft1*pft2;
    pft4 = pft2 * pft2;
    pft5 = pft2 * pft3;
    pft6 = pft3 * pft3;
%     phi_t = alpha_t*pft2 + 2*(2-alpha_t)*pft3 + (alpha_t-3)*pft4;
%     stress = Cmat * ( Bmat * disp - phi_t * gamma0_t * sm_vec);
%     ElemF_fd(block_u) = ElemF_fd(block_u) + c1 * transpose(Bmat) * stress;
%     
%     usm = transpose(sm_vec_2)*Bmat*disp;
%     tmp1 = 2*A_t*(pft1-3*pft2+2*pft3);
%     tmp2 = mu*gamma0_t*( ( alpha_t*pft2 + 2*(2-alpha_t)*pft3 + (alpha_t-3)*pft4 ) * gamma0_t - 2*usm )...
%          * ( 2*alpha_t*pft1 + 6*(2-alpha_t)*pft2 + 4*(alpha_t-3)*pft3 );
%     ElemF_fd(block_t) = ElemF_fd(block_t) + c1 * shl*(tmp1+tmp2) + c1 * 2*k_t*(Qxy*Qxy')*pft_node;
    pfc1 = transpose(shl) * pfc_node;
    g_c = (1-k_c)*(1-pfc1)*(1-pfc1) + k_c;
    gp_c = 2*(1-k_c)*(pfc1-1);
    gpp_c = 2*(1-k_c);
    phi_t = alpha_t*pft2 + 2*(2-alpha_t)*pft3 + (alpha_t-3)*pft4;
    phip_t = 2*alpha_t*pft1 + 6*(2-alpha_t)*pft2 + 4*(alpha_t-3)*pft3;
    phipp_t = 2*alpha_t + 12*(2-alpha_t)*pft1 + 12*(alpha_t-3)*pft2;
    
    % displacement residual
    epse = Bmat * disp - phi_t * gamma0_t * sm_vec;
    stress_eff = Cmat * epse;
    stress = g_c * stress_eff;
    if epse(1)+epse(2)>=0 % volumetric tension
        Hn1_c = 0.5 * transpose(stress_eff) * epse;
    else % volumetric compression
        epsev = epse(1) + epse(2);
        epsed = epse;
        epsed(1) = epsed(1) - 0.5 * epsev;
        epsed(2) = epsed(2) - 0.5 * epsev;
        eded = epsed(1)*epsed(1) + epsed(2)*epsed(2) + 0.5*epsed(3)*epsed(3);
        Hn1_c = mu*eded;
    end
    ElemF_fd(block_u) = ElemF_fd(block_u) + c1 * transpose(Bmat) * stress;
    
    % twinning residual
    usm = transpose(sm_vec_2)*Bmat*disp;
    tmp1 = 2*A_t*(pft1-3*pft2+2*pft3);
    tmp2 = mu*gamma0_t*( ( alpha_t*pft2 + 2*(2-alpha_t)*pft3 + (alpha_t-3)*pft4 ) * gamma0_t - 2*usm )...
         * ( 2*alpha_t*pft1 + 6*(2-alpha_t)*pft2 + 4*(alpha_t-3)*pft3 );
    ElemF_fd(block_t) = ElemF_fd(block_t) + c1 * g_c * shl*(tmp1+tmp2) + c1 * g_c * 2*k_t*(Qxy*Qxy')*pft_node;
    
    % crack residual
    pft_grad = Qxy'*pft_node;
    tmp1 = Gcdl0 * pfc1 + gp_c * (Hn1_c + A_t*(pft4 - 2*pft3 + pft2) + k_t*transpose(pft_grad)*pft_grad);
    tmp2 = Gcml0 * (Qxy*Qxy')*pfc_node;
    ElemF_fd(block_c) = ElemF_fd(block_c) + c1 * shl*tmp1 + c1*tmp2;
    
end %je

ElemK_fd(:,ii) = (ElemF_fd - ElemF)/tol;
end

ElemK_diff = ElemK - ElemK_fd;
ElemK_diff(abs(ElemK_diff)<1.0e-3) = 0;
