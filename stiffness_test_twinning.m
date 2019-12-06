clear
clc

% input material parameters
% @@@@@@@@@@@@@@@@@@@@@ input options @@@@@@@@@@@@@@@@@@@@@
% ul = rand(3,4);
ul = [0.123932277598070,0.873927405861733,0.564979570738201,0.205975515532243;
    0.490357293468018,0.270294332292698,0.640311825162758,0.947933121293169;
    0.852998155340816,0.208461358751314,0.417028951642886,0.0820712070977259];
theta = 38*pi/180;
K = 36.7e3; % MPa
nu = 0.276; % MPa
alpha = 3;
gamma_0 = 0.1295;
gamma_pf = 117e-6; % mJ.mm^{-2}
l_pf = 1.0e-6; % mm
% @@@@@@@@@@@@@@@@@@@@@ input options @@@@@@@@@@@@@@@@@@@@@

% @@@@@@@@@@@@@@@@@@@@@ derived @@@@@@@@@@@@@@@@@@@@@
kappa_pf = 0.75 * gamma_pf * l_pf; % mJ/mm
A_pf = 12 * gamma_pf / l_pf; % MPa
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
% @@@@@@@@@@@@@@@@@@@@@ derived @@@@@@@@@@@@@@@@@@@@@


ndf = 3; % u1, u2, c
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
block_u = [1,2,4,5,7,8,10,11]';
block_c = [3,6,9,12]';

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
    eta = ulres(block_c);
    eta_gp = transpose(shl) * eta;
    pf2 = eta_gp*eta_gp;
    pf3 = eta_gp*pf2;
    pf4 = pf2 * pf2;
    pf5 = pf2 * pf3;
    pf6 = pf3 * pf3;
    
    phi_eta = alpha*pf2 + 2*(2-alpha)*pf3 + (alpha-3)*pf4;
    stress = Cmat * ( Bmat * disp - phi_eta * gamma_0 * sm_vec);
    ElemF(block_u) = ElemF(block_u) + c1 * transpose(Bmat) * stress;
    
    usm = transpose(sm_vec_2)*Bmat*disp;
    tmp1 = 2*A_pf*(eta_gp-3*pf2+2*pf3);
    tmp2 = mu*gamma_0*( ( alpha*pf2 + 2*(2-alpha)*pf3 + (alpha-3)*pf4 ) * gamma_0 - 2*usm )...
         * ( 2*alpha*eta_gp + 6*(2-alpha)*pf2 + 4*(alpha-3)*pf3 );
    ElemF(block_c) = ElemF(block_c) + c1 * shl*(tmp1+tmp2) + c1 * 2*kappa_pf*(Qxy*Qxy')*eta;
    
    % stiffness matrix
    tmp1 = c1 * Bmat' * Cmat * Bmat;
    ElemK(block_u,block_u) = ElemK(block_u,block_u) + tmp1;
    
    tmp1 = 2*c1*A_pf*(1-6*eta_gp+6*pf2)*(shl*shl');
    tmp2 = 2*c1*kappa_pf*(Qxy*Qxy');
    tmp3 = 6*alpha*alpha*pf2*gamma_0 + 40*alpha*(2-alpha)*pf3*gamma_0 + 30*(alpha-1)*(3*alpha-8)*pf4*gamma_0 ...
        +84*(2-alpha)*(alpha-3)*pf5*gamma_0 + 28*(alpha-3)*(alpha-3)*pf6*gamma_0 ...
        - (4*alpha+24*(2-alpha)*eta_gp+24*(alpha-3)*pf2)*usm;
    tmp4 = c1 * mu * gamma_0 * tmp3 * (shl*shl');
    ElemK(block_c,block_c) = ElemK(block_c,block_c) + tmp1 + tmp2 + tmp4;
    
    tmp1 = -2*mu*gamma_0*(2*alpha*eta_gp + 6*(2-alpha)*pf2 + 4*(alpha-3)*pf3);
    tmp2 = Bmat' * sm_vec_2 * shl';
    ElemK(block_u,block_c) = ElemK(block_u,block_c) + c1*tmp1*tmp2;
    ElemK(block_c,block_u) = ElemK(block_c,block_u) + c1*tmp1*tmp2';
end %je

%{
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
    
    disp = ulres_fd(block_u);
    eta = ulres_fd(block_c);
    eta_gp = transpose(shl) * eta;
    pf2 = eta_gp*eta_gp;
    pf3 = eta_gp*pf2;
    pf4 = pf2 * pf2;
    pf5 = pf2 * pf3;
    pf6 = pf3 * pf3;
    
    phi_eta = alpha*pf2 + 2*(2-alpha)*pf3 + (alpha-3)*pf4;
    stress = Cmat * ( Bmat * disp - phi_eta * gamma_0 * sm_vec);
    ElemF_fd(block_u) = ElemF_fd(block_u) + c1 * transpose(Bmat) * stress;
    
    usm = transpose(sm_vec_2)*Bmat*disp;
    tmp1 = 2*A_pf*(eta_gp-3*pf2+2*pf3);
    tmp2 = mu*gamma_0*( ( alpha*pf2 + 2*(2-alpha)*pf3 + (alpha-3)*pf4 ) * gamma_0 - 2*usm )...
         * ( 2*alpha*eta_gp + 6*(2-alpha)*pf2 + 4*(alpha-3)*pf3 );
    ElemF_fd(block_c) = ElemF_fd(block_c) + c1 * shl*(tmp1+tmp2) + c1 * 2*kappa_pf*(Qxy*Qxy')*eta;
    
    
end %je

ElemK_fd(:,ii) = (ElemF_fd - ElemF)/tol;
end

ElemK_diff = ElemK - ElemK_fd;
ElemK_diff(abs(ElemK_diff)<1.0e-3) = 0;
%}