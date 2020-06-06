function [L, E, dt_dual, iter] = sralt_inner_admm(D, J, lambda, alpha, p, tol, maxIter, muBound, rho0)

if nargin < 5
    error('Too few arguments') ;
end

if nargin < 6
    tol = 1e-7;
elseif tol == -1
    tol = 1e-7;
end

if nargin < 7
    maxIter = 1000;
elseif maxIter == -1
    maxIter = 1000;
end

if nargin < 8
    muBound = 1e+10;
end

if nargin < 9
    rho0 = 1.25;
end

DISPLAY_EVERY = 10 ;

[m1, m2, m3] = size(D);

L = D;
Y = zeros(m1,m2,m3);
Q_1 = zeros(m1,m2*m3);
Q_2 = zeros(m2,m3*m1);
Q_3 = zeros(m3,m1*m2);

dt_dual = cell(1,m3); 
dt_dual_matrix = zeros(m1*m2, m3);

rho = rho0;
mu = 1.25/norm(unfold(D,3));


d_norm = norm(unfold(D,3), 'fro');

iter = 0;
converged = false;
while ~converged       
    iter = iter + 1;

    % update M_1,M_2,M_3
    temp = unfold(L,1)+(1/mu)*Q_1;
    [U, S, V] = svd(temp, 'econ');
    diagS = diag(S);
    if p == 1
        M_1 = U * diag(pos(diagS-alpha(1)/mu)) * V';
    else
        M_1 = U * diag(solve_Lp(diagS, alpha(1)/mu, p )) * V';
    end
    
    temp = unfold(L,2)+(1/mu)*Q_2;
    [U, S, V] = svd(temp, 'econ');
    diagS = diag(S);
    if p == 1
        M_2 = U * diag(pos(diagS-alpha(2)/mu)) * V';
    else
        M_2 = U * diag(solve_Lp(diagS, alpha(2)/mu, p )) * V';
    end
    
    temp = unfold(L,3)+(1/mu)*Q_3;
    [U, S, V] = svd(temp, 'econ');
    diagS = diag(S);
    if p == 1
        M_3 = U * diag(pos(diagS-alpha(3)/mu)) * V';
    else
        M_3 = U * diag(solve_Lp(diagS, alpha(3)/mu, p )) * V';
    end
    
    % update E
    temp = D + reshape(dt_dual_matrix,[m1 m2 m3]) - L + (1/mu)*Y;
    if p == 1
        E = sign(temp) .* pos( abs(temp) - lambda/mu );
    else
        E = solve_Lp(temp,  lambda/mu, p );
    end
   
    % update L
    T = D - E;
    dt_dual_tensor = reshape(dt_dual_matrix,[m1 m2 m3]);
    TenTemp1 = refold(M_1-(1/mu)*Q_1,1,[m1 m2 m3]);
    TenTemp2 = refold(M_2-(1/mu)*Q_2,2,[m1 m2 m3]);
    TenTemp3 = refold(M_3-(1/mu)*Q_3,3,[m1 m2 m3]);
    L = (T+dt_dual_tensor+(1/mu)*Y+TenTemp1+TenTemp2+TenTemp3)/4;
    
    % update deltaTau
    H = L-T;
    H_wave = H-(1/mu)*Y;
    H_wave3 = unfold(H_wave,3)';
    
    for i = 1 : m3
        dt_dual{i} =  J{i}'*H_wave3(:,i) ;
        dt_dual_matrix(:, i) = J{i}*dt_dual{i} ;
    end
    
    % update Y
    Z = refold(dt_dual_matrix',3,[m1 m2 m3]) - H;
    Y = Y + mu*Z;
    % update Q_i
    Q_1 = Q_1 + mu*(unfold(L,1)-M_1);
    Q_2 = Q_2 + mu*(unfold(L,2)-M_2);
    Q_3 = Q_3 + mu*(unfold(L,3)-M_3);

    mu = min(muBound, mu*rho);
    stoppingCriterion = norm(unfold(Z,3), 'fro') / d_norm;
    
        
    if stoppingCriterion <= tol
        converged = true ;
    end
    
    if ~converged && iter >= maxIter
        converged = 1 ;       
    end
end
