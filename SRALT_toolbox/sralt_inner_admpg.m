function [L, E, dt_dual, iter] = sralt_inner_admpg(D, J, lambda, alpha, p, tau, tol, maxIter, muBound, rho0)

if nargin < 6
    error('Too few arguments') ;
end

if nargin < 7
    tol = [1e-4 1e-5];
elseif tol == -1
    tol = [1e-4 1e-5];
end

if nargin < 8
    maxIter = 1000;
elseif maxIter == -1
    maxIter = 1000;
end

if nargin < 9
    muBound = 1e-10;
end

if nargin < 10
    rho0 = 1.25;
end

tau1 = tau(1);
tau2 = tau(2);
tol1 = tol(1);
tol2 = tol(2);

[m1, m2, m3] = size(D);

E = zeros(m1,m2,m3);
L = D;
M_1 = unfold(L,1);
M_2 = unfold(L,2);
M_3 = unfold(L,3);

Y = zeros(m1,m2,m3);
Q_1 = zeros(m1,m2*m3);
Q_2 = zeros(m2,m3*m1);
Q_3 = zeros(m3,m1*m2);

dt_dual = cell(1,m3);  
dt_dual_matrix = zeros(m1*m2, m3); 


mu = norm(unfold(D,3))/1.25;
d_norm = norm(unfold(D,3), 'fro'); 

iter = 0;
converged = false;

dt_dual_tensor = reshape(dt_dual_matrix,[m1 m2 m3]);
T = L+E-D-dt_dual_tensor;
while ~converged       
    iter = iter + 1;

    % update M_1,M_2,M_3
    temp = M_1 - tau1*(M_1-unfold(L,1)+mu*Q_1);
    [U, S, V] = svd(temp, 'econ');
    diagS = diag(S);
    if p == 1
        M_1_new = U * diag(pos(diagS-alpha(1)*mu*tau1)) * V';
    else
        M_1_new = U * diag(solve_Lp( diagS, alpha(1)*mu*tau1, p )) * V';
    end
    
    temp = M_2 - tau1*(M_2-unfold(L,2)+mu*Q_2);
    [U, S, V] = svd(temp, 'econ');
    diagS = diag(S);
    if p == 1
        M_2_new = U * diag(pos(diagS-alpha(2)*mu*tau1)) * V';
    else
        M_2_new = U * diag(solve_Lp( diagS, alpha(2)*mu*tau1, p )) * V';
    end
    
    temp = M_3 - tau1*(M_3-unfold(L,3)+mu*Q_3);
    [U, S, V] = svd(temp, 'econ');
    diagS = diag(S);
    if p == 1
        M_3_new = U * diag(pos(diagS-alpha(3)*mu*tau1)) * V';
    else
        M_3_new = U * diag(solve_Lp( diagS, alpha(3)*mu*tau1, p )) * V';
    end
    
     % update L    
    TenTemp1 = L-refold(M_1+mu*Q_1,1,[m1 m2 m3]);
    TenTemp2 = L-refold(M_2+mu*Q_2,2,[m1 m2 m3]);
    TenTemp3 = L-refold(M_3+mu*Q_3,3,[m1 m2 m3]);
    L_new = L-tau1*(TenTemp1+TenTemp2+TenTemp3+T-mu*Y);
    
    % update E
    temp = E-tau1*(T-mu*Y); 
    if p == 1
        E_new = sign(temp) .* pos( abs(temp) - lambda*mu*tau1 );
    else
        E_new = solve_Lp(temp, lambda*mu*tau1, p);
    end
    
    % update deltaTau
    H_wave = dt_dual_tensor-tau2*(dt_dual_tensor-L_new-E_new+D+mu*Y);   
    H_wave3 = unfold(H_wave,3)';
    
    for i = 1 : m3
        dt_dual{i} =  J{i}'*H_wave3(:,i) ; 
        dt_dual_matrix(:, i) = J{i}*dt_dual{i} ;
    end
    
    % update Y
    dt_dual_tensor_new = reshape(dt_dual_matrix,[m1 m2 m3]);
    T = L_new+E_new-D-dt_dual_tensor_new;
    Y = Y-T/mu;

    % update Q_i
    Q_1 = Q_1 - (unfold(L_new,1)-M_1_new)/mu;
    Q_2 = Q_2 - (unfold(L_new,2)-M_2_new)/mu;
    Q_3 = Q_3 - (unfold(L_new,3)-M_3_new)/mu;
    
    % update mu
    M_1_Difnorm = norm(M_1-M_1_new,'fro');
    M_2_Difnorm = norm(M_2-M_2_new,'fro');
    M_3_Difnorm = norm(M_3-M_3_new,'fro');
    L_Difnorm = norm(unfold(L-L_new,3),'fro');
    E_Difnorm = norm(unfold(E-E_new,3),'fro');
    x_norm = M_1_Difnorm+M_2_Difnorm+M_3_Difnorm+L_Difnorm+E_Difnorm;
    deltatau_norm = norm(unfold(dt_dual_tensor-dt_dual_tensor_new,3),'fro');
    stoppingCriterion2 = max(x_norm/sqrt(tau1),deltatau_norm/sqrt(tau2))/(mu*d_norm);
    if stoppingCriterion2 < tol2
        rho = rho0;
    else
        rho = 1;
    end
    mu = max(muBound, mu/ rho);
        
    stoppingCriterion1 = (norm(unfold(T,3), 'fro') + norm(unfold(L_new,1)-M_1_new, 'fro')...
                            + norm(unfold(L_new,2)-M_2_new, 'fro') + norm(unfold(L_new,3)-M_3_new, 'fro'))/ d_norm;
   
     
    M_1 = M_1_new;
    M_2 = M_2_new;
    M_3 = M_3_new;
    L = L_new;
    E = E_new;
    dt_dual_tensor = dt_dual_tensor_new;
    
    if stoppingCriterion1 <= tol1 && stoppingCriterion2 <= tol2
        disp('RASL inner loop is converged at:');
        disp(['#Iteration ' num2str(iter) '  rank(unfold(L,1)) ' num2str(rank(unfold(L,1))) ...
             '  rank(unfold(L,2)) ' num2str(rank(unfold(L,2)))  '  rank(unfold(L,3)) ' num2str(rank(unfold(L,3)))...
            ' ||E||_0 ' num2str(length(find(abs(E)>0)))  '  Stopping Criterion1 ' ...
            num2str(stoppingCriterion1)  '  Stopping Criterion2 ' ...
            num2str(stoppingCriterion2)]) ;
        converged = true ;
    end
    
    if ~converged && iter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
   
end
