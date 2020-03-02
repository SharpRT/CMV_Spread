function [HstarSoln,IrstarSoln,ZrstarSoln] = sim_initSolve(...
    theta,p,sigma,omega,rho_r,lambda_r,alpha_var,tau,gamma_r,P,r,K,rhoRep...
)
    syms Hstar Irstar Zrstar
    H_eqn = theta * sigma + ((1-theta)*sigma*Hstar/(Hstar+(1-p)*Irstar)) - omega*Hstar - (lambda_r*Zrstar*Hstar) + rhoRep * rho_r * Irstar;
    Ir_eqn = ((1-theta)*sigma*(1-p)*Irstar/(Hstar+(1-p)*Irstar)) + (lambda_r*Zrstar*Hstar)-(omega + rho_r) * Irstar;
    Zr_eqn = -(alpha_var+tau)*Zrstar+gamma_r*(P(alpha_var,r,K)-Zrstar)*Irstar;

    S = solve(H_eqn,Ir_eqn,Zr_eqn);

    solutions = [S.Hstar, S.Irstar, S.Zrstar];
    solutions = vpa(solutions);

    viableSolns = all(solutions>0,2);
    numViableSolns = nnz(viableSolns);
    index = find(viableSolns);

    if (numViableSolns < 1) 
        disp('No Viable Solutions for IC s.t. H,I,Z>0; setting to H=sigma/omega,I,Z=0');

        viableSolns = all(solutions>=0,2);
        numViableSolns = nnz(viableSolns);
        if (numViableSolns == 0)
            error('No Viable Solns AT ALL');
        elseif (numViableSolns > 1) 
            error('Multiple WEAKly viable solns')
        else
            index = find(viableSolns);                       
        end
    end

    J = [diff(H_eqn,Hstar) diff(H_eqn,Irstar) diff(H_eqn,Zrstar);...
        diff(Ir_eqn,Hstar) diff(Ir_eqn,Irstar) diff(Ir_eqn,Zrstar);...
        diff(Zr_eqn,Hstar) diff(Zr_eqn,Irstar) diff(Zr_eqn,Zrstar)];

    stableSoln_indices = [];
    for i_soln=1:numViableSolns
        Ji = subs(J,[Hstar,Irstar,Zrstar],solutions(index(i_soln),:));
        eigenvalues = eig(Ji);
        is_stableSoln = all(real(eigenvalues)<0);
        if is_stableSoln
            stableSoln_indices = [stableSoln_indices index(i_soln)];
        end
    end

    if (length(stableSoln_indices)<1) 
        error('No Viable/Stable Solutions for IC');
    elseif (length(stableSoln_indices)>1) 
        error('Multiple viable/stable solutions - must select');
    end    

    for i_soln=1:numViableSolns
        if (~isempty(find(stableSoln_indices == index(i_soln))))
            solutions(stableSoln_indices,:); %vpa(real(S
            HstarSoln = double(vpa(S.Hstar(index(i_soln)))); %vpa(real(S
            IrstarSoln = double(vpa(S.Irstar(index(i_soln)))); %vpa(real(S
            ZrstarSoln = double(vpa(S.Zrstar(index(i_soln)))); %vpa(real(S                                                      
        end
    end           
end
         