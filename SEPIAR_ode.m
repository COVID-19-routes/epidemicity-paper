function dxdt=SEPIAR_ode(~,x,...
    betaP,betaI,betaA,epsilon,chiE,chiP,chiI,chiA,...
    mu,deltaE,deltaP,sigma,gammaI,eta,alpha,gammaA,...
    N,MS,ME,MP,MI,MA,MR,n)

    % some useful matrices
    u=ones(1,n); % all-one row vector of size n
    u1n=1:n; % 1-to-n row vector
    U=sparse(u1n,u1n,u,n,n); % identity matrix of size n
    
    % state variables
    S=x(1:n); 
    E=x(n+1:2*n);
    P=x(2*n+1:3*n); 
    I=x(3*n+1:4*n);
    A=x(4*n+1:5*n); 
    Eq=x(5*n+1:6*n);
    Pq=x(6*n+1:7*n); 
    Ih=x(7*n+1:8*n);
    Aq=x(8*n+1:9*n);
    R=x(9*n+1:10*n);
    
    % total infection rate    
    Sdiag=sparse(u1n,u1n,S,n,n); Ediag=sparse(u1n,u1n,E,n,n);
    Pdiag=sparse(u1n,u1n,P,n,n); Idiag=sparse(u1n,u1n,I,n,n);
    Adiag=sparse(u1n,u1n,A,n,n); Rdiag=sparse(u1n,u1n,R,n,n);
    DeltaInv=sparse(u1n,u1n,1./(u*(Sdiag*MS+Ediag*ME+Pdiag*MP+Idiag*MI+Adiag*MA+Rdiag*MR)),n,n);
    epsilondiag=sparse(u1n,u1n,epsilon,n,n);
    FoI_P=MS*DeltaInv*(U-epsilondiag)*betaP*MP'*P;
    FoI_I=MS*DeltaInv*(U-epsilondiag)*betaI*MI'*I;
    FoI_A=MS*DeltaInv*(U-epsilondiag)*betaA*MA'*A;
    TIR=(FoI_P+FoI_I+FoI_A).*S;
	
    % ODEs
    dSdt=mu*(N-S)-TIR;
    dEdt=TIR-(mu+deltaE+chiE).*E;
    dPdt=deltaE*E-(mu+deltaP+chiP).*P;
    dIdt=sigma*deltaP*P-(mu+alpha+gammaI+eta+chiI).*I;
    dAdt=(1-sigma)*deltaP*P-(mu+gammaA+chiA).*A;
    dEqdt=chiE.*E-(mu+deltaE)*Eq;
    dPqdt=chiP.*P+deltaE*Eq-(mu+deltaP)*Pq;
    dIhdt=(eta+chiI).*I+sigma*deltaP*Pq-(mu+alpha+gammaI)*Ih;
    dAqdt=chiA.*A+(1-sigma)*deltaP*Pq-(mu+gammaA)*Aq;
    dRdt=gammaI*(I+Ih)+gammaA*(A+Aq)-mu*R;
    dCdt=TIR;
    
    % output
    dxdt=[dSdt; dEdt; dPdt; dIdt; dAdt; dEqdt; dPqdt; dIhdt; dAqdt; dRdt; dCdt];
end