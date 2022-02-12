function [ u,G ] = BEAMWARMING( u,C,n,N,dt,dx,NTime,Sigma,Beta )

MTX = zeros(N+1,N+1);
for i = 2:N
    MTX(i,i) = 1.0;
    MTX(i,i+1) = 0.25*Sigma;
    MTX(i,i-1) = -0.25*Sigma;
end
MTX(N+1,N+1) = 1.0;
MTX(N+1,1) = 0.25*Sigma;
MTX(N+1,N) = -0.25*Sigma;
MTX(1,1) = 1.0;
MTX(1,2) = 0.25*Sigma;
MTX(1,N+1) = -0.25*Sigma;

while( n <= NTime )
    for i = 2:N
        RHS(i) = u(n,i)-0.25*Sigma*(u(n,i+1)-u(n,i-1));
    end
    RHS(N+1) = u(n,N+1)-0.25*Sigma*(u(n,1)-u(n,N));
    RHS(1) = u(n,1)-0.25*Sigma*(u(n,1+1)-u(n,N+1));
    U = MTX\RHS';
    u(n+1,:) = U';
    n = n+1;
end

G = (1.0-1i*0.5*Sigma*sin(Beta))/(1.0+1i*0.5*Sigma*sin(Beta));

end

