function [ u ] = DAMPED_BEAMWARMING( u,C,n,N,dt,dx,NTime,Sigma,Beta,Epsilon_e,Epsilon_i )

MTX = zeros(N+1,N+1);
for i = 2:N
    MTX(i,i) = 1.0+2.0*Epsilon_i;
    MTX(i,i+1) = 0.25*Sigma-Epsilon_i;
    MTX(i,i-1) = -0.25*Sigma-Epsilon_i;
end
MTX(N+1,N+1) = 1.0+2.0*Epsilon_i;
MTX(N+1,1) = 0.25*Sigma-Epsilon_i;
MTX(N+1,N) = -0.25*Sigma-Epsilon_i;
MTX(1,1) = 1.0+2.0*Epsilon_i;
MTX(1,2) = 0.25*Sigma-Epsilon_i;
MTX(1,N+1) = -0.25*Sigma-Epsilon_i;

while( n <= NTime )
    for i = 3:N-1
        RHS(i) = u(n,i)-0.25*Sigma*(u(n,i+1)-u(n,i-1))-Epsilon_e*(u(n,i-2)-4.0*u(n,i-1)+6.0*u(n,i)-4.0*u(n,i+1)+u(n,i+2));
    end
    RHS(N+1) = u(n,N+1)-0.25*Sigma*(u(n,1)-u(n,N))-Epsilon_e*(u(n,N-1)-4.0*u(n,N)+6.0*u(n,N+1)-4.0*u(n,1)+u(n,2));
    RHS(N) = u(n,N)-0.25*Sigma*(u(n,N+1)-u(n,N-1))-Epsilon_e*(u(n,N-2)-4.0*u(n,N-1)+6.0*u(n,N)-4.0*u(n,N+1)+u(n,1));
    RHS(1) = u(n,1)-0.25*Sigma*(u(n,2)-u(n,N+1))-Epsilon_e*(u(n,N)-4.0*u(n,N+1)+6.0*u(n,1)-4.0*u(n,2)+u(n,3));
    RHS(2) = u(n,2)-0.25*Sigma*(u(n,3)-u(n,1))-Epsilon_e*(u(n,N+1)-4.0*u(n,1)+6.0*u(n,2)-4.0*u(n,3)+u(n,4));
    U = MTX\RHS';
    u(n+1,:) = U';
    n = n+1;
end

end

