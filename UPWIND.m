function [ u,G ] = UPWIND( u,C,n,N,dt,dx,NTime,Sigma,Beta )

while( n <= NTime )
    for i = 2:N
        u(n+1,i) = u(n,i)-(dt/dx)*((0.5*(C+abs(C)))*(u(n,i)-u(n,i-1))+(0.5*(C-abs(C)))*(u(n,i+1)-u(n,i)));
    end
    u(n+1,1) = u(n,1)-(dt/dx)*((0.5*(C+abs(C)))*(u(n,1)-u(n,N+1))+(0.5*(C-abs(C)))*(u(n,2)-u(n,1)));
    u(n+1,N+1) = u(n,N+1)-(dt/dx)*((0.5*(C+abs(C)))*(u(n,N+1)-u(n,N))+(0.5*(C-abs(C)))*(u(n,1)-u(n,N+1)));
    n = n+1;
end

G = (1.0-Sigma+Sigma*cos(Beta))-1i*Sigma*sin(Beta);

end

