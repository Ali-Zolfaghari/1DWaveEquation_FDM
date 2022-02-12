function [ u,G ] = LAXWENDROFF( u,C,n,N,dt,dx,NTime,Sigma,Beta )

while( n <= NTime )
    for i = 2:N
        u(n+1,i) = u(n,i)-(Sigma/2)*(u(n,i+1)-u(n,i-1))+((Sigma^2)/2)*(u(n,i+1)-2*u(n,i)+u(n,i-1));
    end
    u(n+1,1) = u(n,1)-(Sigma/2)*(u(n,2)-u(n,N+1))+((Sigma^2)/2)*(u(n,2)-2*u(n,1)+u(n,N+1));
    u(n+1,N+1) = u(n,N+1)-(Sigma/2)*(u(n,1)-u(n,N))+((Sigma^2)/2)*(u(n,1)-2*u(n,N+1)+u(n,N));
    n = n+1;
end

G = 1.0+(Sigma^2.0)*(cos(Beta)-1.0)-1i*Sigma*sin(Beta);

end

