%{Uklad wielostanowy probabilistyczny 2}%
x=0;
function [Pr,T] = ProbRing2(N,a,b)
    dt=0.3;
    Pr_curr = 1/N + zeros(N,1);
    diff=1;
    t=0;
    while (diff>1e-10)
      Pr_prev = Pr_curr;
      
      h=0.5* dt*(a(N)*Pr_prev(N) + b*Pr_prev(2) - (b+a(1))*Pr_prev(1));
      Pr_curr(1)+= dt*(a(N)*Pr_prev(N) + b*Pr_prev(2) - (b+a(1))*(Pr_prev(1)+h));
      #disp(h)
      h=0.5* dt*(a(N-1)*Pr_prev(N-1) + b*Pr_prev(1) - (b+a(N))*Pr_prev(N));
      Pr_curr(N)+= dt*(a(N-1)*Pr_prev(N-1) + b*Pr_prev(1) - (b+a(N))*(Pr_prev(N)+h));
      
      f=@(i) a(i-1)*Pr_prev(i-1) + b*Pr_prev(i+1) - (b+a(i))*Pr_prev(i);
      
      for i=2:(N-1)
        h=0.5*dt*f(i);
        Pr_curr(i)+=dt*(a(i-1)*Pr_prev(i-1) + b*Pr_prev(i+1) - (b+a(i))*(Pr_prev(i)+h));
      end    
      diff=max(abs(Pr_prev-Pr_curr)./Pr_prev);
      t++;
    end
    Pr=Pr_curr;
    T=t;
end
N=20;
a=hilb(N)(:,1);
b=0.5;
[Pr,t]=ProbRing2(N,a,b);
disp(Pr)
sum(Pr)
disp(t)
plot(1:N, Pr)
hold

b=0.1;
[Pr,t]=ProbRing2(N,a,b);
disp(Pr)
sum(Pr)
disp(t)
plot(1:N, Pr)
