%{Uklad wielostanowy probabilistyczny 2}%
x=0;
function [Pr,T] = ProbRing2(N,a,b)
    Pr_curr = 1/N + zeros(N,1);
    diff=1;
    t=0;
    while (diff>1e-10)
      Pr_prev = Pr_curr;
      
      h=0.5* (a(N)*Pr_prev(N) + b*Pr_prev(2) - (b+a(1))*Pr_prev(1));
      Pr_curr(1)+= a(N)*Pr_prev(N) + b*Pr_prev(2) - (b+a(1))*(Pr_prev(1)+h);
      
      h=0.5* (a(N-1)*Pr_prev(N-1) + b*Pr_prev(1) - (b+a(N))*Pr_prev(N));
      Pr_curr(N)+= a(N-1)*Pr_prev(N-1) + b*Pr_prev(1) - (b+a(N))*(Pr_prev(N)+h);
      
      f=@(i) a(i-1)*Pr_prev(i-1) + b*Pr_prev(i+1) - (b+a(i))*Pr_prev(i);
      
      for i=2:(N-1)
        h=0.5*f(i);
        Pr_curr(i)+=a(i-1)*Pr_prev(i-1) + b*Pr_prev(i+1) - (b+a(i))*(Pr_prev(i)+h);
      end    
      diff=max(abs(Pr_prev-Pr_curr)./Pr_prev);
      t++;
    end
    Pr=Pr_curr;
    T=t;
end
