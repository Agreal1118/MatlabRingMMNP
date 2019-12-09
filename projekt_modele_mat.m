%{Uklad wielostanowy probabilistyczny 2}%
x=0;
function [Pr,T] = ProbRing2(N,a,b)
    Pr_curr = 1/N + zeros(N,1);
    diff=1;
    t=0;
    while (all(diff>1e-10))
      Pr_prev = Pr_curr;
      Pr_curr(1)=
      Pr_curr(N)=
      f=@(i) a(i-1)*Pr_prev(i-1) + b*Pr_prev(i+1) - (b+a(i-1))*Pr_prev(i);
      for i=2:(N-1)
        h=0.5*
        Pr_curr(i)=
      end    
      diff=abs(Pr_prev-Pr_curr)./Pr_prev;
      t++;
    end
    Pr=Pr_curr;
    T=t;
end