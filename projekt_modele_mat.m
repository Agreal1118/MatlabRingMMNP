%{Uklad wielostanowy probabilistyczny 2}%
x=0;
N=20;

% Pierwszy test
disp("Rozpoczêcie pierwszego testu")
a=hilb(N);                       % (:,1); tworzy macierz hilberta o wielko¶ci N
a=a(:,1);                        % pierwsza kolumna z a
figure(1);                       % Pisanie po wykresie numer 1
hold on
b=0.3 + zeros(N,1);
[Pr,t]=ProbRing2(N,a,b);
disp("Wypisanie ostatecznych prawdopodobieñstw")
disp(Pr);                        
disp("Wypisanie liczby iteracji")
disp(t);                         
plot(1:N, Pr, 'k','LineWidth', 2)                   
title(["Liczba iteracji: ", t]);

% Drugi test
disp("Rozpoczêcie drugiego testu")
figure(2);
a=hilb(N);                       % (:,1); tworzy macierz hilberta o wielko¶ci N
a=a(:,1);                        % pierwsza kolumna z a
b=0.5 + zeros(N,1);
[Pr,t]=ProbRing2(N,a,b);
disp("Wypisanie ostatecznych prawdopodobieñstw")
disp(Pr);
disp("Wypisanie liczby iteracji")
disp(t);
plot(1:N, Pr, 'k','LineWidth', 2)
title(["Liczba iteracji: ", t]);

% Trzeci test
disp("Rozpoczêcie trzeciego testu")
figure(3);
a=0.1 + zeros(N,1);
b=[0.1, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1];
[Pr,t]=ProbRing2(N,a,b);
disp("Wypisanie ostatecznych prawdopodobieñstw")
disp(Pr);
disp("Wypisanie liczby iteracji")
disp(t);
plot(1:N, Pr, 'k','LineWidth', 2)
title(["Liczba iteracji: ", t]);

% Czwarty test 
disp("Rozpoczêcie Czwartego testu")
figure(4);
a=rand(20,1);
b=rand(20,1);
disp(a);
disp(b);
[Pr,t]=ProbRing2(N,a,b);
disp("Wypisanie ostatecznych prawdopodobieñstw")
disp(Pr);
disp("Wypisanie liczby iteracji")
disp(t);
plot(1:N, Pr, 'k','LineWidth', 2)
title(["Liczba iteracji: ", t]);


function [Pr,T] = ProbRing2(N,a,b)
    dt=0.3;
    Pr_curr = 1/N + zeros(N,1);                         % pr_curr startuje jako wektor 20 warto¶ci 1/20
    diff=1;
    t=0;
    while (diff>1e-10)                                  % dopóki diff (1) wiêkszy od 1^-10
      Pr_prev = Pr_curr;                               
      %disp (Pr_curr);
      plot(1:N, Pr_prev);
      hold on
      
      h=0.5* dt*(a(N)*Pr_prev(N) + b(2)*Pr_prev(2) - (b(1)+a(1))*Pr_prev(1));     
      
      Pr_curr(1) = Pr_curr(1) + dt*(a(N)*Pr_prev(N) + b(2)*Pr_prev(2) - (b(1)+a(1))*(Pr_prev(1)+h));
      
      %disp(h)
      h=0.5* dt*(a(N-1)*Pr_prev(N-1) + b(1)*Pr_prev(1) - (b(N)+a(N))*Pr_prev(N));
      
      Pr_curr(N) = Pr_curr(N) + dt*(a(N-1)*Pr_prev(N-1) + b(1)*Pr_prev(1) - (b(N)+a(N))*(Pr_prev(N)+h));
      
      f=@(i) a(i-1)*Pr_prev(i-1) + b(i+1)*Pr_prev(i+1) - (b(i)+a(i))*Pr_prev(i);
      
      for i=2:(N-1)
        h=0.5*dt*f(i);
        Pr_curr(i) = Pr_curr(i) + dt*(a(i-1)*Pr_prev(i-1) + b(i+1)*Pr_prev(i+1) - (b(i)+a(i))*(Pr_prev(i)+h));
      end    
      diff=max(abs(Pr_prev-Pr_curr)./Pr_prev);
      t= t + 1;
    end
    Pr=Pr_curr;
    T=t;
end