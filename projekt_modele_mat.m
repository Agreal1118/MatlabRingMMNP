%{Uklad wielostanowy probabilistyczny 2}%
x=0;
N=20;

% Pierwszy test
disp("Rozpoczêcie pierwszego testu - a = const")
a = 0.4 + zeros(N,1); 
figure(1);                       % Pisanie po wykresie numer 1
hold on
b=0.5 + zeros(N,1);
[Pr,t]=ProbRing2(N,a,b,true);
disp("Wypisanie ostatecznych prawdopodobieñstw")
disp(Pr);                        
disp("Wypisanie liczby iteracji")
disp(t);                         
plot(1:N, Pr, 'k','LineWidth', 2)                   
title(["Test I, a=0.4, liczba N:20, liczba iteracji:", t]);
xlabel('Numer stanu')
ylabel('Prawdopodobieñstwo')

% Drugi test
disp("Rozpoczêcie drugiego testu - a = 1/k")
figure(2);
a=hilb(N);                       % (:,1); tworzy macierz hilberta o wielko¶ci N
a=a(:,1);                        % pierwsza kolumna z a
b=0.5 + zeros(N,1);
[Pr,t]=ProbRing2(N,a,b,true);
disp("Wypisanie ostatecznych prawdopodobieñstw")
disp(Pr);
disp("Wypisanie liczby iteracji")
disp(t);
plot(1:N, Pr, 'k','LineWidth', 2)
title(["Test II, a=1/k, liczba N:20, liczba iteracji:", t]);
xlabel('Numer stanu')
ylabel('Prawdopodobieñstwo')

% Trzeci test
disp("Rozpoczêcie trzeciego testu - a = e^(-a*k)")
figure(3);
a = zeros(N,1);
for i=1:N
    a(i)=exp(-0.5*i);
end
b=0.5 + zeros(N,1);
[Pr,t]=ProbRing2(N,a,b,true);
disp("Wypisanie ostatecznych prawdopodobieñstw")
disp(Pr);
disp("Wypisanie liczby iteracji")
disp(t);
plot(1:N, Pr, 'k','LineWidth', 2)
title(["Test III, a=e^(-a*k), liczba N:20, liczba iteracji:", t]);
xlabel('Numer stanu')
ylabel('Prawdopodobieñstwo')

% test 1
figure(4);
N=2:2:40;
times=zeros(length(N),1);
for i=1:length(N)
  a=0.4 + zeros(N(i),1);
  b=0.5 + zeros(N(i),1);
  [Pr,t]=ProbRing2(N(i),a,b,false);
  times(i)=t;
end
plot(N,times)
title(["Test IV, a=0.4, liczba N:od 2 do 40"]);
xlabel('Liczba stanów')
ylabel('Liczba iteracji')

% test 2
%figure(5);
hold on
N=2:2:40;
times=zeros(length(N),1);
for i=1:length(N)
  a=hilb(N(i));
  a=a(:,1);
  b=0.5 + zeros(N(i),1);
  [Pr,t]=ProbRing2(N(i),a,b,false);
  times(i)=t;
end
plot(N,times)
title(["Test V, a=1/k, liczba N:od 2 do 40"]);
xlabel('Liczba stanów')
ylabel('Liczba iteracji')

% test 3
%figure(6);
hold on
N=2:2:40;
times=zeros(length(N),1);

for i=1:length(N)
    a = zeros(N(i),1);
    for j=1:N(i)
    a(j)=exp(-0.5*j);
    end
  b=0.5 + zeros(N(i),1);
  [Pr,t]=ProbRing2(N(i),a,b,false);
  times(i)=t;
end
plot(N,times)
title(["Wykres zale¿no¶ci iteracji od liczby N, liczba N:od 2 do 40"]);
xlabel('Liczba stanów')
ylabel('Liczba iteracji')

function [Pr,T] = ProbRing2(N,a,b,truefalse)
    dt=1;
    Pr_curr = 1/N + zeros(N,1);                         % pr_curr startuje jako wektor 20 warto¶ci 1/20
    diff=1;
    t=0;
    while (diff>1e-10)                                  % dopóki diff (1) wiêkszy od 1^-10
      Pr_prev = Pr_curr;
      if truefalse==true
        %disp (Pr_curr);
        plot(1:N, Pr_prev);
        hold on
      end
      
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