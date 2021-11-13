%function y = voight(z, y, const)
%f = @(z,y,t) exp(-t * t)/(y^2 + (z - t)^2);
%I = 0;
%f0 = 0;
%x0 = -2.5;
%dx = 2 * abs(x0) / 400;
%for i = 1:400
%  I = I + const * (f0 + f(z, y, i*dx))*dx/2;
%end
%y = I;
%end

function lines(v1, v2, p, p_air, T)
k = 1.38 * 10^(-16);
Na = 6 * 10^23;
c = 3 * 10^10;
M = 18 / Na;
T0 = 296;
n = p*(10^6) / k / T; %концентрация
file = load('input.txt');
v = file(:,2);
S = file(:,3);
g_air = file(:,5);
g = file(:,6);
n_g = file(:,8);
b = file(:,9);
dopler_const = sqrt(2 * 8.31 * (10^7) * T / M * log(2))/c;
i = 1;
while v(i) < v1
    i = i + 1;
end
j = i;
dVc = [];
dVD = [];
while v(i) <= v2
   n_b = 1.5 * n_g(i) + 0.25; % для температурной зависимости
   g(i) = g(i) * (T0/T)^(n_g(i)); % температурная зависимость g(T)
   b(i) = b(i) * (T0/T)^(n_b);
   dVc(i) = p*g(i) + p_air * g_air(i);
   dVD(i) = v(i) * dopler_const;
   i = i + 1;
end
n1 = i - 1;
f = @(x) 0;
for l = j : n1
   Vc = dVc(l);
   VD = dVD(l);
   v_ = v(l) + b(l) * p_air;
   if Vc >= 100*VD % контур Лоренца
   	c1 = Vc * S(l) * n / pi;
   	I = @(x) c1 / ((Vc)^2 + (x-v_)^2);
   end
   if Vc <= 0.01*VD % контур Доплера
        I = @(x) 1/VD/sqrt(pi) * exp(-((x-v_)^2) / VD^2)*S(l)*n;
   end
   if Vc < 100*VD && Vc>0.01*VD
   	y=Vc*sqrt(log(2))/VD;
        const1 = S(l)*n*y/pi;
        const2 = sqrt(log(2))/VD;
        I = @(x) voigt((x-v_)*const2, y, const1); 
    end
        f1 = @(x) f(x)+I(x); % сумма спектральных линий
        f = f1;
    end
    x = linspace(v1, v2, 1000);
    y = linspace(v1, v2, 1000);
    for i = 1:1000
       y(i)=f(x(i));
    end
plot(x,y);
save('C:\graphs\x.txt','x','-ascii');
save('C:\graphs\y.txt','y','-ascii');
end
