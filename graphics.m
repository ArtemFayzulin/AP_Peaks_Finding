clc
[ir, red, flt_ECG, Fs, ap] = readpwdata('Соколова_Евгения_Андреевна_13-04-22_12-04-29_.bin');
[tMin,tMax] = readlab('Соколова_Евгения_Андреевна_13-04-22_12-04-29_.json');

z = length(ap); %считаем длину сигнала
fragment = 15000; %вручную задаем сколько отсчетов мы хотим исследовать
zaderzhka = z - fragment; % переменная обуславшивающая в дальнейшем индексы, а также начало отсчета для фрагмента

%отсортируем необходимые фрагменты
flt_ECG_fg = flt_ECG(zaderzhka:end);
ap_fg = ap(zaderzhka:end);
ir_fg = ir(zaderzhka:end);
red_fg = red(zaderzhka:end);

tMin = tMin';



%переходим от отсчетов ко времени
Dlina = length(ap_fg);
T = 1/Fs;
tmax = Dlina*T;
t = 0:T:tmax - T;


%Выводим необходимые графики
figure ('Name','Выведенные графики')
subplot(2,2,1)
plot(t,flt_ECG_fg)
set(gca,'XLim', [0 tmax])
title('flt_ECG')

subplot(2,2,2)
plot(t,ap_fg)
set(gca,'XLim', [0 tmax])
title('AP')

subplot(2,2,3)
plot(t,-ir_fg)
set(gca,'XLim', [0 tmax])
title('Инфракрасный')

subplot(2,2,4)
plot(t,-red_fg)
set(gca,'XLim', [0 tmax])
title('Красный')
% конец вывода графиков


%1-й шаг, отфильтруем сигнал 
%Фильтрация ФНЧ
Filted_LPF = zeros(1,Dlina);
for n=31:Dlina
Filted_LPF(n)=ap_fg(n)-2*ap_fg(n-15)+ap_fg(n-30)+2*Filted_LPF(n-1)-Filted_LPF(n-2);
end

%Фильтрация ФВЧ
Filted = zeros(1,Dlina);
for n=775:(Dlina)
Filted(n)= Filted(n-1) - (1/774) * Filted_LPF(n) + Filted_LPF(n-387) - Filted_LPF(n-388) + (1/774)*Filted_LPF(n-774);
end



%Построение и сравнение двух графиков
figure ('Name','Сравнение отфильтрованного и неотфильтрованного')

subplot(3,1,1)
plot(t,ap_fg)
set(gca,'XLim', [0 tmax])
title('AP')



%использованеи функции SSF
w = 128; %поскольку 64 для Fd 512, то для 1000 возьмем 128
threshold = 0;
SSF = zeros(1,Dlina);
for k = 1:Dlina-1
   if (Filted(k+1) - Filted(k)) > 0  
        delta_x = Filted(k+1) - Filted(k);
        SSF(k) = SSF(k) + delta_x;
   elseif Filted(k+1) - Filted(k)<=0
        delta_x = 0;
        SSF(k) = SSF(k) + delta_x;
      
        if mod(k,128)==0
          SSF(k) = 0;
          delta_x = 0;
        end
   end
%нахождение максимума в первые 3 секунды
      if t(k)<=3 
            if SSF(k)>=threshold
               threshold = SSF(k);      
            end
     end

end


%нахождение  пиков
tochka = zeros(1,Dlina);
pos_test_min = zeros (1,Dlina);
pos_test_max = zeros (1,Dlina);
for i=1:Dlina-w
    if (SSF(i)<= 0.7*threshold) && (SSF(i+1)>=0.7*threshold) && (SSF(i)~=0)
                         for a=i:-1:i-w
                             if SSF(a)==0 && SSF(a+1)~=0
                                tochka(a) = Filted(a);
                                pos_test_min(a) = a+zaderzhka;
                             end
                         end


                         for a=i:i+w
                             if SSF(a)==0 && SSF(a-1)~=0
                                tochka(a) = Filted(a); 
                                pos_test_max(a) = a+zaderzhka;
                             end
                         end
     end
end 

pos_test_min(pos_test_min==0) = [];
pos_test_max(pos_test_max==0) = [];


%нахождение Se и P
[Se,P,TP,FP,FN] = calcStat(tMin,pos_test_min,300);


%строим график SSF
subplot(3,1,2)
plot(t,SSF)
hold on
yline(0.7*threshold,'--')
set(gca,'XLim', [0 tmax])
title('SSF')

%построим график отфильтрованного сигнала с пиками на нем
subplot(3,1,3)
plot(t,Filted)
hold on
plot(t, tochka,'c*')
set(gca,'XLim', [0 tmax])
title('Filted with peaks')



