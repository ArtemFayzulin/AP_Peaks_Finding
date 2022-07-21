%индексы
refdata = [2,15,32,48,66,92,115];
testdata = [1,10,45,67,90,187];

%вызов функции для дельта 1
[~,~,new_massive_TP_ref,new_massive_TP_test,new_massive_FP_test,new_massive_FN_ref] = calcStat(refdata,testdata);



%построение графиков
subplot(3,1,1)
plot(refdata,ones(length(refdata)),'*k')
hold on
plot(testdata,zeros(length(testdata)),'*b')
hold on
plot(new_massive_TP_ref(2,:),new_massive_TP_ref(1,:),'Og')
hold on
plot(new_massive_TP_test(2,:),new_massive_TP_test(1,:),'Og')
hold on
plot(new_massive_FP_test(2,:),new_massive_FP_test(1,:),'Or')
hold on
plot(new_massive_FN_ref(2,:),new_massive_FN_ref(1,:),'Oy')
ylim([-1 2])
xlim([-10 200])
grid on
title('delta = 1')

%вызов для delta = 3
[~,~,new_massive_TP_ref,new_massive_TP_test,new_massive_FP_test,new_massive_FN_ref] = calcStat(refdata,testdata,3);

subplot(3,1,2)
plot(refdata,ones(length(refdata)),'*k')
hold on
plot(testdata,zeros(length(testdata)),'*b')
hold on
plot(new_massive_TP_ref(2,:),new_massive_TP_ref(1,:),'Og')
hold on
plot(new_massive_TP_test(2,:),new_massive_TP_test(1,:),'Og')
hold on
plot(new_massive_FP_test(2,:),new_massive_FP_test(1,:),'Or')
hold on
plot(new_massive_FN_ref(2,:),new_massive_FN_ref(1,:),'Oy')
ylim([-1 2])
xlim([-10 200])
grid on
title('delta = 3')


%вызов для delta = 5
[~,~,new_massive_TP_ref,new_massive_TP_test,new_massive_FP_test,new_massive_FN_ref] = calcStat(refdata,testdata,5);

subplot(3,1,3)
plot(refdata,ones(length(refdata)),'*k')
hold on
plot(testdata,zeros(length(testdata)),'*b')
hold on
plot(new_massive_TP_ref(2,:),new_massive_TP_ref(1,:),'Og')
hold on
plot(new_massive_TP_test(2,:),new_massive_TP_test(1,:),'Og')
hold on
plot(new_massive_FP_test(2,:),new_massive_FP_test(1,:),'Or')
hold on
plot(new_massive_FN_ref(2,:),new_massive_FN_ref(1,:),'Oy')
ylim([-1 2])
xlim([-10 200])
grid on
title('delta = 5')


function [Se,P,new_massive_TP_ref,new_massive_TP_test,new_massive_FP_test,new_massive_FN_ref] = calcStat(pos_ref,pos_test,delta)

% Se - Sensitivity
% P - Positive Predictive Value

% TP - значения, совпадающие на разметках
% FP - значения лишь на тестовой разметке
% FN - значения лишь на референтной разметке

% pos_ref - позиции референтных значений
% pos_test - позиции полученных (тестовых) значений

% xxx_array - это логические массивы для TP FP FN

%проверка входных значений
if nargin<3
   delta = 1;
end

Dlina  = length (pos_ref);
Dlina2 = length (pos_test);

TP_array =  zeros(1,Dlina);
FN_array = TP_array;
FP_array = ones(1,Dlina2);

for i=1:Dlina
    for j=1:Dlina2
      if abs(pos_ref(i)-pos_test(j))<=delta
        TP_array(i) = 1;
      end
    end
end

for i = 1:Dlina
    if TP_array(i) == 0 
        FN_array(i) = 1 ;
    else 
        FN_array(i) = 0 ;
    end
end


   for i=1:Dlina2
      for j=1:Dlina
          if ((pos_test(i)<=pos_ref(j)+delta) && (pos_test(i)>=pos_ref(j)-delta))
              FP_array(i) = 0;
          end
      end
   end


%создание массивов для построения на референтных индексах
new_massive_TP_ref = zeros(2,length(pos_ref));
new_massive_FN_ref = zeros(2,length(pos_ref));
%создание массивов для построения на тестовых индексах
new_massive_TP_test = zeros(2,length(pos_test));
new_massive_FP_test = zeros(2,length(pos_test));

% строим по типу: х - индекс, а y будет равен 1
for i=1:length(pos_ref)
    if TP_array(i) == 1
        new_massive_TP_ref(1,i) = 1;
        new_massive_TP_ref(2,i) = pos_ref(i);
    end
end
% тоже самое для FN 
for i=1:length(pos_ref)
    if TP_array(i) == 0
        new_massive_FN_ref(1,i) = 1;
        new_massive_FN_ref(2,i) = pos_ref(i);
    end
end

% инвертируем амплитуду в 0 и считаем аналогично предыдущему для TP_test 
for i=1:length(pos_test)
    if FP_array(i) == 0
        new_massive_TP_test(1,i) = 1;
        new_massive_TP_test(2,i) = pos_test(i);
    end
end
% аналогично для FP
for i=1:length(pos_test)
    if FP_array(i) == 1
        new_massive_FP_test(1,i) = 1;
        new_massive_FP_test(2,i) = pos_test(i);
    end
end

%убираем лишние строки
new_massive_TP_ref(:, ~any(new_massive_TP_ref)) = [];
new_massive_FN_ref(:, ~any(new_massive_FN_ref)) = [];

new_massive_TP_test(:, ~any(new_massive_TP_test)) = [];
new_massive_FP_test(:, ~any(new_massive_FP_test)) = [];
new_massive_TP_test(1, :) = 0;
new_massive_FP_test(1, :) = 0;


TP = sum(TP_array);
FN = sum (FN_array);
FP = sum (FP_array);

Se = TP/(TP+FN) * 100;
P = TP/(TP+FP) * 100;
end
