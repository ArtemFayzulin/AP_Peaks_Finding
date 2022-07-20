function [Se,P,TP,FP,FN] = calcStat(pos_ref,pos_test,delta)

% Se - Sensitivity
% P - Positive Predictive Value

% TP - значения, совпадающие на разметках
% FP - значения лишь на тестовой разметке
% FN - значения лишь на референтной разметке

% pos_ref - позиции референтных значений
% pos_test - позиции полученных (тестовых) значений

%проверка входных значений
if nargin<3
   delta = 1;
end

Dlina  = length (pos_ref);
Dlina2 = length (pos_test);

TP =  zeros(1,Dlina);
FN = TP;
FP = zeros(1,Dlina2);

for i=1:Dlina
    for j=1:Dlina2
      if abs(pos_ref(i)-pos_test(j))<=delta
        TP(i) = 1;
      end
    end
end

for i = 1:Dlina
    if TP(i) == 0 
        FN(i) = 1 ;
    else 
        FN(i) = 0 ;
    end
end

for i=1:Dlina2
    for j=1:Dlina
      if abs(pos_test(i)-pos_ref(j))>=delta && TP(i)~=1
        FP(i) = 1;
      end
    end
end
 
TP = sum(TP);
FN = sum (FN);
FP = sum (FP);

Se = TP/(TP+FN) * 100;
P = TP/(TP+FP) * 100;
end
