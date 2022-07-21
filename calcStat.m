function [Se,P,TP,FP,FN,TP_array,FP_array,FN_array] = calcStat(pos_ref,pos_test,delta)

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




TP = sum(TP_array);
FN = sum (FN_array);
FP = sum (FP_array);

Se = TP/(TP+FN) * 100;
P = TP/(TP+FP) * 100;
end
