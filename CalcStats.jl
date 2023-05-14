#function Calculation(pos_test::T,pos_ref::T,delta::T1) where {T::Float64,T1::Int64} ???
#не смог определить типы разные типы для переменных функции
function Calculation(pos_test,pos_ref,delta)
  #pos_test - позиции реферетных значений
  #pos_ref - position of referent marks

   Dlina = length(pos_ref);
   Dlina2 = length(pos_test);

   TP_array = zeros(Int8,1,Dlina);
   FP_array = ones(Int8,1,Dlina2);
   FN_array = zeros(Int8,1,Dlina);

   for i in 1:Dlina
       for j in 1:Dlina2
            if abs(pos_ref[i]-pos_test[j])<=delta
              TP_array[i]=1;
            end
       end
   end
  
   for i in 1:Dlina2
       for j in 1:Dlina
           if (pos_test[i]<=pos_ref[j]+delta) && (pos_test[i]>=pos_ref[j]-delta)
              FP_array[i]=0;
           end
       end
    end

    for i in 1:Dlina
        if TP_array[i] == 0 
            FN_array[i] = 1 ;
        else 
            FN_array[i] = 0 ;
        end
    end
    TP = sum(TP_array);
    FN = sum(FN_array);
    FP = sum(FP_array);

   Se=TP/(TP+FN)*100;
   P = TP/(TP+FP)*100;

   println(Se);
   println(P);
   println(TP_array);
   println(FP_array);
   println(FN_array);
end

test = [1,10,45,67,90,187]
ref  = [2,15,32,48,66,92,115]



Calculation(test,ref,1)
