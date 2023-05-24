#work with AD
using DSP

function AD(ap_fg,fs,w)
   ap_fg = ap_fg[50:end-50]
   Dlina = length(ap_fg)
   T = 1/fs;
   tmax = length(ap_fg)*T;
   t = 0:T:tmax - T;

  Filted_LPF = fill(0.0,length(ap_fg))
  for n in 31:Dlina
    Filted_LPF[n] = ap_fg[n]-(2*ap_fg[n-15])+ap_fg[n-30]+(2*Filted_LPF[n-1])-Filted_LPF[n-2]
  end


 #Фильтрация ФВЧ
  b = fill(0,16)
  b[1],b[16] = -1,1   
  a = [-1,1]
  H = PolynomialRatio(b,a)
  H_LowPass = H^2
  mydata = copy(ap_fg)
  Filted_LPF1 = filt(H_LowPass,mydata)
  Filted_LPF1 = reverse(Filted_LPF1)
  Filted_LPF = reverse(filt(H_LowPass,Filted_LPF1))
 
  Filted = fill(0.0,length(Filted_LPF));
  for n=775:Dlina
     Filted[n]= Filted[n-1] - (1/774) * Filted_LPF[n] + Filted_LPF[n-387] - Filted_LPF[n-388] + (1/774)*Filted_LPF[n-774];
  end


 #использованеи функции SSF
 threshold = 0;
 SSF = fill(0.0,Dlina);
 for k in 1:Dlina-1
    if (Filted[k+1] - Filted[k]) > 0  
         delta_x = Filted[k+1] - Filted[k];
         SSF[k] = SSF[k] + delta_x;
    elseif Filted[k+1] - Filted[k]<=0
         delta_x = 0;
         SSF[k] = SSF[k] + delta_x;
       
         if mod(k,w)==0
           SSF[k] = 0;
           delta_x = 0;
         end
    end
 #нахождение максимума в первые 3 секунды
 if t[k]<=3 
     if SSF[k]>=threshold
        threshold = SSF[k];      
     end
 end

 end


 #нахождение  пиков
 # tochka = fill(0.0,Dlina)
 # pos_test_min = fill(0.0,Dlina)
 # pos_test_max = fill(0.0,Dlina)

 pos_test_min = [];
 pos_test_max = [];
 tochka_Flt = Array{Any}(missing, Dlina);
 tochka = Array{Any}(missing, Dlina);
 zaderzhka_Flt = 387 #задержка вносимая фильтром
 for i=1:Dlina-w
   if (SSF[i]<= 0.8*threshold) && (SSF[i+1]>=0.8*threshold) && (SSF[i]!=0)
                        for a=i:-1:i-w
                            if SSF[a]==0 && SSF[a+1]!=0 
                               (tochka[a-zaderzhka_Flt]=ap_fg[a-zaderzhka_Flt]);
                               (tochka_Flt[a]=Filted[a]);
                               push!(pos_test_min,(a-zaderzhka_Flt)*T);
                            end
                        end


                        for a=i:i+w
                            if SSF[a]==0 && SSF[a-1]!=0                               
                               (tochka[a-zaderzhka_Flt]=ap_fg[a-zaderzhka_Flt]);
                                (tochka_Flt[a]=Filted[a]);
                               push!(pos_test_max,(a-zaderzhka_Flt)*T);
                            end
                        end
    end
end 

  #t_zader = t.- (zaderzhka_Flt*T)
 return tochka,tochka_Flt,Filted,ap_fg,pos_test_max,pos_test_min
 
end

#pos_test_min(pos_test_min==0) = [];
#pos_test_max(pos_test_max==0) = [];


#нахождение Se и P
#[Se,P,TP,FP,FN] = calcStat(tMin,pos_test_min,300);


#строим график SSF
