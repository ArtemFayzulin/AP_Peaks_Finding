#ECG_working
using DSP
using StatsBase
function Pan_Tompkins(ECG_fg,fs)

    """"
    FILTERING ECG SIGNAL
    Аналогично сигналу ФПГ используется каскадный фильтр, который мы
    задаем следующими уравнениями
    Для ФНЧ
    y(nT) = 2y(nT - T) - y(nT - 2T) + 1/32*x(nT) - 1/16*x(nT - 6T) + 1/32*x(nT - 12T)
    Задаваемая задержка фильтром составляет 5 отсчетов или 25 мс при частоте 200 Гц
    Для ФВЧ
    y(nT) = y(nT - T) - 1/32*x(nT ) + x(nT- 16T) - x(nT - 17T) + 1/32*x(nT - 32T)
    Фильтр дает задержку в 80 мс
    Необходимо провести рассчет задержки для другой частоты
    по следующей формуле пропорции
    5/200 = x/1000
    16/200 = x/1000
    
    zaderzhka = (5*fs/200) + (16*fs/200) # рассчет производится в отсчетах (первое слагаемое - для ФНЧ, а второе для ФВЧ)
    
    Однако в данной работе воспользуемся библиотекой DSP, позволяющей оформлять
    цифровую фильтрацию, запишем фильтры со следующими передаточными функциями
    
    H(z) = (1 - z^-6)^2 / (1 - z^-1)^2
    H(z) = (-1/32 + z^-16 - z^-17 + 1/32 * z^-32) / (1 - z^-1)
    """
    b = fill(0,7)
    b[1],b[7] = -1,1
    
    a = [-1,1]
    H_LowPass = PolynomialRatio(b,a)^2
    mydata = copy(ECG_fg)
    
    Filted_LPF = reverse(filt(H_LowPass,reverse(filt(H_LowPass,mydata))))
    
    #numerator_coefs = coefb(H_LowPass)
    #denominator_coefs = coefa(H_LowPass)
    #Filted_LPF = filtfilt(numerator_coefs,mydata)
    
    b = fill(0.0,33)
    b[1],b[17],b[18],b[33] = 1/32,-1,1,-1/32
    H_HighPass = PolynomialRatio(b,a)
    
    Filted_ECG = reverse(filt(H_HighPass,reverse(filt(H_HighPass,Filted_LPF))))

    Filted_ECG  = Filted_ECG[50:end-50] # поскольку выброс
    ECG_fg = ECG_fg[50:end-50] # поскольку не рассматриваем предыдуший фрагмент
    
    Dlina = length(ECG_fg)
    T = 1/fs;
    tmax = length(ECG_fg)*T;
    t = 0:T:tmax - T;
    
    
    ##SQUARING and DERIVATIVE
    # данном блоке мы будем возводить в квадрат данную функцию, чтобь
    #точно учесть то, что мы избавляемся от ложных срабатываний, такие как зубец t
    #По факту мы усиливаем наклон кривой частной характеристики, полученной на шаге производной
    #Фильтр производной выглядит следующим образом: y(n) = 1/8 [2x(n)+x(n-1-x(n-3)-2x(n-4))]
   
    coeff = [1/4,1/8,0,-1/8,-1/4]
    Derivative = filtfilt(coeff,Filted_ECG).^2
    
    
    #MOVING WINDOW INTEGRATOR
    #Получаем информацию о наклоне и ширине комплекса QRS. Размер окна вычисляется, 
    #как 0.15 * частота дискретизации, для более точных результатов 
    #y(n) = [x(n - (N-1)) + x(n - (N-2))+...+ x(n)]/N
    
    window_size = Int(fld(0.15,T))
    Complexes = fill(0.0,length(ECG_fg))
    sum = 0.0
    for n in eachindex(Complexes)
        if n<=window_size
            sum += Derivative[n]/window_size
            Complexes[n] = sum
        elseif n>window_size
            sum += (Derivative[n] - Derivative[n-window_size])/window_size
            Complexes[n] = sum
        end
    end
    
  
    #В данном разделе мы будем отмечать непосредственно точки
    pos_ECG_R = [];
    


    TheSig = maximum(Filted_ECG[1:2*Int(fs)])*(0.25)
    TheNoise = mean(Filted_ECG[1:2*Int(fs)])*(0.5)
    Thereshold = TheNoise + 0.25 * (TheSig - TheNoise)
    
    SignalLine = fill(0.0,length(Filted_ECG))
    NoiseLine = fill(0.0,length(Filted_ECG))
    TheresholdLine = fill(0.0,length(Filted_ECG))
    R_zub = Array{Any}(missing, Dlina);
    R_zub2 = Array{Any}(missing, Dlina);

    limit = 0
    
    for i = 2:length(Filted_ECG)-1
        if (Filted_ECG[i]>Filted_ECG[i+1]) && (Filted_ECG[i]>Filted_ECG[i-1])
            if Filted_ECG[i]<Thereshold && Filted_ECG[i]>TheNoise
                TheNoise = 0.125 * Filted_ECG[i] + 0.850 * TheNoise #was 0.875, we replace in case to reduce noise and find low peaks ECG
            elseif Filted_ECG[i]>=Thereshold 
                TheSig = 0.125 * Filted_ECG[i] + 0.875 * TheSig
                if limit<=0 && Filted_ECG[i]>TheSig
                    R_zub[i] = Filted_ECG[i]
                    temp = 0
                    temp2 = 0
                    for a in -10:10
                        if ECG_fg[i+a]>temp
                            temp = ECG_fg[i+a]
                            temp2 = a
                        end
                    end
                    R_zub2[i+temp2] = temp
                    push!(pos_ECG_R,(i+temp2)*T);
                    limit = 200/(T*1000)
                end
            end
        end
        Thereshold = TheNoise + 0.25 * (TheSig - TheNoise)
        NoiseLine[i] = TheNoise
        SignalLine[i] = TheSig
        TheresholdLine[i] = Thereshold
        limit = limit - 1
    end
    return(Filted_ECG,ECG_fg,R_zub,R_zub2,t,pos_ECG_R)
end