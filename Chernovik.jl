include("readbin.jl")
using DSP
using StatsBase

#filepath = raw"C:\Users\artem\OneDrive\Рабочий стол\Для чего\вуз\Science And Diploma\Програмесы\Сигналы САКР\Мельникова 1\Мельникова_Елизавета_Дмитриевна_21-04-22_11-43-20_.hdr"
#filepath = raw"C:\Users\artem\OneDrive\Рабочий стол\Для чего\вуз\Science And Diploma\Програмесы\Сигналы САКР\Сигналы\Мельникова_Елизавета_Дмитриевна_26-05-22_14-03-23_.hdr"
filepath=raw"C:\Users\artem\OneDrive\Рабочий стол\Для чего\вуз\Science And Diploma\Програмесы\Сигналы САКР\Сигналы\Мельникова_Елизавета_Дмитриевна_26-05-22_14-03-23_.hdr"
named_channels, fs, timestart, units = readbin(filepath)
ap = (named_channels.Ir .- named_channels.Red)./1000# давление
ecg = named_channels.LR ./ 1000 #экг

z = length(ap); #считаем длину сигнала
fragment = 35000; #вручную задаем сколько отсчетов мы хотим исследовать
zaderzhka = z - fragment; # переменная обуславшивающая в дальнейшем индексы, а также начало отсчета для фрагмента

ECG_fg = ecg[zaderzhka:end];
ap_fg = ap[zaderzhka:end];

w = 128
Dlina = length(ap_fg)
T = 1/fs;
tmax = length(ap_fg)*T;
t = 0:T:tmax - T;

plot(t[1:2000],ECG_fg[1:2000])
plot(t[1:2000],ap_fg[1:2000])
################
#ECG SIGNAL working


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
t = t[50:end-50]
#numerator_coefs = coefb(H_HighPass)
#Filted_ECG = filtfilt(numerator_coefs,Filted_LPF)

plot(t, [Filted_ECG,ECG_fg], layout = (2,1))
plot([Filted_ECG[200:2000],ECG_fg[200:2000]],layout = (2,1))
"""""
SQUARING and DERIVATIVE
В данном блоке мы будем возводить в квадрат данную функцию, чтобь
точно учесть то, что мы избавляемся от ложных срабатываний, такие как зубец t
По факту мы усиливаем наклон кривой частной характеристики, полученной на шаге производной
Фильтр производной выглядит следующим образом: y(n) = 1/8 [2x(n)+x(n-1-x(n-3)-2x(n-4))]

"""""
coeff = [1/4,1/8,0,-1/8,-1/4]
Derivative = filtfilt(coeff,Filted_ECG).^2

plot(Derivative[1000:6000])
plot(Filted_ECG[3000:6000])

"""""
MOVING WINDOW INTEGRATOR
Получаем информацию о наклоне и ширине комплекса QRS. Размер окна вычисляется, 
как 0.15 * частота дискретизации, для более точных результатов 


y(n) = [x(n - (N-1)) + x(n - (N-2))+...+ x(n)]/N
"""""
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

plot(Complexes)
plot([Complexes[200:3000],ECG_fg[200:3000]],layout = (2,1))
issleduem = 2000:6000


Derivative

plot(t[issleduem],[ECG_fg[issleduem],Filted_ECG[issleduem],Complexes[issleduem],Derivative[issleduem]],layout = (2,2))
#savefig("4_ECG.png")
plot(t[issleduem],[Complexes[issleduem],Derivative[issleduem]])

""""
В данном разделе мы будем отмечать непосредственно точки
"""""

function indmax(x)
    z = findmax(x)[2]
    return z
end



TheSig = maximum(Filted_ECG[1:2*Int(fs)])*(0.5)
TheNoise = mean(Filted_ECG[1:2*Int(fs)])*(0.5)
Thereshold = TheNoise + 0.25 * (TheSig - TheNoise)

SignalLine = fill(0.0,length(Filted_ECG))
SignalLine[1]


NoiseLine = fill(0.0,length(Filted_ECG))
TheresholdLine = fill(0.0,length(Filted_ECG))
R_zub = Array{Any}(missing, Dlina);
R_zub2 = Array{Any}(missing, Dlina);
limit = 0

for i = 2:length(Filted_ECG)-Int(2*fs)
    if (Filted_ECG[i]>Filted_ECG[i+1]) && (Filted_ECG[i]>Filted_ECG[i-1])
        if Filted_ECG[i]<Thereshold 
            TheNoise = 0.125 * Filted_ECG[i] + 0.875 * TheNoise
        elseif Filted_ECG[i]>Thereshold 
            TheSig = 0.125 * Filted_ECG[i] + 0.875 * TheSig
            Thereshold = TheNoise + 0.25 * (TheSig - TheNoise)
            if limit<=0 && Filted_ECG[i]>=TheSig && (Filted_ECG[i]>=maximum(Filted_ECG[i:i+Int(2*fs)])*(0.5))
                R_zub[i] = Filted_ECG[i]
                temp = 0
                temp2 = 0
                for a in -20:20
                    if ECG_fg[i+a]>temp
                        temp = ECG_fg[i+a]
                        temp2 = a
                    end
                end
                R_zub2[i+temp2] = temp
                limit = 0.5*fs
            end
        end
    end
    NoiseLine[i] = TheNoise
    SignalLine[i] = TheSig
    TheresholdLine[i] = Thereshold
    limit = limit - 1
end
gran  = 1:length(Filted_ECG)
plot(Filted_ECG[gran])
plot!(NoiseLine[gran])
plot!(SignalLine[gran])
plot!(TheresholdLine[gran])
scatter!(R_zub[gran])
plot(ECG_fg[gran])
scatter!(R_zub2[gran])
#savefig("POISK_TOCHEK.png")   

plot(Filted_ECG)
scatter!(R_zub)

plot(t,ECG_fg)
scatter!(t,R_zub2)




plot(Derivative[1000:1200])
plot([Derivative[1000:1110],Filted_ECG[1000:1110]],layout=(2,1))


TheSig = maximum(Derivative[1:2*Int(fs)])*(0.3)
TheNoise = mean(Derivative[1:2*Int(fs)])*(0.5)
Thereshold = TheNoise + 0.25 * (TheSig - TheNoise)
limit = 0

SignalLine = fill(0.0,length(Filted_ECG))
NoiseLine = fill(0.0,length(Filted_ECG))
TheresholdLine = fill(0.0,length(Filted_ECG))
R_zub = Array{Any}(missing, Dlina);
R_zub2 = Array{Any}(missing, Dlina);

for i = 2:length(Derivative)-1
    if (Derivative[i]>Derivative[i+1]) && (Derivative[i]>Derivative[i-1])
        if Derivative[i]<Thereshold && Derivative[i]>TheNoise
            TheNoise = 0.125 * Derivative[i] + 0.875 * TheNoise
        elseif Derivative[i]>Thereshold 
            TheSig = 0.125 * Derivative[i] + 0.875 * TheSig
            Thereshold = TheNoise + 0.25 * (TheSig - TheNoise)
            if limit<=0 && Derivative[i]>TheSig && Filted_ECG[i]>0 && (Derivative[i]>maximum(Derivative[i:i+400])*(0.3))
                R_zub[i] = Filted_ECG[i]
                temp = ECG_fg[i]
                temp2 = 0
                for a in -15:15
                    if ECG_fg[i+a]>temp
                        temp = ECG_fg[i+a]
                        temp2 = a
                    end
                end
                R_zub2[i+temp2] = temp
                limit = 0.6*fs
            end
        end
    end
    NoiseLine[i] = TheNoise
    SignalLine[i] = TheSig
    TheresholdLine[i] = Thereshold
    limit = limit - 1
end

gran  = 1:600
plot([Derivative[gran],Filted_ECG[gran]],layout=(2,1))
plot!(NoiseLine[gran])
plot!(SignalLine[gran])
plot!(TheresholdLine[gran])
scatter!(R_zub[gran])
#savefig("POISK_TOCHEK.png")   
gran1 = 1:600
plot([Filted_ECG[gran1],ECG_fg[gran1]],layout = (2,1)) 
scatter!(R_zub[gran1])  

x = plot(t,Filted_ECG)
    scatter!(t,R_zub)
y = plot(t,ECG_fg)
    scatter!(t,R_zub2)
plot(ECG_fg[gran1])
scatter!(R_zub2[gran1])


plot(t[5000:5200],ECG_fg[5000:5200])
scatter!(t[5000:5200],R_zub2[5000:5200], mc=:red, ms=2, ma=0.5)

plot(x,y,layout = (2,1))
#savefig("PROSMOTR.png")

########

#LowPassFilter 
Filted_LPF = fill(0.0,length(ap_fg))
for n in 31:Dlina
    Filted_LPF[n] = ap_fg[n]-(2*ap_fg[n-15])+ap_fg[n-30]+(2*Filted_LPF[n-1])-Filted_LPF[n-2]
end

#HighPassFilter

Filted = fill(0.0,length(Filted_LPF));
for n=775:Dlina
    Filted[n]= Filted[n-1] - (1/774) * Filted_LPF[n] + Filted_LPF[n-387] - Filted_LPF[n-388] + (1/774)*Filted_LPF[n-774];
end

plot(Filted)

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

plot(Filted[50:length(Filt)])
Dlina = length(Filted)
threshold = maximum(SSF[2*Int(fs):3*Int(fs)])
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
end

plot(t,SSF)
plot(SSF[2000:end])
pos_test_min = [];
pos_test_max = [];
tochka_Flt = Array{Any}(missing, Dlina);
tochka = Array{Any}(missing, Dlina);
zaderzhka_Flt = 387
limit = 1
for i=1:Dlina-w
   if (SSF[i]<= 0.7*threshold) && (SSF[i+1]>=0.7*threshold) && (SSF[i]!=0)
                        for a=i:-1:i-w
                            if SSF[a]==0 && SSF[a+1]!=0 && limit!=0
                               (tochka[a-zaderzhka_Flt]=ap_fg[a-zaderzhka_Flt]);
                               (tochka_Flt[a]=Filted[a]);
                               push!(pos_test_min,a*T);
                               limit = 0
                            end
                        end
                        limit = 1


                        for a=i:i+w
                            if SSF[a]==0 && SSF[a-1]!=0 && limit!=0                              
                               (tochka[a-zaderzhka_Flt]=ap_fg[a-zaderzhka_Flt]);
                                (tochka_Flt[a]=Filted[a]);
                               push!(pos_test_max,a*T);
                               limit = 0
                            end
                        end
                        limit = 1
    end
end 


watch = 2000:8000
plot(t[watch],Filted[watch],label="Filted")
scatter!(t[watch],tochka_Flt[watch], label="data", mc=:red, ms=2, ma=0.5)
#savefig("graphics_Flt.png")

plot(t,ap_fg,label="Original")
scatter!(t,tochka, label="data", mc=:red, ms=2, ma=0.5)


plot(t,[Filted,ap_fg],layout = (2,1))
scatter!(t,tochka_Flt)

print(pos_test_min)








""" WORK MOMENTS
"""

x = [1,2,3,4,5,6,7,8,9,10,11]

for i in 1:length(x)-4
    x[i] = x[i+4]
end


for i in 1:length(Filted)-400
    Filted[i] = Filted[i+400]
end

i=1
while i != 2  
    for n in 31:Dlina
        Filted_LPF[n] = ap_fg[n]-(2*ap_fg[n-15])+ap_fg[n-30]+(2*Filted_LPF[n-1])-Filted_LPF[n-2]
    end
    for n=775:Dlina
        Filted[n]= Filted[n-1] - (1/774) * Filted_LPF[n] + Filted_LPF[n-387] - Filted_LPF[n-388] + (1/774)*Filted_LPF[n-774];
   end
   i+=1
end

function reverse(x::Array)
    temp = 0
    for i in 1:Int(floor(length(x)/2))
        temp =  x[i] 
        x[i] = x[length(x)-i+1]
        x[length(x)-i+1] = temp
    end
    return x   
end


print(reverse(x))

print(x)



Filted_LPF1 = fill(0.0,length(ap_fg))
Filted1 = fill(0.0,length(Filted_LPF1));

for n in 31:Dlina
    Filted_LPF1[n] = ap_fg[n]-(2*ap_fg[n-15])+ap_fg[n-30]+(2*Filted_LPF1[n-1])-Filted_LPF1[n-2]
end
for n=775:Dlina
    Filted1[n]= Filted1[n-1] - (1/774) * Filted_LPF1[n] + Filted_LPF1[n-387] - Filted_LPF1[n-388] + (1/774)*Filted_LPF1[n-774];
end

FiltFilt1 = reverse(Filted1)

Filted_LPF = fill(0.0,length(ap_fg))
Filted = fill(0.0,length(Filted_LPF));

for n in 31:Dlina
    Filted_LPF[n] = FiltFilt1[n]-(2*FiltFilt1[n-15])+FiltFilt1[n-30]+(2*Filted_LPF[n-1])-Filted_LPF[n-2]
end
for n=775:Dlina
    Filted[n]= Filted[n-1] - (1/774) * Filted_LPF[n] + Filted_LPF[n-387] - Filted_LPF[n-388] + (1/774)*Filted_LPF[n-774];
end

Filted = reverse(Filted)

plot(t,[Filted1,ap_fg],layout = (2,1))
plot(t,[Filted,ap_fg],layout = (2,1))

"""Filtering
Реализуем фильтрацию путем ПолиномиальногоRatio
"""
b = fill(0,16)
b[1],b[16] = -1,1

a = [-1,1]
H = PolynomialRatio(b,a)
H_LowPass = H^2
mydata = copy(ap_fg)
Filted_LPF1 = filt(H_LowPass,mydata)
Filted_LPF1 = reverse(Filted_LPF1)
Filted_LPF2 = filt(H_LowPass,Filted_LPF1)

plot(Filted_LPF2)
Filted_LPF2=reverse(Filted_LPF2)

"""
b = fill(0,388)
b[388] = 1
H_HighPass = PolynomialRatio(b,[1])
Filted = filt(H_HighPass,Filted_LPF)
"""

Filted = fill(0.0,length(Filted_LPF));
for n=775:Dlina
    Filted[n]= Filted[n-1] - (1/774) * Filted_LPF[n] + Filted_LPF[n-387] - Filted_LPF[n-388] + (1/774)*Filted_LPF[n-774];
end

plot(t,Filted_LPF)
plot(t,Filted)

plot(mydata)




"
Parametrs solving
"
new_hor = 700:2000
x = plot(t[new_hor],ECG_fg[new_hor])
scatter!(t[new_hor],R_zub[new_hor])
y = plot(t[new_hor],ap_fg[new_hor])
scatter!(t[new_hor],tochka[new_hor])
plot(x,y,layout = (2,1))

print(pos_PPG_min)

print(pos_ECG_R)

struct Complexe 
    name::String
    ECG_R::Float64
    time_PPG_min :: Float64
    time_PPG_max :: Float64
    PPG_min::Float64
    PPG_max::Float64
end

all_Complexes = Vector{Complexe}()

# определим комплексы
for i in eachindex(pos_PPG_min)
    for a in eachindex(pos_PPG_min)
        if ((pos_PPG_min[i] * fs - pos_ECG_R[a]*fs<=700)  && (pos_PPG_min[i]>pos_ECG_R[a]))
            ECG_POSMENT = pos_ECG_R[a]
            push!(all_Complexes,Complexe("Комлпекс $i",ECG_POSMENT,pos_PPG_min[i],pos_PPG_max[i],ap_fg[Int(round(pos_PPG_min[i]*fs))],ap_fg[Int(round(pos_PPG_max[i]*fs))]))
        end
    end
end

print(all_Complexes)
NameofComplex = []
tempRR = []
tempAmplitude1 = []
tempAmplitude2 = []
tempPosPPG_min = []
tempPosPPG_max = []


for point in all_Complexes
    push!(NameofComplex,point.name)
    push!(tempRR,point.ECG_R)
    push!(tempAmplitude1,point.PPG_min)
    push!(tempAmplitude2,point.PPG_max)
    push!(tempPosPPG_min,point.time_PPG_min)
    push!(tempPosPPG_max,point.time_PPG_max)
end


ComplexFrames = DataFrame(Наименование_комплекса = NameofComplex,
                         Время_R_зубца = tempRR,
                         Время_offset = tempPosPPG_min,
                         Время_onset = tempPosPPG_max, 
                         Амплитуда_offset = tempAmplitude1,
                         Амплитуда_onset = tempAmplitude2,
)


print(length(tempAmplitude1)," ", length(tempAmplitude2)," ", length(tempRR))

struct Parametrs
    duration :: Float64
    Amplitude :: Float64
    Weidth :: Float64
    hypotenuse :: Float64
    tangens :: Float64
    time_between_R_and_max :: Float64
    time_between_R_and_min :: Float64
end
ParamsOfSystem = Vector{Parametrs}()

for i in 1:length(tempRR)-1
    push!(ParamsOfSystem,Parametrs(tempPosPPG_max[i]-tempPosPPG_min[i],tempAmplitude2[i]-tempAmplitude1[i],tempPosPPG_min[i+1]-tempPosPPG_min[i],sqrt((tempAmplitude2[i]-tempAmplitude1[i])^2+(tempPosPPG_max[i]-tempPosPPG_min[i])^2),(tempAmplitude2[i]-tempAmplitude1[i])/(tempPosPPG_max[i]-tempPosPPG_min[i]),abs(tempRR[i]-tempPosPPG_max[i]),abs(tempRR[i]-tempPosPPG_min[i])))
end

print(ParamsOfSystem[2])

print((ParamsOfSystem))


tempduration = []
tempAmplitude = []
temphypotenuse = []
temptangens = []
temp_time_between_R_and_max = []
temp_time_between_R_and_min = []

for point in ParamsOfSystem
    push!(tempduration,point.duration)
    push!(tempAmplitude,point.Amplitude)
    push!(temphypotenuse,point.hypotenuse)
    push!(temptangens,point.tangens)
    push!(temp_time_between_R_and_max,point.time_between_R_and_max)
    push!(temp_time_between_R_and_min,point.time_between_R_and_min)
end

