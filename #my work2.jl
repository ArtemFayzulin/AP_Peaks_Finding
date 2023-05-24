#my work
include("readbin.jl")
include("#PPG working.jl")
include("#ECG_working.jl")
include("CalcStats.jl")
using DataFrames
using StatsPlots


filepath = raw"C:\Users\artem\OneDrive\Рабочий стол\Для чего\вуз\Science And Diploma\Програмесы\Сигналы САКР\Мельникова 1\Мельникова_Елизавета_Дмитриевна_21-04-22_11-43-20_.hdr"
named_channels, fs, timestart, units = readbin(filepath)


ap = (named_channels.Ir) #./ 1000 # давление
ecg = named_channels.LR ./ 1000 #экг
ad = named_channels.FPrsNorm./1000  # AD - FPrs1, FPrsNorm
"""
Первый вариант считывания сигнала
"""
#istart = 5000
#iend = 120000

#ECG_fg = ecg[istart:iend];
#ap_fg = ap[istart:iend];
#ad_fg = ad[istart:iend];
"""
Второй вариант рассчитывания сигнала
"""
z = length(ad);#считаем длину сигнала
fragment = 25000; #вручную задаем сколько отсчетов мы хотим исследовать
zaderzhka = z - fragment; # переменная обуславшивающая в дальнейшем индексы, а также начало отсчета для фрагмента

ECG_fg = ecg[zaderzhka:end];
ap_fg = ap[zaderzhka:end];
ad_fg = ad[zaderzhka:end];
""""""


w=128 #поскольку 64 для Fd 512, то для 1000 возьмем 128
tochka,tochka_Flt,Filted,ap_fg,pos_PPG_max,pos_PPG_min = PPG(ap_fg,fs,w) 
tochkaAD,tochka_FltAD,FiltedAD,ad_fg,pos_AD_max,pos_AD_min = PPG(ad_fg,fs,w) 
Filted_ECG,ECG_fg,R_zub_filt,R_zub,t,pos_ECG_R = Pan_Tompkins(ECG_fg,fs)



plot(t,[ap_fg,ad_fg,ECG_fg], layout = (3,1))

plot(t,[ap_fg,Filted],layout = (2,1))
plot(t,[ECG_fg,Filted_ECG],layout = (2,1))
plot(t,[ad_fg,FiltedAD],layout = (2,1))

#нахождение Se и P
#[Se_min,P_min,TP_min,FP_min,FN_min] = calcStat(tMin,pos_test_min,300);


#строим график SSF
#plot(t,SSF)

#yline(0.7*threshold,'--')
#set(gca,'XLim', [0 tmax])
#title('SSF')

#построим график отфильтрованного сигнала с пиками на нем
x = plot(t,Filted,label="Filted")
scatter!(t,tochka_Flt, label="data", mc=:red, ms=2, ma=0.5)

y = plot(t,ap_fg,label="Filted")
scatter!(t,tochka, label="data", mc=:red, ms=2, ma=0.5)

plot(x,y,layout = (2,1))    
#savefig("graphics_PPG.png")

x = plot(t,FiltedAD,label="Filted")
scatter!(t,tochka_FltAD, label="data", mc=:red, ms=2, ma=0.5)

y = plot(t,ad_fg,label="Filted")
scatter!(t,tochkaAD, label="data", mc=:red, ms=2, ma=0.5)

plot(x,y,layout = (2,1))  
#savefig("graphics_AD.png")

x = plot(t,Filted_ECG,label="Filted")
scatter!(t,R_zub_filt, label="data", mc=:red, ms=2, ma=0.5)

y = plot(t,ECG_fg,label="Filted")
scatter!(t,R_zub, label="data", mc=:red, ms=2, ma=0.5)

plot(x,y,layout = (2,1))  
#savefig("graphics_ECG.png")

""""
Отобразим необходимые графики
"""

x = plot(t,ap_fg,label="ФПГ")
scatter!(t,tochka, label="data", mc=:red, ms=2, ma=0.5)

y = plot(t,ECG_fg,label="ЭКГ")
scatter!(t,R_zub, label="data", mc=:red, ms=2, ma=0.5)

z = plot(t,ad_fg,label="АД")
scatter!(t,tochkaAD, label="data", mc=:red, ms=2, ma=0.5)

plot(x,y,z,layout = (3,1))




"
Считаем параметры
"
function Ploshad(istart,iend,signal)
    istart = Int(round(istart*fs))
    iend = Int(round(iend*fs))
    ploshad = 0
    for i in signal[istart:iend]
        ploshad += i 
    end
    return ploshad
end  # рассчитывает площадь фигуры под графиком


Complexes = DataFrame(
    Комлпекс = String[],            #1
    Положение_R_зубца = Float64[],  #2
    Положение_offset_ФПГ = Float64[], #3 
    Положение_onset_ФПГ = Float64[], #4 
    Амплитуда_offset_ФПГ = Float64[], #5
    Амплитуда_onset_ФПГ = Float64[], #6 
    Положение_offset_АД = Float64[], #7
    Положение_onset_АД = Float64[], #8
    Амплитуда_offset_АД = Float64[], #9 
    Амплитуда_onset_АД = Float64[]  #10
)

print(length(pos_AD_max), " ",length(pos_ECG_R), " ",length(pos_PPG_max)," ", length(pos_PPG_min))

for i in eachindex(pos_PPG_min)
    for a in eachindex(pos_PPG_min)
        if ((pos_PPG_min[i] * fs - pos_ECG_R[a]*fs<=700)  && (pos_PPG_min[i]>pos_ECG_R[a]))
            ECG_POSMENT = pos_ECG_R[a]
            push!(Complexes,["Комлпекс $i",
                              ECG_POSMENT,
                              pos_PPG_min[i],
                              pos_PPG_max[i],
                              ap_fg[Int(round(pos_PPG_min[i]*fs))],
                              ap_fg[Int(round(pos_PPG_max[i]*fs))],
                              pos_PPG_min[i],
                              pos_PPG_max[i],
                              ad_fg[Int(round(pos_AD_min[i]*fs))],
                              ad_fg[Int(round(pos_AD_max[i]*fs))]])
        end
    end
end

print(Complexes)


Parametrs = DataFrame(
    ЧСС = Float64[],
    t1_Время_между_onset_и_offset_ФПГ = Float64[],
    t2_Время_между_offset_и_onset_ФПГ = Float64[],
    Амплитуда_между_onset_и_offset_ФПГ = Float64[],
    Отношение_амплитуд_ФПГ = Float64[],
    Длительность_ФПГ_комплекса = Float64[],
    Гиппотенуза_ФПГ= Float64[],
    Тангенс_ФПГ = Float64[],
    Время_между_R_и_offset_ФПГ = Float64[],
    Время_между_R_и_onset_ФПГ = Float64[],
    Площадь_под_ФПГ_t1 = Float64[],
    Площадь_под_ФПГ_t2 = Float64[],
    Время_между_R_и_offset_АД = Float64[],
    Время_между_R_и_onset_АД= Float64[],
    Амплитуда_между_onset_и_offset_АД = Float64[],
    CрАД = Float64[],
    Отношение_амплитуд_АД = Float64[],
    t1_Время_между_onset_и_offset_АД = Float64[],
    t2_Время_между_offset_и_onset_АД = Float64[],
    Площадь_под_АД_t1 = Float64[],
    Площадь_под_АД_t2 = Float64[]
)


for i in 1:length(Complexes[!,1])-1
    push!(Parametrs,[
          Complexes[i+1,2] - Complexes[i,2],
          Complexes[i,4] - Complexes[i,3],
          abs(Complexes[i,4] - Complexes[i+1,3]),
          Complexes[i,6] - Complexes[i,5],
          Complexes[i,6]/Complexes[i,5],
          Complexes[i+1,3]-Complexes[i,3],
          sqrt((Complexes[i,4] - Complexes[i,3])^2+(Complexes[i,6] - Complexes[i,5])^2),
          (Complexes[i,6] - Complexes[i,5])/(Complexes[i,4] - Complexes[i,3]),
          abs(Complexes[i,2] - Complexes[i,3]),
          abs(Complexes[i,2] - Complexes[i,4]),
          Ploshad(Complexes[i,3],Complexes[i,4],ap_fg),
          Ploshad(Complexes[i,4],Complexes[i+1,3],ap_fg),
          abs(Complexes[i,2] - Complexes[i,7]),
          abs(Complexes[i,2] - Complexes[i,8]),
          Complexes[i,10] - Complexes[i,9],
          1/3*Complexes[i,10] + 2/3*Complexes[i,9],
          Complexes[i,10]/Complexes[i,9],
          Complexes[i,8] - Complexes[i,7],
          abs(Complexes[i,8] - Complexes[i+1,7]),
          Ploshad(Complexes[i,7],Complexes[i,8],ad_fg),
          Ploshad(Complexes[i,8],Complexes[i+1,7],ad_fg),
          ])

end



Parametrs

histogram([Parametrs[!,1],Parametrs[!,2],Parametrs[!,3],Parametrs[!,4],Parametrs[!,5],Parametrs[!,6],Parametrs[!,7],Parametrs[!,8],Parametrs[!,9],
           Parametrs[!,10],Parametrs[!,11],Parametrs[!,12],Parametrs[!,13],Parametrs[!,14],Parametrs[!,15],Parametrs[!,16],Parametrs[!,17],Parametrs[!,18],
           Parametrs[!,19],Parametrs[!,20],Parametrs[!,21]],
         #label =["Время между onset и offset" "Амплитуда_между_onset_и_offset" "Гиппотенуза ФПГ" "Тангенс наклона" "Время между R и onset" "Время между R и offset"] ,
         bins = 6,
         layout = (5,5))

#savefig("histogram.png")


     
#cols = [:Время_между_onset_и_offset,:Амплитуда_между_onset_и_offset,:Длительность_ФПГ_комплекса,
 #     :Гиппотенуза_ФПГ,:Тангенс_ФПГ,:Время_между_R_и_offset,:Время_между_R_и_onset]
n = size(Parametrs)[2]
m = size(Parametrs)[2]
M = cor(Matrix(Parametrs))  # correlation matrix
heatmap(M, fc=cgrad([:white,:dodgerblue4]), xrot=90, yflip=true)
annotate!([(j, i, text(round(M[i,j],digits=3), 5,"Computer Modern",:black)) for i in 1:n for j in 1:m])
#savefig("Cor_Matrix.png")
#corrplot([Parametrs[!,1],Parametrs[!,2],Parametrs[!,3],Parametrs[!,4],Parametrs[!,5],Parametrs[!,6],Parametrs[!,6],Parametrs[!,7]])

x = [10,5,3,2,8,4,3,1,10]
y = 0
for i in x[3:5]
    y += 1 * i
end
y


Complexes[1,3]
Complexes[1,4]

print(Ploshad(Complexes[1,3],Complexes[1,4],ap_fg))

watch = 1:1000
x = plot(t[watch],ap_fg[watch],label="ФПГ")


y = plot(t[watch],ECG_fg[watch],label="ЭКГ")


z = plot(t[watch],ad_fg[watch],label="АД")


plot(x,y,z,layout = (3,1))

savefig("orkworkwokr.png")