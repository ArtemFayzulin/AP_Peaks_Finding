#my work
include("readbin.jl")
include("#PPG working.jl")
include("#ECG_working.jl")
include("CalcStats.jl")
using DataFrames
using StatsPlots


filepath = raw"C:\Users\artem\OneDrive\Рабочий стол\Для чего\вуз\Science And Diploma\Програмесы\Сигналы САКР\Мельникова 1\Мельникова_Елизавета_Дмитриевна_21-04-22_11-43-20_.hdr"
named_channels, fs, timestart, units = readbin(filepath)


ap = named_channels.FPrsNorm1 ./ 1000 # давление
ecg = named_channels.LR ./ 1000 #экг

"""
Первый вариант считывания сигнала
"""
istart = 5000
iend = 120000

ECG_fg = ecg[istart:iend];
ap_fg = ap[istart:iend];

"""
Второй вариант рассчитывания сигнала
"""
z = length(ap);#считаем длину сигнала
fragment = 25000; #вручную задаем сколько отсчетов мы хотим исследовать
zaderzhka = z - fragment; # переменная обуславшивающая в дальнейшем индексы, а также начало отсчета для фрагмента

ECG_fg = ecg[zaderzhka:end];
ap_fg = ap[zaderzhka:end];
""""""

w=128 #поскольку 64 для Fd 512, то для 1000 возьмем 128
tochka,tochka_Flt,Filted,ap_fg,pos_PPG_max,pos_PPG_min = PPG(ap_fg,fs,w) 
Filted_ECG,ECG_fg,R_zub_filt,R_zub,t,pos_ECG_R = Pan_Tompkins(ECG_fg,fs)

plot(t,[ap_fg,ECG_fg], layout = (2,1))


plot(t,[ap_fg,Filted],layout = (2,1))
plot(t,[ECG_fg,Filted_ECG],layout = (2,1))
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

x = plot(t,Filted_ECG,label="Filted")
scatter!(t,R_zub_filt, label="data", mc=:red, ms=2, ma=0.5)


y = plot(t,ECG_fg,label="Filted")
scatter!(t,R_zub, label="data", mc=:red, ms=2, ma=0.5)

plot(x,y,layout = (2,1))  
#savefig("graphics_ECG.png")

""""
Отобразим необходимые графики
"""

x = plot(t,ap_fg,label="Filted")
scatter!(t,tochka, label="data", mc=:red, ms=2, ma=0.5)

y = plot(t,ECG_fg,label="Filted")
scatter!(t,R_zub, label="data", mc=:red, ms=2, ma=0.5)

plot(x,y,layout = (2,1))



"
Считаем параметры
"

Complexes = DataFrame(
    Комлпекс = String[],
    Положение_R_зубца = Float64[],
    Положение_offset = Float64[],
    Положение_onset = Float64[],
    Амплитуда_offset = Float64[],
    Амплитуда_onset = Float64[]
)

for i in eachindex(pos_PPG_min)
    for a in eachindex(pos_PPG_min)
        if ((pos_PPG_min[i] * fs - pos_ECG_R[a]*fs<=700)  && (pos_PPG_min[i]>pos_ECG_R[a]))
            ECG_POSMENT = pos_ECG_R[a]
            push!(Complexes,["Комлпекс $i",ECG_POSMENT,pos_PPG_min[i],pos_PPG_max[i],ap_fg[Int(round(pos_PPG_min[i]*fs))],ap_fg[Int(round(pos_PPG_max[i]*fs))]])
        end
    end
end

Complexes


Parametrs = DataFrame(
    Время_между_onset_и_offset = Float64[],
    Амплитуда_между_onset_и_offset = Float64[],
    Длительность_ФПГ_комплекса = Float64[],
    Гиппотенуза_ФПГ= Float64[],
    Тангенс_ФПГ = Float64[],
    Время_между_R_и_offset = Float64[],
    Время_между_R_и_onset = Float64[]
)

for i in 1:length(Complexes[!,1])-1
    push!(Parametrs,[Complexes[i,4] - Complexes[i,3],
          Complexes[i,6] - Complexes[i,5],
          Complexes[i+1,3]-Complexes[i,3],
          sqrt((Complexes[i,4] - Complexes[i,3])^2+(Complexes[i,6] - Complexes[i,5])^2),
          (Complexes[i,6] - Complexes[i,5])/(Complexes[i,4] - Complexes[i,3]),
          abs(Complexes[i,2] - Complexes[i,3]),
          abs(Complexes[i,2] - Complexes[i,4])])
end

Parametrs

histogram([Parametrs[!,1],Parametrs[!,2],Parametrs[!,3],Parametrs[!,4],Parametrs[!,5],Parametrs[!,6]],
         label =["Время между onset и offset" "Амплитуда_между_onset_и_offset" "Гиппотенуза ФПГ" "Тангенс наклона" "Время между R и onset" "Время между R и offset"] ,
         bins = 6,
         layout = (3,2))

#savefig("histogram.png")

     
#cols = [:Время_между_onset_и_offset,:Амплитуда_между_onset_и_offset,:Длительность_ФПГ_комплекса,
 #     :Гиппотенуза_ФПГ,:Тангенс_ФПГ,:Время_между_R_и_offset,:Время_между_R_и_onset]
M = cor(Matrix(Parametrs))  # correlation matrix
heatmap(M, fc=cgrad([:white,:dodgerblue4]), xrot=90, yflip=true)
annotate!([(j, i, text(round(M[i,j],digits=3), 8,"Computer Modern",:black)) for i in 1:n for j in 1:m])
#savefig("Cor_Matrix.png")
#corrplot([Parametrs[!,1],Parametrs[!,2],Parametrs[!,3],Parametrs[!,4],Parametrs[!,5],Parametrs[!,6],Parametrs[!,6],Parametrs[!,7]])