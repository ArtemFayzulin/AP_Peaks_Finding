#pip install -U scikit-learn

import numpy as np
import pandas
import matplotlib.pyplot as plt
import random
#import sklearn
from sklearn.linear_model import LinearRegression #линейная регрессия
from sklearn.preprocessing import PolynomialFeatures # полиномиальная регрессия
from sklearn.model_selection import train_test_split # модуль для деления на тестовое и обучающее множества
from sklearn.metrics import mean_squared_error, mean_absolute_error #модули расчета ошибок


#Исходные данные: dataFr - DataFrame со всеми параметрами, 
# pars_nams - массив с именами параметров, которые используются в качестве иксов (предикторов) 
# например, pars_nams может задаваться как: pars_nams = ['x1', 'x2'] 
# y - значение целевой переменной (ненормированной) в виде массива,
# если y находится в DataFrame dataFr, то может быть определена как: y = dataFr[['DAD']]

dataFr = pandas.read_csv("data.csv")
y = pandas.read_csv("CO.csv")
y = pandas.Series.tolist(y["Величина_сердечного_выброса"])


y_normed = []

for i in range(0,len(y)):
    y_normed.append((y[i] - np.mean(y))/np.std(y))
    
     
#print(y_normed,'',len(y),type(y))
pars_nams = []
for i in dataFr.columns:
    pars_nams.append(i)
 
dictionary = {}     
for i in pars_nams:
    dictionary[i] = (dataFr[i]-np.mean(dataFr[i]))/np.std(dataFr[i])

Frame = pandas.DataFrame.from_dict(dictionary)

 
# dataFr[pars_nams[0]]
# np.corrcoef(dataFr[pars_nams[1]], y)[1,0]
#ищем корреляцию между параметрами
#cor_coffs = []

cor_coffs_new = []
for i in dataFr.columns:
    # if abs(np.corrcoef(Frame[i], y_normed)[1,0])>0.4:
    #    cor_coffs_new.append(i)
       cor_coffs_new.append(np.corrcoef(Frame[i], y_normed)[1,0])

print(cor_coffs_new)
new_dict = {}     
for i in cor_coffs_new:
    new_dict[i] = (dataFr[i]-np.mean(dataFr[i]))/np.std(dataFr[i])
 
Parametrs = pandas.DataFrame.from_dict(new_dict)  
print(type(Parametrs),Parametrs) 
parametrs_nams = []
for i in Parametrs.columns:
    parametrs_nams.append(i)    
    

print(dataFr)
print(pars_nams)

print(Parametrs)
print(parametrs_nams)
# Задаем модель с нужной степенью deg (в данном случае 3)
deg = 1
poly = PolynomialFeatures(degree = deg)
#X_poly = poly.fit_transform(Frame[pars_nams])
X_poly = poly.fit_transform(Parametrs[parametrs_nams])
# делим на тестовое и обучающее множества
x_train, x_test, y_train, y_test = train_test_split(X_poly, y,  test_size=0.3, )


lin_mod = LinearRegression() #создаем модель
lin_mod.fit(x_train, y_train) #обучаем модель на тренировочных данных

yp=lin_mod.predict(x_test) # получаем предсказания на тестовом множестве

# получаем коэффициенты регрессионной кривой (intercept - свободный член получается отдельно)
lin_coeff = lin_mod.coef_
lin_inter = lin_mod.intercept_
print(lin_coeff,len(lin_coeff))
print(lin_inter)

def model1(x1,x2,coef1,coef2,free):
    y = []
    for i in range(1,len(x1)):
        y.append(coef1*x1[i]+coef2*x2[i]+free)
    return y
    

x1 = [random.randrange(78, 91) for i in range(1,10)]
x2 = [random.randrange(11, 14) for i in range(1,10)]
cooo = model1(x1,x2,lin_coeff[1],lin_coeff[2],lin_inter)

print(cooo)




    
# рассчитываем статистики
Rem = y_test - yp #остатки, распределение которых и нужно посмотреть
MAE = mean_absolute_error(y_test, yp) # Средняя ошибка (абсолютная)
RMSE = np.sqrt(mean_squared_error(y_test, yp)) # среднеквадратичная ошибка

print(Rem)
# по MAE и RMSE можно сравнивать модели, по Rem оценивать применимость модели