import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint 


def conveccao(T1,T2,h):
	c=h*A*(T1-T2)
	return c

def conducao(T1,T2,d,k):
	co=k*A*(T1-T2)/d
	return co

def radsol(e):
	r=e*I*A
	return r

def transp(T1,Tp):
	return (cl+((T1+Tp)/2)*(Sg-Sl)+R/Pma*T1*((math.log(1/ur))))

#parametros
#temp incomodo
em=0.90 #coeficiente de emissividade mulher w/m
erp=1 #roupa preta
I=1000 #W/m**2 #indice de insolacao
A=1.64 #m**2 #area 
h= 30 #baixo W/m**2 *K #Coef. de transferência convectiva
hma= 5 #alto #Coef. de transferência convectiva
kra=1/0.047 #mKS/J #Condutância térmica
km=1/0.47  #Condutância térmica
Tp=[295]
Tr=[295]
Ta=320
#Tai=[37]
cl=2.256 #Calor latente de vaporização da agua J/
Pma=18 #Peso molecular da agua g/m
Sl=69.9#entropia agua liquida
Sg=188.8#entropia agua gasosa
R=8.31#contante gasosa
ur=0.15 # umidade relativa
m=61 #massa mulher kg
c=3470
o=5.67e-8
em=0.9
dr=0.001

#Mulher com roupa leve
tempo=np.arange(0,12000,0.1)

def mulhersroupa(Tp,tempo):
    dTmdt=(radsol(em)-(transp(Ta,Tp))*(10**-3)-conveccao(Ta,Tp,hma))/(m*c)
    return dTmdt

temperatura_pele_mulher=odeint(mulhersroupa,Tp[0],tempo)
plt.plot(tempo,temperatura_pele_mulher, 'b')
plt.ylabel('Temperatura da pele mulher(°K)')
plt.title('Temperatura da pele da mulher com roupas leves')
plt.xlabel('tempo(s)')
plt.grid(True)
plt.show()

def mulherroupagrudada(Tp,tempo):
    dTmdt=(radsol(erp)-(transp(Ta,Tp))*(10**-3)-conveccao(Ta,Tp,hma))/(m*c)
    return dTmdt
    
temperatura_mulher_roupagrudada=odeint(mulherroupagrudada,Tp[0],tempo)
plt.plot(tempo,temperatura_mulher_roupagrudada, 'r')
plt.ylabel('Temperatura da pele mulher(°K)')
plt.title('Temperatura da pele da mulher com traje muçulmano')
plt.xlabel('tempo(s)')
plt.grid(True)
plt.show()

T2=[295,295]

t=np.arange(0,12000,0.1)
def mulherroupasolta(T2,t):
    Tai=T2[0]
    Tp=T2[1] 
    dTAdt=(radsol(erp)+(transp(Tai,Tp)*(10**(-3)))-conveccao(Tai,Tp,hma)-conveccao(Ta,Tp,2))/(m*c)
    dTMdt=(conveccao(Ta,Tp,h)-(transp(Tai,Tp)*(10**(-3))))/(m*c)
    return [dTAdt,dTMdt]


temperatura_mulher_roupasolta=odeint(mulherroupasolta,T2,t)
plt.plot(t,temperatura_mulher_roupasolta[:,0],'bo')
plt.plot(t,temperatura_mulher_roupasolta[:,1], 'g')
plt.ylabel('Temperaturas (°K)')
plt.title('Temperatura da pele da mulher e do ar entre o traje beduíno')
plt.xlabel('tempo(s)')
plt.grid(True)
plt.show()

plt.plot(tempo,temperatura_pele_mulher, 'b')
plt.plot(t,temperatura_mulher_roupasolta[:,1], 'g')
plt.plot(tempo,temperatura_mulher_roupagrudada, 'r')
plt.ylabel('Temperaturas (°K)')
plt.title('Temperatura da pele da mulher nas três circunstâncias')
plt.xlabel('tempo(s)')
plt.grid(True)
plt.show()