import streamlit as st
import cantera as ct
from scipy.optimize import fminbound
import numpy as np


with st.sidebar:
    st.write("## Entrada:")
  
    st.write("### Reagentes:")
    mp_ar  =  st.number_input(label="Vazão de ar (g/s):"              , value=210. ,min_value=0., step=1.)
    mp_o2  =  st.number_input(label="Vazão de O2 (g/s):"              , value=40.  ,min_value=0., step=1.)
    mp_ch4 =  st.number_input(label="Vazão de CH4 (g/s):"             , value=20.  ,min_value=0., step=1.)
    Ti     =  st.number_input(label="Temp. inicial dos Reagentes (K):", value=300 ,min_value=0, step=1)    
    

    st.write("### Câmara")
    p0  =  st.number_input(label="Pressão da Câmara(bar):", value=20.0 , min_value=0.0, step=1.0)

    st.write("### Tubeira")
    dA = st.number_input(label="Área da Garganta (mm2): ", value=128.2, min_value=0.0, step=0.1 )
    DA = st.number_input(label="Área da Saída (mm2): "   , value=500.0, min_value=0.0, step=0.1 )
    

## CÀLCULO --------------------------------------------------------------------------------------------------------------------------------
#Camara de Combustao    
n_n2  = 0.767*mp_ar/28
n_o2  = (mp_o2 + 0.233*mp_ar)/32
n_ch4 = mp_ch4/16

X = 'N2:'+str(n_n2)+', O2:'+str(n_o2)+', CH4:'+str(n_ch4)
gas = ct.Solution('gri30.cti')

gas.TPX = Ti, 1e5*p0, X
gas.equilibrate('HP')
h0 = gas.h
g0 = gas.cp/gas.cv  
R0 = ct.gas_constant/gas.mean_molecular_weight   
T0 = gas.T
r0 = gas.density
X0 = gas.X
spec = gas.species_names
RA = (DA/dA)
        

#Vazões de Saída         
M = 1.0
Tg = gas.T*(1 + 0.5*(g0-1)*M**2)**(-1)
pg = gas.P*(1 + 0.5*(g0-1)*M**2)**(-g0/(g0-1))
rg = gas.density*(1 + 0.5*(g0-1)*M**2)**(-1/(g0-1))
gas.TP = Tg,pg
vg = (2*(h0-gas.h))**0.5
mp_calc = 1000*rg*vg*(dA*1e-6)

#Saida
aux = lambda M: np.abs(RA*M - ((0.5*(g0+1.))**(-0.5*(g0+1.)/(g0-1.)))*(1 + 0.5*(g0-1)*M*M)**(0.5*(g0+1.)/(g0-1.)))
Ms = fminbound(aux,1,10,disp=False)
Ts = gas.T*(1 + 0.5*(g0-1)*Ms**2)**(-1)
ps = gas.P*(1 + 0.5*(g0-1)*Ms**2)**(-g0/(g0-1))
rs = gas.density*(1 + 0.5*(g0-1)*Ms**2)**(-1/(g0-1))
gas.TP = Ts,ps
hs = gas.h
vs = (2*(h0-gas.h))**0.5


### Visualizacao de Dados --------------------------------------------------------------------------------------------------------------------------------
st.title("GAV - Gerador de Ar Viciado")

st.write(f"""
##### Entrada:
Pressão de Reserv. Medida (Bar):           {round(p0,2)    }\\
Temperatura Inicial dos Gases (K):         {int(Ti)        }\\
Vazao mássica de ar (g/s):                 {round(mp_ar,3) }\\
Vazao mássica de O2 (g/s):                 {round(mp_o2,3) }\\
Vazao mássica de CH4 (g/s) :               {round(mp_ch4,3)}\\
Razão de Àreas da Tubeira : {round(RA,3)    }
""")

st.write(f""" ##### Vazões
Vazao mássica Total Experimental (g/s) : {round(mp_ch4+mp_o2+mp_ar,1)}\\
Vazao mássica Total Calculada (g/s) :    {round(mp_calc,1)}\\
Vazao Calc/Exp : {round(mp_calc/(mp_ch4+mp_o2+mp_ar),3)}
""")


#with col2: 
st.write("""
##### Câmara de Combustão / Reservatório
Razão de Equiv. :            {round(gas.equivalence_ratio(),2)}\\
Temperatura (K):   {int(T0)}\\
Pressão (Bar):     {round(p0,2)}\\
Densidade (kg/m3): {round(r0,2)}\\
Razão de cal. esp. cp/cv:          {round(gas.cp/gas.cv,3)} 
""")

Texto = "Fração Molar de Espécies > 0,1% :  \n"
for nome in spec:
    if gas[nome].X > 1e-3 :        
        Texto = Texto + str(nome)+": "+str(round(gas[nome].X[0],2))+"  \n"
        
st.write(Texto)
        



st.write(f"""
##### Saída
Temperatura de Saída (K):  {int(Ts)}\\
Pressão de Saída (Bar):  {round(ps/1e5,3)}\\
Densidade de Saída (kg/m3):  {round(rs,3)}\\
Velocidade de Saída (m/s):  {int(vs) }\\
Mach: {round(Ms,3 )} 
""")











