import streamlit as st
import cantera as ct
from scipy.optimize import fminbound
import numpy as np
import pandas as pd

st.set_page_config(
    page_title="BECS - Calculadora",
    layout="wide",
    )

with st.sidebar:
    st.write("## Entrada")
  
    st.write("### Reagentes")
    mp_ar  =  st.number_input(label="Vazão de ar (g/s):"              , value=210. ,min_value=0., step=1.)
    mp_o2  =  st.number_input(label="Vazão de O2 (g/s):"              , value=40.  ,min_value=0., step=1.)
    mp_ch4 =  st.number_input(label="Vazão de CH4 (g/s):"             , value=20.  ,min_value=0., step=1.)
    Ti     =  st.number_input(label="Temp. inicial dos Reagentes (K):", value=300 ,min_value=0, step=1)    
    
    st.write("### Câmara de Combustão")
    p0  =  st.number_input(label="Pressão da Câmara(bar):", value=20.0 , min_value=0.0, step=1.0)
    eta  =  st.number_input(label="Grau de Combustão:", value=1.0, max_value=1.0 , min_value=0.0, step=0.01)

    st.write("### Dimensões da Tubeira")
    dA = st.number_input(label="Área da Garganta (mm2): ", value=128.2, min_value=0.0, step=0.1 )
    DA = st.number_input(label="Área da Saída (mm2): "   , value=500.0, min_value=0.0, step=0.1 )
    

## CÀLCULO --------------------------------------------------------------------------------------------------------------------------------
#Camara de Combustao    
n_n2  = 0.767*mp_ar/28
n_o2  = (mp_o2 + 0.233*mp_ar)/32
n_ch4 = mp_ch4/16

X = 'N2:'+str(n_n2)+', O2:'+str(n_o2)+', CH4:'+str(n_ch4)
gas = ct.Solution('gri30.cti')

gas_brnt = ct.Quantity(gas, constant='HP')
gas_brnt.TPX = Ti, 1e5*p0, X
gas_brnt.equilibrate('HP')

gas_unbrnt = ct.Quantity(gas, constant='HP')
gas_unbrnt.TPX = Ti, 1e5*p0, X

gas_brnt.moles   = eta
gas_unbrnt.moles = 1. - eta

gas = gas_unbrnt + gas_brnt
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
st.title("BECS")

st.write(f"""
#### Parâmetros de Entrada
Pressão Medida da Câmara (Bar):      {round(p0,2)    }\\
Temperatura Inicial dos Gases (K):   {int(Ti)        }\\
Vazao mássica de ar (g/s):           {round(mp_ar,3) }\\
Vazao mássica de O2 (g/s):           {round(mp_o2,3) }\\
Vazao mássica de CH4 (g/s) :         {round(mp_ch4,3)}\\
Grau de Combustão:                   {round(eta,2)}  \\
Razão de Àreas da Tubeira :          {round(RA,3)}
""")

st.write("""Obs: um grau de combustão de $\eta$ significa uma mistura contendo $\eta$ moles dos produtos da combustão e $(1-\eta)$ moles de reagentes não consumidos, ambos à pressão definida pelo usuário.
""")


st.write(f""" #### Vazões Calculadas
Vazao mássica Total Experimental (g/s) : {round(mp_ch4+mp_o2+mp_ar,1)}\\
Vazao mássica Total Estimada para o bocal (g/s) :    {round(mp_calc,1)}\\
Razão entre Vazões Estimada e Experimental : {round(mp_calc/(mp_ch4+mp_o2+mp_ar),3)}
""")


#with col2: 
st.write(f"""
#### Caracterização da Câmara de Combustão
Razão de Equivalência : {round(gas.equivalence_ratio(),2)}\\
Temperatura Total (K):   {int(T0)}\\
Pressão Total (Bar):     {round(p0,2)}\\
Densidade Total (kg/m3): {round(r0,2)}\\
Razão de Calores Específicos cp/cv:  {round(gas.cp/gas.cv,3)} 
""")

st.write(f"###### Espécies com X$_i$ $\ge$ 0,001 ")
Especies = pd.DataFrame()
Especies['Espécie Química'] = spec
Especies['Concentração Molar Xi']   = gas.X
Especies = Especies[  Especies['Concentração Molar Xi'] >= 0.001 ]
Especies.set_index(['Espécie Química'],inplace=True)

st.dataframe(Especies,use_container_width=2)
        

st.write(f"""
##### Escoamento Livre
Temperatura de Saída (K):  {int(Ts)}\\
Pressão de Saída (Bar):  {round(ps/1e5,3)}\\
Densidade de Saída (kg/m3):  {round(rs,3)}\\
Velocidade de Saída (m/s):  {int(vs) }\\
Mach: {round(Ms,3 )} 
""")











