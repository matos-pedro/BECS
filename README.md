# BECS - Calculadora [![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://becs-app.streamlit.app/)
A calculadora utiliza a suíte Cantera para estimação das condições de estagnação de uma câmara de combustão dados como conhecidas as vazões e temperaturas dos reagentes e a pressão dos produtos. Em seguida, o código estima as condições de saída de bocal acoplado ao combustor dadas as áreas da garganta e da saída do bocal.

## Câmara de Combustão - Equilíbrio a HP, solução padrão
A partir das vazões e temperaturas dos reagentes - e portanto a entalpia do sistema - e a pressão dos produtos, especifica-se por completo os produtos da queima, e portanto da câmara de combustão, a partir de um equilíbrio do tipo HP, implementado com auxílio do Cantera como segue:
```
gas.TPX = T_reagentes, P_produtos, X_reagentes
gas.equilibrate('HP')
```
Uma vez equilibrado, o objeto ```gás``` corresponde à estagnação.  


## Câmara de Combustão - Equilíbrio a TP 
A partir das vazões e a pressão dos produtos, especifica-se por completo os produtos da queima, e portanto da câmara de combustão, a partir de um equilíbrio do tipo TP, implementado com auxílio do Cantera como segue:

```
1 - calcula mp_total e identifica frações molares X
    2 - chuta temperatura de estagnação T0
    3 - equilibra o gas: 
        gas.TPX = T_guess, P_produtos, X_reagentes
        gas.equilibrate('TP')
    4 - calcula vazão de saída mp_calc
5 - realiza etapas de 2 a 4 até que mp_calc = mp_total e a temperatura de estagnação seja obtida. 
```
Por fim, equilibra-se a classe ```gas``` na pressão e temperatura de estagnação identificadas.  

### Bocal
O escoamento ao longo do bocal é assumido isentrópico e congelado, portanto pode-se afirmar que:   

- dada a razão de áreas, o número de Mach à saída da tubeira é aquele que satisfaz a seguinte relação:

$$ \frac{A_{Saída}}{A_{Garganta}} = \left( \frac{1}{M} \right)\left( \frac{\gamma+1}{2}\right)^{-\frac{\gamma+1}{2(\gamma-1)}}\left(1 + \frac{\gamma-1}{2}M^2\right)^{\frac{\gamma+1}{2(\gamma-1)}} $$

- a partir do número de Mach, temperatura e pressão podem ser obtidas por:

$$ T_{Saída} = T_{Câmara} \left(1 + \frac{\gamma-1}{2}M^2\right)^{-1} $$

$$ P_{Saída} = P_{Câmara} \left(1 + \frac{\gamma-1}{2}M^2\right)^{-\frac{\gamma}{\gamma-1}} $$

Fim.