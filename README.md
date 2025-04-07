# <center>**ACTIVIDAD #3: Flujo de carga cuasi-estático con series de tiempo**</center>

<div style="text-align: justify;">
En este proyecto implementamos un Flujo de Potencia Cuasi Dinámico en sistemas eléctricos utilizando el Método de Punto Fijo. Su objetivo es determinar el estado del sistema, específicamente las tensiones nodales, considerando la variación temporal de la demanda y la generación. Para ello, se resuelven iterativamente las ecuaciones no lineales del flujo de potencia, permitiendo capturar la evolución dinámica del sistema en un horizonte de tiempo definido. Los resultados obtenidos pueden emplearse en estudios de estabilidad del sistema, planificación operativa y evaluación de contingencias, proporcionando una herramienta eficiente para el análisis de redes eléctricas con alta variabilidad en sus condiciones operativas.

----

# <center>*Marco Teórico*</center>


Este módulo implementa un análisis de flujo de potencia cuasi dinámico en sistemas eléctricos mediante el **método de punto fijo**.  
El objetivo es determinar el estado del sistema, específicamente las **tensiones nodales**, considerando variaciones temporales en la **demanda** y la **generación**.

El algoritmo resuelve iterativamente las ecuaciones no lineales del flujo de potencia y genera resultados que pueden utilizarse en **análisis de estabilidad** y **planificación operativa**.


#### Entradas
- **Datos de líneas**: Resistencias, reactancias, susceptancias y taps de transformadores.
- **Datos de nodos**: Tipo de nodo (Slack, PV, PQ), potencias generadas ($P_g$, $Q_g$) y consumidas ($P_d$, $Q_d$).
- **Demandas**: Perfiles de carga diarios o horarios para cada nodo tipo PQ.
- **Generación**: Perfiles de generación diarios o horarios, como energía solar.

#### Salidas
- **Tensiones en nodos**: Calculadas para cada nodo no-Slack en cada intervalo de tiempo.
- Archivos CSV con las tensiones para análisis posteriores.

#### Clasificación de Nodos
- **Nodos Slack**: Representa la subestación o barra fija (tensión controlada, usualmente 1∠0°). Absorbe las diferencias de potencia del sistema.
- **Nodos PV**: Generadores que controlan su tensión y potencia activa.
- **Nodos PQ**: Nodos donde su potencia activa y reactiva son conocidas y controladas.


#### Método de Punto Fijo
*El método de punto fijo resuelve iterativamente las ecuaciones del flujo de potencia en un espacio de Banach completo (e.g., $\mathbb{C}^n, \mathbb{R}^n$). Se define una transformación $T: \mathbb{C}^n \rightarrow \mathbb{C}^n$* como un mapa de contracción si satisface:
$||T(x)-T(y)||<k||x-y||$, donde $k<1$ es una constante de contracción.

*Esto garantiza un único punto fijo, calculado mediante:*
$x^{(k+1)} = T(x^{(k)})$

En el contexto del flujo de potencia:
1. Los nodos se dividen en:
    - **Nodos S**: Nodo Slack (subestación).
    - **Nodos N**: Nodos restantes (PV y PQ).
2. La ecuación base es: 
    \[
\begin{bmatrix}
\textit{I}_s \\
\textit{I}_n
\end{bmatrix}
=
\begin{bmatrix}
\textit{Y}_{ss} & \textit{Y}_{sn} \\
\textit{Y}_{ns} & \textit{Y}_{nn}
\end{bmatrix}
\begin{bmatrix}
\textit{V}_s \\
\textit{V}_n
\end{bmatrix}
\]

3. Las potencias nodales se relacionan con las tensiones:
$S_n=V_n \cdot I_n^*$

4. La ecuación iterativa para $V_n$ es: $$V_n^{k+1}=Y_{nn}^{-1}\left(\left(\frac{S_n}{V_n^{(k)}}\right)^*-Y_{ns}V_s\right)$$

    *Donde*:
    - $Y_{nn}$ es la matriz de admisión de los nodos no-Slack.
    - $Y_{ns}$ es la matriz de admisión de los nodos no-Slack a los nodos Slack.
    - $S_n$: Vector de potencias netas inyectadas.
    - $V_n$: Vector de tensiones de los nodos no-Slack.
    - $V_s$: Tensión del nodo Slack.

#### Criterio de Convergencia
El proceso iterativo se detiene cuando la norma del error entre dos iteraciones consecutivas es menor que un umbral (e.g., $1*10^{-6}$):
$||V_n^{k+1}-V_n^{k}||<1*10^{-6}$
La norma se define como:
$||V_n|| = \sqrt{V_1^2 + V_2^2 + \cdots + V_n^2}$

---
## **Funciones:**
*Para ejecutar este código, necesitas lo siguiente:*
- *Julia 1.11 o superior*
- *Bibliotecas*:
    - using LinearAlgebra
    - using DataFrames
    - using PrettyTables
    - using Clustering
    - using Statistics
    - using Plots
    - using CSV
    - using StatsPlots


### CALCULAR LA MATRIZ DE ADMITANCIA NODAL
**calcular_ynn_yns:**

*Descripción*

Construye la matriz de admitancias $(Y_{bus})$ del sistema de potencia y extrae las submatrices relevantes para estudios de flujo de potencia. Considera parámetros de líneas como resistencia, reactancia, susceptancia shunt y ajustes de tap en transformadores. La matriz Ybus es esencial para análisis de flujo de potencia y estabilidad.

*Requiere*

    """
        calcular_Ynn_Yns(lines, nodes)

    Calcula la matriz de admitancia nodal del sistema de potencia, separando las submatrices Ynn y Yns.

    # Parámetros
    - `lines::DataFrame`: DataFrame que contiene los datos de las líneas del sistema. 
    Debe incluir las siguientes columnas:
    - `FROM` (nodo de envío)
    - `TO` (nodo de recepción)
    - `R` (resistencia de la línea)
    - `X` (reactancia de la línea)
    - `B` (susceptancia de la línea)
    - `TAP` (relación de transformación del transformador, si aplica)

    - `nodes::DataFrame`: DataFrame que contiene los datos de los nodos del sistema.
    Debe incluir las siguientes columnas:
    - `NAME` (nombre del nodo)
    - `TYPE` (tipo de nodo: 1 = PQ, 2 = PV, 3 = Slack)

    # Retorno
    - `Ynn::Matrix{Complex{Float64}}`: Matriz de admitancia de los nodos distintos al Slack.
    - `Yns::Matrix{Complex{Float64}}`: Matriz de admitancia entre los nodos distintos al Slack y el nodo Slack.

    # Detalles
    1. Se construye la matriz de admitancia nodal `Ybus` a partir de los datos de las líneas.
    2. Se verifica que haya exactamente un nodo Slack en el sistema.
    3. Se llenan las posiciones de `Ybus` considerando:
    - La admitancia de la línea \( y_L = \frac{1}{R + jX} \).
    - La susceptancia de la línea \( B_s = \frac{B}{2} j \).
    - La relación de transformación `TAP`, si está definida.
    4. Se extraen las submatrices `Ynn` y `Yns` eliminando la fila y columna correspondientes al nodo Slack.

    

## FLUJO ESTÁTICO

*Descripción*

Calcula las tensiones en los nodos mediante el método iterativo de punto fijo. Utiliza la matriz $Y_{nn}$ y el vector $Y_{ns}$ para resolver el sistema de ecuaciones no lineales. Incluye un gráfico de convergencia del error en escala logarítmica.

*Requiere*

    """
        Flujo_punto_fijo(lines, nodes)

    Calcula el flujo de carga utilizando el método de punto fijo.

    # Parámetros
    - `lines::DataFrame`: DataFrame con los datos de las líneas del sistema eléctrico.
    - `nodes::DataFrame`: DataFrame con los datos de los nodos del sistema.

    # Retorna
    - `Vn::Vector{Complex{Float64}}`: Vector con los voltajes en los nodos diferentes al slack.

    # Descripción
    Esta función resuelve el flujo de carga en una red eléctrica utilizando el método de punto fijo. Se basa en la formulación de admitancias \( Y \) para calcular los voltajes en los nodos, excluyendo el nodo slack. 

    El método itera hasta que el error entre iteraciones sucesivas sea menor que una tolerancia predefinida (`1e-6`) o hasta alcanzar el número máximo de iteraciones (`100`).  

    Durante la ejecución:
    1. Se calcula la matriz de admitancia nodal (`Ynn`) y la matriz de acoplamiento con el nodo slack (`Yns`).
    2. Se inicializan los voltajes con valores unitarios.
    3. Se define la potencia compleja inyectada en los nodos (`Sn`).
    4. Se actualizan los voltajes en cada iteración utilizando la ecuación del método de punto fijo.
    5. Se almacena la evolución del error y se grafica la convergencia.
    6. Se muestra una tabla con los resultados finales de magnitud y ángulo de voltaje.

    Si la matriz `Ynn` no es invertible, la función lanza un error, ya que esto indica problemas en la topología de la red.


## Ajuste de Dataframes

*Descripción*

En el ajuste de Dataframes se puede modificar el ajuste de tiempo de las bases de datos para que sean de 5, 10, 20, 60 minutos, depende de que tan aproximado se quieran trabajar los datos.

*Requiere*

    """
    Ajusta los datos de generación solar y demanda eléctrica a una resolución de tiempo específica.

    # Parámetros
    - `data::DataFrame`: DataFrame con los datos originales de generación y demanda.
    - `freq_min::Int`: Frecuencia en minutos para la nueva resolución temporal.

    # Retorna
    - `df_ajustado::DataFrame`: DataFrame con los datos interpolados a la nueva resolución de tiempo.

    # Descripción
    Esta función toma un conjunto de datos de generación solar y demanda eléctrica,
    luego ajusta los valores a la frecuencia especificada mediante interpolación.
    El objetivo es obtener una representación más uniforme y adecuada para análisis posteriores.



## Flujo cuasi-dinamico:

*Descripción*

Simula un flujo de potencia cuasi-dinámico para múltiples escenarios (días típicos) e intervalos de tiempo de n minutos. Calcula las tensiones nodales considerando variaciones en generación y demanda, y almacena los resultados en archivos CSV.

*Requiere*

    """
    Flujo_cuasi_din(lines, nodes, Demandas_prom, dias_tipicos)

    Realiza el flujo de potencia cuasi-dinámico en una red eléctrica de distribución.

    # Parámetros
    - `lines::DataFrame`: DataFrame con la información de las líneas de la red.
    - `nodes::DataFrame`: DataFrame con la información de los nodos del sistema.
    - `Demandas_prom::DataFrame`: DataFrame con las demandas promedio de cada nodo en diferentes instantes de tiempo.
    - `dias_tipicos::DataFrame`: DataFrame con perfiles de generación renovable para días típicos.

    # Descripción
    Esta función calcula la evolución temporal de las tensiones en una red de distribución 
    utilizando un método cuasi-dinámico basado en iteraciones de punto fijo.  
    Se realizan los siguientes pasos:

    1. Se calculan las matrices de admitancia de la red (`Ynn` y `Yns`).
    2. Se verifica si la matriz `Ynn` es invertible.
    3. Se inicializan las tensiones en los nodos y la potencia compleja inyectada en la red.
    4. Se actualizan las tensiones nodales mediante un método iterativo de punto fijo.
    5. Se almacenan los resultados en un DataFrame con magnitud y ángulo de tensión en cada instante.
    6. Se exportan los resultados a archivos CSV para cada día típico analizado.

    # Retorna
    - Los resultados se almacenan en archivos CSV con los valores de magnitud y ángulo de tensión para cada día típico.

## Función de Gráficos

    """
    plot_voltage_heatmaps(file_name)

    Genera y muestra mapas de calor para visualizar la magnitud y el ángulo de tensión a lo largo del tiempo en una red eléctrica.

    # Parámetros
    - `file_name::String`: Nombre del archivo CSV que contiene los datos de tensión.

    # Descripción
    La función lee un archivo CSV con datos de magnitud y ángulo de tensión en distintos nodos de una red eléctrica a lo largo del tiempo.
    Se generan dos mapas de calor:
    1. **Mapa de calor de la magnitud de tensión**: Representa la magnitud de la tensión en función del tiempo y el nodo.
    2. **Mapa de calor del ángulo de tensión**: Muestra el ángulo de tensión en función del tiempo y el nodo.

    # Flujo de ejecución
    1. Carga el archivo CSV en un `DataFrame`.
    2. Filtra las columnas de magnitud y ángulo de tensión.
    3. Extrae los valores de tiempo desde los nombres de las columnas.
    4. Convierte los datos en matrices para su representación gráfica.
    5. Ajusta las dimensiones si es necesario para que coincidan con los ejes.
    6. Genera y muestra los mapas de calor utilizando `heatmap`.

    # Notas
    - El archivo CSV debe contener columnas con nombres en el formato `"Mag_HH:MM"` y `"Ang_HH:MM"`, donde `HH:MM` representa la hora del día en formato de 24 horas.
    - Si no se encuentran columnas de magnitud o ángulo, la función lanza un error.
    - Se utilizan las escalas de color `viridis` para la magnitud y `plasma` para el ángulo.
    - El eje X representa el tiempo en minutos desde la medianoche.
    - El eje Y representa los nodos de la red eléctrica.





#### Notas:
- *Puedes instalar las bibliotecas en Julia ejecutando:*
    ```julia
    using Pkg
    Pkg.add(["LinearAlgebra", "DataFrames", "CSV", "Plots", "Dates", "Clustering", "Statistics"])
- flujo_cuasi_dinamico utiliza un enfoque simplificado con iteraciones fijas (sin criterio de convergencia dinámico).

[![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC_BY--NC--ND_4.0-lightgrey)](https://creativecommons.org/licenses/by-nc-nd/4.0/)

---

