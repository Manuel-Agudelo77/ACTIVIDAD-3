using LinearAlgebra
using DataFrames
using PrettyTables
using Clustering
using Statistics
using Plots
using CSV
using StatsPlots

## -------------------------------------------------------------------------------------------------##
## CALCULAR LA MATRIZ DE ADMITANCIA NODAL

function calcular_Ynn_Yns(lines,nodes)
    """
    Entradas:   lines: DataFrame con los datos de las lineas del sistema
                nodes : DataFrame con los datos de los nodos del sistema
    Salida :    Ynn : Matriz de admitancia de los nodos diferentes al Slack.
                Yns : Matriz de admitancia de los nodos diferentes al Slack con respecto al Slack.
    """
    num_nodes = nrow(nodes)
    num_lines = nrow(lines)
    Ybus = zeros(num_nodes, num_nodes)*1im
    nodo_slack = nodes[nodes.TYPE .== 3, "NAME"]
    if length(nodo_slack) != 1
        error("Debe haber exactamente un nodo Slack")
    end
    s = parse(Int64, match(r"\d+", string(nodo_slack[1])).match)
    for k = 1:num_lines
        # Nodo de envío
        n1 = parse(Int64,string(lines.FROM[k])[2:end])
        # Nodo de recibo
        n2 = parse(Int64,string(lines.TO[k])[2:end])
        # Admitancia de la línea
        yL = 1/(lines.R[k]+lines.X[k]*1im)
        # Susceptancia de la línea
        Bs = lines.B[k]*1im/2
        # Valor del TAP
        t = lines.TAP[k]
        if lines.TAP[k] == 0
            Ybus[n1,n1] += yL + Bs   # Dentro de la diagonal
            Ybus[n1,n2] -= yL        # Fuera de la diagonal
            Ybus[n2,n1] -= yL        # Fuera de la diagonal
            Ybus[n2,n2] += yL + Bs   # Dentro de la diagonal
        else
            Ybus[n1,n1] += (t^2 - t)*yL  # Dentro de la diagonal
            Ybus[n1,n2] -= t*yL          # Fuera de la diagonal
            Ybus[n2,n1] -= t*yL          # Fuera de la diagonal
            Ybus[n2,n2] += (1-t)*yL      # Dentro de la diagonal
        end
    end
    # Separación de la matriz Ybus en Ynn y Yns
    Ynn = Ybus[setdiff(1:end,s),setdiff(1:end,s)]
    Yns = Ybus[setdiff(1:end, s), s]
    return Ynn, Yns
end

## -------------------------------------------------------------------------------------------------##
## FLUJO ESTÁTICO

function Flujo_punto_fijo(lines, nodes)
    """
    Entradas:   lines: DataFrame con los datos de las lineas del sistema
                nodes : DataFrame con los datos de los nodos del sistema
    Salida :    Ynn : Matriz de admitancia de los nodos diferentes al Slack.
                Yns : Matriz de admitancia de los nodos diferentes al Slack con respecto al Slack.

    """
    num_nodes = nrow(nodes)
    # Cálculo de Ynn y Yns
    Ynn, Yns = calcular_Ynn_Yns(lines,nodes)
    if rank(Ynn) < size(Ynn, 1)
        error("La matriz Ynn no es invertible. Verifica la topología del sistema.")
    end
    iYnn = inv(Ynn)

    # Se inicializa Vn y Vs.
    Vn = ones(num_nodes-1) + 1im*zeros(num_nodes-1)
    Vs = 1 + 0*1im
    
    # Se define Sn
    @views Sn = (nodes.PGEN[2:end] - nodes.PLOAD[2:end]) + (nodes.QGEN[2:end]-nodes.QLOAD[2:end])*1im
    
    # Se comienza la iteración de punto fijo
    tolerancia = 1e-6  # Error máximo permitido
    max_iter = 100  # Límite de iteraciones
    errores = Float64[]  # Se usa un vector dinámico
    ite = Int[]  # Se usa un vector dinámico para registrar iteraciones

    for i = 1:max_iter
        # Guardar el valor anterior de Vn
        Vn_ant = Vn  
        # Actualización de Vn
        Vn = iYnn * (conj.(Sn ./ Vn) .- Yns * Vs)  
        # Calcular el error
        error = norm(Vn - Vn_ant)  
        push!(errores, error)
        push!(ite, i)
    
        if error < tolerancia
            break  # Se detiene si el error es menor a la tolerancia
        end
    end
    
    theme(:dark)
    if all(errores .> 0)  
        P = plot(ite, errores, yscale=:log10, xlabel="Iteraciones", ylabel="Error", label="Convergencia")
    else  
        P = plot(ite, errores, xlabel="Iteraciones", ylabel="Error", label="Convergencia")
    end
    display(P)
    # Se crea una tabla con los resultados de Vn
    col_nodos = collect(1:num_nodes-1)
    df_resultados = DataFrame(Nodo = 1:num_nodes-1, 
                           Magnitud = abs.(Vn), 
                           Fase = angle.(Vn) * 180 / π)

    pretty_table(df_resultados, 
             header=["Nodo", "Voltaje (|V|)", "Fase (°)"], 
             alignment=:c)
    return Vn    
end

lines = DataFrame(CSV.File("ACTIVIDAD 3/lines.csv"))
nodes = DataFrame(CSV.File("ACTIVIDAD 3/nodes.csv"))
Vn = Flujo_punto_fijo(lines, nodes)

## -------------------------------------------------------------------------------------------------##
## AJUSTE DE DATAFRAMES.

# === Cargar Datos ===
SolarData = DataFrame(CSV.File("ACTIVIDAD 3/SolarData.csv"))
Demandas = DataFrame(CSV.File("ACTIVIDAD 3/Demandas.csv"))

# === Ajustar la base de datos SolarData a una matriz Muestras x Días ===
SolarData_df2 = DataFrame()
for i = 1:288:nrow(SolarData)
    fecha = string(SolarData.Fecha[i])  # Convertir la fecha a string para usarla como nombre de columna
    SolarData_df2[!, fecha] = SolarData.Potencia[i:i+287]  
end

# === Modificar el intervalo de tiempo de ambos DataFrame ===
intervalo = 60  # Puedes cambiarlo a 5, 10, 15, 60 minutos
factor_agrupacion = intervalo ÷ 5

# Convertir el DataFrame a matriz numérica
SolarData_matriz = Matrix(SolarData_df2)

# Agrupar por promedio para ajustar el intervalo de tiempo
SolarData_ajustada = [mean(SolarData_matriz[i:i+factor_agrupacion-1, :], dims=1) |> vec for i in 1:factor_agrupacion:288]
SolarData_ajustada = hcat(SolarData_ajustada...)'  # Convertir la lista en una matriz

# Convertir a DataFrame con las mismas fechas como nombres de columna
SolarData_df = DataFrame(SolarData_ajustada, Symbol.(names(SolarData_df2)))
CSV.write("solardata_matriz.csv", SolarData_df)

# === Aplicar K-means ===
K = 4
result = kmeans(SolarData_ajustada, K)
centroides = result.centers
dias_tipicos = DataFrame(centroides, ["Cluster_$(i)" for i in 1:K])
CSV.write("dias_tipicos.csv", dias_tipicos)

# === Graficar Generación Solar ===
theme(:dark)
p1 = plot(size=(800, 600), legend=false)
plot!(p1, eachcol(SolarData_ajustada), alpha=0.5)
title!("Potencia Generada por día para los 365 días")
xlabel!("Intervalos de $(intervalo) min")
ylabel!("Potencia Generada")

# === Graficar Días Típicos ===
p2 = plot(size=(800, 600), legend=false)
plot!(p2, eachcol(centroides), alpha=0.5)
title!("Potencia Generada para 4 casos típicos")
xlabel!("Intervalos de $(intervalo) min")
ylabel!("Potencia Generada")
display(plot(p1, p2, layout=(2,1)))

## -------------------------------------------------------------------------------------------------##
# === Ajuste de la matriz Demandas ===

# Calcular los máximos antes del promedio
max_values_before = maximum.(eachcol(Demandas))
println("Máximos por columna antes del promedio: ", max_values_before)

# === Promedio de Demanda según el intervalo definido ===
Demandas.grupo = div.(0:(nrow(Demandas)-1), intervalo) .+1
Demandas_prom = combine(groupby(Demandas, :grupo), names(Demandas, Not(:grupo)) .=> mean)
select!(Demandas_prom, Not(:grupo))

# Calcular los máximos después del promedio
max_values_after = maximum.(eachcol(Demandas_prom))
println("Máximos por columna después del promedio: ", max_values_after)

# === Ajustar Nombres de Columnas para nodos PQ ===
nodos_PQ = findall(x -> x != 0, nodes.PLOAD)  # Indices de nodos PQ
encabezados = ["N$(nodos_PQ[i])" for i in 1:length(nodos_PQ)]

# Verificar que el número de columnas coincida con los nodos PQ
if length(encabezados) != size(Demandas_prom, 2)
    error("El número de columnas en Demandas_prom no coincide con la cantidad de nodos PQ.")
end

# Renombrar columnas
rename!(Demandas_prom, Symbol.(encabezados))

# === Guardar el DataFrame de demanda promedio en p.u. ===
CSV.write("Demandas_$(intervalo)min.csv", Demandas_prom)

# === Convertir DataFrame a Matriz (SIN multiplicar por pbase) ===
matriz_dem_prom = Matrix(Demandas_prom)  # No se multiplica aquí

# === Graficar Demanda Promedio Escalada (Multiplicando por pbase solo en la gráfica) ===
pbase = 100
p4 = plot(size=(800, 600), legend=false)
for i in 1:(ncol(Demandas_prom))
    plot!(p4, matriz_dem_prom[:, i] * pbase, label=false, alpha=0.5)  # Multiplicamos aquí
end
title!("Potencia Demandada por día")
xlabel!("Intervalos de $intervalo min")
ylabel!("Potencia Demandada (MW)")
display(p4)

# -------------------------------------------------------------------------------------------------##
# FLUJO CUASI-DINAMICO

function Flujo_cuasi_din(lines, nodes, Demandas_prom, dias_tipicos)
    num_nodes = nrow(nodes)
    Ynn, Yns = calcular_Ynn_Yns(lines, nodes)
    
    if rank(Ynn) < size(Ynn, 1)
        error("La matriz Ynn no es invertible. Verifica la topología del sistema.")
    end
    
    iYnn = inv(Ynn)

    for i in 1:ncol(dias_tipicos)   
        Vns = DataFrame()  

        for j in 1:nrow(dias_tipicos)  
            # Inicialización
            Vn = ones(num_nodes-1) + 1im*zeros(num_nodes-1)
            Vs = 1 + 0*1im
            Sn = zeros(num_nodes-1)

            for k in 1:ncol(Demandas_prom)  
                nodo = parse(Int64, names(Demandas_prom)[k][2:end]) - 1
                if (nodo + 1) == 34  
                    Sn[nodo] = dias_tipicos[j, i] - Demandas_prom[j, k]
                else
                    Sn[nodo] = -Demandas_prom[j, k]
                end
            end

            # Iteración punto fijo
            for iter = 1:4
                Vn = iYnn * (conj.(Sn ./ Vn) .- Yns * Vs)
            end

            # Guardar la magnitud y ángulo
            Vn_mag = abs.(Vn)
            Vn_ang = angle.(Vn) * (180 / π)  # Convierte a grados

            # Crear columnas en la tabla
            nombre_mag = "Mag_$(SolarData.Hora[j])"
            nombre_ang = "Ang_$(SolarData.Hora[j])"

            Vns[!, nombre_mag] = Vn_mag
            Vns[!, nombre_ang] = Vn_ang
        end
        
        # Guardar en CSV
        nombre_archivo = "Vn_caso__$i.csv"
        CSV.write(nombre_archivo, Vns)
    end
end


# Se ejecuta el flujo cuasi dinamico
Flujo_cuasi_din(lines, nodes, Demandas_prom, dias_tipicos)


function plot_voltage_heatmaps(file_name)
    df = CSV.read(file_name, DataFrame)

    # Filtrar columnas de magnitud y ángulo correctamente
    mag_cols = filter(col -> occursin("Mag_", col), names(df))
    ang_cols = filter(col -> occursin("Ang_", col), names(df))

    # Verificar si encontró columnas
    if isempty(mag_cols) || isempty(ang_cols)
        error("No se encontraron columnas de magnitud o ángulo en el archivo.")
    end

    # Crear DataFrames de magnitud y ángulo
    mag_df = df[:, mag_cols]
    ang_df = df[:, ang_cols]

    # Extraer tiempos desde los nombres de columna
    tiempos = replace.(mag_cols, "Mag_" => "")
    tiempos_minutos = [parse(Int, split(t, ":")[1]) * 60 + parse(Int, split(t, ":")[2]) for t in tiempos]

    # Convertir DataFrames a matrices
    Z_mag = Matrix(mag_df)  # No transponemos aún
    Z_ang = Matrix(ang_df)

    # Coordenadas de los ejes
    x_vals = tiempos_minutos  # Tiempo en minutos
    y_vals = 1:size(Z_mag, 1)  # Índices de nodos

    # Ajustar dimensiones si es necesario
    if size(Z_mag) != (length(y_vals), length(x_vals))
        Z_mag = Z_mag'
    end
    if size(Z_ang) != (length(y_vals), length(x_vals))
        Z_ang = Z_ang'
    end

    # Verificación de tamaños
    @show size(Z_mag), size(Z_ang), length(y_vals), length(x_vals)

    # Crear dos mapas de calor
    p1 = heatmap(x_vals, y_vals, Z_mag, xlabel="Tiempo (min)", ylabel="Nodo", title="Mapa de Calor de Magnitud de Tensión", color=:viridis)
    xticks!(0:150:maximum(x_vals))
    
    p2 = heatmap(x_vals, y_vals, Z_ang, xlabel="Tiempo (min)", ylabel="Nodo", title="Mapa de Calor de Ángulo de Tensión", color=:plasma)
    # Personalizar los ticks del eje X
    xticks!(0:150:maximum(x_vals))

    # Mostrar gráficos
    display(plot(p1, p2, layout=(2,1)))
end

# Prueba con un archivo
plot_voltage_heatmaps("Vn_caso__1.csv")