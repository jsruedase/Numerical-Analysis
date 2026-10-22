### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ c43e6158-d1a4-4342-8ba4-2e4c9899e077
md"""
# Tarea opcional Ecuaciones no lineales
Juan Sebastián Rueda Segura
"""

# ╔═╡ c9197620-cdc0-11f1-b0a8-67818a1a5eec
md"""
## Problema 1
**Problema (Kelley)** Considere el método de Shamanskii de orden $m$ definido de la siguiente manera. Dado $x_n$ calcule $x_{n+1}$ usando $m$ pasos del método de la cuerda, esto es,

$z_1=x_n-f(x_n)/f'(x_n)$

$z_{j+1}=z_{j}-f(z_j)/f'(x_n), \hspace{0.5cm}1\leq j\leq m-1$

$x_{n+1}=z_{m}.$

Observe que se avanza m pasos sin actualizar la derivada en el punto. Implemente el método para un order $m$ definido por el usuario y verifique su código con una ecuación escalar. Para la ecuación usada, ¿como depende el tiempo de computo de $m$? """

# ╔═╡ 5b29d15c-4978-4d32-aec0-4791d48245e1
md"""
### Solución:
Vemos que el método propuesto tiene como subrutina el método de Newton pero donde la derivada no se actualiza en m pasos. Así que primero necesitamos la implementación del método de Newton y crearemos una función para calcular los $z_i$, con una función que dé un paso pasándole el valor de la derivada por parámetro para que se fije en m iteraciones el valor de la derivada.

Implementación del método de Newton:
"""

# ╔═╡ 8b982a47-28fe-4557-8eac-04f9f80141fa
function minewton(f,fp,x0,myeps,mydel,max_iter)
    #=f(x)=0 usando newton xn=xa-f(xa)/fp(xa)
    ----------
    f : función
    fp : Derivada
    x0 : Aproximación inicial
    epsilon : Tolerancia residuo.
    M : Número maximo de iteraciones
    ----------
    =#
    x = x0
    n=1
    while n<=max_iter
        fx = f(x)
        if abs(fx) < myeps
            print("residuo f(x(",n,"))=",fx)
            return x
        end
        fpx = fp(x)
        if fpx == 0
            print("Derivada nula. Error!")
            return None
        end
        x = x - fx/fpx
        println( "x(",n,")=",x)
        if abs(x-x0) < mydel
            print("incremento delta(",n,")=",fx)
            return x
        end
        x0=x
        n=n+1
    end
    if n==max_iter
        print("Número maximo de iteraciones")
    end 
    return x0
end

# ╔═╡ 3258697f-b81c-4546-91a8-3a18be1d7bfa
md"""Implementación de la función del paso:"""

# ╔═╡ bc062b1f-1277-4a6d-b7ff-fa1eebea38fb
function paso_shamanskii(f, x, fx, d, m)
    # x: iterado actual, fx = f(x), d = f'(x) 
    z = x - fx/d            # z_1: paso de Newton
    for j in 1:m-1
        z = z - f(z)/d      
    end
    return z
end

# ╔═╡ 65491e8e-b357-4eba-8d2b-c13ff01afb64
md"""Implementación de la función requerida:"""

# ╔═╡ 547192f7-8129-442b-bcd8-41b46170750b
function mishamanskii(f, fp, x0, m, myeps, mydel, max_iter; verbose=true)
    #= f(x)=0 con Shamanskii de orden m
    f, fp : función y derivada
    x0    : aproximación inicial
    m     : número de pasos por iteración (m=1 es Newton)
    Devuelve (x, número de iteraciones, historial de iterados) =#
    x = x0
    hist = [float(x0)]
    for n in 1:max_iter
        fx = f(x)
        if abs(fx) < myeps
            verbose && println("residuo |f(x)| < eps tras ", n-1, " iteraciones")
            return x, n-1, hist
        end
        d = fp(x)
        d == 0 && error("Derivada nula")
        z = paso_shamanskii(f, x, fx, d, m)
        push!(hist, z)
        verbose && println("x(", n, ") = ", z)
        if abs(z - x) < mydel
            verbose && println("incremento < delta en la iteración ", n)
            return z, n, hist
        end
        x = z
    end
    verbose && println("Número máximo de iteraciones")
    return x, max_iter, hist
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.13.0"
manifest_format = "2.1"
project_hash = "71853c6197a6a7f222db0f1978c7cb232b87c5ee"

[deps]
"""

# ╔═╡ Cell order:
# ╟─c43e6158-d1a4-4342-8ba4-2e4c9899e077
# ╟─c9197620-cdc0-11f1-b0a8-67818a1a5eec
# ╟─5b29d15c-4978-4d32-aec0-4791d48245e1
# ╠═8b982a47-28fe-4557-8eac-04f9f80141fa
# ╟─3258697f-b81c-4546-91a8-3a18be1d7bfa
# ╠═bc062b1f-1277-4a6d-b7ff-fa1eebea38fb
# ╟─65491e8e-b357-4eba-8d2b-c13ff01afb64
# ╠═547192f7-8129-442b-bcd8-41b46170750b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
