### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ a1000001-0000-4000-8000-000000000001
using PlutoUI, Plots, Roots, Printf

# ╔═╡ a1000001-0000-4000-8000-000000000002
PlutoUI.TableOfContents(title="Problemas 1 y 3", aside=true)

# ╔═╡ a1000001-0000-4000-8000-000000000003
md"""# Ecuaciones no lineales: Problemas 1 y 3

Se resuelven el Problema 1 (método de Shamanskii) y el Problema 3 (tasa $\omega$ del plano inclinado) usando las herramientas del cuaderno del curso: Newton (`minewton`), punto fijo, bisección y el paquete `Roots`."""

# ╔═╡ a1000001-0000-4000-8000-000000000004
md"""# Problema 1: método de Shamanskii

Dado $x_n$, se hace un paso de Newton y luego $m-1$ pasos de la cuerda **sin actualizar la derivada**:

$z_1=x_n-\frac{f(x_n)}{f'(x_n)},\qquad z_{j+1}=z_j-\frac{f(z_j)}{f'(x_n)},\qquad x_{n+1}=z_m.$

**Justificación.** Cada paso de cuerda es una iteración de punto fijo con $F_n(z)=z-f(z)/f'(x_n)$, que es contractiva cerca de la raíz $r$ porque $F_n'(z)=1-f'(z)/f'(x_n)\approx 0$. Si $e=x_n-r$, por Taylor cada paso de cuerda multiplica el error por $1-f'(r)/f'(x_n)=O(e)$. El paso de Newton deja error $O(e^2)$; tras $m-1$ pasos de cuerda queda $O(e^{m+1})$: **orden $m+1$** con una sola evaluación de $f'$ por iteración. Con $m=1$ es Newton."""

# ╔═╡ a1000001-0000-4000-8000-000000000005
md"""Un paso del método (se separa para reutilizarlo en precisión extendida):"""

# ╔═╡ a1000001-0000-4000-8000-000000000006
function paso_shamanskii(f, x, fx, d, m)
    # x: iterado actual, fx = f(x), d = f'(x) (derivada congelada)
    z = x - fx/d            # z_1: paso de Newton
    for j in 1:m-1
        z = z - f(z)/d      # pasos de la cuerda con la misma derivada
    end
    return z
end

# ╔═╡ a1000001-0000-4000-8000-000000000007
md"""El método completo, con los mismos criterios de parada que `minewton` del cuaderno (residuo $|f(x)|<\epsilon$ o incremento $|x_{n+1}-x_n|<\delta$):"""

# ╔═╡ a1000001-0000-4000-8000-000000000008
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

# ╔═╡ a1000001-0000-4000-8000-000000000009
md"""Función de Newton del cuaderno, para comparar (solo se cambió `None` por `nothing`):"""

# ╔═╡ a1000001-0000-4000-8000-000000000010
function minewton(f,fp,x0,myeps,mydel,max_iter)
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
            return nothing
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
    return x0
end

# ╔═╡ a1000001-0000-4000-8000-000000000011
md"""## Verificación con $f(x)=\exp(x)-\sin(x)$

Misma ecuación y punto inicial del cuaderno ($x_0=-3.5$). La raíz de referencia se obtiene con `fzero` en $[-4,-3]$."""

# ╔═╡ a1000001-0000-4000-8000-000000000012
f(x) = exp(x) - sin(x)

# ╔═╡ a1000001-0000-4000-8000-000000000013
fp(x) = exp(x) - cos(x)

# ╔═╡ a1000001-0000-4000-8000-000000000014
r_ref = fzero(f, -4, -3)

# ╔═╡ a1000001-0000-4000-8000-000000000015
md"""Con $m=1$ el método debe reproducir a `minewton`:"""

# ╔═╡ a1000001-0000-4000-8000-000000000016
minewton(f, fp, -3.5, 1e-15, 1e-14, 50)

# ╔═╡ a1000001-0000-4000-8000-000000000017
mishamanskii(f, fp, -3.5, 1, 1e-15, 1e-14, 50)[1]

# ╔═╡ a1000001-0000-4000-8000-000000000018
md"""Errores $|x_n-r|$ para $m=1,2,3,4$:"""

# ╔═╡ a1000001-0000-4000-8000-000000000019
let
    for m in 1:4
        x, n, hist = mishamanskii(f, fp, -3.5, m, 1e-15, 1e-14, 50; verbose=false)
        println("m = ", m, "   (", n, " iteraciones)")
        for (i, xi) in enumerate(hist)
            @printf("   n = %d   x = %.16f   |x - r| = %.2e\n", i-1, xi, abs(xi - r_ref))
        end
    end
end

# ╔═╡ a1000001-0000-4000-8000-000000000020
let
    p = plot(yscale=:log10, xlabel="iteración n", ylabel="|x_n - r|",
             title="Shamanskii con x0 = -3.5", legend=:bottomleft)
    for m in 1:4
        _, _, hist = mishamanskii(f, fp, -3.5, m, 1e-15, 1e-14, 50; verbose=false)
        errs = max.(abs.(hist .- r_ref), 1e-17)   # evita log(0)
        plot!(p, 0:length(errs)-1, errs, marker=:circle, lw=2, label="m = $m")
    end
    p
end

# ╔═╡ a1000001-0000-4000-8000-000000000021
md"""## Orden de convergencia

En doble precisión se llega al redondeo ($\approx 10^{-16}$) en 2 a 4 pasos, insuficiente para estimar el orden. Se repite con `BigFloat` a 2000 bits ($\approx 600$ dígitos) y se estima

$p_n\approx\frac{\ln(e_{n+1}/e_n)}{\ln(e_n/e_{n-1})}.$

Se espera $p_n\to m+1$."""

# ╔═╡ a1000001-0000-4000-8000-000000000022
ordenes = setprecision(BigFloat, 2000) do
    tolref = BigFloat(10)^-590
    rbig = mishamanskii(f, fp, BigFloat(-3.18), 1, tolref, tolref, 100; verbose=false)[1]
    tol = BigFloat(10)^-500
    resultado = [Float64[] for _ in 1:4]
    for m in 1:4
        x = BigFloat(-3.5)
        errs = [abs(x - rbig)]
        while errs[end] > tol && length(errs) < 15
            x = paso_shamanskii(f, x, f(x), fp(x), m)
            push!(errs, abs(x - rbig))
        end
        errs = filter(e -> e > tol, errs)
        resultado[m] = [Float64(log(errs[i+1]/errs[i]) / log(errs[i]/errs[i-1]))
                        for i in 2:length(errs)-1]
    end
    resultado
end

# ╔═╡ a1000001-0000-4000-8000-000000000023
for m in 1:4
    println("m = ", m, "   órdenes estimados: ", round.(ordenes[m], digits=3))
end

# ╔═╡ a1000001-0000-4000-8000-000000000024
md"""**Conclusión.** Los órdenes se estabilizan en $2,3,4,5$, es decir $m+1$ (el primer valor de cada fila se desvía porque $x_0$ aún no está en la zona asintótica). El costo por iteración es $m$ evaluaciones de $f$ y una de $f'$."""

# ╔═╡ a1000001-0000-4000-8000-000000000025
md"""# Problema 3: tasa $\omega$ del plano inclinado

Que el objeto recorra 1 m en 1 s equivale a resolver

$g(\omega)=s(1,\omega)-1=\frac{9.8}{2\omega^2}\big[\sinh\omega-\sin\omega\big]-1=0.$

**Unicidad.** Por Taylor, $\sinh\omega-\sin\omega=2\left(\frac{\omega^3}{3!}+\frac{\omega^7}{7!}+\cdots\right)$, luego $s(1,\omega)=9.8\left(\frac{\omega}{3!}+\frac{\omega^5}{7!}+\cdots\right)$: todos los coeficientes son positivos, así que $g$ es estrictamente creciente en $\omega>0$, con $g(0^+)=-1$. Hay una única raíz positiva. La serie también da una estimación inicial: $9.8\,\omega/6=1\Rightarrow\omega\approx 0.612$."""

# ╔═╡ a1000001-0000-4000-8000-000000000026
s(t, ω) = 9.8/(2ω^2) * (sinh(ω*t) - sin(ω*t))

# ╔═╡ a1000001-0000-4000-8000-000000000027
gw(ω) = s(1, ω) - 1

# ╔═╡ a1000001-0000-4000-8000-000000000028
md"""Derivada (regla del cociente), necesaria para Newton:

$g'(\omega)=\frac{9.8}{2}\cdot\frac{\omega(\cosh\omega-\cos\omega)-2(\sinh\omega-\sin\omega)}{\omega^3}.$"""

# ╔═╡ a1000001-0000-4000-8000-000000000029
gwp(ω) = 9.8/2 * (ω*(cosh(ω) - cos(ω)) - 2*(sinh(ω) - sin(ω))) / ω^3

# ╔═╡ a1000001-0000-4000-8000-000000000030
md"""## Localización"""

# ╔═╡ a1000001-0000-4000-8000-000000000031
let
    ω = 0.01:0.001:2
    plot(ω, 0*ω, ls=:dash, label="y=0", lw=3)
    plot!(ω, gw.(ω), label="g(ω)", lw=4, xlabel="ω")
end

# ╔═╡ a1000001-0000-4000-8000-000000000032
[(ω, gw(ω)) for ω in (0.5, 0.6, 0.7, 1.0)]

# ╔═╡ a1000001-0000-4000-8000-000000000033
md"""Hay cambio de signo en $[0.5,1]$ (teorema del valor intermedio).

## Bisección: garantía de 12 decimales

En el algoritmo de bisección del cuaderno, $c=a+e$ es el punto medio de un intervalo de longitud $2e$ que contiene la raíz, así que $|\omega^*-c|\le e=(b-a)/2^k$. Para 12 decimales exactos basta $e<0.5\times10^{-12}$, o sea $k\ge 40$. Se usa $\delta=10^{-14}$ y $\epsilon=10^{-16}$ para que pare por la longitud del intervalo."""

# ╔═╡ a1000001-0000-4000-8000-000000000034
function mibiseccion(fun, a, b, M, delta, epsilon)
    # Algoritmo de bisección del cuaderno, en forma de función
    u, v = fun(a), fun(b)
    sign(u) == sign(v) && error("f(a) y f(b) tienen el mismo signo")
    e = b - a
    c, w = a, u
    for k in 1:M
        e = e/2
        c = a + e
        w = fun(c)
        if abs(e) < delta || abs(w) < epsilon
            return c, e, w, k
        end
        if sign(w) != sign(u)
            b, v = c, w
        else
            a, u = c, w
        end
    end
    return c, e, w, M
end

# ╔═╡ a1000001-0000-4000-8000-000000000035
ω_bis, e_bis, w_bis, k_bis = mibiseccion(gw, 0.5, 1.0, 100, 1e-14, 1e-16)

# ╔═╡ a1000001-0000-4000-8000-000000000036
@printf("k = %d,  ω ≈ %.16f,  cota del error e = %.2e,  g(ω) = %.2e", k_bis, ω_bis, e_bis, w_bis)

# ╔═╡ a1000001-0000-4000-8000-000000000037
md"""## Newton, Shamanskii y `Roots`

Newton (`minewton`) desde $\omega_0=1$: los incrementos se elevan aproximadamente al cuadrado en cada paso."""

# ╔═╡ a1000001-0000-4000-8000-000000000038
ω_newton = minewton(gw, gwp, 1.0, 1e-16, 1e-14, 50)

# ╔═╡ a1000001-0000-4000-8000-000000000039
ω_sham = mishamanskii(gw, gwp, 1.0, 2, 1e-16, 1e-14, 50)[1]

# ╔═╡ a1000001-0000-4000-8000-000000000040
ω_roots = fzero(gw, 0.5, 1.0)

# ╔═╡ a1000001-0000-4000-8000-000000000041
md"""## Comparación y resultado"""

# ╔═╡ a1000001-0000-4000-8000-000000000042
let
    for (nombre, ω) in (("Bisección ", ω_bis), ("Newton    ", ω_newton),
                        ("Shamanskii", ω_sham), ("fzero     ", ω_roots))
        @printf("%s  ω = %.15f   |ω - ω_bis| = %.1e   s(1,ω) = %.15f\n",
                nombre, ω, abs(ω - ω_bis), s(1, ω))
    end
    println("\nω con 12 decimales: ", round(ω_newton, digits=12))
end

# ╔═╡ a1000001-0000-4000-8000-000000000043
md"""**Conclusión.** Todos los métodos coinciden en $\omega\approx 0.612142570703094$. La cota de bisección ($e\approx7\times10^{-15}$) y la diferencia entre métodos ($<10^{-15}$) son mucho menores que $0.5\times10^{-12}$, y los dígitos 13–15 (094) no están cerca de alterar el redondeo. Por tanto

$\boxed{\omega = 0.612142570703\ \text{s}^{-1}}$

con 12 decimales exactos. La cancelación en $\sinh\omega-\sin\omega$ ($\approx0.65-0.57$) cuesta menos de un dígito y no afecta este resultado."""

# ╔═╡ Cell order:
# ╠═a1000001-0000-4000-8000-000000000001
# ╟─a1000001-0000-4000-8000-000000000002
# ╟─a1000001-0000-4000-8000-000000000003
# ╟─a1000001-0000-4000-8000-000000000004
# ╟─a1000001-0000-4000-8000-000000000005
# ╠═a1000001-0000-4000-8000-000000000006
# ╟─a1000001-0000-4000-8000-000000000007
# ╠═a1000001-0000-4000-8000-000000000008
# ╟─a1000001-0000-4000-8000-000000000009
# ╠═a1000001-0000-4000-8000-000000000010
# ╟─a1000001-0000-4000-8000-000000000011
# ╠═a1000001-0000-4000-8000-000000000012
# ╠═a1000001-0000-4000-8000-000000000013
# ╠═a1000001-0000-4000-8000-000000000014
# ╟─a1000001-0000-4000-8000-000000000015
# ╠═a1000001-0000-4000-8000-000000000016
# ╠═a1000001-0000-4000-8000-000000000017
# ╟─a1000001-0000-4000-8000-000000000018
# ╠═a1000001-0000-4000-8000-000000000019
# ╠═a1000001-0000-4000-8000-000000000020
# ╟─a1000001-0000-4000-8000-000000000021
# ╠═a1000001-0000-4000-8000-000000000022
# ╠═a1000001-0000-4000-8000-000000000023
# ╟─a1000001-0000-4000-8000-000000000024
# ╟─a1000001-0000-4000-8000-000000000025
# ╠═a1000001-0000-4000-8000-000000000026
# ╠═a1000001-0000-4000-8000-000000000027
# ╟─a1000001-0000-4000-8000-000000000028
# ╠═a1000001-0000-4000-8000-000000000029
# ╟─a1000001-0000-4000-8000-000000000030
# ╠═a1000001-0000-4000-8000-000000000031
# ╠═a1000001-0000-4000-8000-000000000032
# ╟─a1000001-0000-4000-8000-000000000033
# ╠═a1000001-0000-4000-8000-000000000034
# ╠═a1000001-0000-4000-8000-000000000035
# ╠═a1000001-0000-4000-8000-000000000036
# ╟─a1000001-0000-4000-8000-000000000037
# ╠═a1000001-0000-4000-8000-000000000038
# ╠═a1000001-0000-4000-8000-000000000039
# ╠═a1000001-0000-4000-8000-000000000040
# ╟─a1000001-0000-4000-8000-000000000041
# ╠═a1000001-0000-4000-8000-000000000042
# ╟─a1000001-0000-4000-8000-000000000043
