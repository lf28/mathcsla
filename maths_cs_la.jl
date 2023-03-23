### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 9edaaf3c-2973-11ed-193f-d500cbe5d239
begin
	using LinearAlgebra
	using Statistics
	using PlutoUI
	using PlutoTeachingTools
	# import PlutoTeachingTools as PT
	# import CalculusWithJulia as CJ
	using Plots
	# using PlotlyBase
	using LaTeXStrings
	using Latexify
	# using PlutoPlotly
end

# ╔═╡ 38967916-87c6-4bb8-8ca6-f843b482e2ca
TableOfContents(; indent=true, depth=2, aside=true)

# ╔═╡ cb9259ce-beab-4e94-88a9-66057a88b633
present_button() 

# ╔═╡ 7f825237-767e-4b3e-927f-e3b0d6779316
md"# Linear Algebra"

# ╔═╡ 457b21c1-6ea9-4983-b724-fb5cbb69739d
md"""


$(Resource("https://www.st-andrews.ac.uk/assets/university/brand/logos/standard-vertical-black.png", :width=>130, :align=>"right"))

Lei Fang (lf28@st-andrews.ac.uk)

*School of Computer Science*

*University of St Andrews, UK*

*Jan 2023*

"""

# ╔═╡ 57945aa2-b984-407a-89b1-57b0b0a2cd75
md"""

## Linear algebra

Some basic concepts covered here


* Vectors
* Vector operations, Norms 
* Linear independence, span
* Matrix
* Matrix operations
* Rank

"""

# ╔═╡ b18cf2b6-0e45-4158-9742-2d4e48ec3765
md"""

## A few words on `Julia` and `Pluto`

This note is written in [`Julia`](https://julialang.org) and [`Pluto`](https://plutojl.org)

* `Julia` is a very popular programming language for numerical computations, such as machine learning, Bayesian inference *etc*

* `Julia` can be as fast as `Fortran` and `C` but with simple syntax such as `Python/R` and `Matlab`


`Pluto` is a notebook-like computational environment for `Julia`
* easy to use 
* interactive 


You are recommended to download the notebook (top right of this page) and run them locally
* *maths* and *computing* should really go hand-in-hand
* *programming* should help you understand the maths greatly



Here is a very quick and short [introduction](https://lf28.github.io/BayesianModelling/section0_julia.html) to `Julia` 

* covers some essential and useful features of the language
"""

# ╔═╡ 9ca74f71-e1f6-43c7-b08a-7afb7d3dfe0d
md"""
## Vectors


A vector is an ordered list of ``n`` *scalars*, e.g.

```math

\mathbf{a} = \begin{bmatrix} 2 \\ 1\end{bmatrix}, \mathbf{b} =\begin{bmatrix} -1\\2 \end{bmatrix}
```

Visually, a vector has
* a direction 
* also a length, measured by **Euclidean distance** between the point and origin ``\mathbf{0} =[0,0]^\top``
  * due to *Pythagorean theorem*
```math
\text{dist}(\mathbf{a}, \mathbf{0}) = \sqrt{2^2 + 1^2} =\sqrt{5}
```


$(begin

aside(tip(md"All vectors by default are **column vectors**. You may see notations like

```math 
\mathbf{a} = [2,1]^\top
```
* it means a *column* vector
* ``\top`` means transpose
* to save space in text") )
end)


"""

# ╔═╡ bd2498df-069b-4ea0-9e44-738142a3080e
a = [2, 1]; b = [-1,2]

# ╔═╡ 5d4e270a-b480-40d1-97ac-d7eaa4700766
let
	gr()
 	plot(xlim=[-4,4], ylim=[-1,4],ratio=1, framestyle=:origin)
	# quiver([0,0,0],[0,0,0],quiver=([1,1,1],[1,2,3]))
	oo = [0,0]
	quiver!([0], [0], quiver=([a[1]], [a[2]]), lw=2)
	quiver!([0], [0],  quiver=([b[1]], [b[2]]), lw=2)
	plot!([a[1],a[1]], [a[2],0], ls=:dash, lc=:gray, lw=1, label="")
	plot!([a[1],0], [a[2],a[2]], ls=:dash, lc=:gray, lw=1, label="")
	plot!([b[1],b[1]], [b[2],0], ls=:dash, lc=:gray, lw=1, label="")
	plot!([b[1],0], [b[2],b[2]], ls=:dash, lc=:gray, lw=1, label="")
	annotate!(a[1],a[2], text(L"a", :top))
	annotate!(b[1],b[2], text(L"b", :bottom))
	θ = acos(dot(a, [1,0])/ norm(a)) * 180/π
	annotate!(0.5*a[1], 0.5*a[2], text(L"d=\sqrt{%$(a[1])^2+%$(a[2])^2}=\sqrt{%$(sum(a.^2))}", 8, :dark, :top, rotation = θ ))
	θb = acos(dot(b, [1,0])/ norm(b)) * 180/π
	annotate!(0.5*b[1], 0.5*b[2], text(L"d=\sqrt{(%$(b[1]))^2+%$(b[2])^2}=\sqrt{%$(sum(a.^2))}", 8, :dark, :bottom, rotation = θb ))
end

# ╔═╡ 6d4861b6-4a0d-4aac-8f08-e6861e2ecf70
md"""
## Vectors and scalars

Recall how we were taught numbers (scalars) in school:  **Real Axis** 

$(Resource("https://mathworld.wolfram.com/images/eps-svg/RealLine_1201.svg", :width=>500, :align=>"middle"))

* **Vectors** visualisation is a generalisation:  ``R \Rightarrow R^n``

"""

# ╔═╡ 3ea11e49-9edd-43ac-9908-01d766c331a1
md"""

## Vectors in R³

Vectors in ``R^n`` where ``n=3``



```math

 \mathbf{c} =\begin{bmatrix} 2 \\ 0\\0\end{bmatrix}, \mathbf{d} = \begin{bmatrix}2 \\ 4\\0 \end{bmatrix}, \mathbf{f} = \begin{bmatrix} 1 \\ 1\\4\end{bmatrix},
```

* ``\textcolor{red}{\mathbf{c}}, \textcolor{green}{\mathbf{d}}`` live in the ``xy`` -- plane
* ``\textcolor{blue}{\mathbf{f}}``, however, "sticks out"
"""

# ╔═╡ 4c1e7871-a931-4781-a200-304a2ef253e1
cv=[2,0,0]; dv=[2,4,0]; fv = [1,1,4]; # create vectors in Julia

# ╔═╡ fd80776b-1bc2-466a-b29c-28df81cbee3f
md"""

## Higher dimensional vectors


Vectors in ``R^n`` 

```math

\mathbf{a} =\begin{bmatrix} a_1\\ a_2\\ \vdots\\ a_n\end{bmatrix} 

```

* cannot be plotted exactly, 
* but can still be "visualised" in your mind



**Example:** a short document can be summarised into a vector based on a dictionary

!!! infor ""
	Orientation **Week** 2022. **Welcome** to the School of **Computer Science** orientation **week**. Please visit the page on information for current **computer science** students to access the orientation **week** timetable and other resources relevant to your year of **study**.


* a small dictionary (left) and word count vector (right)

```math

\begin{array}{l}
\texttt{welcome}\\
\texttt{study}\\
\texttt{dog}\\
\texttt{computer}\\
\texttt{science}\\
\texttt{week}\\
\end{array}

\;\;\;\;\;\;\;\;

\begin{bmatrix}
1 \\
1\\
0 \\
2 \\
2 \\
3\\
\end{bmatrix}

```
* a document then can be viewed as a vector in ``R^n`` 
* real world vectors are much longer (larger dictionary!)

"""

# ╔═╡ c56ec750-0514-4796-a594-325b930bf4d2
md"""

## Vector operations -- addition

Vectors of the same length can be added together:

```math
\mathbf{a}+\mathbf{b} = \begin{bmatrix}
           a_1 \\
           a_2 \\
           \vdots\\
           a_n
         \end{bmatrix} +  \begin{bmatrix}
           b_1 \\
           b_2 \\
           \vdots\\
           b_n
         \end{bmatrix} = \begin{bmatrix}
           a_1+b_1 \\
           a_2+b_2 \\
           \vdots\\
           a_n+b_n
         \end{bmatrix}

```


Example: ``\mathbf{b} = [2,0,0]^\top, \mathbf{c} =[3,2,0]^\top``

```math

\mathbf{b}+\mathbf{c}=\begin{bmatrix} 2 \\ 0\\0\end{bmatrix}+ \begin{bmatrix}3 \\ 2\\0 \end{bmatrix} = \begin{bmatrix} 2+3 \\ 0+2\\ 0+0\end{bmatrix} = \begin{bmatrix} 5 \\ 2\\ 0\end{bmatrix}
```

"""

# ╔═╡ 1dd53149-6749-4c86-b63a-eb801a827808
md"""

Addition visualisation: parallelogram rule (or "tip to tail")
* shift ``\mathbf{b}``'s *origin* to ``\mathbf{a}``'s *end* point (dashed vector) then connect from ``\mathbf{0}``
* or ``\mathbf{a}``'s tip connect to ``\mathbf{b}``'s tail

"""

# ╔═╡ 057908d5-7f6b-4976-84f4-eabe8ff88c33
let
	gr()
 	plot(xlim=[-1,4], ylim=[-1, 4.5], ratio=1, framestyle=:origin)
	# quiver([0,0,0],[0,0,0],quiver=([1,1,1],[1,2,3]))
	oo = [0,0]
	ab₊ = a +b 
	quiver!([0], [0], quiver=([a[1]], [a[2]]), lw=2)
	quiver!([0], [0],  quiver=([b[1]], [b[2]]), lw=2)
	quiver!([a[1]], [a[2]],  quiver=([b[1]], [b[2]]), lw=1, lc=2, ls=:dash)
	quiver!([0], [0],  quiver=([ab₊[1]], [ab₊[2]]), lw=2)
	annotate!(a[1],a[2], text(L"a", :top))
	annotate!(b[1],b[2], text(L"b", :bottom))
	annotate!(ab₊[1],ab₊[2], text(L"a+b", :bottom))
end

# ╔═╡ 8539909b-0f49-4aa3-a950-89bf1f09e696
md"""
## Zero, ones, and unit vectors

The **zero** vector

```math
\mathbf{0} =\begin{bmatrix} 0, 0, \ldots, 0\end{bmatrix}^\top
```

* ``\mathbf{0} + \mathbf{a} =\mathbf{a}``

The **one** vector

```math
\mathbf{1} =\begin{bmatrix} 1, 1, \ldots, 1\end{bmatrix}^\top
```


Standard basis vectors, *a.k.a.* unit vectors, denoted as ``\mathbf{e}_i``, *e.g.* ``\mathbf{e}_i\in R^3``:


```math

\mathbf{e}_1 = \begin{bmatrix}1 \\0 \\0 \end{bmatrix}, \mathbf{e}_2 = \begin{bmatrix}0 \\1 \\0 \end{bmatrix}; \mathbf{e}_3 = \begin{bmatrix}0 \\0 \\1 \end{bmatrix}

```


"""

# ╔═╡ b5199746-42a4-4462-81ef-d329482f5a3c
md"""
More generally, for ``R^n``

```math
\mathbf{e}_i=[e_{i1}, e_{i2}, \ldots, e_{in}]^\top,
```

where ``e_{ii}=1, e_{i\neq j} =0`` for ``i,j= 1,2, \ldots, n``


"""

# ╔═╡ 94027fb4-cbb0-46d3-abd2-795af9019cd7
md"""

## Vector operations -- scaling



$$k\cdot \mathbf{a} = k \cdot \begin{bmatrix}
           a_1 \\
           a_2 \\
           \vdots\\
           a_n
         \end{bmatrix} = \begin{bmatrix}
           k\times a_1 \\
           k\times a_2 \\
           \vdots\\
           k\times a_n
         \end{bmatrix},\;\; k\in R \text{ or a scalar}$$

* *geometrically*, scaling means shrinking or stretching a vector 
  * ``k>0``, the same direction
  * ``k<0``, the opposite direction
  * ``k=0``, we get ``\mathbf{0}``
* *arithmetically*, adding ``k`` copies of ``\mathbf{a}`` together

```math

k\cdot \mathbf{a} = \underbrace{\begin{bmatrix}
           a_1 \\
           a_2 \\
           \vdots\\
           a_n
         \end{bmatrix} + \begin{bmatrix}
           a_1 \\
           a_2 \\
           \vdots\\
         a_n
         \end{bmatrix} + \ldots + \begin{bmatrix}
           a_1 \\
           a_2 \\
           \vdots\\
         a_n
         \end{bmatrix}}_{k}
```

* ``k\cdot \mathbf{a} = \mathbf{a} \cdot k``, *i.e.* scaling is associative

"""

# ╔═╡ 46a8bf2d-5fea-46c3-ba65-4cb02e2279f2
aside(tip(md"If you run the notebook locally, try to move the radio botton to see the effect."))

# ╔═╡ a54dbf58-d082-440f-bc3d-ceebdfabbda6
md"""

**Example**
"""

# ╔═╡ 8c5201e7-a0a0-497b-9459-4e518add61de
@bind ka Slider(-2.2:0.1:2.3, default=1)

# ╔═╡ a37be5c2-2c38-4625-8ca4-3c37ffbc3bac
md"``k``=$(ka)"

# ╔═╡ d4d8448c-1148-4c56-a460-acc6279013ba
let
	gr()
 	plot(xlim=[-4.5,4.5], ylim=[-4.5, 4.5], ratio=1, framestyle=:origin)
	# quiver([0,0,0],[0,0,0],quiver=([1,1,1],[1,2,3]))
	oo = [0,0]
	# ab₊ = a +b 
	ka_ = ka *a
	# kb_ = kb * b
	quiver!([0], [0], quiver=([a[1]], [a[2]]), lw=2)
	quiver!([0], [0],  quiver=([b[1]], [b[2]]), lw=2)
	quiver!([0], [0],  quiver=([ka_[1]], [ka_[2]]), lw=2)
	# quiver!([0], [0],  quiver=([kb_[1]], [kb_[2]]), lw=2)
	# quiver!([a[1]], [a[2]],  quiver=([b[1]], [b[2]]), lw=2, lc=2, ls=:dash)
	# quiver!([0], [0],  quiver=([ab₊[1]], [ab₊[2]]), lw=2)
	annotate!(a[1],a[2], text(L"a", :top))
	annotate!(b[1],b[2], text(L"b", :bottom))
	annotate!(ka_[1],ka_[2], text(L"%$(ka)a", :bottom))
	# annotate!(kb_[1],kb_[2], text(L"%$(kb)b", :bottom))
	# annotate!(ab₊[1],ab₊[2], text(L"a+b", :bottom))
end

# ╔═╡ fc9fb8cd-6557-4eff-abab-30be01366c91
plt1, plt2, plt3=let
	gr()
 	plt1 = plot(xlim=[-.1,5], ylim=[-.1, 3], ratio=1, framestyle=:origin)
	plt2 = plot(xlim=[-.1,5], ylim=[-.1, 3], ratio=1, framestyle=:origin)
	plt3 = plot(xlim=[-.1,5], ylim=[-.1, 3], ratio=1, framestyle=:origin)
	oo = [0,0]
	a = [1,0]
	b = [1,1]
	c = [0, 1]
	# ab₊ = a +b 
	k₁ = 2.5
	k₂ = 1.5
	k₃ = 1
	ka = k₁ * a
	kb = k₂ * b
	kc = k₃ * c
	kab= ka + kb
	kabc = ka + kb + kc
	quiver!(plt1, [0], [0], quiver=([a[1]], [a[2]]), lc=1, lw=2)
	quiver!(plt1, [0], [0],  quiver=([b[1]], [b[2]]), lc=1, lw=2)
	quiver!(plt1, [0], [0],  quiver=([c[1]], [c[2]]), lc=1, lw=2)
	quiver!(plt2, [0], [0],  quiver=([ka[1]], [ka[2]]), lc=1, lw=2)
	quiver!(plt2, [0], [0],  quiver=([kb[1]], [kb[2]]), lc=1, lw=2)
	quiver!(plt3, [0], [0],  quiver=([kc[1]], [kc[2]]), lc=1, lw=2)
	quiver!(plt2, [0], [0],  quiver=([kab[1]], [kab[2]]), lc=2, lw=2)
	quiver!(plt3, [0], [0],  quiver=([kab[1]], [kab[2]]), lc=2, lw=2)
	quiver!(plt2, [ka[1]], [ka[2]],  quiver=([kb[1]], [kb[2]]), lw=2, lc=1, ls=:dash)

	quiver!(plt3, [kab[1]], [kab[2]],  quiver=([kc[1]], [kc[2]]), lw=2, lc=1, ls=:dash)
	quiver!(plt3, [0], [0],  quiver=([kabc[1]], [kabc[2]]), lc=3, lw=2)
	# quiver!([0], [0],  quiver=([ab₊[1]], [ab₊[2]]), lw=2)
	annotate!(plt1, a[1],a[2], text(L"a", :bottom))
	annotate!(plt1, b[1],b[2], text(L"b", :bottom))
	annotate!(plt1, c[1],c[2], text(L"c", :bottom))
	annotate!(plt2, ka[1],ka[2], text(L"%$(k₁)a", :bottom))
	θb = acos(dot(kb, [1,0])/ norm(kb)) * 180/π
	annotate!(plt2, (ka[1]+kab[1])/2,(ka[2]+kab[2])/2, text(L"%$(k₂)b", :top, rotation= θb))
	annotate!(plt3, (kab[1]+kabc[1])/2,(kab[2]+kabc[2])/2, text(L"%$(k₃)c", :left))
	θab = acos(dot(kab, [1,0])/ norm(kab)) * 180/π
	annotate!(plt2, kab[1]/2,kab[2]/2, text(L"%$(k₁) a+ %$(k₂)b", :red, :bottom, rotation=θab))
	annotate!(plt3, kab[1]/2,kab[2]/2, text(L"%$(k₁) a+ %$(k₂)b", :red, :bottom, rotation=θab))
	θabc = acos(dot(kabc, [1,0])/ norm(kabc)) * 180/π
	annotate!(plt3, kabc[1]/2,kabc[2]/2, text(L"%$(k₁) a+ %$(k₂)b + %$(k₃)c", :green, :bottom, rotation=θabc))
	plt1, plt2, plt3
end;

# ╔═╡ dba57307-80ce-48e7-97fd-e354da79f8fb
md"""

## Vector operation: linear combination

**Linear combination**:  combine **scaling** and **addition** together

* *scale* each vector first
* then add the scaled vectors


```math

k_1 \mathbf{a}_1 + k_2 \mathbf{a}_2 = k_1\begin{bmatrix}
           a_1\\
           a_2\\
           \vdots\\
           a_n
         \end{bmatrix} + k_2\begin{bmatrix}
           a_1\\
           a_2\\
           \vdots\\
           a_n
         \end{bmatrix}= \begin{bmatrix}
           k_1 a_1 +k_2  a_1\\
           k_1 a_2 +k_2  a_2\\
           \vdots\\
           k_1 a_n+k_2  a_n
         \end{bmatrix}

```

And the idea can be generalised to multiple vectors ``\{\mathbf{a}_1, \ldots, \mathbf{a}_m\}``, *i.e.* a set of ``m`` vectors of length ``n``


```math

k_1 \mathbf{a}_1 + k_2 \mathbf{a}_2+\ldots+ k_m \mathbf{a}_m = k_1\begin{bmatrix}
\vert \\
\mathbf{a}_{1}  \\
\vert 
\end{bmatrix} + k_2\begin{bmatrix}
\vert \\
\mathbf{a}_{2}  \\
\vert 
\end{bmatrix}+\ldots +  k_m\begin{bmatrix}
\vert \\
\mathbf{a}_{m}  \\
\vert 
\end{bmatrix} = \sum_{j=1}^m k_j \mathbf{a}_j

```
"""

# ╔═╡ 93457317-b283-4305-b426-b6d058d388f5
md"""

For example, ``
2.5\cdot \mathbf{a} + 1.5 \cdot \mathbf{b} + 1\cdot \mathbf{c}
``


```math
2.5 \begin{bmatrix} 1\\ 0\end{bmatrix} + 1.5 \begin{bmatrix} 1\\ 1\end{bmatrix} + 1 \begin{bmatrix} 0\\ 1\end{bmatrix}=\begin{bmatrix}2.5 \times 1+ 1.5 \times 1 + 1 \times 0\\ 2.5 \times 0 + 1.5\times 1 + 1\times 1\end{bmatrix}=\begin{bmatrix}4 \\ 2.5 \end{bmatrix}
```
"""

# ╔═╡ b347db05-0f86-4139-b477-f6c34e58826b
md"Geometrically, just the same \"tip-tail\" operation of the scaled vectors"

# ╔═╡ 532cd8ac-2f07-43a4-ae0a-76825d9731b7
plot(plt1, plt2,  plt3, size=(800, 500))

# ╔═╡ bdd4fa08-d6e7-42d5-89d4-2a5d0ee4f791
md"""
## Span


**Span** is the set of all possible linear combinations, *e.g.*


$$\text{span}( \{\mathbf{a}_{1}, \mathbf{a}_2, \ldots, \mathbf{a}_m\}) = \left\{\sum_{j=1}^m k_j \mathbf{a}_j ;\;k_j \in R\;  \text{for } j= 1,2, \ldots, m\right \}$$


!!! question "Question" 
	What is the span of the basis vectors in ``R^2``, *i.e.*

	$\text{span}(\{[1,0]^\top, [0,1]^\top\})?$ 


"""

# ╔═╡ db4fad25-2a6e-4cff-bcfa-9f9787bdcf11
md"""

**Answer**: ``R^2``!

Since all ``\mathbf{v} = [x,y]^\top \in R^2`` can be represented as a linear combination of the two vectors:

$$\begin{bmatrix} x\\ y\end{bmatrix} = x \begin{bmatrix} 1 \\ 0 \end{bmatrix}+ y \begin{bmatrix} 0 \\ 1 \end{bmatrix}$$

On the other hand, all linear combinations of the two vectors are in ``R^2``: 

```math
x \begin{bmatrix} 1 \\ 0 \end{bmatrix}+ y \begin{bmatrix} 0 \\ 1 \end{bmatrix} = \begin{bmatrix} x\\ y\end{bmatrix} \in R^2
```


Therefore, the span is exactly ``R^2``.


Now consider another span

!!! question "Question" 
	What is the span of the basis vectors in ``R^2``, *i.e.*

	$\text{span}(\{[2,1]^\top, [1,2]^\top\})?$ 


**Answer**: still ``R^2`` (all points in ``R^2`` can be represented by a linear combination of the two vectors; you can visually convince yourself. Remember the linear coefficients can be negative); therefore,

$$R^2 = \text{span}(\{[1,0]^\top, [0,1]^\top\}) =  \text{span}\{[2,1]^\top, [1,2]^\top\}$$

* to be more specific and formal, any ``\mathbf{v}^\top \in R^2`` can be represented as 

  $\mathbf{v} = k_1 \begin{bmatrix} 2 \\ 1 \end{bmatrix} + k_2 \begin{bmatrix} 1 \\ 2 \end{bmatrix} = \begin{bmatrix} 2 & 1 \\ 1 & 2\end{bmatrix} \begin{bmatrix} k_1 \\ k_2\end{bmatrix}$

  * and the coefficients ``k_1, k_2`` can be found by inverting the matrix (with the columns being the basis)
  $\begin{bmatrix} k_1 \\ k_2 \end{bmatrix} = \begin{bmatrix} 2 & 1 \\ 1 & 2\end{bmatrix}^{-1} \mathbf{v}$
  * and the other direction, all linear combinations of the two vectors are still in ``R^2``
* just a different set of basis (or languages) to represent the same space!

"""

# ╔═╡ 01c453a2-5e89-40f8-9a45-2f2e03b91bff
md"""

!!! question "Question"

	How about 

	$\text{span}( \{[2,1]^\top, [4,2]^\top\})?$

"""

# ╔═╡ 387eb588-7d35-45cd-8b39-33cbe611e5be
Foldable("Answer", md"""Not ``R^2``! A thin line along ``k\cdot [2,1]^\top``

$(begin
	gr()
 	plot(xlim=[-5,5], ylim=[-5, 5], ratio=1, framestyle=:origin)
	# quiver([0,0,0],[0,0,0],quiver=([1,1,1],[1,2,3]))
	oo = [0,0]
	quiver!([0], [0], quiver=([2], [1]), lw=2)
	quiver!([0], [0],  quiver=([4], [2]), lw=2)

	annotate!(2,1, text(L"v_1", :top))
	annotate!(4,2, text(L"v_2", :bottom))
	plot!(-5:1:5, (x)-> 0.5x, lw=2, ls=:dash, label=L"\texttt{span}")

end)

""")

# ╔═╡ ed1c8460-3fee-41a0-9883-372ee873aee0
md"""
!!! question "Question"
	How about 

	$\text{span}(\{[2,0,0]^\top, [2,4,0]^\top\})?$ 
"""

# ╔═╡ 0852156f-40f1-49ae-abf7-862079bb3664
md"""

## Linear independence

!!! question "Question"

	Whether adding ``\mathbf{a}_3 = \begin{bmatrix}2\\ 2 \end{bmatrix}`` to 
	```math
	\left \{\mathbf{a}_1 = \begin{bmatrix}2\\ 1 \end{bmatrix}, \mathbf{a}_2 = \begin{bmatrix}1\\ 2 \end{bmatrix}\right \}\cup \left \{\mathbf{a}_3 = \begin{bmatrix}2\\ 2 \end{bmatrix}\right \}
	``` increase the span ? *i.e.*
	
	
	```math
	
	\left \{k_1 \begin{bmatrix} 2 \\1\end{bmatrix} + k_2 \begin{bmatrix} 1 \\2\end{bmatrix}, k_1,k_2\in R \right \} \overset{?}{=} \left \{k_1 \begin{bmatrix} 2 \\1\end{bmatrix} + k_2 \begin{bmatrix} 1 \\2\end{bmatrix} + k_3 \begin{bmatrix} 2 \\2\end{bmatrix}, k_1,k_2, k_3 \in R \right \}
	```

"""

# ╔═╡ d78ec021-ae5a-4ecb-b7e8-8e9b16bb202e
Foldable("Answer", md"""

Answer: **No**

The third vector (or any one of the three) is **redundant**. Since 


```math
\mathbf{a}_3 = \frac{2}{3}\mathbf{a}_1 + \frac{2}{3}\mathbf{a}_2 =[2, 2]^\top
```

Therefore 

```math
\begin{align}
k_1 \begin{bmatrix} 2 \\1\end{bmatrix} + k_2 \begin{bmatrix} 1 \\2\end{bmatrix} + k_3 \begin{bmatrix} 2 \\2\end{bmatrix} &= k_1 \mathbf{a}_1 + k_2 \mathbf{a}_2+ k_3 \left (\frac{2}{3}\mathbf{a}_1+ \frac{2}{3}\mathbf{a}_2\right ) \\
&= (k_1+ \frac{2}{3}k_3)\mathbf{a}_1 + (k_2 + \frac{2}{3}k_3) \mathbf{a}_2\\
&= k_1' \mathbf{a}_1 + k_2' \mathbf{a}_2
\end{align}
```

That is: it is reduced to a linear combo involving the first two vectors only. We do not need the third to represent all linear combinations in the span.
""")

# ╔═╡ 58eae636-9daa-4e81-ab60-47e4c7fcace9
md"""
## Linear independence

"""

# ╔═╡ d4396c3a-d41b-41d0-aee5-0777ab6c8cd1
md"""

!!! note "Linear dependence"
	Given a set of vectors 
	```math
		\{\mathbf{a}_1, \mathbf{a}_2, \ldots, \mathbf{a}_m\}
	```

	When a vector ``\mathbf{a}_j\in \{\mathbf{a}_1, \mathbf{a}_2, \ldots, \mathbf{a}_m\}`` can be represented as a linear combination of the rest, 

	```math
		\mathbf{a}_j = \sum_{i\neq j} \lambda_i \mathbf{a}_i
	```
	then we call them **linearly dependent**.
"""

# ╔═╡ 1dce3371-3815-44ed-a35b-5cbf67e157df
md"""


**Examples** the following vectors

```math
	\left \{\mathbf{a}_1 = \begin{bmatrix}1\\ 1 \end{bmatrix}, \mathbf{a}_2 = \begin{bmatrix}2\\ 2 \end{bmatrix}\right \}
```

are linearly dependent. So is 

```math
	\left \{\mathbf{a}_1 = \begin{bmatrix}2\\ 1 \end{bmatrix}, \mathbf{a}_2 = \begin{bmatrix}1\\ 2 \end{bmatrix}, \mathbf{a}_3 = \begin{bmatrix}2\\ 2 \end{bmatrix}\right \}
```
"""

# ╔═╡ 22ac98e6-eb0e-42e0-815e-4d354235ef72
md"""
**Linear independence**: a set of vectors are **NOT linearly dependent**. *Intuitively* when none of them can be represented by the rest.

"""

# ╔═╡ 5a8dab5b-b05d-40dd-be9f-0af3b2dbfe0e
Foldable("Formal definition of linear independence", md"""

!!! note "Linear independence"
	A set of vectors 
	```math
		\{\mathbf{a}_1, \mathbf{a}_2, \ldots, \mathbf{a}_m\}
	```
	is called linearly independent if there exists no ``\lambda_1, \lambda_2, \ldots, \lambda_m`` (except the trivial solution ``\lambda_1=\lambda_2=\ldots=0`` ) such that

	```math
		\lambda_1 \mathbf{a}_1 + \lambda_2\mathbf{a}_2 +\ldots + \lambda_m \mathbf{a}_m =\mathbf{0}
	```
""")

# ╔═╡ 5d3d7bfd-3d6e-4e5c-932f-d0f5e7326737
md"""
## Vector operation: inner product 


Inner product (dot product) for real vectors 


$$\mathbf{a}^\top \mathbf{b} = [a_1, a_2\ldots, a_n] \cdot \begin{bmatrix}
            b_1 \\
            b_2 \\
           \vdots\\
            b_n
         \end{bmatrix} = \sum_{i=1}^n a_ib_i$$

**Key fact:** an inner product is a scalar rather than a vector!




"""

# ╔═╡ f5851fd7-c7bf-4f80-9436-8dba8f78f0c6
md"""



$(begin

aside(keep_working(md"I expect you to verify them or at least convince yourself!") )
end)

Some results can be proved from the definition:

  - ``\mathbf{a}^\top \mathbf{b} = \mathbf{b}^\top\mathbf{a}`` 
  - ``\mathbf{a}^\top(\mathbf{b}+\mathbf{c}) = \mathbf{a}^\top\mathbf{b}+ \mathbf{a}^\top\mathbf{c}``
  - ``(k\mathbf{a})^\top\mathbf{b} = \mathbf{a}^\top(k\mathbf{b})= k(\mathbf{a}^\top\mathbf{b})``

"""

# ╔═╡ f6fec48f-b7d9-442c-94d5-5fa938e20526
md"""
## aᵀa: inner product with itself


Inner product of a vector with itself

$$\mathbf{a}^\top \mathbf{a} = [a_1, a_2\ldots, a_n] \cdot \begin{bmatrix}
            a_1 \\
            a_2 \\
           \vdots\\
            a_n
         \end{bmatrix} = \sum_{i=1}^n a_i^2 \triangleq \|\mathbf{a}\|^2$$

- ``\mathbf{a}^\top \mathbf{a} = \sum_{i=1}^n (a_i-0)^2`` is *squared* **Euclidean distance** between $\mathbf{a}$ and $\mathbf{0}$
  * the length of vector ``\mathbf{a}``


"""

# ╔═╡ a3eeb0df-37fb-4d13-963b-e0491bc2bf7a
md"""

$(begin

aside(tip(md"``L_2`` norm measures the length of a vector from the origin."))
end)
- ``\sqrt{\mathbf{a}^\top \mathbf{a}}`` is also called **``L_2`` norm** of ``\mathbf{a}``
- ``\mathbf{a}^\top \mathbf{a} \geq 0``, sum of squares is positive
  * distance has to be positive!

- ``\mathbf{a}=\mathbf{0} \Leftrightarrow \mathbf{a}^\top\mathbf{a} =0``
  * **zero** vector has a length of zero

"""

# ╔═╡ 6a1784fe-1069-498c-bd1c-606aa3526657
md"""

## Unitify a vector


Sometimes we want to reduce a vector's norm to 1, *a.k.a.* unit directional vector; 

Given a vector ``\mathbf{a}``, the scaled vector

```math 
\frac{\mathbf{a}}{\sqrt{\mathbf{a}^\top\mathbf{a}}} 
``` 

has a unit norm (or length of 1).

*Proof:* easy to check the norm of the new vector is
```math
\frac{\mathbf{a}^\top}{\sqrt{\mathbf{a}^\top\mathbf{a}}} \frac{\mathbf{a}}{\sqrt{\mathbf{a}^\top\mathbf{a}}} = \frac{1}{(\sqrt{\mathbf{a}^\top\mathbf{a}})^2}\mathbf{a}^\top\mathbf{a}=1
```


"""

# ╔═╡ 16ec5cea-0722-4ce8-97f6-447d7a6b8acb
md"""
## Some inner products


Inner product with ``\mathbf{0}``

```math
\mathbf{0}^\top \mathbf{a} = \begin{bmatrix} 0 & \ldots & 0\end{bmatrix} \begin{bmatrix} a_1 \\ \vdots \\ a_n\end{bmatrix} = 0
```
* zero vector returns zero inner product

Inner product with ``\mathbf{1}``
```math
\mathbf{1}^\top \mathbf{a} = \begin{bmatrix} 1 & \ldots & 1\end{bmatrix} \begin{bmatrix} a_1 \\ \vdots \\ a_n\end{bmatrix} = \sum_{i=1}^n a_i 
```

* return the sum of the other vector, 
* obviously ``\mathbf{1}`` here is ``\mathbf{1}_n``, *i.e.* a length-``n`` ones vector
* therefore, ``\frac{1}{n} \mathbf{1}^\top \mathbf{a}`` is the mean 
 

!!! question "Question"
	Verify ``\frac{1}{\mathbf{1}^\top \mathbf{1}} \mathbf{1}^\top \mathbf{a}`` is also the sample mean of ``\mathbf{a}``'s entries.

"""

# ╔═╡ 9c0f7a54-690c-447e-9b20-04bb586e36ec
Foldable("Answer", md"

Note that
```math
\mathbf{1}^\top\mathbf{1} = \sum_{i=1}^n 1= \underbrace{1+1+\ldots+1}_{\text{n of 1s}} =n.
```
")

# ╔═╡ 77af9dac-fae3-4ab2-9d72-2367fcfe22e0
md"""

## Inner product relates to angle


Another interpretation: 

$$\mathbf{a}^\top \mathbf{b} = \|\mathbf{a}\| \|\mathbf{b}\|\cos \theta$$

* ``\theta``: the angle between the two vectors
* product of the two norms and the cosine of the angle between the two vectors



Measures **angular similarity** between ``\mathbf{a}`` and ``\mathbf{b}``; Remember: ``\cos\theta \in [-1, 1]``
  * maximum ``\cos(\theta)=1`` when ``\theta = 0`` (``\mathbf{a}, \mathbf{b}`` point at the same direction)
  * minimum ``\cos(\theta)=-1`` when ``\theta = \pi`` (``\mathbf{a}, \mathbf{b}`` point at opposite directions) 

"""

# ╔═╡ 41782dde-3f35-4879-aa89-1ebadd4bf8af
md"""

**Example**


``\mathbf{a} = [3,0]^\top, \mathbf{b} =[2,2]^\top``

Normal way
```math
\mathbf{a}^\top \mathbf{b} = 3\times 2 + 0 \times 2 =6
```

With the other definition

```math
\mathbf{a}^\top \mathbf{b} =  \sqrt{3^2+0^2} \sqrt{2^2+2^2} \cdot \cos 45^\circ=6
```
"""

# ╔═╡ d1f7eba1-0c09-4dcd-af6b-087699869f31
let
	gr()
 	plot(xlim=[-1,3], ylim=[-1, 3], ratio=1, framestyle=:origin)
	# quiver([0,0,0],[0,0,0],quiver=([1,1,1],[1,2,3]))
	oo = [0,0]
	a = [3,0]
	b=[2,2]
	# bp = dot(a,b)/dot(a,a)*a
	quiver!([0], [0], quiver=([3], [0]), lw=2)
	quiver!([0], [0],  quiver=([2], [2]), lw=2)
	# plot!([2,2], [2,0], ls=:dash, label="")
	annotate!(0+0.6, 0+0.3, text(L"\theta=45^\circ", :top))
	annotate!(a[1],a[2], text(L"a", :bottom))
	annotate!(b[1],b[2], text(L"b", :bottom))
	# annotate!(bp[1],bp[2]-0.1, text(L"b_{\texttt{proj}}", :top))
end

# ╔═╡ e53de410-3bad-4e07-b94f-2285c9ed8c61
md"""
## Orthogonal vectors

It can be proved that (``\mathbf{a} \perp \mathbf{b}`` means ``\mathbf{a}, \mathbf{b}`` are **orthogonal** to each other):

```math

\mathbf{a} \perp \mathbf{b} \Leftrightarrow \mathbf{a}^\top \mathbf{b}=0
```


Why ? ``\mathbf{a} \perp \mathbf{b} \Leftrightarrow  \theta=90^\circ\Leftrightarrow\cos\theta =0 \Leftrightarrow \mathbf{a}^\top\mathbf{b} =0.``

*Remarks. Maths makes so much sense, you hardly need to remember anything!*
"""

# ╔═╡ 12269268-3fde-4e4e-a954-1e9e3f27fb71
md"""

**Examples**


* the standard basis are all orthogonal to each other, *e.g.* ``[1,0] \begin{bmatrix}0\\ 1\end{bmatrix}=0``
* another example: ``[2,1] \begin{bmatrix}-1\\ 2\end{bmatrix}= -2+2=0``
* also ``\mathbf{0}`` ``\perp`` to all vectors (by definition but it also makes sense, just picture it!)

"""

# ╔═╡ 28af146b-9eb2-4490-b89e-fd9cd2965d37
md"""

## Projection by inner product


The projection of ``\mathbf{b}`` on ``\mathbf{a}`` is

```math
\mathbf{\mathbf{b}}_{\text{proj}} = \|\mathbf{b}\| \cos \theta \times \frac{\mathbf{a}}{\|\mathbf{a}\|}
```

* ``\|\mathbf{b}\| \cos\theta `` is the projected length
* ``\mathbf{a}/\|\mathbf{a}\|`` is a unit directional vector (length of one) in ``\mathbf{a}``'s direction


And the projection can also be computed with inner products only

```math

\mathbf{b}_{\text{proj}} = \frac{\mathbf{a}^\top\mathbf{b}}{\mathbf{a}^\top\mathbf{a}} \mathbf{a}
```

"""

# ╔═╡ f9011bcf-d9fd-41dd-a7b5-e86a461ef390
Foldable("Proof", md"

First

```math
\Vert \mathbf{b}\Vert \cos \theta \times \frac{\mathbf{a}}{\|\mathbf{a}\|} = \|\mathbf{a}\|\|\mathbf{b}\| \cos \theta \times \frac{\mathbf{a}}{\|\mathbf{a}\|^2} =\mathbf{a}^\top\mathbf{b} \frac{\mathbf{a}}{\mathbf{a}^\top\mathbf{a}}
```


```math
\mathbf{b}_{\text{proj}} = \frac{\mathbf{a}^\top \mathbf{b}}{\mathbf{a}^\top\mathbf{a}}\mathbf{a} = \mathbf{a}\frac{\mathbf{a}^\top \mathbf{b}}{\mathbf{a}^\top\mathbf{a}}= \frac{\mathbf{a} \mathbf{a}^\top}{\mathbf{a}^\top\mathbf{a}}\mathbf{b}
```
")

# ╔═╡ ee6f5d82-6b11-4d09-8299-2f9b3b97e809
# md"""

# Note that
# # 
# \frac{\mathbf{a}\mathbf{a}^\top}{\mathbf{a}^\top\mathbf{a}}
# `` is a ``d\times d`` **matrix**, it is called projection matrix

# * it projects ``\mathbf{b}`` to ``\mathbf{b}_{\text{proj}}``
# """

# ╔═╡ e9be3b56-87ac-4ec1-9c8b-dfd7cb12c762
proj(x::Vector{T}, a::Vector{T}) where T <: Real = dot(a,x)/dot(a,a) * a # project vector x to vector a in Julia

# ╔═╡ 3b72ff7b-2d04-44fc-8271-3509b726005d
proj([2, 2], [3, 0])

# ╔═╡ 0a108e03-08a7-4a73-be44-cb75b2943058
md"""
## Sample mean as vector inner product

We use sample mean a lot to summarise a bunch of number, *e.g.* ``\mathbf{d} = \{d_1, d_2\ldots, d_n\}``:


```math

\mu = \frac{1}{n} \sum_i d_i
```
It *compresses* a bunch of number into one scalar.


We already know the inner product is the sum

```math
\mathbf{1}^\top\mathbf{a} = \sum_{i=1}^n a_i
```


And the ratio is the sample average:

```math
\frac{\mathbf{1}^\top \mathbf{a}}{\mathbf{1}^\top\mathbf{1}} = \frac{\sum_i a_i}{n} \triangleq \bar{{a}}
```

Multiply ``\mathbf{1}`` on both side:

```math
\frac{\mathbf{1}^\top \mathbf{a}}{\mathbf{1}^\top\mathbf{1}}\mathbf{1} =\bar{{a}}\mathbf{1} =\begin{bmatrix} \bar{a} \\\bar{a} \\ \vdots\\\bar{a}\end{bmatrix}
```

* it is ``\mathbf{a}``'s projection on ``\mathbf{1}``! 

* again, the original vector is compressed to this projected vector.

"""

# ╔═╡ 82b1689d-e968-4434-9100-48c18a295080
data = [2.0,1.0]; 

# ╔═╡ 846bf22d-6ee5-4ed5-9165-7428bf348a89
mean(data)

# ╔═╡ af3a712b-5d10-4ec5-8f74-786f15169af5
proj(data, ones(length(data)))

# ╔═╡ 3f185276-9f70-499e-a1be-4d1eed1fec2d
# md"""

# ## Matrix 



# A rectangular collection of ordered numbers $\mathbf{A}\in R^{n\times m}$



# ```math
# \mathbf{A} = \begin{bmatrix}
# a_{11} & a_{12} & \ldots & a_{1m}\\
# a_{21} & a_{22} & \ldots & a_{2m} \\
# \vdots & \vdots &\ddots & \vdots \\
# a_{n1} & a_{n2} & \ldots & a_{nm}
# \end{bmatrix}
# ```
# - sometimes written as $\mathbf{A} = (a_{ij})$ $i= 1,\ldots n; j= 1,\ldots, m$
# - vectors are specific matricies
#   - row vector: ``n=1``, ``R^{1\times m}``
#   - column vector: ``m=1``, ``R^{n}``

# """

# ╔═╡ 18ed2e4c-f512-4662-b23b-317862a483bc
md"""
## Matrix 

Matrix is an extension of vectors: a matrix is a collection of ``m`` column **vectors** of length ``n``:

```math
\{\mathbf{a}_1, \mathbf{a}_2, \ldots \mathbf{a}_m\}
```


```math
\mathbf{A} = \begin{pmatrix}
a_{11} & \columncolor{lightsalmon}a_{12} & \ldots & a_{1m}\\
a_{21} & a_{22} & \ldots & a_{2m} \\
\vdots & \vdots &\ddots & \vdots \\
a_{n1} & a_{n2} & \ldots & a_{nm}
\end{pmatrix} = \begin{pmatrix}
\vert & \columncolor{lightsalmon} \vert &  & \vert\\
\mathbf{a}_1 & \mathbf{a}_2 & \ldots & \mathbf{a}_{m} \\
\vert & \vert & & \vert 
\end{pmatrix} 
```

where for ``j = 1,\ldots, m``
```math
\mathbf{a}_j = \begin{bmatrix} a_{1j}\\ a_{2j} \\ \vdots \\ a_{nj} \end{bmatrix}
```

Also written in a short hand notation as 

```math
\mathbf{A} = (a_{ij})\in R^{n\times m}\;\; \text{for}\; i= 1,\ldots n; j= 1,\ldots, m

```
* ``\mathbf{A}`` is of order (size) ``n \times m``
* ``n`` rows, ``m`` columns


**Vectors** are specific kinds of matrix
  - row vector: ``n=1``, ``R^{1\times m}``
  - column vector: ``m=1``, ``R^{n\times 1}``


"""

# ╔═╡ fe9f3c1c-1241-40da-9a9b-6c317c78cd96
md"""
## Matrix -- row view

A matrix ``\mathbf{A}\in R^{n\times m}`` can also be viewed as a collection of ``n`` **row vectors** of length ``m``:
``
\{\boldsymbol{\alpha}_1^\top, \boldsymbol{\alpha}_2^\top, \ldots \boldsymbol{\alpha}_n^\top\}
``


```math
\mathbf{A} = \begin{pmatrix}
a_{11} & a_{12} & \ldots & a_{1m}\\
\rowcolor{lightsalmon}  a_{21} & a_{22} & \ldots & a_{2m} \\
\vdots & \vdots &\ddots & \vdots \\
a_{n1} & a_{n2} & \ldots & a_{nm}
\end{pmatrix} = \begin{pmatrix}
  \rule[.5ex]{2.5ex}{0.5pt} & \boldsymbol{\alpha}_{1}^\top & \rule[.5ex]{2.5ex}{0.5pt} \\
 \rowcolor{lightsalmon}  \rule[.5ex]{2.5ex}{0.5pt} & \boldsymbol{\alpha}_{2}^\top  &  \rule[.5ex]{2.5ex}{0.5pt} \\
 & \vdots &  \\
\rule[.5ex]{2.5ex}{0.5pt} & \boldsymbol{\alpha}_{n}^\top  &  \rule[.5ex]{2.5ex}{0.5pt} \\
\end{pmatrix} 
```

where for ``i = 1,\ldots, n``
```math
\boldsymbol{\alpha}_i^\top = \begin{bmatrix} a_{i1}& a_{i2} & \ldots & a_{im} \end{bmatrix}
```
"""

# ╔═╡ eeed2564-43a8-4287-aace-6446093b07e0
md"""

## Matrix shapes


An ``n\times m`` matrix ``\mathbf{A}`` is


* **tall** if ``n > m``
* **wide** or **fat** if ``n < m``
* **square** if ``n=m``


"""

# ╔═╡ 0ee314e0-5be0-4d92-9a11-7ba10e8c920d
md"""

## Matrix addition and scaling

Matrix **addition** is in the same vain as vectors


```math
\mathbf{A} +\mathbf{B} = \begin{pmatrix}
a_{11} + b_{11}& a_{12}+b_{12} & \ldots & a_{1m}+b_{1m}\\
a_{21}+b_{21} & a_{22}+b_{22} & \ldots & a_{2m}+b_{2m} \\
\vdots & \vdots &\ddots & \vdots \\
a_{n1}+b_{n1} & a_{n2}+b_{n2} & \ldots & a_{nm}+b_{nm}
\end{pmatrix}
```

> Short-hand notation: ``(\mathbf{A}+ \mathbf{B})_{ij} = a_{ij} +b_{ij}`` for ``i = 1,\ldots, n, j = 1,\ldots, m``.

Matrix **scaling** is also the same 


```math
k\mathbf{A} = \begin{pmatrix}
k \cdot a_{11} & k\cdot a_{12} & \ldots & k\cdot a_{1m}\\
k \cdot a_{21} & k \cdot a_{22} & \ldots &k \cdot a_{2m} \\
\vdots & \vdots &\ddots & \vdots \\
k \cdot a_{n1} & k \cdot a_{n2} & \ldots & k \cdot a_{nm}
\end{pmatrix} 
```

> Short-hand notation: ``(k\mathbf{A})_{ij} = k \cdot a_{ij}`` for ``i = 1,\ldots, n, j = 1,\ldots, m``.

And the combination of the two operations:


```math
(k_1 \mathbf{A} + k_2 \mathbf{B})_{ij} = k_1\cdot a_{ij}+ k_2\cdot b_{ij}
```
for ``i = 1,\ldots, n, j = 1,\ldots, m``.

"""

# ╔═╡ e5c74165-cde4-44a4-8387-71e06894e785
md"""
## Matrix multiplication -- inner product view


Matrices of **conforming orders** can be multiplied together, *i.e.* ``\mathbf{A}\in R^{n\times m}`` and ``\mathbf{B}\in R^{m\times l}``

```math

\mathbf{A} = \begin{pmatrix}
  \rule[.5ex]{2.5ex}{0.5pt} & \boldsymbol{\alpha}_{1}^\top & \rule[.5ex]{2.5ex}{0.5pt} \\
 & \vdots &  \\
 \rowcolor{lightgreen}\rule[.5ex]{2.5ex}{0.5pt} & \boldsymbol{\alpha}_{i}^\top  &  \rule[.5ex]{2.5ex}{0.5pt} \\
 & \vdots &  \\
\rule[.5ex]{2.5ex}{0.5pt} & \boldsymbol{\alpha}_{n}^\top  &  \rule[.5ex]{2.5ex}{0.5pt} \\
\end{pmatrix} ,\;\; 

\mathbf{B} = \begin{pmatrix}
\vert &  & \columncolor{lightsalmon} \vert &  & \vert\\
\mathbf{b}_1& \ldots & \mathbf{b}_j & \ldots & \mathbf{b}_{l} \\
\vert&  & \vert & & \vert 
\end{pmatrix} 
```
```math

\mathbf{AB} =\begin{pmatrix}\boldsymbol{\alpha}_1^{\top} \mathbf{b}_1 &  \boldsymbol{\alpha}_1^{\top} \mathbf{b}_2 & \ldots & \boldsymbol{\alpha}_1^{\top} \mathbf{b}_l \\

\boldsymbol{\alpha}_2^{\top} \mathbf{b}_1 &  \boldsymbol{\alpha}_2^{\top} \mathbf{b}_2 & \ldots & \boldsymbol{\alpha}_2^{\top} \mathbf{b}_l \\
\vdots & \vdots & \ddots & \vdots\\
\boldsymbol{\alpha}_n^{\top} \mathbf{b}_1 &  \boldsymbol{\alpha}_n^{\top} \mathbf{b}_2 & \ldots & \boldsymbol{\alpha}_n^{\top} \mathbf{b}_l 
\end{pmatrix}
```


> shorthand notation: 
> ```math
> (\mathbf{AB})_{\textcolor{green}{i}\textcolor{salmon}{j}}  = \textcolor{green}{\boldsymbol{\alpha}_i}^\top \textcolor{salmon}{\mathbf{b}_j} = \sum_{r=1}^m \alpha_{ir} b_{rj}
> ``` for ``i= 1,\ldots,n``, ``j=1,\ldots, l``.
- remember we use greek letters to denote row vectors of a matrix
"""

# ╔═╡ c9c5716a-74f3-4365-98bc-daaabe47dd88
md"""

**Example**

``16`` (with ``i=2, j=2``) is the inner product between 
* ``i=2``th row of ``\mathbf{A}`` and 
* ``j=2``th column of ``\mathbf{B}``

```math
  \left(\begin{array}{cc}
    2  & 3  \\
    \rowcolor{lightgreen}6   & 4   \\
    1   & 1    \\
  \end{array}\right)

  \left(\begin{array}{ccc}
    5  & \columncolor{lightsalmon}2  & 2  \\
    1   & 1  & 2 
  \end{array}\right) = \begin{pmatrix} 13 & 7 & 10 \\
34 & \cellcolor{lightblue}{16} & 20\\
6 & 3 & 4\end{pmatrix}
```


"""

# ╔═╡ 084998bc-14f9-4e54-9704-14e02295577c
md"""

## Matrix operations 


> **Transpose**: ``(\mathbf{A}^\top)_{ij} = (\mathbf{A})_{ji}``; and ``(\mathbf{A}^\top)^\top=\mathbf{A}``



> Matrix multiplication is **not** commmutative: ``\mathbf{A}\mathbf{B} \neq \mathbf{BA}``


> Matrix multiplication is **associative**: ``\mathbf{A}\mathbf{B}\mathbf{C} = \mathbf{A}(\mathbf{BC})= (\mathbf{AB})\mathbf{C}``


> Matrix multiplication is **distributive**: 
> ``\mathbf{A}(\mathbf{B}+\mathbf{C}) = \mathbf{A}\mathbf{B}+\mathbf{AC}``,  and ``(\mathbf{A}+\mathbf{B})\mathbf{C} = \mathbf{AC}+\mathbf{BC}``

> Matrix multi+transpose: ``(\mathbf{AB})^\top = \mathbf{B}^\top\mathbf{A}^\top``


All the above can be recursively applied together, *e.g.* 

> ``(\mathbf{ABC})^\top = \mathbf{C}^\top \mathbf{B}^\top \mathbf{A}^\top``

*Proof:*
``(\mathbf{ABC})^\top = ((\mathbf{AB})\mathbf{C})^\top =\mathbf{C}^\top (\mathbf{AB})^\top= \mathbf{C}^\top \mathbf{B}^\top \mathbf{A}^\top``; 
*You will get the same result if associate the other way around:* *i.e.* ``(\mathbf{ABC})^\top = (\mathbf{A}(\mathbf{B}\mathbf{C}))^\top``*, which is left as an exercise.*
"""

# ╔═╡ a571a9c3-7eb2-4fcb-b5db-052fb01b811a
md"""

## Special matrices


> **Identity matrix**: a square matrix with ones on its diagonal
> ```math
> \mathbf{I}_{n} = \begin{pmatrix} 1 & 0 & \ldots & 0 \\
> 0 & 1 & \ldots & 0 \\
> \vdots & \vdots & \ddots & \vdots\\
> 0 & 0 & \ldots & 1 \\
> \end{pmatrix}
> ```
> * ``\mathbf{I}_{m}\mathbf{A} = \mathbf{AI}_{n} = \mathbf{A}`` (you should verify yourself!)
> * serve the same role of ``1``, note that ``1\times x = x``

> **Diagonal matrix**: a square matrix with non-zero entries on the diagonal entries,
> *e.g.* 
> ```math
>  \text{diag}(\{1,2,3\}) = \begin{pmatrix} 1 & 0 & 0\\
> 0 & 2 & 0\\
> 0 & 0 & 3 \end{pmatrix}
> ```


> **Symmetric matrix**: a matrix such that ``\mathbf{A}^\top =\mathbf{A}``


**Example**

All diagonal matrices are symmetric, such as ``\mathbf{I}``


Gram matrix ``\mathbf{X}^\top\mathbf{X}`` is always symmetric (again you should verify it!)
"""

# ╔═╡ e9dc184d-253f-45d7-ac44-0a1cd02ee746
md"""

## Matrix inversion

The generalisation of the division operation between scalars 
* but a lot more restricted.

What is the number ``a^{-1}``?
> for all ``a\neq 0``, ``a^{-1}\in R`` is defined as a scalar *s.t.*
> ```math 
> a a^{-1} =1\;\; \text{or}\;\; a^{-1}a =1
> ```

Now matrix inversion:

> Given square ``\mathbf{A} \in R^{n\times n}``, its **inversion matrix**, if exists, ``\mathbf{A}^{-1}`` is defined as a ``n\times n`` matrix such that 
> ```math
> \mathbf{A}\mathbf{A}^{-1} =\mathbf{I}_{n}
> ```
> or equivalently ``\mathbf{A}^{-1}\mathbf{A}=\mathbf{I}_{n}``.
> * not all number has their inverse, such as ``0``, so are matrices;
> *  only some square matrices have inverses
"""

# ╔═╡ f7183a91-8cf2-4b99-8cd5-dd250a9d07a1
md"""
## More matrix identities


You should try to verify them.

> ``(\mathbf{A}^\top)^{-1} = (\mathbf{A}^{-1})^\top \triangleq \mathbf{A}^{-\top}``

> ``(\mathbf{A}^{-1})^{-1} = \mathbf{A}``

> ``(\mathbf{A}\mathbf{B})^{-1} = \mathbf{B}^{-1}\mathbf{A}^{-1}`` if the individual inverses exist

> ``\text{diag}(\{x_1, x_2, \ldots, x_n\})^{-1} = \text{diag}(\{x_1^{-1}, x_2^{-1}, \ldots, x_n^{-1}\})`` 
> if ``x_1, x_2, \ldots, x_n \neq 0``
**Example**

```math
\begin{pmatrix}
1 & 0 &0\\
0 & 2 & 0\\
0 & 0 & 3
\end{pmatrix}^{-1} = \begin{pmatrix}
1/1 & 0 &0\\
0 & 1/2 & 0\\
0 & 0 & 1/3
\end{pmatrix}

```

"""

# ╔═╡ 9b898455-35b4-40db-b585-994f06ebc496
md"""

## Matrix operations in `Julia`\*
"""

# ╔═╡ deabaa69-aae4-4630-8e34-804979c391c2
md"""

**Matrix transpose** ``\top`` with '

```julia
A' # Matrix A's transpose
```
"""

# ╔═╡ eb0ba4b6-4df9-487f-9451-30872181c523
let
	X = rand(3,2)
	X, X'
end

# ╔═╡ 2bb2f9f0-9e97-4ff4-97a8-acd62afe1d8d
md"Gram matrix is **symmetric**"

# ╔═╡ 146f5e6d-9f7d-439f-8bd8-5f18ba206e80
let
	# randomly generate a matrix
	X = rand(5,3)
	# the gram matrix is symmetric (which should be obvious based on the multiplication and transpose rule)
	X' * X, issymmetric(X'*X)
end

# ╔═╡ 255e1650-8b9f-47f3-b92e-93e10df238bd
md"**Diagonal matrix**"

# ╔═╡ 62570da5-f5db-45f8-bbcf-584ed486b9fe
let
	# create a diagonal matrix
	diagm([1,2,3])
end

# ╔═╡ a46b4c33-0d8d-4a85-a61b-5fe11edc76e0
md"**Matrix inversion**"

# ╔═╡ 0d911b10-a700-4ddb-88bf-8a3201995f3c
let
	A = diagm([1,2,3])
	A^(-1), inv(A), A * A^(-1)
end

# ╔═╡ f3a56c3e-52d2-4e33-9220-c6b6ff413973
md"""

## Why some A not invertible?

I have claimed some matrices, even square, *i.e.* ``\mathbf{A}\in R^{n\times n}`` are not invertible (*a.k.a.* **singular**)
* but why ?


To see it, let's view matrix operation as a *linear* transformation (or just transformation)

```math
\phi(\mathbf{v})\triangleq\mathbf{Av}

```
* a function ``\phi: R^n \rightarrow R^n``
* input vector ``\mathbf{v}\in R^{n}`` 
* it maps ``\mathbf{v}`` to an output vector ``\phi(\mathbf{v}) = \mathbf{Av}``

And we consider whether such a mapping can be reversed: ``\phi^{-1}: R^n \rightarrow R^n``
* some ``\mathbf{A}``'s transformation can be reversed
* some cannot!


Formally, ``\mathbf{A}`` is only invertible ``\Leftrightarrow`` the mapping ``\phi`` is one to one (or **bijective**)

**Example** consider matrix ``\mathbf{A}=\mathbf{I}``, ``\phi`` maps a vector to itself:

```math
\phi(\mathbf{v}) =\mathbf{v};
```
it is one-to-one mapping; therefore, such an operation can be reversed easily

```math
\phi^{-1}(\mathbf{v}) =\mathbf{v}
```

Therefore, the identity matrix ``\mathbf{I}`` is invertible.
"""

# ╔═╡ 8d3f8d08-353f-4896-83af-cf37c0bd2f06
begin
	Rmat(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]
end

# ╔═╡ cac10de9-50f4-4dd9-9238-338c0e8dc37b
md"""

## Invertible example -- rotation


A rotation matrix in ``R^2`` 


```math
\mathbf{R} =\begin{bmatrix} \cos\theta & -\sin\theta \\
  \sin\theta & \cos\theta\end{bmatrix} 

```
- it *rotates* $\mathbf{v}$ anti-clockwise by $\theta$

  $$\mathbf{R}\mathbf{v} = \begin{bmatrix} \cos\theta & -\sin\theta \\
  \sin\theta & \cos\theta\end{bmatrix} \begin{bmatrix} v_1\\ v_2\end{bmatrix} =  \begin{bmatrix} \cos\theta v_1 - \sin\theta v_2\\ \sin\theta v_1 + \cos\theta v_2\end{bmatrix}$$ *e.g.* $\theta = \pi/4$

  ``\mathbf{R}=``$(begin latexify_md(Rmat(π/4)) end)
"""

# ╔═╡ f7dabe72-edb7-4188-91c7-a3285e82d408
md"You can also rotate more than once: ``\mathbf{R}\mathbf{R}\ldots\mathbf{R}\mathbf{v}``"

# ╔═╡ d5198c78-e073-448f-8e54-c5b854daa773
let
	gr()
 	plot(xlim=[-3,3], ylim=[-.5, 3], ratio=1, framestyle=:origin)
	# quiver([0,0,0],[0,0,0],quiver=([1,1,1],[1,2,3]))
	a = [2,1]
	R = Rmat(π/4)
	Ra= R * a
	RRa= R * Ra
	# bp = dot(a,b)/dot(a,a)*a
	quiver!([0], [0], quiver=([a[1]], [a[2]]), lw=2)
	quiver!([0], [0],  quiver=([Ra[1]], [Ra[2]]), lw=2)
	quiver!([0], [0],  quiver=([RRa[1]], [RRa[2]]), lw=2)
	# plot!([b[1], bp[1]], [b[2],bp[2]], ls=:dash, lc=:gray, lw=2, label="")
	# annotate!(0+0.3, 0+0.3, text(L"\theta", :top))
	annotate!(a[1],a[2], text(L"\mathbf{v}", :bottom))
	annotate!(Ra[1],Ra[2], text(L"\mathbf{Rv}", :bottom))
	annotate!(RRa[1],RRa[2], text(L"\mathbf{RRv}", :bottom))
	# plot!(Shape([0, aunit[1], abunit[1], bunit[1]], [0, aunit[2], abunit[2], bunit[2]]), lw=1, fillcolor=false)
	# plot!(perp_square([0,0], a, b; δ=0.1), lw=1, fillcolor=false, label="")
	# annotate!(bp[1],bp[2]-0.1, text(L"b_{\texttt{proj}}", :top))

end

# ╔═╡ 32fc8230-82fd-461c-b219-65c3ba717dd1
md"""
## Is rotation matrix R invertible ?

**Inversion of ``\mathbf{R}``:**  *intuitively,* rotation operations can be reversed/cancelled
* the inversion ``\mathbf{R}^{-1}\mathbf{w}`` just *rotate back*
* no information is lost in the rotation operation


Based on trigonometry, ``\mathbf{R}^{-1}`` can be found as 

```math
\mathbf{R}^{-1} = \begin{bmatrix} \cos\theta & \sin\theta \\
  -\sin\theta & \cos\theta\end{bmatrix} 

```
* you should verify ``\mathbf{R}\mathbf{R}^{-1} =\mathbf{I}``

"""

# ╔═╡ 03e45498-c6b5-4d3d-934a-60cc0a4ac620
let
	gr()
 	plot(xlim=[-3,3], ylim=[-.5, 3], ratio=1, framestyle=:origin)
	# quiver([0,0,0],[0,0,0],quiver=([1,1,1],[1,2,3]))
	a = [2,1]
	R = Rmat(π/4)
	Ra, RRa= R * a, R * R*a
	# Rinva, RRinva = Rinv*Ra, 
	# bp = dot(a,b)/dot(a,a)*a
	quiver!([0], [0], quiver=([a[1]], [a[2]]), lw=2)
	quiver!([0], [0],  quiver=([Ra[1]], [Ra[2]]), lw=2)
	quiver!([0], [0],  quiver=([RRa[1]], [RRa[2]]), lw=2)
	# plot!([b[1], bp[1]], [b[2],bp[2]], ls=:dash, lc=:gray, lw=2, label="")
	# annotate!(0+0.3, 0+0.3, text(L"\theta", :top))
	annotate!(a[1],a[2], text(L"\mathbf{R}^{-1}\mathbf{R}^{-1}\mathbf{w}", :bottom))
	annotate!(Ra[1],Ra[2], text(L"\mathbf{R}^{-1}\mathbf{w}", :bottom))
	annotate!(RRa[1],RRa[2], text(L"\mathbf{w}", :bottom))
	# plot!(Shape([0, aunit[1], abunit[1], bunit[1]], [0, aunit[2], abunit[2], bunit[2]]), lw=1, fillcolor=false)
	# plot!(perp_square([0,0], a, b; δ=0.1), lw=1, fillcolor=false, label="")
	# annotate!(bp[1],bp[2]-0.1, text(L"b_{\texttt{proj}}", :top))

end

# ╔═╡ 27ee0554-8df8-4449-bcc2-352dc4556bd8
md"""

## Non-invertible example -- projection


Projection: vector ``\mathbf{b}`` to vector ``\mathbf{a}`` is

```math

\mathbf{b}_{\text{proj}} = \frac{\mathbf{a}^\top\mathbf{b}}{\mathbf{a}^\top\mathbf{a}} \mathbf{a}
```

which can be rearranged as a matrix operation


```math
\mathbf{b}_{\text{proj}} = \frac{\mathbf{aa}^\top}{\mathbf{a}^\top\mathbf{a}} \mathbf{b} = \mathbf{P}_{a} \mathbf{b}
```

* note that ``\mathbf{P}_a = \frac{\mathbf{aa}^\top}{\mathbf{a}^\top\mathbf{a}}`` is a ``n\times n`` matrix
* it projects ``\mathbf{b}`` to ``\mathbf{a}``

Projection operations **NOT invertible**, therefore ``\mathbf{P}_a^{-1}`` **does not** exist

* projection is a many-to-one operation
  * information is lost in the projection operation
* in other words, projection, once done, cannot be reversed!

"""

# ╔═╡ d445d308-66ec-468c-836e-d2a15127ba2b
md"**Example** when you use sample average to summarise datasets with the same size
* information is lost in this operation
* as long as the sum is the same, it returns the same vector
* geometrically, the projection is the same
"

# ╔═╡ 78dd7f6c-568f-40bf-96af-0c74bde4e98d
md"""

## Rank


Recall a matrix can be viwed as  a collection of vectors ``\{\mathbf{a}_1, \mathbf{a}_2, \ldots, \mathbf{a}_m\}``:



```math
\mathbf{A} =  \begin{pmatrix}
\vert & \columncolor{lightsalmon} \vert &  & \vert\\
\mathbf{a}_1 & \mathbf{a}_2 & \ldots & \mathbf{a}_{m} \\
\vert & \vert & & \vert 
\end{pmatrix} 
```


!!! info "Rank"
	**Rank** of matrix ``\mathbf{A}`` can be defined as the maximum number of linearly independent column vectors among the set ``\{\mathbf{a}_1, \ldots, \mathbf{a}_m\}.``



Rank and invertibility:

> ``\mathbf{A}`` is invertible ``\Leftrightarrow`` ``\texttt{rank}(\mathbf{A})=m``, *i.e.*  ``\mathbf{A}`` is full rank
"""

# ╔═╡ 9bd6dc7f-f449-4ef1-bd16-ca2d47f983eb
md"""

**Example** consider vector set ``
	\left \{\mathbf{a}_1 = \begin{bmatrix}1\\ 1 \end{bmatrix}, \mathbf{a}_2 = \begin{bmatrix}2\\ 2 \end{bmatrix}\right \}
``, the associated column matrix's rank is 1.

```math
\texttt{rank}\left (\begin{bmatrix} 1 & 2 \\ 1 & 2\end{bmatrix}\right ) = 1
```

Given another set

```math
	\left \{\mathbf{a}_1 = \begin{bmatrix}2\\ 1 \end{bmatrix}, \mathbf{a}_2 = \begin{bmatrix}1\\ 2 \end{bmatrix}, \mathbf{a}_3 = \begin{bmatrix}2\\ 2 \end{bmatrix}\right \},
```its matrix's rank is 2:


```math
\texttt{rank}\left (\begin{bmatrix} 2 & 1 & 2\\ 1 & 2 & 2\end{bmatrix}\right ) = 2.
```
"""

# ╔═╡ fcdf7707-c675-4104-9a4c-03b12ffd1024
md"**Demonstration:** we can use `rank()` method to find the rank of matrix in Julia. 

"

# ╔═╡ 87e5238d-01c6-422a-b847-16fe1dc9971f
rank([1 2; 1 2])

# ╔═╡ cb2dab5c-a719-47e3-8df2-b6a191f8dccb
rank([2 1 2; 1 2 2])

# ╔═╡ 62658932-19b1-4103-a0bc-653825f8bfe3
md"For example, rotation matrix is always full rank:"

# ╔═╡ e40ede69-afa7-4d0a-b2b4-67c91d1f5faa
rank(Rmat(π/4))

# ╔═╡ 67f12dea-67db-49a6-816f-678f88577369
md"A projection matrix to ``\mathbf{1}`` is rank deficit: ``\texttt{rank}(\mathbf{P}_{\mathbf{1}}) = 1 < 2``; therefore, not invertible or singular"

# ╔═╡ 1a2cac05-52ab-4483-9b97-d4453173f626
let
	ones₂ = ones(2)
	P = ones₂*ones₂'/dot(ones₂, ones₂) 
	rank(P)
end

# ╔═╡ 77b40da6-4577-4968-a48d-f4ba7c6d1bca
md"""

# Appendix
"""

# ╔═╡ be040c96-da49-44e6-9a73-7e26a1960261
# as: arrow head size 0-1 (fraction of arrow length)
# la: arrow alpha transparency 0-1
function arrow3d!(x, y, z,  u, v, w; as=0.1, lc=:black, la=1, lw=0.4, scale=:identity)
    (as < 0) && (nv0 = -maximum(norm.(eachrow([u v w]))))
    for (x,y,z, u,v,w) in zip(x,y,z, u,v,w)
        nv = sqrt(u^2 + v^2 + w^2)
        v1, v2 = -[u,v,w]/nv, nullspace(adjoint([u,v,w]))[:,1]
        v4 = (3*v1 + v2)/3.1623  # sqrt(10) to get unit vector
        v5 = v4 - 2*(v4'*v2)*v2
        (as < 0) && (nv = nv0) 
        v4, v5 = -as*nv*v4, -as*nv*v5
        plot!([x,x+u], [y,y+v], [z,z+w], lc=lc, la=la, lw=lw, scale=scale, label=false)
        plot!([x+u,x+u-v5[1]], [y+v,y+v-v5[2]], [z+w,z+w-v5[3]], lc=lc, la=la, lw=lw, label=false)
        plot!([x+u,x+u-v4[1]], [y+v,y+v-v4[2]], [z+w,z+w-v4[3]], lc=lc, la=la, lw=lw, label=false)
    end
end

# ╔═╡ 85843ba6-f480-4c0e-a80b-2d5742c0cd72
let
	gr()
 	plt=plot(xlim=[-1,5], ylim=[-1, 5], zlim =[-2,5], framestyle=:zerolines, camera=(20,15), size=(800,600))
	# quiver([0,0,0],[0,0,0],quiver=([1,1,1],[1,2,3]))
	# quiver!([0], [0], [0], quiver=([c[1]], [c[2]], [c[3]]), lw=2)
	# quiver!([0], [0],  quiver=([b[1]], [b[2]]), lw=2)
	arrow3d!([0], [0], [0], [cv[1]], [cv[2]], [cv[3]]; as=0.1, lc=2, la=1, lw=2, scale=:identity)
	# annotate!(c[1], c[2], c[3], text("c"))
	arrow3d!([0], [0], [0], [dv[1]], [dv[2]], [dv[3]]; as=0.1, lc=3, la=1, lw=2, scale=:identity)
	scatter!([0], [0], [0], mc=:black, label="")

	arrow3d!([0], [0], [0], [fv[1]], [fv[2]], [fv[3]]; as=0.1, lc=1, la=1, lw=2, scale=:identity)

	surface!(-1:0.1:5, -1:0.1:5, (x,y) -> 0, colorbar=false, alpha=0.1)
	plt
	# annotate!(a[1],a[2], text(L"a", :top))
	# annotate!(b[1],b[2], text(L"b", :top))
end

# ╔═╡ d91f4692-fce4-4cbc-9f41-2d25cfab93cb
Foldable("Answer", md"""

Not ``R^3`` but only the ``xy``-plane!

$(let
	c, d, e= fv, cv, dv
 	plt=plot(xlim=[-1,5], ylim=[-1, 5], zlim =[-2,5], framestyle=:zerolines, camera=(20,15), size=(700,500))
	arrow3d!([0], [0], [0], [c[1]], [c[2]], [c[3]]; as=0.1, lc=1, la=1, lw=2, scale=:identity)
	# annotate!(c[1], c[2], c[3], text("c"))
	arrow3d!([0], [0], [0], [d[1]], [d[2]], [d[3]]; as=0.1, lc=2, la=1, lw=2, scale=:identity)
	scatter!([0], [0], [0], mc=:black, label="")

	arrow3d!([0], [0], [0], [e[1]], [e[2]], [e[3]]; as=0.1, lc=3, la=1, lw=2, scale=:identity)

	surface!(-1:0.1:5, -1:0.1:5, (x,y) -> 0, colorbar=false, alpha=0.1)
	plt
end)

""")

# ╔═╡ 6863be47-bced-4a2b-8578-ff4fc7bc1070
function perp_square(origin, vx, vy; δ=0.1) 
	x = δ * vx/sqrt(norm(vx))
	y = δ * vy/sqrt(norm(vy))
	xyunit = origin+ x + y
	xunit = origin + x
	yunit = origin +y
	Shape([origin[1], xunit[1], xyunit[1], yunit[1]], [origin[2], xunit[2], xyunit[2], yunit[2]])
end

# ╔═╡ faf91165-0026-48ce-8228-7877426401ab
let
	gr()
 	plot(xlim=[-1.5,3], ylim=[-1, 2.5], ratio=1, framestyle=:origin)
	# quiver([0,0,0],[0,0,0],quiver=([1,1,1],[1,2,3]))
	oo = [0,0]
	a = [2,1]
	b= [-1,2]
	# bp = dot(a,b)/dot(a,a)*a
	quiver!([0], [0], quiver=([a[1]], [a[2]]), lw=2)
	quiver!([0], [0],  quiver=([b[1]], [b[2]]), lw=2)
	# plot!([b[1], bp[1]], [b[2],bp[2]], ls=:dash, lc=:gray, lw=2, label="")
	# annotate!(0+0.3, 0+0.3, text(L"\theta", :top))
	annotate!(a[1],a[2], text(L"a=[2,1]^\top", :bottom))
	annotate!(b[1],b[2], text(L"b=[-1,2]^\top", :bottom))

	# plot!(Shape([0, aunit[1], abunit[1], bunit[1]], [0, aunit[2], abunit[2], bunit[2]]), lw=1, fillcolor=false)
	plot!(perp_square([0,0], a, b; δ=0.1), lw=1, fillcolor=false, label="")
	# annotate!(bp[1],bp[2]-0.1, text(L"b_{\texttt{proj}}", :top))
end

# ╔═╡ 79ff49e9-4c0d-4e9e-bd6b-080abdf86020
let
	gr()
 	plot(xlim=[-1,3], ylim=[-1, 3], ratio=1, framestyle=:origin)
	# quiver([0,0,0],[0,0,0],quiver=([1,1,1],[1,2,3]))
	oo = [0,0]
	a = [3,0]
	b= [2,2]
	bp = dot(a,b)/dot(a,a)*a
	quiver!([0], [0], quiver=([3], [0]), lw=2)
	quiver!([0], [0],  quiver=([2], [2]), lw=2)
	plot!([b[1], bp[1]], [b[2],bp[2]], ls=:dash, lc=:gray, lw=2, label="")
	annotate!(0+0.3, 0+0.3, text(L"\theta", :top))
	annotate!(a[1],a[2], text(L"a", :bottom))
	annotate!(b[1],b[2], text(L"b", :bottom))
	annotate!(bp[1],bp[2]-0.1, text(L"b_{\texttt{proj}}", :top))
	plot!(perp_square(bp, [1,0], [0,1]; δ=0.15), lw=1, label="", fillcolor=false)
end

# ╔═╡ 4110293d-7f40-4d9c-aefa-268042c8db07
let

	gr()
 	plot( ratio=1, framestyle=:origin)
	# quiver([0,0,0],[0,0,0],quiver=([1,1,1],[1,2,3]))
	oo = [0,0]
	a = [1,1]
	b = data
	bp = dot(a,b)/dot(a,a)*a
	quiver!([0], [0], quiver=([a[1]], [a[2]]), lw=2)
	quiver!([0], [0],  quiver=([b[1]], [b[2]]), lw=2)
	plot!([b[1], bp[1]], [b[2],bp[2]], ls=:dash, lc=:gray, lw=2, label="")

	quiver!([0], [0],  quiver=([bp[1]], [bp[2]]), ls=:dash, lw=2)
	annotate!(a[1],a[2], text(L"\mathbf{1}", :top))
	annotate!(b[1],b[2], text(L"\mathbf{d}", :bottom))
	# annotate!(bp[1]+0.2,bp[2], text(L"b_{\texttt{proj}} =latexify(:(x = $t))", :left))
	annotate!(bp[1]+0.2,bp[2], text(L"\mathbf{d}_{\texttt{proj}} =%$(bp)", 10,:left))
	plot!(perp_square(bp, a, b-bp; δ=0.1), lw=1, label="", fillcolor=false)
end

# ╔═╡ 7c00704e-1272-4b15-a84b-8110627274ad
let

	gr()
 	plot(xlim =[-2.5, 4.9], ylim=[-1.5, 5],  ratio=1, framestyle=:origin)
	# quiver([0,0,0],[0,0,0],quiver=([1,1,1],[1,2,3]))
	oo = [0,0]
	a = [1,1]
	b = data
	bp = dot(a,b)/dot(a,a)*a
	quiver!([0], [0], quiver=([a[1]], [a[2]]), lw=2)
	# quiver!([0], [0],  quiver=([b[1]], [b[2]]), lw=2)
	# plot!([b[1], bp[1]], [b[2],bp[2]], ls=:dash, lc=:gray, lw=2, label="")
	plot!(-2:0.5:4.5, x -> -x +3, ls=:dash, lc=:gray, lw=2, label="")
	quiver!([0], [0],  quiver=([bp[1]], [bp[2]]), ls=:dash, lw=2)
	annotate!(a[1],a[2], text(L"\mathbf{1}", :top))
	# annotate!(b[1],b[2], text(L"\mathbf{d}_1", :top))

	data2 = [[4,-1], [3,0], [2,1], [1,2], [0,3], [-1, 4]]

	for i in 1:length(data2)
		quiver!([0], [0],  quiver=([data2[i][1]], [data2[i][2]]), lw=1)
		annotate!(data2[i][1], data2[i][2], text(L"\mathbf{d}", :left))
	end
	# data3 = [0,3]
	# data4 = [-1,4]
	# quiver!([0], [0],  quiver=([data2[1]], [data2[2]]), lw=2)
	# annotate!(data2[1],data2[2], text(L"\mathbf{d}_2", :left))
	# quiver!([0], [0],  quiver=([data3[1]], [data3[2]]), lw=2)
	# annotate!(data3[1],data3[2], text(L"\mathbf{d}_3", :left))
	# quiver!([0], [0],  quiver=([data4[1]], [data4[2]]), lw=2)
	# annotate!(data4[1],data4[2], text(L"\mathbf{d}_4", :left))
	# annotate!(bp[1]+0.2,bp[2], text(L"b_{\texttt{proj}} =latexify(:(x = $t))", :left))
	annotate!(bp[1]+0.2,bp[2], text(L"\mathbf{P}_{1}\mathbf{d} =%$(bp)", 10,:left))
	plot!(perp_square(bp, a, b-bp; δ=0.1), lw=1, label="", fillcolor=false)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
LaTeXStrings = "~1.3.0"
Latexify = "~0.15.16"
Plots = "~1.31.7"
PlutoTeachingTools = "~0.1.7"
PlutoUI = "~0.7.39"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "d83bcd3c1142c9f8e9bfe5b6e2122c62c0f72af5"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "195c5505521008abea5aee4f96930717958eac6f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.4.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "80ca332f6dcb2508adba68f22f551adb2d00a624"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.3"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "1833bda4a027f4b2a1c984baddcf755d77266818"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.1.0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "1fd869cc3875b57347f7027521f561cf46d1fcd8"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.19.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "78bee250c6826e1cf805a88b7f1e86025275d208"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.46.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "5158c2b41018c5f7eb1470d558127ac274eca0c9"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.1"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "cf0a9940f250dc3cb6cc6c6821b4bf8a4286cf9c"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.66.2"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "2d908286d120c584abbe7621756c341707096ba4"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.66.2+0"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "fb28b5dc239d0174d7297310ef7b84a11804dfab"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.0.1"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "a7a97895780dab1085a97769316aa348830dc991"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.3"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "59ba44e0aa49b87a8c7a8920ec76f8afe87ed502"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.3.3"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "b3364212fb5d870f724876ffcd34dd8ec6d98918"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.7"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "0f960b1404abb0b244c1ece579a0ec78d056a5d1"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.15"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "1a43be956d433b5d0321197150c2f94e16c0aaa0"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.16"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "94d9c52ca447e23eac0c0f074effbcd38830deb5"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.18"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "5d4d2d9904227b8bd66386c1138cf4d5ffa826bf"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "0.4.9"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "dedbebe234e06e1ddad435f5c6f4b85cd8ce55f7"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.2.2"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "ae6676d5f576ccd21b6789c2cbe2ba24fcc8075d"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.5"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e60321e3f2616584ff98f0a4f18d98ae6f89bbb3"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.17+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "3d5bf43e3e8b412656404ed9466f1dcbf7c50269"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.4.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "8162b2f8547bc23876edd0c5181b27702ae58dce"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.0.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "9888e59493658e476d3073f1ce24348bdc086660"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.0"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "a19652399f43938413340b2068e11e55caa46b65"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.31.7"

[[deps.PlutoHooks]]
deps = ["InteractiveUtils", "Markdown", "UUIDs"]
git-tree-sha1 = "072cdf20c9b0507fdd977d7d246d90030609674b"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0774"
version = "0.0.5"

[[deps.PlutoLinks]]
deps = ["FileWatching", "InteractiveUtils", "Markdown", "PlutoHooks", "Revise", "UUIDs"]
git-tree-sha1 = "0e8bcc235ec8367a8e9648d48325ff00e4b0a545"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0420"
version = "0.1.5"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "LaTeXStrings", "Latexify", "Markdown", "PlutoLinks", "PlutoUI", "Random"]
git-tree-sha1 = "67c917d383c783aeadd25babad6625b834294b30"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.1.7"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8d1f54886b9037091edf146b517989fc4a09efec"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.39"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "e7eac76a958f8664f2718508435d058168c7953d"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.3"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "22c5201127d7b243b9ee1de3b43c408879dff60f"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.3.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "dad726963ecea2d8a81e26286f625aee09a91b7c"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.4.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "dfec37b90740e3b9aa5dc2613892a3fc155c3b42"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.6"

[[deps.StaticArraysCore]]
git-tree-sha1 = "ec2bd695e905a3c755b33026954b119ea17f2d22"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.3.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArraysCore", "Tables"]
git-tree-sha1 = "8c6ac65ec9ab781af05b08ff305ddc727c25f680"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.12"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "8a75929dcd3c38611db2f8d08546decb514fcadf"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.9"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIs]]
git-tree-sha1 = "e59ecc5a41b000fa94423a578d29290c7266fc10"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "58443b63fb7e465a8a7210828c91c08b92132dff"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.14+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╟─9edaaf3c-2973-11ed-193f-d500cbe5d239
# ╟─38967916-87c6-4bb8-8ca6-f843b482e2ca
# ╟─cb9259ce-beab-4e94-88a9-66057a88b633
# ╟─7f825237-767e-4b3e-927f-e3b0d6779316
# ╟─457b21c1-6ea9-4983-b724-fb5cbb69739d
# ╟─57945aa2-b984-407a-89b1-57b0b0a2cd75
# ╟─b18cf2b6-0e45-4158-9742-2d4e48ec3765
# ╟─9ca74f71-e1f6-43c7-b08a-7afb7d3dfe0d
# ╠═bd2498df-069b-4ea0-9e44-738142a3080e
# ╟─5d4e270a-b480-40d1-97ac-d7eaa4700766
# ╟─6d4861b6-4a0d-4aac-8f08-e6861e2ecf70
# ╟─3ea11e49-9edd-43ac-9908-01d766c331a1
# ╠═4c1e7871-a931-4781-a200-304a2ef253e1
# ╟─85843ba6-f480-4c0e-a80b-2d5742c0cd72
# ╟─fd80776b-1bc2-466a-b29c-28df81cbee3f
# ╟─c56ec750-0514-4796-a594-325b930bf4d2
# ╟─1dd53149-6749-4c86-b63a-eb801a827808
# ╟─057908d5-7f6b-4976-84f4-eabe8ff88c33
# ╟─8539909b-0f49-4aa3-a950-89bf1f09e696
# ╟─b5199746-42a4-4462-81ef-d329482f5a3c
# ╟─94027fb4-cbb0-46d3-abd2-795af9019cd7
# ╟─46a8bf2d-5fea-46c3-ba65-4cb02e2279f2
# ╟─a54dbf58-d082-440f-bc3d-ceebdfabbda6
# ╟─a37be5c2-2c38-4625-8ca4-3c37ffbc3bac
# ╟─8c5201e7-a0a0-497b-9459-4e518add61de
# ╟─d4d8448c-1148-4c56-a460-acc6279013ba
# ╟─fc9fb8cd-6557-4eff-abab-30be01366c91
# ╟─dba57307-80ce-48e7-97fd-e354da79f8fb
# ╟─93457317-b283-4305-b426-b6d058d388f5
# ╟─b347db05-0f86-4139-b477-f6c34e58826b
# ╟─532cd8ac-2f07-43a4-ae0a-76825d9731b7
# ╟─bdd4fa08-d6e7-42d5-89d4-2a5d0ee4f791
# ╟─db4fad25-2a6e-4cff-bcfa-9f9787bdcf11
# ╟─01c453a2-5e89-40f8-9a45-2f2e03b91bff
# ╟─387eb588-7d35-45cd-8b39-33cbe611e5be
# ╟─ed1c8460-3fee-41a0-9883-372ee873aee0
# ╟─d91f4692-fce4-4cbc-9f41-2d25cfab93cb
# ╟─0852156f-40f1-49ae-abf7-862079bb3664
# ╟─d78ec021-ae5a-4ecb-b7e8-8e9b16bb202e
# ╟─58eae636-9daa-4e81-ab60-47e4c7fcace9
# ╟─d4396c3a-d41b-41d0-aee5-0777ab6c8cd1
# ╟─1dce3371-3815-44ed-a35b-5cbf67e157df
# ╟─22ac98e6-eb0e-42e0-815e-4d354235ef72
# ╟─5a8dab5b-b05d-40dd-be9f-0af3b2dbfe0e
# ╟─5d3d7bfd-3d6e-4e5c-932f-d0f5e7326737
# ╟─f5851fd7-c7bf-4f80-9436-8dba8f78f0c6
# ╟─f6fec48f-b7d9-442c-94d5-5fa938e20526
# ╟─a3eeb0df-37fb-4d13-963b-e0491bc2bf7a
# ╟─6a1784fe-1069-498c-bd1c-606aa3526657
# ╟─16ec5cea-0722-4ce8-97f6-447d7a6b8acb
# ╟─9c0f7a54-690c-447e-9b20-04bb586e36ec
# ╟─77af9dac-fae3-4ab2-9d72-2367fcfe22e0
# ╟─41782dde-3f35-4879-aa89-1ebadd4bf8af
# ╟─d1f7eba1-0c09-4dcd-af6b-087699869f31
# ╟─e53de410-3bad-4e07-b94f-2285c9ed8c61
# ╟─12269268-3fde-4e4e-a954-1e9e3f27fb71
# ╟─faf91165-0026-48ce-8228-7877426401ab
# ╟─28af146b-9eb2-4490-b89e-fd9cd2965d37
# ╟─f9011bcf-d9fd-41dd-a7b5-e86a461ef390
# ╟─79ff49e9-4c0d-4e9e-bd6b-080abdf86020
# ╟─ee6f5d82-6b11-4d09-8299-2f9b3b97e809
# ╠═e9be3b56-87ac-4ec1-9c8b-dfd7cb12c762
# ╠═3b72ff7b-2d04-44fc-8271-3509b726005d
# ╟─0a108e03-08a7-4a73-be44-cb75b2943058
# ╠═82b1689d-e968-4434-9100-48c18a295080
# ╠═846bf22d-6ee5-4ed5-9165-7428bf348a89
# ╠═af3a712b-5d10-4ec5-8f74-786f15169af5
# ╟─4110293d-7f40-4d9c-aefa-268042c8db07
# ╟─3f185276-9f70-499e-a1be-4d1eed1fec2d
# ╟─18ed2e4c-f512-4662-b23b-317862a483bc
# ╟─fe9f3c1c-1241-40da-9a9b-6c317c78cd96
# ╟─eeed2564-43a8-4287-aace-6446093b07e0
# ╟─0ee314e0-5be0-4d92-9a11-7ba10e8c920d
# ╟─e5c74165-cde4-44a4-8387-71e06894e785
# ╟─c9c5716a-74f3-4365-98bc-daaabe47dd88
# ╟─084998bc-14f9-4e54-9704-14e02295577c
# ╟─a571a9c3-7eb2-4fcb-b5db-052fb01b811a
# ╟─e9dc184d-253f-45d7-ac44-0a1cd02ee746
# ╟─f7183a91-8cf2-4b99-8cd5-dd250a9d07a1
# ╟─9b898455-35b4-40db-b585-994f06ebc496
# ╟─deabaa69-aae4-4630-8e34-804979c391c2
# ╠═eb0ba4b6-4df9-487f-9451-30872181c523
# ╟─2bb2f9f0-9e97-4ff4-97a8-acd62afe1d8d
# ╠═146f5e6d-9f7d-439f-8bd8-5f18ba206e80
# ╟─255e1650-8b9f-47f3-b92e-93e10df238bd
# ╠═62570da5-f5db-45f8-bbcf-584ed486b9fe
# ╟─a46b4c33-0d8d-4a85-a61b-5fe11edc76e0
# ╠═0d911b10-a700-4ddb-88bf-8a3201995f3c
# ╟─f3a56c3e-52d2-4e33-9220-c6b6ff413973
# ╟─cac10de9-50f4-4dd9-9238-338c0e8dc37b
# ╠═8d3f8d08-353f-4896-83af-cf37c0bd2f06
# ╟─f7dabe72-edb7-4188-91c7-a3285e82d408
# ╟─d5198c78-e073-448f-8e54-c5b854daa773
# ╟─32fc8230-82fd-461c-b219-65c3ba717dd1
# ╟─03e45498-c6b5-4d3d-934a-60cc0a4ac620
# ╟─27ee0554-8df8-4449-bcc2-352dc4556bd8
# ╟─d445d308-66ec-468c-836e-d2a15127ba2b
# ╟─7c00704e-1272-4b15-a84b-8110627274ad
# ╟─78dd7f6c-568f-40bf-96af-0c74bde4e98d
# ╟─9bd6dc7f-f449-4ef1-bd16-ca2d47f983eb
# ╟─fcdf7707-c675-4104-9a4c-03b12ffd1024
# ╠═87e5238d-01c6-422a-b847-16fe1dc9971f
# ╠═cb2dab5c-a719-47e3-8df2-b6a191f8dccb
# ╟─62658932-19b1-4103-a0bc-653825f8bfe3
# ╠═e40ede69-afa7-4d0a-b2b4-67c91d1f5faa
# ╟─67f12dea-67db-49a6-816f-678f88577369
# ╠═1a2cac05-52ab-4483-9b97-d4453173f626
# ╟─77b40da6-4577-4968-a48d-f4ba7c6d1bca
# ╟─be040c96-da49-44e6-9a73-7e26a1960261
# ╟─6863be47-bced-4a2b-8578-ff4fc7bc1070
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
