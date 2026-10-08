# [10 · History & breakthroughs](@id ch-history)

```@setup ch10
using Riemannian, CairoMakie
set_theme!(riemann_theme())
```

## Timeline

```@example ch10
plot_timeline()
```

The data behind the plot is available as a vector of named tuples:

```@example ch10
first(breakthroughs(), 3)
```

```@example ch10
timeline(category = :zeros)
```

```@example ch10
timeline(category = :computation)
```

## Riemann's 1859 paper

*Ueber die Anzahl der Primzahlen unter einer gegebenen Grösse* runs to eight pages. In it Riemann:

1. continued ``\zeta(s)`` to the whole complex plane;
2. proved the functional equation, in two ways;
3. introduced ``\xi(s)`` and conjectured its product formula (proved by Hadamard, 1893);
4. stated the zero-counting formula ``N(T)`` (proved by von Mangoldt, 1905);
5. gave the explicit formula for ``J(x)`` (proved by von Mangoldt, 1895);
6. remarked that it is "very likely" that all zeros of ``\Xi(t)`` are real, which is the Hypothesis.

His unpublished notes contained the Riemann–Siegel formula and numerical values of the first zeros.

## Landmarks by theme

**Zeros on the line.** Hardy (1914) proved that there are infinitely many. Selberg (1942) proved a
positive proportion, Levinson (1974) more than a third, Conrey (1989) more than two fifths, and
Pratt–Robles–Zaharescu–Zeindler (2020) more than five twelfths.

**Computation.** Gram (1903) computed 15 zeros. Titchmarsh (1936) reached 1041 and Turing (1953) went
further. Odlyzko's work near zero ``10^{20}`` (from 1987) revealed the GUE statistics. Gourdon (2004) checked
the first ``10^{13}`` zeros, and Platt–Trudgian (2021) verified RH up to height ``3\cdot10^{12}``.

**Equivalents.** von Koch (1901), Robin (1984), Li (1997), Lagarias (2002), and de Bruijn–Newman, where
Rodgers–Tao (2020) proved ``\Lambda \ge 0``.

**Bounds.** Zero-free regions (de la Vallée Poussin; Vinogradov–Korobov). Zero-density estimates
(Ingham 1940, Huxley 1972, **Guth–Maynard 2024**).

**Analogues.** For curves over finite fields the analogue of RH was proved by Weil (1948), and for
general varieties by Deligne (1974). Its success rests on a *geometric* interpretation (a Frobenius
operator acting on cohomology) that is still missing for ``\zeta(s)``.

## Read in the REPL

```@example ch10
explain(:breakthroughs)
```

## Functions in this chapter

[`breakthroughs`](@ref), [`timeline`](@ref), [`explain`](@ref).
