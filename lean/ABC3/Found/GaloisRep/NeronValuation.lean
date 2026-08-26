import ABC3.Found.GaloisRep.NeronFinite

/-!
# Galois (G7) 第 315 ブロック —— **★★★★★★★素点の付値で測った Néron 指数**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★到達点

> Néron 指数を **`p` の adic 付値**で測り直した `neronExp p W` を作り、
> **一意性・変数変換則・有限性**をすべて移した。

★★★これで **mathlib の分数イデアルの API(`p.valuation` で書かれている)に載る**。

## ★★★★★★★★2 つの付値は「同値」で足りる

`minimalExp`(第 312)は**局所化 `R_p` の極大イデアルの付値**で測っていた。
一方 mathlib の分数イデアル論は **`p` の adic 付値**で書かれている。
★この 2 つが**等しい**ことを示すのは(正規化の一致まで要るので)手間だが、
★★★★**同値(`Valuation.IsEquiv`)で十分**である:

* 一意性には `IsEquiv.eq_iff`(`v₁ x = v₁ y ↔ v₂ x = v₂ y`)
* 有限性には `IsEquiv.eq_one_iff_eq_one`(`v₁ x = 1 ↔ v₂ x = 1`)

★★★★★同値は**付値環が一致すること**から出る(`isEquiv_iff_valuationSubring`)。
そして付値環はどちらも `p.valuationSubringAtPrime L`——
片方は第 301 の `dvr_mem_of_nonneg`、もう片方は mathlib の
`valuationSubringAtPrime_eq_valuationSubring` から。

## ★★★★「等しい」を避けて「同値」で通す

★★等式が要るのは**値そのものを使う**ときだけで、本ブロックで要るのは
**「等しいかどうか」「1 かどうか」**だけである。
★★★★**必要な強さを見極めると、正規化の比較を丸ごと避けられた。**

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `valuationSubring_primeSubring_eq` | ★★★★局所化の付値環は `primeSubring` |
| `valuation_isEquiv` | ★★★★★★**2 つの付値は同値** |
| `valAdd`・`valAdd_mul` | ★★`p` の加法付値 |
| `neronExp`・`neronExp_eq` | ★★★★★**素点の付値で測った Néron 指数** |
| `neronExp_variableChange`・`finite_bad_primes'` | ★★★★★★変数変換則と有限性 |
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField WeierstrassCurve

/-! ## ★★★★局所化の付値環 -/

section Equiv

variable {L : Type} [Field L] [NumberField L]

/-- ★★★★局所化の極大イデアルの付値環は `primeSubring p` そのもの。 -/
theorem valuationSubring_primeSubring_eq (p : HeightOneSpectrum (𝓞 L)) :
    ((IsDiscreteValuationRing.maximalIdeal (primeSubring p)).valuation L).valuationSubring
      = (primeSubring p) := by
  ext x
  rw [Valuation.mem_valuationSubring_iff]
  constructor
  · intro h
    rcases eq_or_ne x 0 with rfl | hx
    · exact zero_mem _
    · have h1 : 0 ≤ vAdd (tateDvrVal (primeSubring p) L) (Units.mk0 x hx) :=
        (tateDvrVal_nonneg_iff (Units.mk0 x hx)).2 h
      obtain ⟨y, hy⟩ := dvr_mem_of_nonneg (Units.mk0 x hx) h1
      have hy' : ((y : L)) = x := hy
      rw [← hy']
      exact y.2
  · intro h
    exact IsDedekindDomain.HeightOneSpectrum.valuation_le_one _ (⟨x, h⟩ : primeSubring p)

/-- ★★★★★★**2 つの付値は同値**——付値環が一致するから。 -/
theorem valuation_isEquiv (p : HeightOneSpectrum (𝓞 L)) :
    ((IsDiscreteValuationRing.maximalIdeal (primeSubring p)).valuation L).IsEquiv
      (p.valuation L) := by
  rw [Valuation.isEquiv_iff_valuationSubring, valuationSubring_primeSubring_eq]
  exact (IsDedekindDomain.HeightOneSpectrum.valuationSubringAtPrime_eq_valuationSubring p)

end Equiv

/-! ## ★★`p` の加法付値 -/

section ValAdd

variable {L : Type} [Field L] [NumberField L]

theorem valuationP_ne_zero (p : HeightOneSpectrum (𝓞 L)) (x : Lˣ) :
    (p.valuation L) (x : L) ≠ 0 := (Valuation.ne_zero_iff _).2 x.ne_zero

/-- ★★素点 `p` での加法付値(`Lˣ` 上)。 -/
noncomputable def valAdd (p : HeightOneSpectrum (𝓞 L)) (x : Lˣ) : ℤ :=
  -Multiplicative.toAdd (WithZero.unzero (valuationP_ne_zero p x))

theorem valAdd_eq_of_valuation_eq (p : HeightOneSpectrum (𝓞 L)) (x y : Lˣ)
    (h : (p.valuation L) (x : L) = (p.valuation L) (y : L)) : valAdd p x = valAdd p y := by
  rw [valAdd, valAdd]
  congr 1
  refine congrArg Multiplicative.toAdd ?_
  apply WithZero.coe_injective
  rw [WithZero.coe_unzero, WithZero.coe_unzero, h]

theorem valAdd_mul (p : HeightOneSpectrum (𝓞 L)) (x y : Lˣ) :
    valAdd p (x * y) = valAdd p x + valAdd p y := by
  have h : WithZero.unzero (valuationP_ne_zero p (x * y))
      = WithZero.unzero (valuationP_ne_zero p x) * WithZero.unzero (valuationP_ne_zero p y) := by
    apply WithZero.coe_injective
    rw [WithZero.coe_unzero, WithZero.coe_mul, WithZero.coe_unzero, WithZero.coe_unzero]
    push_cast
    exact Valuation.map_mul _ _ _
  rw [valAdd, valAdd, valAdd, h]
  simp
  omega

theorem valAdd_eq_zero_iff (p : HeightOneSpectrum (𝓞 L)) (x : Lˣ) :
    valAdd p x = 0 ↔ (p.valuation L) (x : L) = 1 := by
  rw [valAdd, neg_eq_zero]
  constructor
  · intro h
    have h2 : WithZero.unzero (valuationP_ne_zero p x) = 1 := by
      have h3 := congrArg Multiplicative.ofAdd h
      simpa using h3
    rw [← WithZero.coe_unzero (valuationP_ne_zero p x), h2]
    rfl
  · intro h
    have h2 : WithZero.unzero (valuationP_ne_zero p x) = 1 := by
      apply WithZero.coe_injective
      rw [WithZero.coe_unzero, h]
      rfl
    rw [h2]
    rfl

end ValAdd

/-! ## ★★★★★素点の付値で測った Néron 指数 -/

section NeronExp

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

theorem valuation_eq_of_vAdd_eq (x y : Kˣ)
    (h : vAdd (tateDvrVal R K) x = vAdd (tateDvrVal R K) y) :
    (IsDiscreteValuationRing.maximalIdeal R).valuation K (x : K)
      = (IsDiscreteValuationRing.maximalIdeal R).valuation K (y : K) := by
  rw [valuation_eq_ofAdd_neg_vAdd, valuation_eq_ofAdd_neg_vAdd, h]

theorem vAdd_eq_zero_iff (x : Kˣ) :
    vAdd (tateDvrVal R K) x = 0
      ↔ (IsDiscreteValuationRing.maximalIdeal R).valuation K (x : K) = 1 := by
  rw [valuation_eq_ofAdd_neg_vAdd]
  constructor
  · intro h
    rw [h]
    rfl
  · intro h
    have h2 : (Multiplicative.ofAdd (-(vAdd (tateDvrVal R K) x)) : Multiplicative ℤ) = 1 :=
      WithZero.coe_injective h
    have h3 := congrArg Multiplicative.toAdd h2
    simpa using h3

end NeronExp

section Neron

variable {L : Type} [Field L] [NumberField L]

/-- ★★★★★**Néron 指数**(素点 `p` の付値で測ったもの)。 -/
noncomputable def neronExp (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L) : ℤ :=
  valAdd p (WeierstrassCurve.exists_isMinimal (primeSubring p) W).choose.u

/-- ★★★★どの極小化変換で測っても同じ。 -/
theorem neronExp_eq (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L) (hΔ : W.Δ ≠ 0)
    (C : VariableChange L) (h : IsMinimal (primeSubring p) (C • W)) :
    neronExp p W = valAdd p C.u := by
  refine valAdd_eq_of_valuation_eq p _ _ ?_
  refine ((valuation_isEquiv p).eq_iff).1 ?_
  refine valuation_eq_of_vAdd_eq _ _ ?_
  exact minimal_u_vAdd_eq W hΔ _ C
    (WeierstrassCurve.exists_isMinimal (primeSubring p) W).choose_spec h

/-- ★★★★2 つの指数は同時に `0` になる。 -/
theorem neronExp_eq_zero_iff (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L) :
    neronExp p W = 0 ↔ minimalExp (primeSubring p) W = 0 := by
  rw [neronExp, minimalExp, valAdd_eq_zero_iff, vAdd_eq_zero_iff]
  exact ((valuation_isEquiv p).eq_one_iff_eq_one).symm

/-- ★★★★★★**変数変換則**——界面 `omegaFrac_variableChange` の局所版。 -/
theorem neronExp_variableChange (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    (hΔ : W.Δ ≠ 0) (C₀ : VariableChange L) :
    neronExp p (C₀ • W) = neronExp p W - valAdd p C₀.u := by
  obtain ⟨C, hC⟩ := WeierstrassCurve.exists_isMinimal (primeSubring p) (C₀ • W)
  have h1 : neronExp p (C₀ • W) = valAdd p C.u :=
    neronExp_eq p _ (variableChange_Delta_ne_zero W hΔ C₀) C hC
  have h2 : IsMinimal (primeSubring p) ((C * C₀) • W) := by
    rw [mul_smul]
    exact hC
  have h3 : neronExp p W = valAdd p (C * C₀).u := neronExp_eq p W hΔ _ h2
  rw [h1, h3, show (C * C₀).u = C.u * C₀.u from rfl, valAdd_mul]
  omega

/-- ★★★★★★**指数が `0` でない素点は有限個**(第 314 から移した)。 -/
theorem finite_bad_primes' (W : WeierstrassCurve L) (hΔ : W.Δ ≠ 0) :
    {p : HeightOneSpectrum (𝓞 L) | neronExp p W ≠ 0}.Finite := by
  refine (finite_bad_primes W hΔ).subset ?_
  intro p hp
  simp only [Set.mem_setOf_eq] at hp ⊢
  intro hc
  exact hp ((neronExp_eq_zero_iff p W).2 hc)

end Neron

/-! ## ★出典の紐付け(`.src`) -/

def valuation_isEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(局所化の付値と p の adic 付値の同値性)",
    sectionId := "genell-def-3-3" }

def neronExp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(素点の付値で測った Néron 指数)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
