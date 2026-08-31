/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ArithDiv
import Mathlib.RingTheory.DedekindDomain.Factorization
import ABC3.Meta.Claim

/-!
# **分数イデアルの群 ≃+ 有限素点の自由アーベル群**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where cv ∈ Z if v ∈ V(F )non and cv ∈ R if v ∈ V(F )arc [cf. [Szp], §1.1]. Here, if

## ★★★★★★★★台帳 `arakelov-degF-finite-places` の第 2 段

`ADiv(F) = (FinitePlace F →₀ ℤ) × (InfinitePlace F →₀ ℝ)`（`ArithDiv.lean`）の
**有限素点部分**を、分数イデアルの群と繋ぐ。

    `Additive (FractionalIdeal (𝓞_F)⁰ F)ˣ ≃+ (HeightOneSpectrum (𝓞_F) →₀ ℤ)`

★これは Dedekind 環の**一意分解**そのものである。

## ★★★mathlib の在庫は「部品」であって「群同型」ではない

2026-08-28 の実測:

| mathlib | 内容 |
|---|---|
| `FractionalIdeal.count K v I` | ★`v` での指数 |
| `FractionalIdeal.count_mul` / `count_one` | ★★加法性（非零のとき） |
| `FractionalIdeal.finite_factors` | ★台が有限であること |
| `FractionalIdeal.count_finsuppProd` | ★★`count (∏ v^{e v}) = e v` |
| `FractionalIdeal.finprod_heightOneSpectrum_factorization'` | ★★★`∏ v^{count v I} = I` |

★★**群同型として束ねられてはいない**。本ファイルがそれを束ねる。

## ★★★★★摩擦 —— `finprod` と `Finsupp.prod`

`finprod_heightOneSpectrum_factorization'` は `∏ᶠ`（`finprod`）で述べられ、
`Finsupp.prod` は `∏ x ∈ support` である。
★橋は `finprod_eq_finsetProd_of_mulSupport_subset` で、
`mulSupport ⊆ support` は「`count v I = 0` なら `v.asIdeal ^ 0 = 1`」から出る
（mathlib の `count_finprod` の証明と同じ形）。

## ★残っている段（明示）

★★これで `arakelov-degF-finite-places` は
**第 1 段（`Pic(Spec 𝓞_F) ≃* Cl(F)`）と第 2 段（本ファイル）**が揃った。
残るのは第 3 段 —— `ADiv(F)/APrc(F) ≃* APicOf (Spec 𝓞_F)` と `deg_F` の転送である。
-/

namespace ABC3.Found.GenEll

open IsDedekindDomain FractionalIdeal Function NumberField

variable {R : Type} [CommRing R] (K : Type) [Field K] [Algebra R K] [IsFractionRing R K]
  [IsDedekindDomain R]

/-! ## ★分数イデアルの単元の因子 -/

/-- ★**分数イデアルの単元の因子**——`count` を `Finsupp` に束ねたもの。 -/
noncomputable def unitCountFinsupp (I : (FractionalIdeal (nonZeroDivisors R) K)ˣ) :
    HeightOneSpectrum R →₀ ℤ :=
  Finsupp.mk
    (FractionalIdeal.finite_factors (K := K)
      (I : FractionalIdeal (nonZeroDivisors R) K)).toFinset
    (fun v => FractionalIdeal.count K v (I : FractionalIdeal (nonZeroDivisors R) K))
    (fun _ => (FractionalIdeal.finite_factors _).mem_toFinset)

@[simp] theorem unitCountFinsupp_apply (I : (FractionalIdeal (nonZeroDivisors R) K)ˣ)
    (v : HeightOneSpectrum R) :
    unitCountFinsupp K I v
      = FractionalIdeal.count K v (I : FractionalIdeal (nonZeroDivisors R) K) :=
  rfl

theorem unit_ne_zero (I : (FractionalIdeal (nonZeroDivisors R) K)ˣ) :
    (I : FractionalIdeal (nonZeroDivisors R) K) ≠ 0 := I.ne_zero

/-- ★★**因子は加法的**（`count_mul`）。 -/
noncomputable def unitCountHom :
    Additive (FractionalIdeal (nonZeroDivisors R) K)ˣ →+ (HeightOneSpectrum R →₀ ℤ) where
  toFun I := unitCountFinsupp K (Additive.toMul I)
  map_zero' := by
    ext v
    show FractionalIdeal.count K v ((1 : (FractionalIdeal (nonZeroDivisors R) K)ˣ) :
      FractionalIdeal (nonZeroDivisors R) K) = 0
    rw [Units.val_one, FractionalIdeal.count_one]
  map_add' I J := by
    ext v
    show FractionalIdeal.count K v
        (((Additive.toMul I * Additive.toMul J : (FractionalIdeal (nonZeroDivisors R) K)ˣ)) :
          FractionalIdeal (nonZeroDivisors R) K) = _
    rw [Units.val_mul, FractionalIdeal.count_mul K v (unit_ne_zero K _) (unit_ne_zero K _)]
    rfl

/-! ## ★★逆向き —— 因子から分数イデアルへ -/

theorem finsuppProd_ne_zero (exps : HeightOneSpectrum R →₀ ℤ) :
    (exps.prod (fun v n => (v.asIdeal : FractionalIdeal (nonZeroDivisors R) K) ^ n)) ≠ 0 := by
  rw [Finsupp.prod]
  exact Finset.prod_ne_zero_iff.2 (fun v _ => zpow_ne_zero _ (coeIdeal_ne_zero.mpr v.ne_bot))

/-- ★★**因子から分数イデアルの単元へ**——Dedekind 環では非零な分数イデアルは可逆。 -/
noncomputable def finsuppToUnit (exps : HeightOneSpectrum R →₀ ℤ) :
    (FractionalIdeal (nonZeroDivisors R) K)ˣ :=
  (Ne.isUnit (finsuppProd_ne_zero K exps)).unit

@[simp] theorem finsuppToUnit_val (exps : HeightOneSpectrum R →₀ ℤ) :
    ((finsuppToUnit K exps : (FractionalIdeal (nonZeroDivisors R) K)ˣ)
      : FractionalIdeal (nonZeroDivisors R) K)
      = exps.prod (fun v n => (v.asIdeal : FractionalIdeal (nonZeroDivisors R) K) ^ n) :=
  IsUnit.unit_spec _

theorem unitCountFinsupp_finsuppToUnit (exps : HeightOneSpectrum R →₀ ℤ) :
    unitCountFinsupp K (finsuppToUnit K exps) = exps := by
  ext v
  show FractionalIdeal.count K v
    ((finsuppToUnit K exps : (FractionalIdeal (nonZeroDivisors R) K)ˣ)
      : FractionalIdeal (nonZeroDivisors R) K) = _
  rw [finsuppToUnit_val, FractionalIdeal.count_finsuppProd]

/-- ★★★**一意分解**——`∏ v^{count v I} = I`。

★`finprod` と `Finsupp.prod` の橋は `finprod_eq_finsetProd_of_mulSupport_subset` である。 -/
theorem finsuppProd_count (I : (FractionalIdeal (nonZeroDivisors R) K)ˣ) :
    (unitCountFinsupp K I).prod
        (fun v n => (v.asIdeal : FractionalIdeal (nonZeroDivisors R) K) ^ n)
      = (I : FractionalIdeal (nonZeroDivisors R) K) := by
  have hsub : mulSupport
      (fun v : HeightOneSpectrum R =>
        (v.asIdeal : FractionalIdeal (nonZeroDivisors R) K)
          ^ FractionalIdeal.count K v (I : FractionalIdeal (nonZeroDivisors R) K))
      ⊆ ((unitCountFinsupp K I).support : Set (HeightOneSpectrum R)) := by
    intro v hv
    simp only [Finsupp.mem_support_iff, Finset.mem_coe, unitCountFinsupp_apply]
    intro h
    rw [mem_mulSupport, h, zpow_zero] at hv
    exact hv rfl
  have hfin := finprod_eq_finsetProd_of_mulSupport_subset
    (fun v : HeightOneSpectrum R =>
      (v.asIdeal : FractionalIdeal (nonZeroDivisors R) K)
        ^ FractionalIdeal.count K v (I : FractionalIdeal (nonZeroDivisors R) K)) hsub
  rw [Finsupp.prod]
  show ∏ a ∈ (unitCountFinsupp K I).support,
      (a.asIdeal : FractionalIdeal (nonZeroDivisors R) K)
        ^ FractionalIdeal.count K a (I : FractionalIdeal (nonZeroDivisors R) K) = _
  rw [← hfin]
  exact FractionalIdeal.finprod_heightOneSpectrum_factorization' (hI := unit_ne_zero K I)

/-! ## ★★★★★★★★群同型 -/

/-- ★★★★★★★★**分数イデアルの群 ≃+ 有限素点の自由アーベル群**（Dedekind の一意分解）。

原文 (GenEll p.4):
> — where cv ∈ Z if v ∈ V(F )non and cv ∈ R if v ∈ V(F )arc [cf. [Szp], §1.1]. Here, if

★これが `ADiv(F)` の**有限素点部分**の意味である。 -/
noncomputable def fractionalIdealCountEquiv :
    Additive (FractionalIdeal (nonZeroDivisors R) K)ˣ ≃+ (HeightOneSpectrum R →₀ ℤ) where
  toFun := unitCountHom K
  invFun exps := Additive.ofMul (finsuppToUnit K exps)
  left_inv I := by
    apply Additive.toMul.injective
    apply Units.ext
    show ((finsuppToUnit K (unitCountFinsupp K (Additive.toMul I))
      : (FractionalIdeal (nonZeroDivisors R) K)ˣ) : FractionalIdeal (nonZeroDivisors R) K) = _
    rw [finsuppToUnit_val, finsuppProd_count]
  right_inv exps := unitCountFinsupp_finsuppToUnit K exps
  map_add' := (unitCountHom K).map_add

/-- ★★★★★★★★★**`ADiv(F)` の有限素点部分は分数イデアルの群である**。

原文 (GenEll p.4):
> — where cv ∈ Z if v ∈ V(F )non and cv ∈ R if v ∈ V(F )arc [cf. [Szp], §1.1]. Here, if -/
noncomputable def adivFinEquiv (F : Type) [Field F] [NumberField F] :
    Additive (FractionalIdeal (nonZeroDivisors (𝓞 F)) F)ˣ ≃+ (FinitePlace F →₀ ℤ) :=
  fractionalIdealCountEquiv F

/-! ### ★出典の紐付け(`.src`) -/

def fractionalIdealCountEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "§1 地の文(算術因子の有限素点部分——分数イデアルの群 ≃+ 自由アーベル群)",
    sectionId := "genell-deg" }

def adivFinEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.2, (i)(層 D——ADiv(F) の有限素点部分は分数イデアルの群)",
    sectionId := "genell-def-1-2-i" }

def fractionalIdealCountEquiv.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "FractionalIdeal.count(素点での指数)"
      (.inMathlib "FractionalIdeal.count") 4,
    .citation "[mathlib]" "FractionalIdeal.count_mul / count_one(加法性)"
      (.inMathlib "FractionalIdeal.count_mul") 4,
    .citation "[mathlib]" "FractionalIdeal.finite_factors(台の有限性)"
      (.inMathlib "FractionalIdeal.finite_factors") 4,
    .citation "[mathlib]" "FractionalIdeal.finprod_heightOneSpectrum_factorization'(一意分解)"
      (.inMathlib "FractionalIdeal.finprod_heightOneSpectrum_factorization'") 4,
    .implicitStep
      ("★mathlib の在庫は**部品**であって群同型ではない(2026-08-28 実測)。" ++
       "本ファイルがそれを束ねる") 4,
    .implicitStep
      ("★★これで台帳 arakelov-degF-finite-places は第 1 段" ++
       "(Pic(Spec 𝓞_F) ≃* Cl(F))と第 2 段(本ファイル)が揃った。" ++
       "残るのは第 3 段——ADiv(F)/APrc(F) ≃* APicOf (Spec 𝓞_F) と deg_F の転送") 4 ]

end ABC3.Found.GenEll
