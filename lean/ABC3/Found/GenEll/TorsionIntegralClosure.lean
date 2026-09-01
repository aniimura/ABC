/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TorsionIntegralGood
import ABC3.Found.GenEll.VeluDescent
import ABC3.Meta.Claim

/-!
# 第 1155 ブロック —— **捩れ点の座標は `L̄` でも整**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——節点 2 の道の第 1 歩

`Skeleton/GenEll/LCyclicReading.lean` の節点 2 は
「`Lemma 3.5` を安定直線の側で述べ直す」であった。

★★★★**2026-09-01（第 1155）の測定——道が 1 本短くなる**

`Found/GaloisRep/TorsionIntegralGood.lean` は捩れ点の座標が
**`L` の付値環 `primeSubring p` に属する**ことを、付値の言葉で示している
（`v(x) < 0` なら深さ `m` が取れて… という議論）。

☆安定直線の側では点の座標は `L̄` にあり、`L̄` に `p` の付値は（一意には）伸びない。
★**だが `Lemma 3.5` が実際に要るのは Vélu の和 `v`・`w` が整であることだけ**である。
★★そして第 1154 で **`v`・`w` は `L` の元**だと分かった。

☆したがって次の 2 段でよい:

| 段 | 内容 |
|---|---|
| 1 | `L̄` の捩れ点の座標は `primeSubring p` 上**整**である（付値は要らない） |
| 2 | `v`・`w` はその多項式なので整、かつ `L` の元。`primeSubring p` は `L` で整閉なので**属する** |

★**付値の議論がまるごと要らなくなる**——本ファイルは第 1 段を取る。

## ★機構

`ΨSq_l` の主係数は `l²`（`leadingCoeff_ΨSq`、偶奇不問）で、`hlu` からそれは単元。
☆したがって `x` は `primeSubring p` 上整である（`isIntegral_of_isUnit_leadingCoeff`）。
★`y` は Weierstrass 方程式が `y` について**モニック**なので、`x` が整なら整である。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep WeierstrassCurve IsDedekindDomain NumberField
open scoped Classical

variable {L : Type} [Field L] [NumberField L]

local notation "Lbar" => AlgebraicClosure L

/-! ## ★★★★★★★★捩れ点の `x` は整 -/

/-- ★★★★★★★★★★★★★★**`L̄` の位数 `l` の点の `x` は `primeSubring p` 上整**
——★**偶奇を問わない**（第 1155）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`Found/GaloisRep/TorsionIntegralGood.lean` の `mem_primeSubring_x_of_addOrderOf_prime'`
（第 1148）は `L` の点で「付値環に属する」を出す。
★本定理は `L̄` の点で「整である」を出す——**付値を使わない**ので `L̄` でもそのまま通る。 -/
theorem isIntegral_x_of_addOrderOf_prime (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [WeierstrassCurve.IsIntegral (primeSubring p) W]
    {l : ℕ} (hl : Nat.Prime l) (hlu : IsUnit ((l : primeSubring p)))
    {x y : Lbar}
    (h : (W.map (algebraMap L Lbar)).toAffine.Nonsingular x y)
    (hQ : addOrderOf (WeierstrassCurve.Affine.Point.some x y h) = l) :
    IsIntegral (primeSubring p) x := by
  have hroot : ((W.map (algebraMap L Lbar)).ΨSq (l : ℤ)).eval x = 0 :=
    ΨSq_eval_eq_zero_of_addOrderOf_prime _ h hl hQ
  set R := primeSubring p with hR
  set Wi := WeierstrassCurve.integralModel R W with hWi
  have hbc : Wi.baseChange L = W := WeierstrassCurve.baseChange_integralModel_eq R W
  -- ☆`W⁄L̄` は `Wi` を `R → L̄` で送ったものである
  have hcomp : (algebraMap L Lbar).comp (algebraMap R L) = algebraMap R Lbar := by
    ext r
    show algebraMap L Lbar (algebraMap R L r) = algebraMap R Lbar r
    rw [IsScalarTower.algebraMap_apply R L Lbar]
  have hWmap : W.map (algebraMap L Lbar) = Wi.map (algebraMap R Lbar) := by
    conv_lhs => rw [← hbc]
    rw [WeierstrassCurve.baseChange, WeierstrassCurve.map_map, hcomp]
  -- ☆分点多項式も同じように送られる
  have hmap : (W.map (algebraMap L Lbar)).ΨSq (l : ℤ)
      = (Wi.ΨSq (l : ℤ)).map (algebraMap R Lbar) := by
    rw [hWmap]
    exact Wi.map_ΨSq (algebraMap R Lbar) (l : ℤ)
  have hne : (((l : ℤ)) : R) ≠ 0 := by
    have hc : (((l : ℤ)) : R) = ((l : ℕ) : R) := by push_cast; ring
    rw [hc]
    exact hlu.ne_zero
  have hlc : (Wi.ΨSq (l : ℤ)).leadingCoeff = (((l : ℤ)) : R) ^ 2 :=
    Wi.leadingCoeff_ΨSq hne
  have hu : IsUnit (Wi.ΨSq (l : ℤ)).leadingCoeff := by
    rw [hlc]
    have hc : (((l : ℤ)) : R) = ((l : ℕ) : R) := by push_cast; ring
    rw [hc]
    exact hlu.pow 2
  have haeval : Polynomial.aeval x (Wi.ΨSq (l : ℤ)) = 0 := by
    rw [Polynomial.aeval_def, ← Polynomial.eval_map, ← hmap]
    exact hroot
  exact isIntegral_of_isUnit_leadingCoeff haeval hu

/-! ## ☆付値環の元の像は整 -/

/-- ☆`primeSubring p` の元の `L̄` での像は `primeSubring p` 上整である。 -/
theorem isIntegral_algebraMap_of_mem (p : HeightOneSpectrum (𝓞 L)) {z : L}
    (hz : z ∈ primeSubring p) : IsIntegral (primeSubring p) (algebraMap L Lbar z) := by
  have h1 : algebraMap (primeSubring p) Lbar ⟨z, hz⟩ = algebraMap L Lbar z := by
    rw [IsScalarTower.algebraMap_apply (primeSubring p) L Lbar]
    rfl
  rw [← h1]
  exact isIntegral_algebraMap

/-- ☆数は整である。 -/
theorem isIntegral_ofNat (p : HeightOneSpectrum (𝓞 L)) (n : ℕ) :
    IsIntegral (primeSubring p) ((n : Lbar)) := by
  have h1 : ((n : Lbar)) = algebraMap (primeSubring p) Lbar (n : primeSubring p) := by
    rw [map_natCast]
  rw [h1]
  exact isIntegral_algebraMap

/-! ## ★★★★★★★★`x` が整なら `y` も整 -/

/-- ★★★★★★★★★★**`x` が整なら `y` も整**——★**無条件**（第 1156）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆Weierstrass 方程式 `y² + a₁xy + a₃y = x³ + a₂x² + a₄x + a₆` は
`y` について**モニックな 2 次式**である。★係数は `x` と `aᵢ` の多項式で、
どれも `R = primeSubring p` 上整なので `A ≔ integralClosure R L̄` の元である。
☆したがって `y` は `A` 上整であり、`isIntegral_trans` で `R` 上整である。 -/
theorem isIntegral_y_of_isIntegral_x (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [WeierstrassCurve.IsIntegral (primeSubring p) W]
    {x y : Lbar}
    (heq : (W.map (algebraMap L Lbar)).toAffine.Equation x y)
    (hx : IsIntegral (primeSubring p) x) :
    IsIntegral (primeSubring p) y := by
  set R := primeSubring p with hR
  obtain ⟨ha1, ha2, ha3, ha4, ha6⟩ := mem_primeSubring_of_isIntegral p W
  -- ☆`R` の元の像は整
  have hint : ∀ z : L, z ∈ R → IsIntegral R (algebraMap L Lbar z) := fun _ hz =>
    isIntegral_algebraMap_of_mem p hz
  have h1 : IsIntegral R (algebraMap L Lbar W.a₁) := hint _ ha1
  have h2 : IsIntegral R (algebraMap L Lbar W.a₂) := hint _ ha2
  have h3 : IsIntegral R (algebraMap L Lbar W.a₃) := hint _ ha3
  have h4 : IsIntegral R (algebraMap L Lbar W.a₄) := hint _ ha4
  have h6 : IsIntegral R (algebraMap L Lbar W.a₆) := hint _ ha6
  -- ★係数を `A = integralClosure R L̄` の元として取る
  set A := integralClosure R Lbar with hA
  have hbmem : (algebraMap L Lbar W.a₁ * x + algebraMap L Lbar W.a₃) ∈ A :=
    (h1.mul hx).add h3
  have hcmem : (x ^ 3 + algebraMap L Lbar W.a₂ * x ^ 2
      + algebraMap L Lbar W.a₄ * x + algebraMap L Lbar W.a₆) ∈ A :=
    (((hx.pow 3).add (h2.mul (hx.pow 2))).add (h4.mul hx)).add h6
  set b : A := ⟨_, hbmem⟩ with hb
  set c : A := ⟨_, hcmem⟩ with hc
  refine isIntegral_trans (A := A) y ?_
  refine ⟨Polynomial.X ^ 2 + (Polynomial.C b * Polynomial.X - Polynomial.C c), ?_, ?_⟩
  · monicity!
  · rw [WeierstrassCurve.Affine.equation_iff] at heq
    simp only [Polynomial.eval₂_add, Polynomial.eval₂_sub, Polynomial.eval₂_mul,
      Polynomial.eval₂_pow, Polynomial.eval₂_X, Polynomial.eval₂_C]
    show y ^ 2 + ((algebraMap L Lbar W.a₁ * x + algebraMap L Lbar W.a₃) * y
      - (x ^ 3 + algebraMap L Lbar W.a₂ * x ^ 2
          + algebraMap L Lbar W.a₄ * x + algebraMap L Lbar W.a₆)) = 0
    simp only [WeierstrassCurve.map_a₁, WeierstrassCurve.map_a₂, WeierstrassCurve.map_a₃,
      WeierstrassCurve.map_a₄, WeierstrassCurve.map_a₆] at heq
    linear_combination heq

/-! ## ★★★★★★★★★★★★`v` の側は整 -/

/-- ★★★★★★★★**座標が整なら `veluV2` も整**——★**無条件**（第 1157）。 -/
theorem isIntegral_veluV2 (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [WeierstrassCurve.IsIntegral (primeSubring p) W]
    {x y : Lbar} (hx : IsIntegral (primeSubring p) x) (hy : IsIntegral (primeSubring p) y) :
    IsIntegral (primeSubring p) (veluV2 (W.map (algebraMap L Lbar)) x y) := by
  obtain ⟨ha1, ha2, ha3, ha4, ha6⟩ := mem_primeSubring_of_isIntegral p W
  have h3 : IsIntegral (primeSubring p) ((3 : Lbar)) := by
    have := isIntegral_ofNat (L := L) p 3
    rwa [Nat.cast_ofNat] at this
  have h2 : IsIntegral (primeSubring p) ((2 : Lbar)) := by
    have := isIntegral_ofNat (L := L) p 2
    rwa [Nat.cast_ofNat] at this
  show IsIntegral (primeSubring p) (veluGx (W.map (algebraMap L Lbar)) x y)
  rw [veluGx]
  simp only [WeierstrassCurve.map_a₁, WeierstrassCurve.map_a₂, WeierstrassCurve.map_a₄]
  exact (((h3.mul (hx.pow 2)).add
      ((h2.mul (isIntegral_algebraMap_of_mem p ha2)).mul hx)).add
      (isIntegral_algebraMap_of_mem p ha4)).sub
      ((isIntegral_algebraMap_of_mem p ha1).mul hy)

/-- ★★★★★★★★★★**座標がすべて整なら `veluVFull` も整**——★**無条件**（第 1157）。 -/
theorem isIntegral_veluVFull (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [WeierstrassCurve.IsIntegral (primeSubring p) W]
    (S : Finset (Lbar × Lbar))
    (hS : ∀ z ∈ S, IsIntegral (primeSubring p) z.1 ∧ IsIntegral (primeSubring p) z.2) :
    IsIntegral (primeSubring p) (veluVFull (W.map (algebraMap L Lbar)) S) := by
  rw [veluVFull]
  refine Subalgebra.sum_mem (integralClosure (primeSubring p) Lbar) (fun z hz => ?_)
  exact isIntegral_veluV2 p W (hS z hz).1 (hS z hz).2

/-! ## ★★★★★★★★★★★★★★★★整で `L` の元なら付値環に属する -/

/-- ★★★★★★★★★★★★★★★★
**`L` の元がその像で整なら付値環に属する**——★**無条件**（第 1157）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★これが節点 2 の要である——`Gal`-安定な部分群の Vélu の和は
第 1154 で **`L` の元**、本ファイルで **`R` 上整**、そして `R` は `L` で整閉なので
**`R` に属する**。☆付値の議論はどこにも要らない。 -/
theorem mem_primeSubring_of_isIntegral_image (p : HeightOneSpectrum (𝓞 L)) {z : L}
    (hz : IsIntegral (primeSubring p) (algebraMap L Lbar z)) :
    z ∈ primeSubring p := by
  have hinj : Function.Injective (algebraMap L Lbar) := (algebraMap L Lbar).injective
  have hz' : IsIntegral (primeSubring p) z := (isIntegral_algebraMap_iff hinj).mp hz
  obtain ⟨w, hw⟩ := IsIntegrallyClosed.isIntegral_iff.1 hz'
  exact hw ▸ w.2

/-! ## ★出典の紐付け(`.src`) -/

def isIntegral_x_of_addOrderOf_prime.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(L̄ の位数 l の点の x は primeSubring p 上整。★偶奇不問、付値を使わない)",
    sectionId := "genell-lemma-3-5" }

def isIntegral_veluVFull.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(座標がすべて整なら veluVFull も整。★無条件)",
    sectionId := "genell-lemma-3-5" }

def mem_primeSubring_of_isIntegral_image.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(L の元がその像で整なら付値環に属する——整閉だから。★無条件)",
    sectionId := "genell-lemma-3-5" }

def isIntegral_y_of_isIntegral_x.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(x が整なら y も整——Weierstrass 方程式は y についてモニック。★無条件)",
    sectionId := "genell-lemma-3-5" }

def isIntegral_y_of_isIntegral_x.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "isIntegral_trans(A 上整で A が R 上整なら R 上整)"
      (.inMathlib "isIntegral_trans") 1 ]

def isIntegral_x_of_addOrderOf_prime.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "ΨSq_eval_eq_zero_of_addOrderOf_prime(第 1148、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.ΨSq_eval_eq_zero_of_addOrderOf_prime") 1,
    .citation "[mathlib]" "WeierstrassCurve.leadingCoeff_ΨSq(主係数は n²、偶奇不問)"
      (.inMathlib "WeierstrassCurve.leadingCoeff_ΨSq") 1,
    .implicitStep
      ("★★★★**2026-09-01（第 1155）の測定**——`TorsionIntegralGood.lean` は" ++
       "捩れ点の座標が付値環に**属する**ことを付値の言葉で示しているが、" ++
       "安定直線の側では点の座標は `L̄` にあり `p` の付値は一意には伸びない。" ++
       "☆しかし `Lemma 3.5` が要るのは Vélu の和 `v`・`w` が整であることだけで、" ++
       "第 1154 でそれらは `L` の元だと分かった。" ++
       "★したがって「`L̄` で整」→「`L` の元で整」→「整閉なので属する」の 3 段でよく、" ++
       "**付値の議論がまるごと要らなくなる**。") 6 ]

end ABC3.Found.GenEll
