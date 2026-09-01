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

/-! ## ☆次の葉（本ブロックでは取らない）

☆`x` が整なら `y` も整である——Weierstrass 方程式は `y` についてモニックだからである。
★正確には `(2y + a₁x + a₃)² = 4x³ + b₂x² + 2b₄x + b₆` の右辺が整なので、
`t ≔ 2y + a₁x + a₃` は `Z² − u`（`u` は整）の根であり、整の推移律で `t` も整、
`2` が単元なら `y` も整である。

☆mathlib の `IsIntegral.trans` 周辺の navigation が要るので別の葉として置く。 -/

/-! ## ★出典の紐付け(`.src`) -/

def isIntegral_x_of_addOrderOf_prime.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(L̄ の位数 l の点の x は primeSubring p 上整。★偶奇不問、付値を使わない)",
    sectionId := "genell-lemma-3-5" }

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
