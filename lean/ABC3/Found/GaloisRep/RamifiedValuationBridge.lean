/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.DegInfLocal
import ABC3.Meta.Claim

/-!
# 第 1370 ブロック —— **分岐した拡大でも付値指数は `e` 倍でしか崩れない**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★これは何か——`hp` の**分岐版**

在庫の局所定理群は「不分岐の仮定」

    hp : ∀ x : L, v_{𝔪_R}(algebraMap L Lv x) = v_p(x)

を受ける形で書かれている（`vAdd_algebraMap_eq_valAdd`、第 321）。
☆ところが `ζ_l` を足す拡大 `L_p(ζ_l)/L_p` は剰余標数が `l` のとき**分岐する**。

★★★測定（2026-09-02、第 1370）——`exists_h2_h1_of_bad_prime`（第 1320）の証明で
`hp` が使われるのは**ただ 1 箇所**である:

    exists_h2_h1_of_bad_prime
      → local_inputs_of_bad_prime
        → not_dvd_vAdd_tateParam_of_not_dvd_jExp
          → vAdd_tateParam_eq_neg_jExp
            → jExp_eq_neg_vAdd_of_j_tateCurveAt
              → vAdd_algebraMap_eq_valAdd   ← ★ここだけ

☆他（`IsMinimal`・`HasSplitMultiplicativeReduction`・`IsElliptic`）は
インスタンス引数なので呼び出し側が別途与える。

★したがって**この 1 本を `e` 倍版に一般化すれば道が通る**。

    hpe : ∀ x : L, v_{𝔪_R}(algebraMap L Lv x) = (v_p(x))^e

☆`e` は分岐指数であり、第 1369 で `1 ≤ e ≤ [L_v′ : L_v]` を与えた。
★`L_p(ζ_l)/L_p` は `≤ l−1` 次なので `l ∤ e` である。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField

open scoped Classical

/-- ★★★★**`exp` の冪は係数倍**（第 1370）。 -/
theorem withZero_exp_pow (s : ℤ) (n : ℕ) :
    (WithZero.exp s) ^ n = WithZero.exp ((n : ℤ) * s) := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [pow_succ, ih, ← WithZero.exp_add]
    congr 1
    push_cast
    ring

variable {L : Type} [Field L] [NumberField L]
  {Lv : Type} [Field Lv] [Algebra L Lv]
  {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [Algebra R Lv] [IsFractionRing R Lv]

/-- ★★★★★★★★**`p` の付値は `exp (−valAdd)`**（第 1370）。 -/
theorem valuation_eq_exp_neg_valAdd (p : HeightOneSpectrum (𝓞 L)) (x : Lˣ) :
    (HeightOneSpectrum.valuation L p) (x : L) = WithZero.exp (-(valAdd p x)) := by
  rw [valAdd_eq_neg_log, neg_neg]
  exact (WithZero.exp_log (valuationP_ne_zero p x)).symm

/-- ★★★★★★★★★★★★★★★★
**分岐した拡大でも付値指数は `e` 倍**——★**無条件**（第 1370）。

☆`e = 1` に取れば在庫の `vAdd_algebraMap_eq_valAdd`（第 321）に戻る。 -/
theorem vAdd_algebraMap_eq_mul_valAdd {e : ℕ} (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (x : L) (hx : x ≠ 0) (hx' : algebraMap L Lv x ≠ 0) :
    vAdd (tateDvrVal R Lv) (Units.mk0 (algebraMap L Lv x) hx')
      = (e : ℤ) * valAdd p (Units.mk0 x hx) := by
  have h1 : (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
      (algebraMap L Lv x)
      = WithZero.exp (-(vAdd (tateDvrVal R Lv) (Units.mk0 (algebraMap L Lv x) hx'))) :=
    valuation_eq_ofAdd_neg_vAdd (R := R) (K := Lv) (Units.mk0 (algebraMap L Lv x) hx')
  have h2 : (HeightOneSpectrum.valuation L p) x
      = WithZero.exp (-(valAdd p (Units.mk0 x hx))) :=
    valuation_eq_exp_neg_valAdd p (Units.mk0 x hx)
  have h3 := hpe x
  rw [h1, h2, withZero_exp_pow] at h3
  have h4 := WithZero.exp_inj.1 h3
  rw [mul_neg] at h4
  exact neg_inj.mp h4

/-! ## ★出典の紐付け(`.src`) -/

def vAdd_algebraMap_eq_mul_valAdd.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(分岐した拡大でも付値指数は e 倍。★無条件)",
    sectionId := "genell-thm-3-8" }

def vAdd_algebraMap_eq_mul_valAdd.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "vAdd_algebraMap_eq_valAdd(第 321、証明済み。e = 1 の場合)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.vAdd_algebraMap_eq_valAdd") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1370）**——`exists_h2_h1_of_bad_prime`（第 1320）で " ++
       "`hp` が使われるのは `vAdd_algebraMap_eq_valAdd` の**ただ 1 箇所**である。" ++
       "☆したがってこの 1 本の `e` 倍版で分岐した拡大まで道が通る。") 19 ]

end ABC3.Found.GaloisRep
