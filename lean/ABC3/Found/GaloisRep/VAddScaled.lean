/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.VAddC4Velu
import ABC3.Found.GaloisRep.DualTate
import ABC3.Found.GaloisRep.RamifiedValuationBridge
import ABC3.Meta.Claim

/-!
# 第 1426 ブロック —— **`l` 倍ずれた `c₄`・`c₆` の付値**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★これは何か——第 1425 の入力を作る

第 1425（`semistableAt_of_valAdd_c4_c6_scaled`）が要求するのは

    `valAdd p (c₄ E′) = 4·valAdd p (n)`、`valAdd p (c₆ E′) = 6·valAdd p (n)`

である。★`μ_l` 型の核では `K` の水準で
`c₄(veluCurve) = l⁴·c₄(E_{q^l})`・`c₆(veluCurve) = l⁶·c₆(E_{q^l})`
（第 1129-1138、`hlu` なし）が既にあるので、あとは

1. `c₆(E_q)` が単元であること（`c₆ = −1 + 504·s₅` の定数項が `−1`）
2. `c₆` の付値も変数変換で `−6·v(u)` だけ動くこと（`c₄` 版の第 1058 系の相棒）
3. 曲線の等式から `v(c₄) = 4·v(l)`・`v(c₆) = 6·v(l)` を読むこと
4. `Lv` の付値から `p` の付値へ降ろすこと（`e` で割る）

を並べればよい。

| 定理 | 内容 |
|---|---|
| `tateC6_isUnit` / `tateCurveAt_c6_isUnit` | ★`c₆(E_q)` は単元 |
| `vAdd_c6_variableChange` | ★`c₆` の付値は変数変換で `−6·v(u)` 動く |
| `vAdd_c4_scaled_of_eq` / `vAdd_c6_scaled_of_eq` | ★★★★曲線の等式から付値へ |
| `valAdd_scaled_of_vAdd` | ★★★★`Lv` から `p` へ降ろす |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField PowerSeries ABC3.Meta

open scoped Classical

/-! ## ★`c₆(E_q)` は単元 -/

/-- ★`c₆` の定数項は `−1`（`c₆ = −1 + 504·s₅`）。 -/
theorem constantCoeff_tateC6 : PowerSeries.constantCoeff tateCurve.c₆ = -1 := by
  rw [tateCurve_c6]
  have h : PowerSeries.constantCoeff (sigmaSeries 5) = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff]
    exact coeff_zero_sigmaSeries 5
  simp [h]

/-- ★★**`c₆` は単元である**。 -/
theorem tateC6_isUnit : IsUnit tateCurve.c₆ := by
  rw [PowerSeries.isUnit_iff_constantCoeff, constantCoeff_tateC6]
  exact isUnit_one.neg

/-- ★★**`E_q` の `c₆` は単元**（定数項が `−1` だから）。 -/
theorem tateCurveAt_c6_isUnit {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]
    (q : R) (hq : q ∈ I) : IsUnit (tateCurveAt q hq).c₆ := by
  rw [tateCurveAt_c6]
  exact tateC6_isUnit.map (evalAdicHom q hq)

/-! ## ★変数変換と付値 -/

section VC

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

/-- ★**`c₆` の付値は変数変換で `−6·v(u)` だけ動く**——`c₄` 版の相棒。 -/
theorem vAdd_c6_variableChange (W : WeierstrassCurve K) (hc6 : W.c₆ ≠ 0)
    (C : WeierstrassCurve.VariableChange K) (hc6' : (C • W).c₆ ≠ 0) :
    vAdd (tateDvrVal R K) (Units.mk0 ((C • W).c₆) hc6')
      = vAdd (tateDvrVal R K) (Units.mk0 W.c₆ hc6) - 6 * vAdd (tateDvrVal R K) C.u := by
  have hu : (Units.mk0 ((C • W).c₆) hc6') = C.u⁻¹ ^ 6 * Units.mk0 W.c₆ hc6 := by
    refine Units.ext ?_
    show (C • W).c₆ = _
    rw [WeierstrassCurve.variableChange_c₆]
    push_cast
    simp
  rw [hu, vAdd_mul, vAdd_pow, vAdd_inv]
  omega

/-- ☆`Units.mk0` の積の付値。 -/
theorem vAdd_mk0_mul {a b : K} (ha : a ≠ 0) (hb : b ≠ 0) :
    vAdd (tateDvrVal R K) (Units.mk0 (a * b) (mul_ne_zero ha hb))
      = vAdd (tateDvrVal R K) (Units.mk0 a ha) + vAdd (tateDvrVal R K) (Units.mk0 b hb) := by
  have h : Units.mk0 (a * b) (mul_ne_zero ha hb) = Units.mk0 a ha * Units.mk0 b hb :=
    Units.ext rfl
  rw [h, vAdd_mul]

/-- ☆`Units.mk0` の冪の付値。 -/
theorem vAdd_mk0_pow {a : K} (ha : a ≠ 0) (k : ℕ) :
    vAdd (tateDvrVal R K) (Units.mk0 (a ^ k) (pow_ne_zero k ha))
      = (k : ℤ) * vAdd (tateDvrVal R K) (Units.mk0 a ha) := by
  have h : Units.mk0 (a ^ k) (pow_ne_zero k ha) = (Units.mk0 a ha) ^ k :=
    Units.ext (by push_cast; rfl)
  rw [h, vAdd_pow]

/-- ★★★★**曲線の等式から `v(c₄) = 4·v(n)`**（第 1426）。

☆`C₀ • X = V`・`v(u(C₀)) = 0`・`V.c₄ = n⁴·z`（`z` は付値 `0`）から出る。 -/
theorem vAdd_c4_scaled_of_eq (X : WeierstrassCurve K) (hXc4 : X.c₄ ≠ 0)
    (C₀ : WeierstrassCurve.VariableChange K) (V : WeierstrassCurve K)
    (hEq : C₀ • X = V)
    (hu : vAdd (tateDvrVal R K) C₀.u = 0)
    {n z : K} (hn : n ≠ 0) (hz : z ≠ 0) (hz0 : vAdd (tateDvrVal R K) (Units.mk0 z hz) = 0)
    (hV : V.c₄ = n ^ 4 * z) :
    vAdd (tateDvrVal R K) (Units.mk0 X.c₄ hXc4)
      = 4 * vAdd (tateDvrVal R K) (Units.mk0 n hn) := by
  have hVne : V.c₄ ≠ 0 := by rw [hV]; exact mul_ne_zero (pow_ne_zero 4 hn) hz
  have hc4' : (C₀ • X).c₄ ≠ 0 := by rw [hEq]; exact hVne
  have hchg := vAdd_c4_variableChange (R := R) X hXc4 C₀ hc4'
  rw [hu, mul_zero, sub_zero] at hchg
  have hEqU : Units.mk0 ((C₀ • X).c₄) hc4' = Units.mk0 (n ^ 4 * z)
      (mul_ne_zero (pow_ne_zero 4 hn) hz) := by
    refine Units.ext ?_
    show (C₀ • X).c₄ = n ^ 4 * z
    rw [hEq, hV]
  rw [hEqU, vAdd_mk0_mul (pow_ne_zero 4 hn) hz, vAdd_mk0_pow hn 4, hz0, add_zero] at hchg
  exact hchg.symm

/-- ★★★★**曲線の等式から `v(c₆) = 6·v(n)`**（第 1426）。 -/
theorem vAdd_c6_scaled_of_eq (X : WeierstrassCurve K) (hXc6 : X.c₆ ≠ 0)
    (C₀ : WeierstrassCurve.VariableChange K) (V : WeierstrassCurve K)
    (hEq : C₀ • X = V)
    (hu : vAdd (tateDvrVal R K) C₀.u = 0)
    {n z : K} (hn : n ≠ 0) (hz : z ≠ 0) (hz0 : vAdd (tateDvrVal R K) (Units.mk0 z hz) = 0)
    (hV : V.c₆ = n ^ 6 * z) :
    vAdd (tateDvrVal R K) (Units.mk0 X.c₆ hXc6)
      = 6 * vAdd (tateDvrVal R K) (Units.mk0 n hn) := by
  have hVne : V.c₆ ≠ 0 := by rw [hV]; exact mul_ne_zero (pow_ne_zero 6 hn) hz
  have hc6' : (C₀ • X).c₆ ≠ 0 := by rw [hEq]; exact hVne
  have hchg := vAdd_c6_variableChange (R := R) X hXc6 C₀ hc6'
  rw [hu, mul_zero, sub_zero] at hchg
  have hEqU : Units.mk0 ((C₀ • X).c₆) hc6' = Units.mk0 (n ^ 6 * z)
      (mul_ne_zero (pow_ne_zero 6 hn) hz) := by
    refine Units.ext ?_
    show (C₀ • X).c₆ = n ^ 6 * z
    rw [hEq, hV]
  rw [hEqU, vAdd_mk0_mul (pow_ne_zero 6 hn) hz, vAdd_mk0_pow hn 6, hz0, add_zero] at hchg
  exact hchg.symm

end VC

/-! ## ★★★★`Lv` の付値から `p` の付値へ -/

section Descend

variable {L : Type} [Field L] [NumberField L]
  {Lv : Type} [Field Lv] [Algebra L Lv]
  {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [Algebra R Lv] [IsFractionRing R Lv]

/-- ★★★★**`v_{Lv}(x) = k·v_{Lv}(n)` なら `v_p(x) = k·v_p(n)`**（第 1426）。

☆`hpe` から両辺が `e` 倍になるので、`e ≥ 1` で割ればよい。 -/
theorem valAdd_scaled_of_vAdd {e : ℕ} (he : 1 ≤ e) (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (x n : L) (hx : x ≠ 0) (hn : n ≠ 0) (k : ℕ)
    (hx' : algebraMap L Lv x ≠ 0) (hn' : algebraMap L Lv n ≠ 0)
    (h : vAdd (tateDvrVal R Lv) (Units.mk0 (algebraMap L Lv x) hx')
       = (k : ℤ) * vAdd (tateDvrVal R Lv) (Units.mk0 (algebraMap L Lv n) hn')) :
    valAdd p (Units.mk0 x hx) = (k : ℤ) * valAdd p (Units.mk0 n hn) := by
  rw [vAdd_algebraMap_eq_mul_valAdd (R := R) p hpe x hx hx',
    vAdd_algebraMap_eq_mul_valAdd (R := R) p hpe n hn hn'] at h
  have hepos : (0 : ℤ) < (e : ℤ) := by exact_mod_cast he
  have hcan : (e : ℤ) * valAdd p (Units.mk0 x hx)
      = (e : ℤ) * ((k : ℤ) * valAdd p (Units.mk0 n hn)) := by linarith [h]
  exact mul_left_cancel₀ (by omega) hcan

end Descend

/-! ## ★出典の紐付け(`.src`) -/

def tateC6_isUnit.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 曲線の c₆ は単元——定数項が −1)",
    sectionId := "genell-def-3-3" }

def vAdd_c6_variableChange.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(c₆ の付値は変数変換で −6·v(u) 動く。★無条件)",
    sectionId := "genell-lemma-3-5" }

def vAdd_c4_scaled_of_eq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(曲線の等式から v(c₄) = 4·v(n)。★無条件)",
    sectionId := "genell-lemma-3-5" }

def vAdd_c6_scaled_of_eq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(曲線の等式から v(c₆) = 6·v(n)。★無条件)",
    sectionId := "genell-lemma-3-5" }

def valAdd_scaled_of_vAdd.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Lv の付値から p の付値へ——分岐版。★無条件)",
    sectionId := "genell-lemma-3-5" }

def valAdd_scaled_of_vAdd.needs : List ProofObligation :=
  [ .citation "[ABC3]" "vAdd_algebraMap_eq_mul_valAdd(第 1370、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.vAdd_algebraMap_eq_mul_valAdd") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1426）**——第 1425 の入力" ++
       "（`valAdd p (c₄) = 4·valAdd p (n)`・`valAdd p (c₆) = 6·valAdd p (n)`）を" ++
       "作るための 5 本を並べた。" ++
       "☆`μ_l` 型の核では `K` の水準で `c₄(veluCurve) = l⁴·c₄(E_{q^l})` が既にあるので" ++
       "（第 1129-1138、`hlu` なし）、あとは変数変換と `e` 倍の付値を通すだけである。") 17 ]

end ABC3.Found.GaloisRep
