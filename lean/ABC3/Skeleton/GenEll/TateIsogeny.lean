/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.DegInfLocal
import ABC3.Found.GaloisRep.Lemma35Concrete
import ABC3.Found.GaloisRep.TateParamJ
import ABC3.Found.GaloisRep.TateDXNeZero
import ABC3.Found.GaloisRep.TateVeluMu
import ABC3.Found.GaloisRep.TateSetupDvr
import ABC3.Found.GaloisRep.AdicCompleteIntegers
import ABC3.Found.GaloisRep.DegInfTateParam
import ABC3.Found.GenEll.MuPrimitiveRoot
import ABC3.Found.GenEll.CyclotomicUnits
import ABC3.Found.GaloisRep.TateModelPoint
import ABC3.Found.GaloisRep.BadPrimeData
import ABC3.Found.GaloisRep.CompletionValuationBridge
import ABC3.Found.GaloisRep.TateMuInvolution
import ABC3.Found.GaloisRep.VeluTateDelta
import Mathlib.NumberTheory.NumberField.Completion.FinitePlace
import ABC3.Found.GaloisRep.VeluMuSum
import ABC3.Found.GenEll.JScale
import ABC3.Meta.Claim
import ABC3.Skeleton.GenEll.TateODE

/-!
# `Lemma 3.2, (ii)` の曲線の水準 —— **`E_q/μ_l` は母数 `q^l` の Tate 曲線**（`Skeleton`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★★★これは何か

`Found/GaloisRep/DegInfLocal.lean` の `minDeltaExp_eq_mul_of_tateModel`（`§9-1153`、第 716）は

> `E` の Tate モデルの母数が `q`、`E′` の母数が `q^l` なら
> `v_p(Δ_min(E′)) = l·v_p(Δ_min(E))`

を**無条件で**証明している。★したがって `Lemma 3.5`（および `Lemma 3.7` の第 3 の主張）に
残っている入力は、**本ファイルが固定する 1 つの命題だけ**である:

> **`E′ = E/μ_l` は母数 `q^l` の Tate モデルを持つ。**

## ★★★なぜ `Skeleton` に置くか

`Lemma 3.2` は本プロジェクトで**完了扱い**の項目だが、`Found/GenEll/Lemma32.lean` の
`lemma_3_2` が取っているのは**群論の核**

    `(Lˣ/q^ℤ) ⧸ (x ↦ x^l の核の像) ≃* Lˣ/(q^l)^ℤ`

であって、**曲線の水準の主張**（`E_q/μ_l` の Weierstrass モデルが `tateCurveAt (q^l)` と
変数変換で一致すること）ではない。★★この差が `Lemma 3.5` の最後の穴である。

## ★★★★★証明の筋（`q`-展開の恒等式）

`tateCurveAt q` は `⟨1, 0, 0, tateA4(q), tateA6(q)⟩`（`tateA4 = −5·σ₃`）である。
`veluCurve` は `a₁, a₂, a₃` を変えないので、**変数変換は要らず係数の一致を見ればよい**:

    `tateA4(q^l) = tateA4(q) − 5·v`,
    `tateA6(q^l) = tateA6(q) − b₂·v − 7·w`   （`b₂ = a₁² + 4a₂ = 1`）

ここで `v = Σ_{ζ ∈ μ_l∖{1}} g^x(P_ζ)`、`w = Σ (u_{P_ζ}/2 + g^x(P_ζ)·x_{P_ζ})` であり、
`P_ζ` の座標は `tateXpair ζ (qζ⁻¹) q`・`tateYpair ζ (qζ⁻¹) q`（`Found/GaloisRep/TateOrigin.lean`）。

★★`ζ` について足すと `ζ` の指数が `l` の倍数の項だけが残る——これが
`σₖ(q) → σₖ(q^l)` を生む機構である。

☆本プロジェクトの Tate 冪級数の在庫（`Found/GaloisRep/Tate*.lean` が 97 ファイル）で
計算できる。**新しい数学は要らない。**

## ★次に何をすれば終わりか

1. `μ_l` の点の座標を `q`-冪級数として書く（`tateXpair`・`tateYpair` の特殊化）
2. `Σ_{ζ ∈ μ_l}` を取ると `ζ` の指数が `l` の倍数の項だけ残ることを示す
3. 上の 2 式（`a₄`・`a₆`）を確かめる
4. `minDeltaExp_eq_mul_of_tateModel`（第 716）へ渡す
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Found.GaloisRep WeierstrassCurve IsDedekindDomain NumberField ABC3.Found.GenEll
open scoped Classical

/-- **[GenEll] `Lemma 3.2, (ii)` の曲線の水準**——`E_q/μ_l` は母数 `q^l` の Tate 曲線。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★★★★★★★☆**これが `Lemma 3.5` に残る唯一の穴である。**
`minDeltaExp_eq_mul_of_tateModel`（`§9-1153`、第 716）がこの結論を消費して
`v_p(Δ_min(E′)) = l·v_p(Δ_min(E))` を出す。

☆仮説 `hquot` は「`E′` は `E` を位数 `l` の巡回部分群で割ったもので、
その部分群は各悪い素点で `μ_l` に対応する」という `Lemma 3.5` の設定
（原文の *global rank one subgroup*）を型にしたものである。 -/
theorem tateModel_of_quot_mu {R : Type} [CommRing R] [IsDomain R]
    [IsDiscreteValuationRing R] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]
    (W W' : WeierstrassCurve K) [W.IsElliptic] [W.IsMinimal R]
    [W'.IsElliptic] [W'.IsMinimal R]
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (l : ℕ) (hl : 0 < l)
    (hql : q ^ l ∈ IsLocalRing.maximalIdeal R)
    (D : VariableChange R) (hD : D • integralModel R W = tateCurveAt q hq)
    (hsplit : W'.HasSplitMultiplicativeReduction R)
    (hparam : tateParamR W' hsplit = q ^ l) :
    ∃ D' : VariableChange R,
      D' • integralModel R W' = tateCurveAt (q ^ l) hql := by
  obtain ⟨hq', C, hne, hC⟩ := tateParamR_spec W' hsplit
  refine ⟨C, ?_⟩
  rw [hC]
  congr 1

def tateModel_of_quot_mu.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tateParamR_spec(Tate モデルの存在、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateParamR_spec") 1 ]

def tateParam_quot_mu.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(E′ の Tate 母数は q^l である)",
    sectionId := "genell-lemma-3-2" }

def tateParam_quot_mu.needs : List ProofObligation :=
  [ .citation "[ABC3]" "j_velu_tate_mu_map(j(E_q/μ_l) = j(E_{q^l})、K の中、第 881、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.j_velu_tate_mu_map") 1,
    .citation "[ABC3]" "tateParamR_eq_of_j_tateCurveAt(j が E_{q₀} の j なら母数は q₀、第 882、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateParamR_eq_of_j_tateCurveAt") 1,
    .implicitStep
      ("★★★到達点(2026-08-31、第 883): hquot はもはや True ではなく" ++
       "「W′ の j が Vélu の商の j に等しい」という**型のついた仮説**であり、" ++
       "定理自体は**証明済み**である。" ++
       "☆残るのは「H が各悪い素点で μ_l に対応する」から hquot を導く段、" ++
       "すなわち**大域の Vélu の商と局所の Vélu の商を繋ぐ配管**だけである") 5 ]

def tateModel_of_quot_mu.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(曲線の水準——E_q/μ_l は母数 q^l の Tate 曲線)",
    sectionId := "genell-lemma-3-2" }

/-! ## ★★★★★★★★★★★★★★★★★★★★訂正後の目標（`c₄`・`c₆`） -/

open Finset in
/-- **[GenEll] 葉 1 の訂正後の目標 (1)**——`c₄` は `l⁴` 倍。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★★★第 834-835 の数値測定により、**係数 `a₄` の恒等式は偽**であり、
`c₄` で書いたこの形が真である（`l = 5, 7` で `q^21` まで全係数一致）。

☆`c₄(veluCurve W v w) = c₄ W + 240v`（`Found/GenEll/Velu.lean` の `veluCurve_c₄`）なので
左辺は Vélu の商の `c₄` そのものである。 -/
theorem c4_velu_tate {R : Type} [CommRing R] [IsDomain R] [CharZero R] {I : Ideal R}
    [IsAdicComplete I R] {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hlu : IsUnit ((l : R))) (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (q : R) (hq : q ∈ I) (hql : q ^ l ∈ I) (h2 : (2 : R) ≠ 0)
    (hDX : ∀ i ∈ (range l).erase 0,
      tateDXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq ≠ 0) :
    (tateCurveAt q hq).c₄
        + 240 * (∑ i ∈ (range l).erase 0,
            veluV2 (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
      = (l : R) ^ 4 * (tateCurveAt (q ^ l) hql).c₄ := by
  have hζl : ζ ^ l = 1 := hζ.pow_eq_one
  have hsum1 := sum_veluV2_eq_sum_tateDYpair hl.pos hζl hu q hq hDX
  have hstep : ∑ i ∈ (range l).erase 0, tateD2Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
      = 2 * (∑ i ∈ (range l).erase 0, tateDYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        + ∑ i ∈ (range l).erase 0, tateDXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq := by
    rw [Finset.mul_sum, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl (fun i hi => ?_)
    exact tateD2Xpair_eq _ _ q hq (hu i hi)
      (isUnit_one_sub (I := I) (Ideal.mul_mem_right _ _ hq))
  have hzero := sum_mu_dxpair_zero hl hζ hu q hq h2
  have hd2 := sum_mu_d2xpair hl hζ hu q hq hql
  rw [hstep, hzero, add_zero] at hd2
  rw [hsum1, tateCurveAt_c4_eq, tateCurveAt_c4_eq]
  linear_combination hd2

open Finset in
/-- **[GenEll] 葉 1 の訂正後の目標 (2)**——`c₆` は `l⁶` 倍。

★★★★**2026-08-31 の再訂正（第 867）**——以前の statement
`2c₆ + 1008∑v + 3024∑(u+2vX) = −2l⁶c₆(q^l)` は**偽**である（数値で落ちた）。
正しいのは

    `c₆ + 504∑v + 3024∑(u + 2vX) = l⁶·c₆(q^l)`

である（`l = 5, 7` で数値確認）。☆Vélu の和は代表元の集合 `R` の上で取るので
`v_Vélu = ∑_{i≠0} g^x`、`w_Vélu = (1/2)∑_{i≠0} u + ∑_{i≠0} x·g^x` となる。
☆`g^x_Q + g^x_{-Q} = v_Q` が `∑_{i≠0} veluV2 = v_Vélu` の理由である。 -/
theorem c6_velu_tate {R : Type} [CommRing R] [IsDomain R] [CharZero R] {I : Ideal R}
    [IsAdicComplete I R] {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (q : R) (hq : q ∈ I) (hql : q ^ l ∈ I) (h2 : (2 : R) ≠ 0)
    (hDX : ∀ i ∈ (range l).erase 0,
      tateDXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq ≠ 0) :
    (tateCurveAt q hq).c₆
        + 504 * (∑ i ∈ (range l).erase 0,
            veluV2 (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
        + 3024 * (∑ i ∈ (range l).erase 0,
            (veluU (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              + 2 * (veluV2 (tateCurveAt q hq)
                      (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                      (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                    * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)))
      = (l : R) ^ 6 * (tateCurveAt (q ^ l) hql).c₆ := by
  have hζl : ζ ^ l = 1 := hζ.pow_eq_one
  have hawi : ∀ i ∈ (range l).erase 0,
      (ζ ^ i) * (q * (ζ ^ i) ^ (l - 1)) = q := by
    intro i hi
    have hpow : (ζ ^ i) * (ζ ^ i) ^ (l - 1) = 1 := by
      rw [← pow_succ']
      rw [Nat.sub_add_cancel hl.pos, ← pow_mul, mul_comm, pow_mul, hζl, one_pow]
    calc (ζ ^ i) * (q * (ζ ^ i) ^ (l - 1)) = q * ((ζ ^ i) * (ζ ^ i) ^ (l - 1)) := by ring
      _ = q := by rw [hpow, mul_one]
  have hwu : ∀ i : ℕ, IsUnit (1 - q * (ζ ^ i) ^ (l - 1)) := fun i =>
    isUnit_one_sub (I := I) (Ideal.mul_mem_right _ _ hq)
  -- 項ごとの恒等式
  have hterm : ∀ i ∈ (range l).erase 0,
      12 * (veluU (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            + 2 * (veluV2 (tateCurveAt q hq)
                    (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                    (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                  * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
        = tateD4Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
          - tateD3Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
          - tateD2Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
          + tateDXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq := by
    intro i hi
    have hd4 := tate_d4x (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq (hawi i hi) (hu i hi)
      (hwu i) (hDX i hi)
    have hd3 := tate_d3x (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq (hawi i hi) (hu i hi)
      (hwu i) (hDX i hi)
    have hd2 := tate_d2x (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq (hawi i hi) (hu i hi)
      (hwu i) (hDX i hi)
    have hdx := tateDXpair_eq (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq (hu i hi) (hwu i)
    rw [veluU_tateCurveAt, veluV2_tateCurveAt, hd4, hd3, hd2, hdx]
    ring
  -- 和に直す
  have hsum12 : 12 * (∑ i ∈ (range l).erase 0,
      (veluU (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            + 2 * (veluV2 (tateCurveAt q hq)
                    (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                    (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                  * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)))
      = (∑ i ∈ (range l).erase 0, tateD4Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        - (∑ i ∈ (range l).erase 0, tateD3Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        - (∑ i ∈ (range l).erase 0, tateD2Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        + (∑ i ∈ (range l).erase 0, tateDXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq) := by
    rw [Finset.mul_sum, Finset.sum_congr rfl hterm, Finset.sum_add_distrib,
      Finset.sum_sub_distrib, Finset.sum_sub_distrib]
  -- ∑DX = 0、∑D³X = 0
  have hz1 := sum_mu_dxpair_zero hl hζ hu q hq h2
  have hz3 := sum_mu_d3xpair_zero hl hζ hu q hq h2
  -- ∑v = ∑DY と ∑D²X = 2∑DY + ∑DX
  have hsum1 := sum_veluV2_eq_sum_tateDYpair hl.pos hζl hu q hq hDX
  have hstep : ∑ i ∈ (range l).erase 0, tateD2Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
      = 2 * (∑ i ∈ (range l).erase 0, tateDYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        + ∑ i ∈ (range l).erase 0, tateDXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq := by
    rw [Finset.mul_sum, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl (fun i hi => ?_)
    exact tateD2Xpair_eq _ _ q hq (hu i hi) (hwu i)
  have hd4sum := sum_mu_d4xpair hl hζ hu q hq hql
  rw [hz1, add_zero] at hstep
  rw [hz1, hz3, sub_zero, add_zero] at hsum12
  rw [hsum1, tateCurveAt_c6_eq, tateCurveAt_c6_eq]
  linear_combination (252 : R) * hsum12 + hd4sum - (252 : R) * hstep

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★**Vélu の商の `j` は `j(E_{q^l})` に等しい**。

★葉 1 の目標 (1)(2) を `j` の言葉にしたもの。`v`・`w` は
`veluVFull`・`veluWFull`（`w = ∑(u/2 + g^x·x)`、すなわち `2w = ∑(u + 2vx)`）の形で与える。

☆`c₄ + 240v = l⁴c₄(q^l)` と `c₆ + 504v + 6048w = l⁶c₆(q^l)` から
`x ↦ l²x + r` の変換で `j` が保たれる（`j_eq_of_c4_c6_scale_pos`）。 -/
theorem j_velu_tate_mu {R : Type} [CommRing R] [IsDomain R] [CharZero R] {I : Ideal R}
    [IsAdicComplete I R] {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hlu : IsUnit ((l : R)))
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (q : R) (hq : q ∈ I) (hql : q ^ l ∈ I) (h2 : (2 : R) ≠ 0)
    (hDX : ∀ i ∈ (range l).erase 0,
      tateDXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq ≠ 0)
    (v w : R)
    (hv : v = ∑ i ∈ (range l).erase 0,
      veluV2 (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
    (hw : 2 * w = ∑ i ∈ (range l).erase 0,
      (veluU (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
          (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        + 2 * (veluV2 (tateCurveAt q hq)
                (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)))
    [(veluCurve (tateCurveAt q hq) v w).IsElliptic]
    [(tateCurveAt (q ^ l) hql).IsElliptic] :
    (veluCurve (tateCurveAt q hq) v w).j = (tateCurveAt (q ^ l) hql).j := by
  have h4 := c4_velu_tate hl hζ hlu hu q hq hql h2 hDX
  have h6 := c6_velu_tate hl hζ hu q hq hql h2 hDX
  refine j_velu_tate_eq q hq l hql v w ?_ ?_
  · rw [hv]
    exact h4
  · rw [hv]
    linear_combination h6 + (3024 : R) * hw

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★**`K` に落とした形**の `j_velu_tate_mu`。

★★★★**2026-08-31 の測定（第 881）**——`j_velu_tate_mu` は `R` の上の曲線に
`[IsElliptic]` を仮定しているが、Tate 曲線は `Δ = q·(単元)` なので
**`R` の上では `IsElliptic` ではない**——つまりその形は空虚である。
★本節点が**実際に使える形**である。 -/
theorem j_velu_tate_mu_map {R : Type} [CommRing R] [IsDomain R] [CharZero R] {I : Ideal R}
    [IsAdicComplete I R] {K : Type} [Field K] [CharZero K] [Algebra R K]
    {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hlu : IsUnit ((l : R)))
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (q : R) (hq : q ∈ I) (hql : q ^ l ∈ I) (h2 : (2 : R) ≠ 0)
    (hDX : ∀ i ∈ (range l).erase 0,
      tateDXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq ≠ 0)
    (v w : R)
    (hv : v = ∑ i ∈ (range l).erase 0,
      veluV2 (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
    (hw : 2 * w = ∑ i ∈ (range l).erase 0,
      (veluU (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
          (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        + 2 * (veluV2 (tateCurveAt q hq)
                (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)))
    [((veluCurve (tateCurveAt q hq) v w).map (algebraMap R K)).IsElliptic]
    [((tateCurveAt (q ^ l) hql).map (algebraMap R K)).IsElliptic] :
    ((veluCurve (tateCurveAt q hq) v w).map (algebraMap R K)).j
      = ((tateCurveAt (q ^ l) hql).map (algebraMap R K)).j := by
  have h4 := c4_velu_tate hl hζ hlu hu q hq hql h2 hDX
  have h6 := c6_velu_tate hl hζ hu q hq hql h2 hDX
  refine j_velu_tate_eq_map q hq l hql v w ?_ ?_
  · rw [hv]
    exact h4
  · rw [hv]
    linear_combination h6 + (3024 : R) * hw

def j_velu_tate_mu_map.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(K に落とした形の Vélu の商の j)",
    sectionId := "genell-lemma-3-2" }

def j_velu_tate_mu_map.needs : List ProofObligation :=
  [ .citation "[ABC3]" "c4_velu_tate(第 853、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.c4_velu_tate") 1,
    .citation "[ABC3]" "c6_velu_tate(第 867、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.c6_velu_tate") 1,
    .citation "[ABC3]" "j_velu_tate_eq_map(第 881、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.j_velu_tate_eq_map") 1 ]

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★**[GenEll] 葉 1 の残る中身**——
`E′` の Tate 母数は `q^l` である。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★★★**2026-08-31 の空欄の埋め戻し（第 883）**——第 873 では
`hquot : True` という**空欄**だった。★本ブロックでそれを

    `hquot : W′.j = j((E_q/μ_l) ⊗ K)`

という**型のついた仮説**に置き換え、定理を**証明した**。
☆これが原文の「H が各悪い素点で μ_l に対応する」
（global rank one subgroup）の、曲線の水準での内容である。

☆道は 2 段だけである:

1. `j_velu_tate_mu_map`（第 881）で `j(E_q/μ_l) = j(E_{q^l})`
2. `tateParamR_eq_of_j_tateCurveAt`（第 882）で `q_{E′} = q^l` -/
theorem tateParam_quot_mu {R : Type} [CommRing R] [IsDomain R] [CharZero R]
    [IsDiscreteValuationRing R] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [CharZero K] [Algebra R K] [IsFractionRing R K]
    {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hlu : IsUnit ((l : R)))
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R)
    (hql : q ^ l ∈ IsLocalRing.maximalIdeal R) (h2 : (2 : R) ≠ 0)
    (hDX : ∀ i ∈ (range l).erase 0,
      tateDXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq ≠ 0)
    (v w : R)
    (hv : v = ∑ i ∈ (range l).erase 0,
      veluV2 (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
    (hw : 2 * w = ∑ i ∈ (range l).erase 0,
      (veluU (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
          (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        + 2 * (veluV2 (tateCurveAt q hq)
                (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)))
    (W' : WeierstrassCurve K) [W'.IsElliptic] [W'.IsMinimal R]
    (hsplit : W'.HasSplitMultiplicativeReduction R)
    [((veluCurve (tateCurveAt q hq) v w).map (algebraMap R K)).IsElliptic]
    [((tateCurveAt (q ^ l) hql).map (algebraMap R K)).IsElliptic]
    (hquot : W'.j = ((veluCurve (tateCurveAt q hq) v w).map (algebraMap R K)).j) :
    tateParamR W' hsplit = q ^ l := by
  refine tateParamR_eq_of_j_tateCurveAt W' hsplit (q ^ l) hql ?_
  rw [hquot]
  exact j_velu_tate_mu_map hl hζ hlu hu q hq hql h2 hDX v w hv hw

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★**[GenEll] 葉 1 —— `hquot` を
`⟨Φ(ζ)⟩` で書いた形**。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★★★**2026-08-31（第 891）**——第 883 で置いた `hquot`（`j` の一致）をさらに
**`W′` が `E_q` の `⟨Φ(ζ)⟩` による Vélu の商である**という形に退けた。
☆これが原文の「`H` は乗法還元の素点で `μ_l` に対応する」の、
曲線の水準での直訳である。

☆道は 2 段だけである:

1. `veluQuotientFull_tate_mu`（第 890）で商を `veluCurve (E_q) v w` に直す
2. `tateParam_quot_mu`（第 883） -/
theorem tateParam_quot_velu {R : Type} [CommRing R] [IsDomain R] [CharZero R]
    [IsDiscreteValuationRing R] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [CharZero K] [Algebra R K] [IsFractionRing R K]
    (S : TateSetup R (IsLocalRing.maximalIdeal R) K)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hlu : IsUnit ((l : R)))
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (uζ : Kˣ) (hζu : algebraMap R K ζ = (uζ : K)) (hζl : uζ ^ l = 1)
    (hord : ∀ n : ℕ, 0 < n → n < l → uζ ^ n ≠ 1)
    (hql : S.q ^ l ∈ IsLocalRing.maximalIdeal R)
    (h2 : (2 : R) ≠ 0) (h2K : (2 : K) ≠ 0) (hodd : l ≠ 2)
    (v w : R)
    (hv : v = ∑ i ∈ (range l).erase 0,
      veluV2 (tateCurveAt S.q S.hq)
        (tateXpair (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
        (tateYpair (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq))
    (hw : 2 * w = ∑ i ∈ (range l).erase 0,
      (veluU (tateCurveAt S.q S.hq)
          (tateXpair (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
          (tateYpair (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
        + 2 * (veluV2 (tateCurveAt S.q S.hq)
                (tateXpair (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
                (tateYpair (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
              * tateXpair (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)))
    (W' : WeierstrassCurve K) [W'.IsElliptic] [W'.IsMinimal R]
    (hsplit : W'.HasSplitMultiplicativeReduction R)
    [((veluCurve (tateCurveAt S.q S.hq) v w).map (algebraMap R K)).IsElliptic]
    [((tateCurveAt (S.q ^ l) hql).map (algebraMap R K)).IsElliptic]
    (hW' : W' = veluQuotientFull ((tateCurveAt S.q S.hq).map (algebraMap R K))
      (((range l).erase 0).image
        (fun k : ℕ => pointCoords (k • tatePhi S hΔ (QuotientGroup.mk uζ))))) :
    tateParamR W' hsplit = S.q ^ l := by
  -- ★★★★**2026-08-31（第 895）**——侧条件 `hDX` はここでは**定理**である
  have hDX : ∀ i ∈ (range l).erase 0,
      tateDXpair (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq ≠ 0 := by
    intro i hi
    exact tateDXpair_ne_zero_of_mu S hΔ Φ hΦ hl hodd ζ uζ hζu hζl hord i
      (Finset.mem_erase.1 hi).1 (Finset.mem_range.1 (Finset.mem_erase.1 hi).2) (hu i hi)
  have hquot := veluQuotientFull_tate_mu S hΔ Φ hΦ hl.pos ζ uζ hζu hζl hord hu v w h2K hv hw
  have hWW : W' = (veluCurve (tateCurveAt S.q S.hq) v w).map (algebraMap R K) :=
    hW'.trans hquot
  subst hWW
  exact tateParam_quot_mu hl hζ hlu hu S.q S.hq hql h2 hDX v w hv hw _ hsplit rfl

def tateParam_quot_velu.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(hquot を ⟨Φ(ζ)⟩ の Vélu の商で書いた形)",
    sectionId := "genell-lemma-3-2" }

def tateParam_quot_velu.needs : List ProofObligation :=
  [ .citation "[ABC3]" "veluQuotientFull_tate_mu(⟨Φ(ζ)⟩ の Vélu の商、第 890、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.veluQuotientFull_tate_mu") 1,
    .citation "[ABC3]" "tateParam_quot_mu(j から母数へ、第 883、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.tateParam_quot_mu") 1,
    .citation "[ABC3]" "tateDXpair_ne_zero_of_mu(侧条件 hDX は定理、第 894、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateDXpair_ne_zero_of_mu") 1,
    .implicitStep
      ("☆残るのは大域の `E′ = E/H` を各悪い素点で完備化に落とし、" ++
       "H の像が `⟨Φ(ζ)⟩` になることを言う段（Lemma 3.2, (i) の帰結）だけである") 4 ]

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**[GenEll] 葉 1 ——
`Φ` を受けない形**。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★★★**2026-08-31（第 900）**——第 899 の `dvrTatePhiAddEquiv` により、
完備な DVR なら Tate 一意化 `Φ` は**無条件に存在する**。
☆したがって `tateParam_quot_velu`（第 891）から `Φ`・`hΦ` を**落とせる**。

★残る仮説はすべて**幾何のデータ**である——Tate 母数 `q`、1 の原始 `l` 乗根 `ζ`、
そして「`W′` は `E_q` の `⟨Φ(ζ)⟩` による Vélu の商である」。 -/
theorem tateParam_quot_velu_dvr {R : Type} [CommRing R] [IsDomain R] [CharZero R]
    [IsDiscreteValuationRing R] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [CharZero K] [Algebra R K] [IsFractionRing R K]
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (hq0 : q ≠ 0)
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hlu : IsUnit ((l : R)))
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (uζ : Kˣ) (hζu : algebraMap R K ζ = (uζ : K)) (hζl : uζ ^ l = 1)
    (hord : ∀ n : ℕ, 0 < n → n < l → uζ ^ n ≠ 1)
    (hql : q ^ l ∈ IsLocalRing.maximalIdeal R)
    (h2 : (2 : R) ≠ 0) (h2K : (2 : K) ≠ 0) (hodd : l ≠ 2)
    (v w : R)
    (hv : v = ∑ i ∈ (range l).erase 0,
      veluV2 (tateCurveAt q hq)
        (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
    (hw : 2 * w = ∑ i ∈ (range l).erase 0,
      (veluU (tateCurveAt q hq)
          (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
          (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        + 2 * (veluV2 (tateCurveAt q hq)
                (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)))
    (W' : WeierstrassCurve K) [W'.IsElliptic] [W'.IsMinimal R]
    (hsplit : W'.HasSplitMultiplicativeReduction R)
    [((veluCurve (tateCurveAt q hq) v w).map (algebraMap R K)).IsElliptic]
    [((tateCurveAt (q ^ l) hql).map (algebraMap R K)).IsElliptic]
    (hW' : W' = veluQuotientFull ((tateCurveAt q hq).map (algebraMap R K))
      (((range l).erase 0).image
        (fun k : ℕ => pointCoords
          (k • tatePhi (mkTateSetup (K := K) q hq hq0) hΔ (QuotientGroup.mk uζ))))) :
    tateParamR W' hsplit = q ^ l := by
  -- ★(mkTateSetup q hq hq0).q は q と定義上等しいが、
  --   インスタンス探索は構文的なので手で渡す
  haveI i1 : ((veluCurve (tateCurveAt (mkTateSetup (K := K) q hq hq0).q
      (mkTateSetup (K := K) q hq hq0).hq) v w).map (algebraMap R K)).IsElliptic :=
    inferInstanceAs (((veluCurve (tateCurveAt q hq) v w).map (algebraMap R K)).IsElliptic)
  haveI i2 : ((tateCurveAt ((mkTateSetup (K := K) q hq hq0).q ^ l) hql).map
      (algebraMap R K)).IsElliptic :=
    inferInstanceAs (((tateCurveAt (q ^ l) hql).map (algebraMap R K)).IsElliptic)
  exact tateParam_quot_velu (mkTateSetup q hq hq0) hΔ
    (dvrTatePhiAddEquiv q hq hq0 hΔ) (fun _ => rfl)
    hl hζ hlu hu uζ hζu hζl hord hql h2 h2K hodd v w hv hw W' hsplit hW'

def tateParam_quot_velu_dvr.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(Φ を受けない形——完備 DVR なら Φ は無条件に存在する)",
    sectionId := "genell-lemma-3-2" }

def tateParam_quot_velu_dvr.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tateParam_quot_velu(第 891、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.tateParam_quot_velu") 1,
    .citation "[ABC3]" "dvrTatePhiAddEquiv(完備 DVR の Tate 一意化、第 899、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.dvrTatePhiAddEquiv") 1,
    .implicitStep
      ("☆残るのは 2 つだけである: (1) E ⊗ Lv が分裂乗法還元をもつこと、" ++
       "(2) H の像が ⟨Φ(ζ)⟩ になること（Lemma 3.2, (i) の対偶）") 4 ]

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] 葉 1 の終点——数体の素点での `Δ_min` の関係**。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★★★**2026-08-31（第 904）**——局所の連鎖を**数体の素点の完備化に当てた**形である。

    `R := p.adicCompletionIntegers L`   `Lv := p.adicCompletion L`

☆道は 2 段だけ:

1. `tateParam_quot_velu_dvr`（第 900）で `q_{E′} = q_E^l`
2. `minDeltaExp_eq_mul_of_tateParamR`（第 892）で `v_p(Δ_min(E′)) = l·v_p(Δ_min(E))`

★これが `lemma_3_5_velu_bad_delta`（第 903）の入力そのものである。 -/
theorem minDeltaExp_eq_mul_of_veluMu {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E E' : WeierstrassCurve L)
    [(E.baseChange (p.adicCompletion L)).IsElliptic]
    [(E.baseChange (p.adicCompletion L)).IsMinimal (p.adicCompletionIntegers L)]
    [(E'.baseChange (p.adicCompletion L)).IsElliptic]
    [(E'.baseChange (p.adicCompletion L)).IsMinimal (p.adicCompletionIntegers L)]
    (h : (E.baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    (h' : (E'.baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation (p.adicCompletion L)
        (IsDiscreteValuationRing.maximalIdeal (p.adicCompletionIntegers L)))
        (algebraMap L (p.adicCompletion L) x) = (HeightOneSpectrum.valuation L p) x)
    (C : VariableChange L) (hC : IsMinimal (primeSubring p) (C • E))
    (hc4ne : (C • E).c₄ ≠ 0) (hc4 : valAdd p (Units.mk0 ((C • E).c₄) hc4ne) = 0)
    (C' : VariableChange L) (hC' : IsMinimal (primeSubring p) (C' • E'))
    (hc4ne' : (C' • E').c₄ ≠ 0) (hc4' : valAdd p (Units.mk0 ((C' • E').c₄) hc4ne') = 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hΔ : ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).toAffine.Δ ≠ 0)
    {ζ : (p.adicCompletionIntegers L)} (hζ : IsPrimitiveRoot ζ l)
    (hlu : IsUnit ((l : (p.adicCompletionIntegers L))))
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (uζ : (p.adicCompletion L)ˣ)
    (hζu : algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L) ζ = (uζ : _))
    (hζl : uζ ^ l = 1)
    (hord : ∀ n : ℕ, 0 < n → n < l → uζ ^ n ≠ 1)
    (hql : (tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l
      ∈ IsLocalRing.maximalIdeal (p.adicCompletionIntegers L))
    (h2 : (2 : (p.adicCompletionIntegers L)) ≠ 0)
    (h2K : (2 : (p.adicCompletion L)) ≠ 0)
    (v w : (p.adicCompletionIntegers L))
    (hv : v = ∑ i ∈ (range l).erase 0,
      veluV2 (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
          (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
        (tateXpair (ζ ^ i)
          ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
          (tateParamR (E.baseChange (p.adicCompletion L)) h)
          (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
        (tateYpair (ζ ^ i)
          ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
          (tateParamR (E.baseChange (p.adicCompletion L)) h)
          (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)))
    (hw : 2 * w = ∑ i ∈ (range l).erase 0,
      (veluU (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
          (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
        (tateXpair (ζ ^ i)
          ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
          (tateParamR (E.baseChange (p.adicCompletion L)) h)
          (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
        (tateYpair (ζ ^ i)
          ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
          (tateParamR (E.baseChange (p.adicCompletion L)) h)
          (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
        + 2 * (veluV2 (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
                (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
              (tateXpair (ζ ^ i)
                ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                (tateParamR (E.baseChange (p.adicCompletion L)) h)
                (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
              (tateYpair (ζ ^ i)
                ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                (tateParamR (E.baseChange (p.adicCompletion L)) h)
                (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
            * tateXpair (ζ ^ i)
                ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                (tateParamR (E.baseChange (p.adicCompletion L)) h)
                (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))))
    [((veluCurve (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)) v w).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).IsElliptic]
    [((tateCurveAt ((tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) hql).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).IsElliptic]
    (hW' : E'.baseChange (p.adicCompletion L)
      = veluQuotientFull ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
          (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
            (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
        (((range l).erase 0).image
          (fun k : ℕ => pointCoords (k • tatePhi
            (mkTateSetup (K := p.adicCompletion L)
              (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h))
            hΔ (QuotientGroup.mk uζ))))) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  have hqpow := tateParam_quot_velu_dvr
    (tateParamR (E.baseChange (p.adicCompletion L)) h)
    (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
    (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h)
    hΔ hl hζ hlu hu uζ hζu hζl hord hql h2 h2K hodd v w hv hw
    (E'.baseChange (p.adicCompletion L)) h' hW'
  exact minDeltaExp_eq_mul_of_tateParamR (R := p.adicCompletionIntegers L) E E' l h h' p hp
    C hC hc4ne hc4 C' hC' hc4ne' hc4' hqpow

def minDeltaExp_eq_mul_of_veluMu.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(数体の素点での Δ_min の関係——局所の連鎖の終点)",
    sectionId := "genell-lemma-3-2" }

def minDeltaExp_eq_mul_of_veluMu.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tateParam_quot_velu_dvr(q_{E′} = q_E^l、第 900、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.tateParam_quot_velu_dvr") 1,
    .citation "[ABC3]" "minDeltaExp_eq_mul_of_tateParamR(Δ_min へ、第 892、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.minDeltaExp_eq_mul_of_tateParamR") 1,
    .implicitStep
      ("☆残るのは局所の幾何データを大域の仮説から作る段だけである: " ++
       "分裂乗法還元、極小モデルの完備化への移行、そして H の像が ⟨Φ(ζ)⟩ になること") 4 ]

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**[GenEll] 葉 1 ——
`j` で受ける形**。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★★★**2026-09-01（第 914）**——`tateParam_quot_velu_dvr`（第 900）は
`W′` が Vélu の商に**等しい**ことを要求していたが、
実際に得られるのは変数変換を挑んだ `C • W′ = …` である。
☆`j` は変数変換で不変なので、**`j` で受ければその隔たりが消える**。 -/
theorem tateParam_quot_velu_j_dvr {R : Type} [CommRing R] [IsDomain R] [CharZero R]
    [IsDiscreteValuationRing R] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [CharZero K] [Algebra R K] [IsFractionRing R K]
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (hq0 : q ≠ 0)
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hlu : IsUnit ((l : R)))
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (uζ : Kˣ) (hζu : algebraMap R K ζ = (uζ : K)) (hζl : uζ ^ l = 1)
    (hord : ∀ n : ℕ, 0 < n → n < l → uζ ^ n ≠ 1)
    (hql : q ^ l ∈ IsLocalRing.maximalIdeal R)
    (h2 : (2 : R) ≠ 0) (h2K : (2 : K) ≠ 0) (hodd : l ≠ 2)
    (v w : R)
    (hv : v = ∑ i ∈ (range l).erase 0,
      veluV2 (tateCurveAt q hq)
        (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
    (hw : 2 * w = ∑ i ∈ (range l).erase 0,
      (veluU (tateCurveAt q hq)
          (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
          (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        + 2 * (veluV2 (tateCurveAt q hq)
                (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)))
    (W' : WeierstrassCurve K) [W'.IsElliptic] [W'.IsMinimal R]
    (hsplit : W'.HasSplitMultiplicativeReduction R)
    [((veluCurve (tateCurveAt q hq) v w).map (algebraMap R K)).IsElliptic]
    [((tateCurveAt (q ^ l) hql).map (algebraMap R K)).IsElliptic]
    [(veluQuotientFull ((tateCurveAt q hq).map (algebraMap R K))
        (((range l).erase 0).image
          (fun k : ℕ => pointCoords
            (k • tatePhi (mkTateSetup (K := K) q hq hq0) hΔ
              (QuotientGroup.mk uζ))))).IsElliptic]
    (hW'j : W'.j = (veluQuotientFull ((tateCurveAt q hq).map (algebraMap R K))
      (((range l).erase 0).image
        (fun k : ℕ => pointCoords
          (k • tatePhi (mkTateSetup (K := K) q hq hq0) hΔ (QuotientGroup.mk uζ))))).j) :
    tateParamR W' hsplit = q ^ l := by
  haveI i1 : ((veluCurve (tateCurveAt (mkTateSetup (K := K) q hq hq0).q
      (mkTateSetup (K := K) q hq hq0).hq) v w).map (algebraMap R K)).IsElliptic :=
    inferInstanceAs (((veluCurve (tateCurveAt q hq) v w).map (algebraMap R K)).IsElliptic)
  haveI i2 : ((tateCurveAt ((mkTateSetup (K := K) q hq hq0).q ^ l) hql).map
      (algebraMap R K)).IsElliptic :=
    inferInstanceAs (((tateCurveAt (q ^ l) hql).map (algebraMap R K)).IsElliptic)
  have hquot := veluQuotientFull_tate_mu (mkTateSetup q hq hq0) hΔ
    (dvrTatePhiAddEquiv q hq hq0 hΔ) (fun _ => rfl) hl.pos ζ uζ hζu hζl hord hu v w h2K hv hw
  exact tateParam_quot_mu hl hζ hlu hu q hq hql h2
    (fun i hi => tateDXpair_ne_zero_of_mu (mkTateSetup q hq hq0) hΔ
      (dvrTatePhiAddEquiv q hq hq0 hΔ) (fun _ => rfl) hl hodd ζ uζ hζu hζl hord i
      (Finset.mem_erase.1 hi).1 (Finset.mem_range.1 (Finset.mem_erase.1 hi).2) (hu i hi))
    v w hv hw W' hsplit (hW'j.trans (ABC3.Found.GenEll.j_congr_curve hquot))

def tateParam_quot_velu_j_dvr.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(j で受ける形——変数変換の隔たりが消える)",
    sectionId := "genell-lemma-3-2" }

def tateParam_quot_velu_j_dvr.needs : List ProofObligation :=
  [ .citation "[ABC3]" "veluQuotientFull_tate_mu(第 890、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.veluQuotientFull_tate_mu") 1,
    .citation "[ABC3]" "tateParam_quot_mu(第 883、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.tateParam_quot_mu") 1,
    .citation "[ABC3]" "j_congr_curve(曲線が等しければ j も等しい、第 913、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.j_congr_curve") 1 ]

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★**[GenEll] 葉 1 ——
`veluCurve` の形で受ける形**。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★★★**2026-09-01（第 927）**——捧りを使うと商は
`veluQuotientFull`（点集合）ではなく **`veluCurve`**（`v`・`w` の形）で出てくる
（`quadTwist_veluQuotientFull`、第 925）。
☆本節点はその形を直接受け、侧条件 `hDX` だけを内部で消す。
★**非分裂の降下に使える形**である。 -/
theorem tateParam_quot_veluCurve_dvr {R : Type} [CommRing R] [IsDomain R] [CharZero R]
    [IsDiscreteValuationRing R] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [CharZero K] [Algebra R K] [IsFractionRing R K]
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (hq0 : q ≠ 0)
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hlu : IsUnit ((l : R)))
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (uζ : Kˣ) (hζu : algebraMap R K ζ = (uζ : K)) (hζl : uζ ^ l = 1)
    (hord : ∀ n : ℕ, 0 < n → n < l → uζ ^ n ≠ 1)
    (hql : q ^ l ∈ IsLocalRing.maximalIdeal R)
    (h2 : (2 : R) ≠ 0) (hodd : l ≠ 2)
    (v w : R)
    (hv : v = ∑ i ∈ (range l).erase 0,
      veluV2 (tateCurveAt q hq)
        (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
    (hw : 2 * w = ∑ i ∈ (range l).erase 0,
      (veluU (tateCurveAt q hq)
          (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
          (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        + 2 * (veluV2 (tateCurveAt q hq)
                (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)))
    (W' : WeierstrassCurve K) [W'.IsElliptic] [W'.IsMinimal R]
    (hsplit : W'.HasSplitMultiplicativeReduction R)
    [((veluCurve (tateCurveAt q hq) v w).map (algebraMap R K)).IsElliptic]
    [((tateCurveAt (q ^ l) hql).map (algebraMap R K)).IsElliptic]
    (hquot : W'.j = ((veluCurve (tateCurveAt q hq) v w).map (algebraMap R K)).j) :
    tateParamR W' hsplit = q ^ l :=
  tateParam_quot_mu hl hζ hlu hu q hq hql h2
    (fun i hi => tateDXpair_ne_zero_of_mu (mkTateSetup q hq hq0) hΔ
      (dvrTatePhiAddEquiv q hq hq0 hΔ) (fun _ => rfl) hl hodd ζ uζ hζu hζl hord i
      (Finset.mem_erase.1 hi).1 (Finset.mem_range.1 (Finset.mem_erase.1 hi).2) (hu i hi))
    v w hv hw W' hsplit hquot

def tateParam_quot_veluCurve_dvr.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(veluCurve の形で受ける形——非分裂の降下に使える)",
    sectionId := "genell-lemma-3-2" }

def tateParam_quot_veluCurve_dvr.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tateParam_quot_mu(第 883、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.tateParam_quot_mu") 1,
    .citation "[ABC3]" "tateDXpair_ne_zero_of_mu(侧条件は定理、第 894、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateDXpair_ne_zero_of_mu") 1 ]

def j_velu_tate_mu.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(Vélu の商の j は j(E_{q^l}) に等しい)",
    sectionId := "genell-lemma-3-2" }

def j_velu_tate_mu.needs : List ProofObligation :=
  [ .citation "[ABC3]" "c4_velu_tate(目標 (1)、第 853、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.c4_velu_tate") 1,
    .citation "[ABC3]" "c6_velu_tate(目標 (2)、第 867、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.c6_velu_tate") 1,
    .citation "[ABC3]" "j_velu_tate_eq(c₄・c₆ から j へ、第 838、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.j_velu_tate_eq") 1 ]

def c4_velu_tate.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(訂正後の目標 (1)——c₄ は l⁴ 倍)",
    sectionId := "genell-lemma-3-2" }

def c4_velu_tate.needs : List ProofObligation :=
  [ .citation "[ABC3]" "sum_veluV2_eq_sum_tateDYpair(v = ∑ DY、第 846、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.sum_veluV2_eq_sum_tateDYpair") 1,
    .citation "[ABC3]" "tateD2Xpair_eq(D²X = 2DY + DX、第 852、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateD2Xpair_eq") 1,
    .citation "[ABC3]" "sum_mu_dxpair_zero(∑ DX = 0、第 853)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.sum_mu_dxpair_zero") 1,
    .citation "[ABC3]" "sum_mu_d2xpair(∑ D²X の閉じた式、第 853)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.sum_mu_d2xpair") 1,
    .citation "[ABC3]" "tateCurveAt_c4_eq(c₄ = 1 + 240 s₃、第 853、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateCurveAt_c4_eq") 1 ]

def c6_velu_tate.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(訂正後の目標 (2)——c₆ は −l⁶ 倍)",
    sectionId := "genell-lemma-3-2" }

def c6_velu_tate.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tate_d4x(D⁴X = 12X·D²X + 12(DX)² + D²X、第 866、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.tate_d4x") 1,
    .citation "[ABC3]" "sum_mu_d4xpair(252·∑D⁴X の閉じた式、第 866、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.sum_mu_d4xpair") 1,
    .citation "[ABC3]" "tateCurveAt_c6_eq(c₆ = −1 + 504s₅、第 867、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateCurveAt_c6_eq") 1 ]

/-! ## ★★★★★★★★★★★★★★★★`Lemma 3.5` が直接消費する形（`j` の付値） -/

/-- **[GenEll] `Lemma 3.5` の残る入力 (1)**——悪い還元の素点での `v_p(j′) = l·v_p(j)`。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★★`Found/GaloisRep/Lemma35Concrete.lean` の `lemma_3_5_velu_j`（`§9-1177`、第 751）が
**そのまま消費する形**である。

☆機構: 悪い還元の素点では `E` は Tate 曲線 `E_q` であり `v_p(j) = −v_p(q)`。
`H` が `μ_l` に対応するとき `E′ = E_q/μ_l = E_{q^l}` なので
`v_p(j′) = −l·v_p(q) = l·v_p(j)`。★これが `tateModel_of_quot_mu` の内容である。

☆仮説 `hmu` は「`H` は各悪い素点で `μ_l` に対応する」という `Lemma 3.5` の設定
（原文の *global rank one subgroup*）を型にしたものである。 -/
theorem jExp_velu_bad {L : Type} [Field L] [NumberField L]
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    (l : ℕ) (hl : Nat.Prime l) (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E (((Finset.range l).erase 0).image
      (fun k : ℕ => pointCoords (k • Q))))
    (hmu : True)
    (p : HeightOneSpectrum (𝓞 L)) (hss : SemistableAt p E) (hbad : jExp p E < 0) :
    jExp p E' = (l : ℤ) * jExp p E := by
  sorry

/-- **[GenEll] `Lemma 3.5` の残る入力 (2)**——良い還元は同種写像で保たれる。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`v_p(j) ≥ 0` かつ半安定なら `E` は `p` で良い還元をもつ。★★同種な曲線は
同じ還元型をもつ（Néron–Ogg–Shafarevich、あるいは Tate 加群の不分岐性）ので
`E′` も良い還元をもち、`v_p(j′) ≥ 0` である。

☆こちらは `v_p(j′) = l·v_p(j)` **より弱く**てよい——半安定なら
`v_p(Δ_min) = max(0, −v_p(j))` なので、良い還元の素点では両辺とも `0` になる。 -/
theorem jExp_velu_good {L : Type} [Field L] [NumberField L]
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    (l : ℕ) (hl : Nat.Prime l) (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E (((Finset.range l).erase 0).image
      (fun k : ℕ => pointCoords (k • Q))))
    (p : HeightOneSpectrum (𝓞 L)) (hss : SemistableAt p E) (hss' : SemistableAt p E')
    (hgood : 0 ≤ jExp p E) :
    0 ≤ jExp p E' := by
  sorry

def jExp_velu_bad.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(残る入力 (1)——悪い素点での v_p(j′) = l·v_p(j))",
    sectionId := "genell-lemma-3-5" }

def jExp_velu_bad.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tateModel_of_quot_mu(E_q/μ_l は母数 q^l の Tate 曲線)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.tateModel_of_quot_mu") 15,
    .implicitStep
      ("★Tate 曲線では v_p(j) = −v_p(q)(TateJ.lean)。母数が q^l なら " ++
       "v_p(j′) = −l·v_p(q) = l·v_p(j)") 4 ]

def jExp_velu_good.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(残る入力 (2)——良い還元は同種写像で保たれる)",
    sectionId := "genell-lemma-3-5" }

def jExp_velu_good.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★★★★★**後日談（2026-08-31、第 902）——この節点はもう要らない**。" ++
       "`Lemma 3.5` の証明が使うのは `l·deg∞(E) ≤ deg∞(E′)` の片側だけであり、" ++
       "良い素点では左辺が 0 なので `minDeltaExp_nonneg` だけで自動的に成り立つ" ++
       "（`Found/GaloisRep/Lemma35Ineq.lean` の `minDeltaExp_le_of_jExp_bad`）。" ++
       "★★★★したがって Néron–Ogg–Shafarevich の義務は**消えた**") 1,
    .implicitStep
      ("★★同種な楼円曲線は同じ還元型をもつ（Neron-Ogg-Shafarevich、"
       ++ "あるいは Tate 加群の不分岐性）。mathlib には無い") 8,
    .implicitStep
      ("★★★★**測定（2026-08-31、第 827）——葉 1 に帰着する**。"
       ++ "対偶同種写像 `E′ → E`（核は位数 l）を使うと、"
       ++ "対偶を取って `jExp p E′ < 0`（= E′ が p で乗法還元）と仮定すれば、"
       ++ "葉 1（jExp_velu_bad）を E′ に適用して `jExp p E = l·jExp p E′ < 0` となり、"
       ++ "仮説 `0 ≤ jExp p E` に矛盾する。★★ゆえに本命題は"
       ++ "**葉 1 ＋「E は E′ の Vélu の商でもある」の 2 つに分解される**。"
       ++ "☆Néron-Ogg-Shafarevich も Néron 模型も要らない") 4,
    .implicitStep
      ("☆(a) p ∤ l のときの別道（第 785）: l-捧じれ点の x 座標は p で整だから"
       ++ "Vélu の和も整で、Vélu の式は還元と可換する。"
       ++ "☆(b) p ∣ l のときは有限平坦群スキームか Tate 曲線が要る") 4,
    .citation "[ABC3]" "jExp_velu_bad（対偶側へ適用する、第 827 の帰着）"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.jExp_velu_bad") 1 ]

/-! ## ★★★★★★★★★★★★★★★★★★★★第 947 と 927 を繋ぐ——`ζ` を消す -/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★**[GenEll] 有理な `l`-捉れ点だけで
`q_{E′} = q_E^l`**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 948）**——`tateParam_quot_velu_dvr`（第 927）は
`ζ`・`uζ`・`hζ`・`hζu`・`hζl`・`hord` の 6 つを受けていた。
☆本定理はそれを `exists_primitiveRoot_of_torsion_point`（第 947）で埋め、

    **`P` が位数 `l` の点で、`l ∤ v(q)`**

だけに置き換える。★`ζ` はもはや引数に現れない——
Vélu の帳簿（`hu`・`hv`・`hw`）だけが `ζ` について全称で残る。

☆これが `isMuAtBadPrimes_of_veluQuotient` の局所の段そのものである。 -/
theorem tateParam_quot_velu_of_torsion {R : Type} [CommRing R] [IsDomain R] [CharZero R]
    [IsDiscreteValuationRing R] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [CharZero K] [Algebra R K] [IsFractionRing R K]
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (hq0 : q ≠ 0)
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hcop : ¬ ((l : ℤ) ∣ vAdd (mkTateSetup (K := K) q hq hq0).v
      (mkTateSetup (K := K) q hq hq0).Q))
    (hlu : IsUnit ((l : R)))
    (hql : q ^ l ∈ IsLocalRing.maximalIdeal R)
    (h2 : (2 : R) ≠ 0) (h2K : (2 : K) ≠ 0)
    (hvw : ∀ ζ : R, IsPrimitiveRoot ζ l → ∃ v w : R,
      v = ∑ i ∈ (range l).erase 0,
          veluV2 (tateCurveAt q hq)
            (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        ∧ 2 * w = ∑ i ∈ (range l).erase 0,
          (veluU (tateCurveAt q hq)
              (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            + 2 * (veluV2 (tateCurveAt q hq)
                    (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                    (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                  * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
        ∧ ((veluCurve (tateCurveAt q hq) v w).map (algebraMap R K)).IsElliptic)
    [((tateCurveAt (q ^ l) hql).map (algebraMap R K)).IsElliptic]
    (P : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Point)
    (hP : l • P = 0) (hP0 : P ≠ 0)
    (W' : WeierstrassCurve K) [W'.IsElliptic] [W'.IsMinimal R]
    (hsplit : W'.HasSplitMultiplicativeReduction R)
    (hW' : W' = veluQuotientFull ((tateCurveAt q hq).map (algebraMap R K))
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))) :
    tateParamR W' hsplit = q ^ l := by
  obtain ⟨ζ, uζ, hζ, hζu, hζl, hord, hPz⟩ :=
    exists_primitiveRoot_of_torsion_point q hq hq0 hΔ hl hcop P hP hP0
  obtain ⟨v, w, hv, hw, hell⟩ := hvw ζ hζ
  haveI := hell
  refine tateParam_quot_velu_dvr q hq hq0 hΔ hl hζ hlu
    (ABC3.Found.GenEll.isUnit_one_sub_pow_of_isUnit_natCast hl.pos hζ hlu)
    uζ hζu hζl hord hql h2 h2K hodd v w hv hw W' hsplit ?_
  rw [hW', hPz]
  rfl

def tateParam_quot_velu_of_torsion.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(有理な l-捉れ点だけで q_{E′} = q_E^l)",
    sectionId := "genell-lemma-3-5" }

def tateParam_quot_velu_of_torsion.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_primitiveRoot_of_torsion_point(ζ を作る、第 947、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_primitiveRoot_of_torsion_point") 1,
    .citation "[ABC3]" "tateParam_quot_velu_dvr(Vélu の商の Tate 母数、第 927、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.tateParam_quot_velu_dvr") 1 ]

/-! ## ★★★★★★★★★★★★★★★★★★★★`j` で受ける形の捉れ点版 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★**[GenEll] 有理な `l`-捉れ点と `j` の一致だけで
`q_{E′} = q_E^l`**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 955）**——第 948 は `W′` が Vélu の商に**等しい**ことを
要求していたが、実際に得られるのは変数変換を除いた `j` の一致である。
☆`j_veluQuotientFull_nsmul_variableChange`（第 950）がその `j` を与えるので、
本定理が**そのままの形で受けられる**。

★`ζ` は引数に現れない——`exists_primitiveRoot_of_torsion_point`（第 947）が作る。 -/
theorem tateParam_quot_velu_j_of_torsion {R : Type} [CommRing R] [IsDomain R] [CharZero R]
    [IsDiscreteValuationRing R] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [CharZero K] [Algebra R K] [IsFractionRing R K]
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (hq0 : q ≠ 0)
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hcop : ¬ ((l : ℤ) ∣ vAdd (mkTateSetup (K := K) q hq hq0).v
      (mkTateSetup (K := K) q hq hq0).Q))
    (hlu : IsUnit ((l : R)))
    (hql : q ^ l ∈ IsLocalRing.maximalIdeal R)
    (h2 : (2 : R) ≠ 0) (h2K : (2 : K) ≠ 0)
    (hvw : ∀ ζ : R, IsPrimitiveRoot ζ l → ∃ v w : R,
      v = ∑ i ∈ (range l).erase 0,
          veluV2 (tateCurveAt q hq)
            (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        ∧ 2 * w = ∑ i ∈ (range l).erase 0,
          (veluU (tateCurveAt q hq)
              (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            + 2 * (veluV2 (tateCurveAt q hq)
                    (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                    (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                  * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
        ∧ ((veluCurve (tateCurveAt q hq) v w).map (algebraMap R K)).IsElliptic)
    [((tateCurveAt (q ^ l) hql).map (algebraMap R K)).IsElliptic]
    (P : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Point)
    (hP : l • P = 0) (hP0 : P ≠ 0)
    (hellQ : (veluQuotientFull ((tateCurveAt q hq).map (algebraMap R K))
        (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))).IsElliptic)
    (W' : WeierstrassCurve K) [W'.IsElliptic] [W'.IsMinimal R]
    (hsplit : W'.HasSplitMultiplicativeReduction R)
    (hW'j : W'.j = (veluQuotientFull ((tateCurveAt q hq).map (algebraMap R K))
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))).j) :
    tateParamR W' hsplit = q ^ l := by
  obtain ⟨ζ, uζ, hζ, hζu, hζl, hord, hPz⟩ :=
    exists_primitiveRoot_of_torsion_point q hq hq0 hΔ hl hcop P hP hP0
  subst hPz
  haveI hQell : (veluQuotientFull ((tateCurveAt q hq).map (algebraMap R K))
      (((range l).erase 0).image (fun k : ℕ => pointCoords
        (k • tatePhi (mkTateSetup (K := K) q hq hq0) hΔ
          (QuotientGroup.mk uζ))))).IsElliptic := hellQ
  obtain ⟨v, w, hv, hw, hell⟩ := hvw ζ hζ
  haveI := hell
  exact tateParam_quot_velu_j_dvr q hq hq0 hΔ hl hζ hlu
    (ABC3.Found.GenEll.isUnit_one_sub_pow_of_isUnit_natCast hl.pos hζ hlu)
    uζ hζu hζl hord hql h2 h2K hodd v w hv hw W' hsplit hW'j

def tateParam_quot_velu_j_of_torsion.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(有理な l-捉れ点と j の一致だけで q_{E′} = q_E^l)",
    sectionId := "genell-lemma-3-5" }

def tateParam_quot_velu_j_of_torsion.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_primitiveRoot_of_torsion_point(ζ を作る、第 947、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_primitiveRoot_of_torsion_point") 1,
    .citation "[ABC3]" "tateParam_quot_velu_j_dvr(j で受ける形、第 914、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.tateParam_quot_velu_j_dvr") 1,
    .citation "[ABC3]" "isUnit_one_sub_pow_of_isUnit_natCast(hu は hlu から、第 951、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.isUnit_one_sub_pow_of_isUnit_natCast") 1 ]

/-! ## ★★★★★★★★★★★★★★★★★★★★捉れ点で受ける形の局所の終点 -/

open ABC3.Skeleton.GenEll ABC3.Found.GaloisRep ABC3.Found.GenEll IsDedekindDomain
  NumberField WeierstrassCurve Finset in
open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★**[GenEll] 数体の素点での `Δ_min` の関係——
捉れ点と `j` で受ける形**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 965）**——第 904 は `ζ`・`uζ`・`hu`・`v`・`w` を受け、
さらに `W′` が Vélu の商に**等しい**ことを要求していた。
☆本定理はそれを

* 位数 `l` の点 `P`（`P ≠ 0`）と `l ∤ v(q)`（第 947・946）
* `j` の一致（第 950 が与える）
* Vélu の帳簿 `hvw`（第 961・962 が与える）

に置き換える。★これが `isMuAtBadPrimes_of_veluQuotient` の
**各悪い素点での終点**である。 -/
theorem minDeltaExp_eq_mul_of_torsion {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E E' : WeierstrassCurve L)
    [(E.baseChange (p.adicCompletion L)).IsElliptic]
    [(E.baseChange (p.adicCompletion L)).IsMinimal (p.adicCompletionIntegers L)]
    [(E'.baseChange (p.adicCompletion L)).IsElliptic]
    [(E'.baseChange (p.adicCompletion L)).IsMinimal (p.adicCompletionIntegers L)]
    (h : (E.baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    (h' : (E'.baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation (p.adicCompletion L)
        (IsDiscreteValuationRing.maximalIdeal (p.adicCompletionIntegers L)))
        (algebraMap L (p.adicCompletion L) x) = (HeightOneSpectrum.valuation L p) x)
    (C : VariableChange L) (hC : IsMinimal (primeSubring p) (C • E))
    (hc4ne : (C • E).c₄ ≠ 0) (hc4 : valAdd p (Units.mk0 ((C • E).c₄) hc4ne) = 0)
    (C' : VariableChange L) (hC' : IsMinimal (primeSubring p) (C' • E'))
    (hc4ne' : (C' • E').c₄ ≠ 0) (hc4' : valAdd p (Units.mk0 ((C' • E').c₄) hc4ne') = 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hΔ : ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).toAffine.Δ ≠ 0)
    (hcop : ¬ ((l : ℤ) ∣ vAdd (mkTateSetup (K := p.adicCompletion L)
        (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h)).v
      (mkTateSetup (K := p.adicCompletion L)
        (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h)).Q))
    (hlu : IsUnit ((l : (p.adicCompletionIntegers L))))
    (hql : (tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l
      ∈ IsLocalRing.maximalIdeal (p.adicCompletionIntegers L))
    (h2 : (2 : (p.adicCompletionIntegers L)) ≠ 0)
    (h2K : (2 : (p.adicCompletion L)) ≠ 0)
    (hvw : ∀ ζ : (p.adicCompletionIntegers L), IsPrimitiveRoot ζ l →
      ∃ v w : (p.adicCompletionIntegers L),
      v = ∑ i ∈ (range l).erase 0,
          veluV2 (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
            (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
            (tateXpair (ζ ^ i)
              ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
              (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
            (tateYpair (ζ ^ i)
              ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
              (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
        ∧ 2 * w = ∑ i ∈ (range l).erase 0,
          (veluU (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
              (tateXpair (ζ ^ i)
                ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                (tateParamR (E.baseChange (p.adicCompletion L)) h)
                (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
              (tateYpair (ζ ^ i)
                ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                (tateParamR (E.baseChange (p.adicCompletion L)) h)
                (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
            + 2 * (veluV2 (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
                    (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
                    (tateXpair (ζ ^ i)
                      ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR (E.baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
                    (tateYpair (ζ ^ i)
                      ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR (E.baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
                  * tateXpair (ζ ^ i)
                      ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR (E.baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)))
        ∧ ((veluCurve (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)) v w).map
            (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).IsElliptic)
    [((tateCurveAt ((tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) hql).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).IsElliptic]
    (P : ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).toAffine.Point)
    (hP : l • P = 0) (hP0 : P ≠ 0)
    (hellQ : (veluQuotientFull
        ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
          (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
          (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
        (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))).IsElliptic)
    (hW'j : (E'.baseChange (p.adicCompletion L)).j
      = (veluQuotientFull
        ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
          (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
          (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
        (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))).j) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  have hqpow := tateParam_quot_velu_j_of_torsion
    (tateParamR (E.baseChange (p.adicCompletion L)) h)
    (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
    (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h)
    hΔ hl hodd hcop hlu hql h2 h2K hvw P hP hP0 hellQ
    (E'.baseChange (p.adicCompletion L)) h' hW'j
  exact minDeltaExp_eq_mul_of_tateParamR (R := p.adicCompletionIntegers L) E E' l h h' p hp
    C hC hc4ne hc4 C' hC' hc4ne' hc4' hqpow

def minDeltaExp_eq_mul_of_torsion.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(数体の素点での Δ_min の関係——捉れ点と j で受ける形)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_mul_of_torsion.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tateParam_quot_velu_j_of_torsion(q_{E′} = q_E^l、第 955、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.tateParam_quot_velu_j_of_torsion") 1,
    .citation "[ABC3]" "minDeltaExp_eq_mul_of_tateParamR(Δ_min へ、第 892、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.minDeltaExp_eq_mul_of_tateParamR") 1 ]

/-! ## ★★★★★★★★★★★★★★★★★★★★大域の Vélu の商で受ける形 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★**[GenEll] 各悪い素点での `Δ_min` の関係——
大域の Vélu の商で受ける形**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 972）**——第 965 は Tate モデルの上の点 `P` と
`j` の一致 `hW′j`、それに商の楕円性 `hellQ` を受けていた。
☆本定理はその 4 つを**大域のデータ `Q`・`hQ`・`hE′` から作る**:

* `P`・`hP`・`hP0`・`hW′j` ← `exists_point_j_tateModel`（第 970）
* `ζ`・`uζ` ← `exists_primitiveRoot_of_torsion_point`（第 947）
* `hellQ` ← `isElliptic_veluQuotient_tate_mu`（第 971）

★残るのは**各悪い素点で局所データを供給する段**だけである。 -/
theorem minDeltaExp_eq_mul_of_globalVelu {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange (p.adicCompletion L)).IsElliptic]
    [(E.baseChange (p.adicCompletion L)).IsMinimal (p.adicCompletionIntegers L)]
    [(E'.baseChange (p.adicCompletion L)).IsElliptic]
    [(E'.baseChange (p.adicCompletion L)).IsMinimal (p.adicCompletionIntegers L)]
    (h : (E.baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    (h' : (E'.baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation (p.adicCompletion L)
        (IsDiscreteValuationRing.maximalIdeal (p.adicCompletionIntegers L)))
        (algebraMap L (p.adicCompletion L) x) = (HeightOneSpectrum.valuation L p) x)
    (C : VariableChange L) (hC : IsMinimal (primeSubring p) (C • E))
    (hc4ne : (C • E).c₄ ≠ 0) (hc4 : valAdd p (Units.mk0 ((C • E).c₄) hc4ne) = 0)
    (C' : VariableChange L) (hC' : IsMinimal (primeSubring p) (C' • E'))
    (hc4ne' : (C' • E').c₄ ≠ 0) (hc4' : valAdd p (Units.mk0 ((C' • E').c₄) hc4ne') = 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hΔ : ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).toAffine.Δ ≠ 0)
    (hcop : ¬ ((l : ℤ) ∣ vAdd (mkTateSetup (K := p.adicCompletion L)
        (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h)).v
      (mkTateSetup (K := p.adicCompletion L)
        (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h)).Q))
    (hlu : IsUnit ((l : (p.adicCompletionIntegers L))))
    (hql : (tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l
      ∈ IsLocalRing.maximalIdeal (p.adicCompletionIntegers L))
    (h2 : (2 : (p.adicCompletionIntegers L)) ≠ 0)
    (h2K : (2 : (p.adicCompletion L)) ≠ 0)
    (hvw : ∀ ζ : (p.adicCompletionIntegers L), IsPrimitiveRoot ζ l →
      ∃ v w : (p.adicCompletionIntegers L),
      v = ∑ i ∈ (range l).erase 0,
          veluV2 (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
            (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
            (tateXpair (ζ ^ i)
              ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
              (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
            (tateYpair (ζ ^ i)
              ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
              (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
        ∧ 2 * w = ∑ i ∈ (range l).erase 0,
          (veluU (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
              (tateXpair (ζ ^ i)
                ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                (tateParamR (E.baseChange (p.adicCompletion L)) h)
                (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
              (tateYpair (ζ ^ i)
                ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                (tateParamR (E.baseChange (p.adicCompletion L)) h)
                (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
            + 2 * (veluV2 (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
                    (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
                    (tateXpair (ζ ^ i)
                      ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR (E.baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
                    (tateYpair (ζ ^ i)
                      ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR (E.baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
                  * tateXpair (ζ ^ i)
                      ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR (E.baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)))
        ∧ ((veluCurve (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)) v w).map
            (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).IsElliptic)
    [((tateCurveAt ((tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) hql).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).IsElliptic]
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  obtain ⟨P, hP, hP0, hj⟩ := exists_point_j_tateModel p E E' h hl hQ h2K hE'
  obtain ⟨ζ, uζ, hζ, hζu, hζl, hord, hPz⟩ :=
    exists_primitiveRoot_of_torsion_point
      (tateParamR (E.baseChange (p.adicCompletion L)) h)
      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
      (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h) hΔ hl hcop P hP hP0
  obtain ⟨v, w, hv, hw, hell⟩ := hvw ζ hζ
  have hellQ := isElliptic_veluQuotient_tate_mu
      (tateParamR (E.baseChange (p.adicCompletion L)) h)
      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
      (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h) hΔ hl hlu h2K
      ζ uζ hζ hζu hζl hord v w hv hw hell
  rw [← hPz] at hellQ
  exact minDeltaExp_eq_mul_of_torsion p E E' h h' hp C hC hc4ne hc4 C' hC' hc4ne' hc4'
    hl hodd hΔ hcop hlu hql h2 h2K hvw P hP hP0 hellQ (hj hellQ)

def minDeltaExp_eq_mul_of_globalVelu.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(各悪い素点での Δ_min の関係——大域の Vélu の商で受ける形)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_mul_of_globalVelu.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_point_j_tateModel(P と hW′j、第 970、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_point_j_tateModel") 1,
    .citation "[ABC3]" "isElliptic_veluQuotient_tate_mu(hellQ、第 971、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isElliptic_veluQuotient_tate_mu") 1,
    .citation "[ABC3]" "minDeltaExp_eq_mul_of_torsion(終点、第 965、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_of_torsion") 1 ]

/-! ## ★★★★★★★★★★★★★★★★★★★★第 996 —— `j` の段で止める

★第 995 により `E′` の側に要るのは **`j` の一致だけ**になった。
☆本ブロックは第 914（`tateParam_quot_velu_j_dvr`）の**前半**を取り出す:

    `j(μ_l による Vélu の商) = j(E_{q^l})`

★第 914 はこの後 `tateParamR_eq_of_j_tateCurveAt` に流して `q_{E′} = q_E^l` を出すが、
そこで `E′` の分裂性が要る。☆`j` の段で止めれば要らない。 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★**[GenEll] `μ_l` による Vélu の商の `j` は
`E_{q^l}` の `j`**。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★★★**2026-09-01（第 996）**——第 914 の前半である。
☆`veluQuotientFull_tate_mu`（第 890）で商を `veluCurve` の形にし、
`j_velu_tate_mu_map`（第 886）で `E_{q^l}` に繋ぐ。
★`E′` の分裂性も極小モデルも要らない。 -/
theorem j_veluQuot_eq_j_tate_pow {R : Type} [CommRing R] [IsDomain R] [CharZero R]
    [IsDiscreteValuationRing R] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [CharZero K] [Algebra R K] [IsFractionRing R K]
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (hq0 : q ≠ 0)
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hlu : IsUnit ((l : R)))
    (uζ : Kˣ) (hζu : algebraMap R K ζ = (uζ : K)) (hζl : uζ ^ l = 1)
    (hord : ∀ n : ℕ, 0 < n → n < l → uζ ^ n ≠ 1)
    (hql : q ^ l ∈ IsLocalRing.maximalIdeal R)
    (h2 : (2 : R) ≠ 0) (h2K : (2 : K) ≠ 0) (hodd : l ≠ 2)
    (v w : R)
    (hv : v = ∑ i ∈ (range l).erase 0,
      veluV2 (tateCurveAt q hq)
        (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
    (hw : 2 * w = ∑ i ∈ (range l).erase 0,
      (veluU (tateCurveAt q hq)
          (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
          (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        + 2 * (veluV2 (tateCurveAt q hq)
                (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)))
    [((veluCurve (tateCurveAt q hq) v w).map (algebraMap R K)).IsElliptic]
    [((tateCurveAt (q ^ l) hql).map (algebraMap R K)).IsElliptic]
    [(veluQuotientFull ((tateCurveAt q hq).map (algebraMap R K))
        (((range l).erase 0).image
          (fun k : ℕ => pointCoords
            (k • tatePhi (mkTateSetup (K := K) q hq hq0) hΔ
              (QuotientGroup.mk uζ))))).IsElliptic] :
    (veluQuotientFull ((tateCurveAt q hq).map (algebraMap R K))
        (((range l).erase 0).image
          (fun k : ℕ => pointCoords
            (k • tatePhi (mkTateSetup (K := K) q hq hq0) hΔ (QuotientGroup.mk uζ))))).j
      = ((tateCurveAt (q ^ l) hql).map (algebraMap R K)).j := by
  haveI i1 : ((veluCurve (tateCurveAt (mkTateSetup (K := K) q hq hq0).q
      (mkTateSetup (K := K) q hq hq0).hq) v w).map (algebraMap R K)).IsElliptic :=
    inferInstanceAs (((veluCurve (tateCurveAt q hq) v w).map (algebraMap R K)).IsElliptic)
  haveI i2 : (veluQuotientFull ((tateCurveAt (mkTateSetup (K := K) q hq hq0).q
      (mkTateSetup (K := K) q hq hq0).hq).map (algebraMap R K))
      (((range l).erase 0).image (fun k : ℕ => pointCoords
        (k • tatePhi (mkTateSetup (K := K) q hq hq0) hΔ
          (QuotientGroup.mk uζ))))).IsElliptic :=
    inferInstanceAs ((veluQuotientFull ((tateCurveAt q hq).map (algebraMap R K))
      (((range l).erase 0).image (fun k : ℕ => pointCoords
        (k • tatePhi (mkTateSetup (K := K) q hq hq0) hΔ
          (QuotientGroup.mk uζ))))).IsElliptic)
  have hu := ABC3.Found.GenEll.isUnit_one_sub_pow_of_isUnit_natCast hl.pos hζ hlu
  have hquot := veluQuotientFull_tate_mu (mkTateSetup q hq hq0) hΔ
    (dvrTatePhiAddEquiv q hq hq0 hΔ) (fun _ => rfl) hl.pos ζ uζ hζu hζl hord hu v w h2K hv hw
  -- ☆`(mkTateSetup q hq hq0).q` と `q` は defeq だが構文的には違うので、
  -- `have` でゴールの形に言い直してから `rw` する（第 971 と同じ穴）。
  have hquot' : veluQuotientFull ((tateCurveAt q hq).map (algebraMap R K))
      (((range l).erase 0).image (fun k : ℕ => pointCoords
        (k • tatePhi (mkTateSetup (K := K) q hq hq0) hΔ (QuotientGroup.mk uζ))))
      = (veluCurve (tateCurveAt q hq) v w).map (algebraMap R K) := hquot
  rw [ABC3.Found.GenEll.j_congr_curve hquot']
  exact j_velu_tate_mu_map hl hζ hlu hu q hq hql h2
    (fun i hi => tateDXpair_ne_zero_of_mu (mkTateSetup q hq hq0) hΔ
      (dvrTatePhiAddEquiv q hq hq0 hΔ) (fun _ => rfl) hl hodd ζ uζ hζu hζl hord i
      (Finset.mem_erase.1 hi).1 (Finset.mem_range.1 (Finset.mem_erase.1 hi).2) (hu i hi))
    v w hv hw

def j_veluQuot_eq_j_tate_pow.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(μ_l による Vélu の商の j は E_{q^l} の j)",
    sectionId := "genell-lemma-3-2" }

def j_veluQuot_eq_j_tate_pow.needs : List ProofObligation :=
  [ .citation "[ABC3]" "veluQuotientFull_tate_mu(商を veluCurve の形に、第 890、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.veluQuotientFull_tate_mu") 1,
    .citation "[ABC3]" "j_velu_tate_mu_map(veluCurve の j は E_{q^l} の j、第 886、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.j_velu_tate_mu_map") 1 ]

/-! ## ★★★★★★★★★★★★★★★★★★★★第 999 —— 大域の Vélu の商で受ける形（軽い版）

★第 972 は `E′` について
分裂乗法還元 `h′`・極小モデル `[IsMinimal (E′ ⊗ Lv)]`・`C′`・`hC′`・`hc4ne′`・`hc4′`、
さらに `E` について `C`・`hC`・`hc4ne`・`hc4` を要求していた。

☆第 996（`j` の段で止める）＋第 997（`E′` の仮説は 2 本）＋第 998（`E′.j ≠ 0` は自動）で、
それらは**すべて半安定性 2 本に置き換わる**。★残るのは `E` の側の局所データだけである。 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★**[GenEll] 各悪い素点での `Δ_min` の関係——
大域の Vélu の商で受ける形（`E′` の局所データを要らない版）**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 999）**——第 972 の軽量版。
☆`E′` について要るのは `SemistableAt p E′` だけになった。 -/
theorem minDeltaExp_eq_mul_of_globalVelu' {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange (p.adicCompletion L)).IsElliptic]
    [(E.baseChange (p.adicCompletion L)).IsMinimal (p.adicCompletionIntegers L)]
    [(E'.baseChange (p.adicCompletion L)).IsElliptic]
    (h : (E.baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation (p.adicCompletion L)
        (IsDiscreteValuationRing.maximalIdeal (p.adicCompletionIntegers L)))
        (algebraMap L (p.adicCompletion L) x) = (HeightOneSpectrum.valuation L p) x)
    (hssE : SemistableAt p E) (hssE' : SemistableAt p E') (hjneg : jExp p E < 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hΔ : ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).toAffine.Δ ≠ 0)
    (hcop : ¬ ((l : ℤ) ∣ vAdd (mkTateSetup (K := p.adicCompletion L)
        (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h)).v
      (mkTateSetup (K := p.adicCompletion L)
        (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h)).Q))
    (hlu : IsUnit ((l : (p.adicCompletionIntegers L))))
    (hql : (tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l
      ∈ IsLocalRing.maximalIdeal (p.adicCompletionIntegers L))
    (h2 : (2 : (p.adicCompletionIntegers L)) ≠ 0)
    (h2K : (2 : (p.adicCompletion L)) ≠ 0)
    (hvw : ∀ ζ : (p.adicCompletionIntegers L), IsPrimitiveRoot ζ l →
      ∃ v w : (p.adicCompletionIntegers L),
      v = ∑ i ∈ (range l).erase 0,
          veluV2 (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
            (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
            (tateXpair (ζ ^ i)
              ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
              (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
            (tateYpair (ζ ^ i)
              ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
              (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
        ∧ 2 * w = ∑ i ∈ (range l).erase 0,
          (veluU (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
              (tateXpair (ζ ^ i)
                ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                (tateParamR (E.baseChange (p.adicCompletion L)) h)
                (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
              (tateYpair (ζ ^ i)
                ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                (tateParamR (E.baseChange (p.adicCompletion L)) h)
                (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
            + 2 * (veluV2 (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
                    (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
                    (tateXpair (ζ ^ i)
                      ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR (E.baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
                    (tateYpair (ζ ^ i)
                      ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR (E.baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
                  * tateXpair (ζ ^ i)
                      ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR (E.baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)))
        ∧ ((veluCurve (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)) v w).map
            (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).IsElliptic)
    [((tateCurveAt ((tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) hql).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).IsElliptic]
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  obtain ⟨P, hP, hP0, hj⟩ := exists_point_j_tateModel p E E' h hl hQ h2K hE'
  obtain ⟨ζ, uζ, hζ, hζu, hζl, hord, hPz⟩ :=
    exists_primitiveRoot_of_torsion_point
      (tateParamR (E.baseChange (p.adicCompletion L)) h)
      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
      (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h) hΔ hl hcop P hP hP0
  obtain ⟨v, w, hv, hw, hell⟩ := hvw ζ hζ
  -- ★`μ_l` の形での楕円性（第 971）
  haveI hellMu := isElliptic_veluQuotient_tate_mu
      (tateParamR (E.baseChange (p.adicCompletion L)) h)
      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
      (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h) hΔ hl hlu h2K
      ζ uζ hζ hζu hζl hord v w hv hw hell
  -- ☆`P` の形での楕円性
  haveI hellP : (veluQuotientFull
      ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
        (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))).IsElliptic := by
    rw [hPz]; exact hellMu
  -- ★曲線の水準で `P` と `μ_l` を繋ぐ（`j` を跨がないので motive が壊れない）
  have hcurveEq : veluQuotientFull
      ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
        (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))
      = veluQuotientFull
      ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
        (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
      (((range l).erase 0).image (fun k : ℕ => pointCoords
        (k • tatePhi (mkTateSetup (K := p.adicCompletion L)
          (tateParamR (E.baseChange (p.adicCompletion L)) h)
          (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
          (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h)) hΔ
          (QuotientGroup.mk uζ)))) := by
    rw [hPz]
    rfl
  -- ☆第 996 で `E_{q^l}` に繋ぐ
  have hjtate := j_veluQuot_eq_j_tate_pow
    (tateParamR (E.baseChange (p.adicCompletion L)) h)
    (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
    (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h) hΔ hl hζ hlu
    uζ hζu hζl hord hql h2 h2K hodd v w hv hw
  exact minDeltaExp_eq_mul_of_j_tate_pow p hp E E' h hssE hssE' hjneg hql
    ((hj hellP).trans ((ABC3.Found.GenEll.j_congr_curve hcurveEq).trans hjtate))

def minDeltaExp_eq_mul_of_globalVelu'.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(各悪い素点での Δ_min の関係——E′ の局所データを要らない版)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_mul_of_globalVelu'.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_point_j_tateModel(P と hW′j、第 970、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_point_j_tateModel") 1,
    .citation "[ABC3]" "j_veluQuot_eq_j_tate_pow(j の段で止める、第 996、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.j_veluQuot_eq_j_tate_pow") 1,
    .citation "[ABC3]" "minDeltaExp_eq_mul_of_j_tate_pow(終点、第 997、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.minDeltaExp_eq_mul_of_j_tate_pow") 1 ]

/-! ## ★★★★★★★★★★★★★★★★★★★★第 1000 —— 悪い素点で局所データを供給する

★第 999 が受ける局所データのうち、**自前で作れるものはすべて作る**:

| 入力 | 出どころ |
|---|---|
| `hp`（付値の橋） | 第 964 `valuation_algebraMap_adicCompletion` |
| `hΔ` | 第 977 `tateModel_map_Delta_ne_zero` |
| `hcop`（付値の言葉） | 第 978 `not_dvd_vAdd_tateParam_of_not_dvd_jExp` |
| `hql` | 第 977 `pow_mem_of_mem_ideal` |
| `h2`・`h2K` | 第 977 |
| `E_{q^l} ⊗ Lv` の楕円性 | `tateCurveAt_map_isElliptic` ＋ `tateCurveAt_c4_isUnit` |

☆残るのは **3 本だけ**である:

* `hmin`（完備化で極小）——第 973
* `h`（完備化で分裂乗法還元）——第 976＋993（`p ∣ 2` の非分裂が残件）
* `hlu`（`l` が単元）・`hvw`（Vélu の係数が整）
-/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★**[GenEll] 悪い素点での `Δ_min` の関係——
局所データを自前で作る形**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1000）**——第 999 の 10 本の局所入力を 3 本に絞る。 -/
theorem minDeltaExp_eq_mul_at_bad_prime {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange (p.adicCompletion L)).IsElliptic]
    [(E'.baseChange (p.adicCompletion L)).IsElliptic]
    (hssE : SemistableAt p E) (hssE' : SemistableAt p E') (hjneg : jExp p E < 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hcop : ¬ ((l : ℤ) ∣ jExp p E))
    (hmin : WeierstrassCurve.IsMinimal (p.adicCompletionIntegers L)
      (E.baseChange (p.adicCompletion L)))
    (h : (E.baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    (hlu : IsUnit ((l : (p.adicCompletionIntegers L))))
    (hvw : ∀ ζ : (p.adicCompletionIntegers L), IsPrimitiveRoot ζ l →
      ∃ v w : (p.adicCompletionIntegers L),
      v = ∑ i ∈ (range l).erase 0,
          veluV2 (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
            (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
            (tateXpair (ζ ^ i)
              ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
              (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
            (tateYpair (ζ ^ i)
              ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
              (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
        ∧ 2 * w = ∑ i ∈ (range l).erase 0,
          (veluU (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
              (tateXpair (ζ ^ i)
                ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                (tateParamR (E.baseChange (p.adicCompletion L)) h)
                (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
              (tateYpair (ζ ^ i)
                ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                (tateParamR (E.baseChange (p.adicCompletion L)) h)
                (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
            + 2 * (veluV2 (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
                    (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
                    (tateXpair (ζ ^ i)
                      ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR (E.baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
                    (tateYpair (ζ ^ i)
                      ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR (E.baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
                  * tateXpair (ζ ^ i)
                      ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR (E.baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)))
        ∧ ((veluCurve (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)) v w).map
            (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).IsElliptic)
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  haveI := hmin
  have hp := ABC3.Found.GaloisRep.valuation_algebraMap_adicCompletion L p
  have hql : (tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l
      ∈ IsLocalRing.maximalIdeal (p.adicCompletionIntegers L) :=
    pow_mem_of_mem_ideal (tateParamR_mem (E.baseChange (p.adicCompletion L)) h) hl.pos
  -- ★`E_{q^l} ⊗ Lv` の楕円性——`c₄` は定数項 `1`、`1/j` は `q^l · 単元`
  have hqlne : algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)
      ((tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) ≠ 0 :=
    (map_ne_zero_iff _ (IsFractionRing.injective _ _)).2
      (pow_ne_zero l (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h))
  have hc4T' : algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)
      (tateCurveAt ((tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) hql).c₄ ≠ 0 :=
    ((tateCurveAt_c4_isUnit ((tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) hql).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).ne_zero
  obtain ⟨u', hu', hueq'⟩ := evalAdic_tateJinvSeries_eq_mul_unit
    (I := IsLocalRing.maximalIdeal (p.adicCompletionIntegers L))
    ((tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) hql
  have hev' : algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)
      (evalAdic tateJinvSeries
        ((tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) hql) ≠ 0 := by
    rw [hueq', map_mul]
    exact mul_ne_zero hqlne
      ((hu'.map (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).ne_zero)
  haveI : ((tateCurveAt ((tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) hql).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).IsElliptic :=
    tateCurveAt_map_isElliptic _ hql hev' hc4T'
  exact minDeltaExp_eq_mul_of_globalVelu' p E E' h hp hssE hssE' hjneg hl hodd
    (tateModel_map_Delta_ne_zero (E.baseChange (p.adicCompletion L)) h)
    (not_dvd_vAdd_tateParam_of_not_dvd_jExp p hp E h hjneg hcop)
    hlu hql (two_ne_zero_adicCompletionIntegers L p) (two_ne_zero_adicCompletion L p)
    hvw hQ hE'

def minDeltaExp_eq_mul_at_bad_prime.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点での Δ_min の関係——局所データを自前で作る形)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_mul_at_bad_prime.needs : List ProofObligation :=
  [ .citation "[ABC3]" "minDeltaExp_eq_mul_of_globalVelu′(第 999、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_of_globalVelu'") 1,
    .citation "[ABC3]" "not_dvd_vAdd_tateParam_of_not_dvd_jExp(hcop の言い換え、第 978、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.not_dvd_vAdd_tateParam_of_not_dvd_jExp") 1 ]

/-! ## ★★★★★★★★★★★★★★★★★★★★第 1001-1002 —— 極小化した対に当てる

★第 973／976 が与えるのは `E` ではなく **`C • E`** についての極小性・乗法還元である
（`C` は `p` で極小にする変数変換、第 954 が `SemistableAt` から取り出す）。
☆したがって第 1000 は `C • E` に当てることになる。

★そのとき Vélu の商も `C • E′` に移す必要があるが、
第 969（`veluQuotientFull_vcPoint_eq`）が**曲線の等式として**それを与える:

    `C • E′ = veluQuotientFull (C • E) (vcPoint C E Q の生成する集合)`

☆`Δ_min` も `jExp` も変数変換で不変（`minDeltaExp_variableChange`・`jExp_variableChange`）、
`SemistableAt` も不変（第 1001）、点の位数も不変（`addOrderOf_vcPoint`）。
★★したがって**そのまま `E`・`E′` の主張に戻る**。 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★**[GenEll] 悪い素点での `Δ_min` の関係——
極小化した対に当てる形**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1002）**——第 1000 を `C • E` に当て、
第 969 で Vélu の商を `C • E′` に移し、変数変換の不変性で `E`・`E′` に戻す。 -/
theorem minDeltaExp_eq_mul_at_bad_prime_vc {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic]
    (hssE : SemistableAt p E) (hssE' : SemistableAt p E') (hjneg : jExp p E < 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hcop : ¬ ((l : ℤ) ∣ jExp p E))
    (C : WeierstrassCurve.VariableChange L)
    [((C • E).baseChange (p.adicCompletion L)).IsElliptic]
    [((C • E').baseChange (p.adicCompletion L)).IsElliptic]
    (hmin : WeierstrassCurve.IsMinimal (p.adicCompletionIntegers L)
      ((C • E).baseChange (p.adicCompletion L)))
    (h : ((C • E).baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    (hlu : IsUnit ((l : (p.adicCompletionIntegers L))))
    (hvw : ∀ ζ : (p.adicCompletionIntegers L), IsPrimitiveRoot ζ l →
      ∃ v w : (p.adicCompletionIntegers L),
      v = ∑ i ∈ (range l).erase 0,
          veluV2 (tateCurveAt (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
            (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h))
            (tateXpair (ζ ^ i)
              ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
              (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
              (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h))
            (tateYpair (ζ ^ i)
              ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
              (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
              (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h))
        ∧ 2 * w = ∑ i ∈ (range l).erase 0,
          (veluU (tateCurveAt (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
              (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h))
              (tateXpair (ζ ^ i)
                ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
                (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h))
              (tateYpair (ζ ^ i)
                ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
                (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h))
            + 2 * (veluV2 (tateCurveAt (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
                    (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h))
                    (tateXpair (ζ ^ i)
                      ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h))
                    (tateYpair (ζ ^ i)
                      ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h))
                  * tateXpair (ζ ^ i)
                      ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h)))
        ∧ ((veluCurve (tateCurveAt (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
              (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h)) v w).map
            (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).IsElliptic)
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  haveI hCE : (C • E).IsElliptic := by
    rw [WeierstrassCurve.isElliptic_iff, WeierstrassCurve.variableChange_Δ]
    exact ((C.u⁻¹).isUnit.pow 12).mul E.isUnit_Δ
  haveI hCE' : (C • E').IsElliptic := by
    rw [WeierstrassCurve.isElliptic_iff, WeierstrassCurve.variableChange_Δ]
    exact ((C.u⁻¹).isUnit.pow 12).mul E'.isUnit_Δ
  have h2L : (2 : L) ≠ 0 := two_ne_zero
  have hQ' : addOrderOf (ABC3.Found.GenEll.vcPoint C E Q) = l := by
    rw [ABC3.Found.GenEll.addOrderOf_vcPoint C E Q]; exact hQ
  have hEq := ABC3.Found.GenEll.veluQuotientFull_vcPoint_eq C E E' hQ h2L hE'
  have hjC : jExp p (C • E) < 0 := by rw [jExp_variableChange p E C]; exact hjneg
  have hcopC : ¬ ((l : ℤ) ∣ jExp p (C • E)) := by
    rw [jExp_variableChange p E C]; exact hcop
  have hkey := minDeltaExp_eq_mul_at_bad_prime p (C • E) (C • E')
    (semistableAt_variableChange p E C hssE) (semistableAt_variableChange p E' C hssE')
    hjC hl hodd hcopC hmin h hlu hvw hQ' hEq
  rwa [minDeltaExp_variableChange p E' C, minDeltaExp_variableChange p E C] at hkey

def minDeltaExp_eq_mul_at_bad_prime_vc.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点での Δ_min の関係——極小化した対に当てる形)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_mul_at_bad_prime_vc.needs : List ProofObligation :=
  [ .citation "[ABC3]" "minDeltaExp_eq_mul_at_bad_prime(第 1000、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_at_bad_prime") 1,
    .citation "[ABC3]" "veluQuotientFull_vcPoint_eq(Vélu の商を変数変換で移す、第 969、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.veluQuotientFull_vcPoint_eq") 1,
    .citation "[ABC3]" "semistableAt_variableChange(半安定性の不変性、第 1001、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.semistableAt_variableChange") 1 ]

/-! ## ★★★★★★★★★★★★★★★★第 1003 —— `hvw` を作る

★第 1000／1002 の `hvw` は「原始 `l` 乗根 `ζ` ごとに Vélu の係数 `v`・`w` が `R` に取れ、
その `veluCurve` が楕円である」ことを要求する。

☆`v` は和そのものだから自明、`w` は第 961（`exists_veluW_mu`）——
和が **`i ↦ l−i` の対で偶数になる**（第 956-960）ことから来る。
★楕円性は第 962（`isElliptic_veluCurve_tate_map`）と `c₄`・`c₆` の関係式である。 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★**[GenEll] `μ_l` に対する Vélu の係数は `R` に取れて、
その商は楕円である**。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★★★**2026-09-01（第 1003）**——第 1000／1002 の `hvw` の中身である。 -/
theorem exists_vw_tate_mu {R : Type} [CommRing R] [IsDomain R] [CharZero R]
    [IsDiscreteValuationRing R] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [CharZero K] [Algebra R K] [IsFractionRing R K]
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (hq0 : q ≠ 0)
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2) (hlu : IsUnit ((l : R)))
    (hql : q ^ l ∈ IsLocalRing.maximalIdeal R) (h2 : (2 : R) ≠ 0)
    [((tateCurveAt (q ^ l) hql).map (algebraMap R K)).IsElliptic]
    (ζ : R) (hζ : IsPrimitiveRoot ζ l) :
    ∃ v w : R,
      v = ∑ i ∈ (range l).erase 0,
          veluV2 (tateCurveAt q hq)
            (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
      ∧ 2 * w = ∑ i ∈ (range l).erase 0,
          (veluU (tateCurveAt q hq)
              (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            + 2 * (veluV2 (tateCurveAt q hq)
                    (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                    (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                  * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
      ∧ ((veluCurve (tateCurveAt q hq) v w).map (algebraMap R K)).IsElliptic := by
  have hu := ABC3.Found.GenEll.isUnit_one_sub_pow_of_isUnit_natCast hl.pos hζ hlu
  -- ★`ζ` を `Kˣ` に上げる
  have hζ0 : ζ ≠ 0 := by
    intro hc
    have := hζ.pow_eq_one
    rw [hc, zero_pow hl.pos.ne'] at this
    exact zero_ne_one this
  have hne : algebraMap R K ζ ≠ 0 :=
    (map_ne_zero_iff _ (IsFractionRing.injective R K)).2 hζ0
  have hζu : algebraMap R K ζ = ((Units.mk0 (algebraMap R K ζ) hne : Kˣ) : K) := rfl
  have hζl : (Units.mk0 (algebraMap R K ζ) hne) ^ l = 1 := by
    ext
    rw [Units.val_pow_eq_pow_val, Units.val_mk0, ← map_pow, hζ.pow_eq_one, map_one,
      Units.val_one]
  have hord : ∀ n : ℕ, 0 < n → n < l →
      (Units.mk0 (algebraMap R K ζ) hne) ^ n ≠ 1 := by
    intro n hn hnl hcon
    have h1 : (algebraMap R K ζ) ^ n = 1 := by
      rw [← Units.val_mk0 hne, ← Units.val_pow_eq_pow_val, hcon, Units.val_one]
    have hz : ζ ^ n = 1 := by
      refine IsFractionRing.injective R K ?_
      rw [map_pow, h1, map_one]
    exact absurd (Nat.le_of_dvd hn ((hζ.pow_eq_one_iff_dvd n).1 hz)) (by omega)
  -- ☆`tateDXpair ≠ 0`
  have hDX : ∀ i ∈ (range l).erase 0,
      tateDXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq ≠ 0 := fun i hi =>
    tateDXpair_ne_zero_of_mu (mkTateSetup q hq hq0) hΔ
      (dvrTatePhiAddEquiv q hq hq0 hΔ) (fun _ => rfl) hl hodd ζ _ hζu hζl hord i
      (Finset.mem_erase.1 hi).1 (Finset.mem_range.1 (Finset.mem_erase.1 hi).2) (hu i hi)
  have h4 := c4_velu_tate hl hζ hlu hu q hq hql h2 hDX
  have h6 := c6_velu_tate hl hζ hu q hq hql h2 hDX
  -- ★`l` は奇素数なので `l = 2m+1`
  obtain ⟨m, rfl⟩ : ∃ m, l = 2 * m + 1 := hl.odd_of_ne_two hodd
  obtain ⟨w, hw0⟩ := ABC3.Found.GaloisRep.exists_veluW_mu (mkTateSetup q hq hq0) hΔ
    (dvrTatePhiAddEquiv q hq hq0 hΔ) (fun _ => rfl) m ζ _ hζu hζl hord hu
  -- ☆`(mkTateSetup q hq hq0).q` を `q` に言い直す
  have hw : 2 * w = ∑ i ∈ (range (2 * m + 1)).erase 0,
      (veluU (tateCurveAt q hq)
          (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (2 * m + 1 - 1)) q hq)
          (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (2 * m + 1 - 1)) q hq)
        + 2 * (veluV2 (tateCurveAt q hq)
                (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (2 * m + 1 - 1)) q hq)
                (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (2 * m + 1 - 1)) q hq)
              * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (2 * m + 1 - 1)) q hq)) := hw0
  refine ⟨_, w, rfl, hw, ?_⟩
  refine ABC3.Found.GaloisRep.isElliptic_veluCurve_tate_map q hq (2 * m + 1) hql _ w hlu
    h4 ?_ inferInstance
  rw [← h6]
  linear_combination (3024 : R) * hw

def exists_vw_tate_mu.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(μ_l に対する Vélu の係数は R に取れ、商は楕円である)",
    sectionId := "genell-lemma-3-2" }

def exists_vw_tate_mu.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_veluW_mu(w が取れる、第 961、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_veluW_mu") 1,
    .citation "[ABC3]" "isElliptic_veluCurve_tate_map(商は楕円、第 962、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isElliptic_veluCurve_tate_map") 1 ]

/-! ## ★★★★★★★★★★★★★★★★★★★★第 1004 —— `hvw` を内側で作る

★第 1003 が `hvw` の中身を与えるので、第 1002 から `hvw` を落とせる。
☆残る局所入力は **`hmin`・`h`・`hlu` の 3 本だけ**である。 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★**[GenEll] 悪い素点での `Δ_min` の関係——
残る入力は極小性・分裂性・`l` の単元性の 3 本**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1004）**——第 1003 で `hvw` を内側から供給する。 -/
theorem minDeltaExp_eq_mul_at_bad_prime_full {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic]
    (hssE : SemistableAt p E) (hssE' : SemistableAt p E') (hjneg : jExp p E < 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hcop : ¬ ((l : ℤ) ∣ jExp p E))
    (C : WeierstrassCurve.VariableChange L)
    [((C • E).baseChange (p.adicCompletion L)).IsElliptic]
    [((C • E').baseChange (p.adicCompletion L)).IsElliptic]
    (hmin : WeierstrassCurve.IsMinimal (p.adicCompletionIntegers L)
      ((C • E).baseChange (p.adicCompletion L)))
    (h : ((C • E).baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    (hlu : IsUnit ((l : (p.adicCompletionIntegers L))))
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  haveI := hmin
  have hql : (tateParamR ((C • E).baseChange (p.adicCompletion L)) h) ^ l
      ∈ IsLocalRing.maximalIdeal (p.adicCompletionIntegers L) :=
    pow_mem_of_mem_ideal
      (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h) hl.pos
  -- ★`E_{q^l} ⊗ Lv` の楕円性（第 1000 と同じ作り方）
  have hqlne : algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)
      ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) ^ l) ≠ 0 :=
    (map_ne_zero_iff _ (IsFractionRing.injective _ _)).2
      (pow_ne_zero l (tateParamR_ne_zero ((C • E).baseChange (p.adicCompletion L)) h))
  have hc4T' : algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)
      (tateCurveAt ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) ^ l)
        hql).c₄ ≠ 0 :=
    ((tateCurveAt_c4_isUnit
        ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) ^ l) hql).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).ne_zero
  obtain ⟨u', hu', hueq'⟩ := evalAdic_tateJinvSeries_eq_mul_unit
    (I := IsLocalRing.maximalIdeal (p.adicCompletionIntegers L))
    ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) ^ l) hql
  have hev' : algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)
      (evalAdic tateJinvSeries
        ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) ^ l) hql) ≠ 0 := by
    rw [hueq', map_mul]
    exact mul_ne_zero hqlne
      ((hu'.map (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).ne_zero)
  haveI : ((tateCurveAt ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) ^ l)
      hql).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).IsElliptic :=
    tateCurveAt_map_isElliptic _ hql hev' hc4T'
  exact minDeltaExp_eq_mul_at_bad_prime_vc p E E' hssE hssE' hjneg hl hodd hcop C hmin h hlu
    (fun ζ hζ => exists_vw_tate_mu _
      (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h)
      (tateParamR_ne_zero ((C • E).baseChange (p.adicCompletion L)) h)
      (tateModel_map_Delta_ne_zero ((C • E).baseChange (p.adicCompletion L)) h)
      hl hodd hlu hql (two_ne_zero_adicCompletionIntegers L p) ζ hζ)
    hQ hE'

def minDeltaExp_eq_mul_at_bad_prime_full.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点での Δ_min の関係——残る入力は極小性・分裂性・l の単元性)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_mul_at_bad_prime_full.needs : List ProofObligation :=
  [ .citation "[ABC3]" "minDeltaExp_eq_mul_at_bad_prime_vc(第 1002、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_at_bad_prime_vc") 1,
    .citation "[ABC3]" "exists_vw_tate_mu(hvw の中身、第 1003、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.exists_vw_tate_mu") 1 ]

end ABC3.Skeleton.GenEll
