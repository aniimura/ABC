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

end ABC3.Skeleton.GenEll
