/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.C6DenomFree
import ABC3.Found.GaloisRep.TateVelu
import ABC3.Meta.Claim

/-!
# VeluTateMuField —— `[GenEll] Lemma 3.5` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GaloisRep

open ABC3.Found.GenEll Finset

/-! ## ★`c₄`・`c₆` が `l⁴`・`l⁶` 倍なら `Δ` は `l¹²` 倍 -/

/-- ★★★★**`c₄ ↦ l⁴c₄`・`c₆ ↦ l⁶c₆` なら `Δ ↦ l¹²Δ`**（`1728` が可逆な体で）。

☆`1728·Δ = c₄³ − c₆²`（mathlib の `WeierstrassCurve.c_relation`）から出る。 -/
theorem Delta_of_c4_c6 {K : Type} [Field K] [CharZero K] (A B : WeierstrassCurve K) (l : K)
    (h4 : A.c₄ = l ^ 4 * B.c₄) (h6 : A.c₆ = l ^ 6 * B.c₆) :
    A.Δ = l ^ 12 * B.Δ := by
  have hA := WeierstrassCurve.c_relation A
  have hB := WeierstrassCurve.c_relation B
  have h1728 : (1728 : K) ≠ 0 := by norm_num
  refine mul_left_cancel₀ h1728 ?_
  rw [hA, h4, h6]
  linear_combination (-(l ^ 12)) * hB

/-! ## ★★★★★★★★★★★★`p ∣ l` でも通る Vélu の係数 -/

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**[GenEll] `μ_l` に対する Vélu の係数は商体に取れて、`c₄`・`c₆` は `l⁴`・`l⁶` 倍になる**
（第 1131）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★**`IsUnit ((l : R))` を仮説に置いていない**——`p ∣ l` でも成り立つ。
☆代わりに `v`・`w` は `R` ではなく**商体 `K` に取る**（`l⁶`・`l⁸` で割る）。 -/
theorem exists_vw_tate_mu_field {R : Type} [CommRing R] [IsDomain R] [CharZero R]
    {I : Ideal R} [IsAdicComplete I R]
    {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (q : R) (hq : q ∈ I) (hql : q ^ l ∈ I) :
    ∃ v w : K,
      (veluCurve ((tateCurveAt q hq).map (algebraMap R K)) v w).c₄
          = (l : K) ^ 4 * ((tateCurveAt (q ^ l) hql).map (algebraMap R K)).c₄
        ∧ (veluCurve ((tateCurveAt q hq).map (algebraMap R K)) v w).c₆
          = (l : K) ^ 6 * ((tateCurveAt (q ^ l) hql).map (algebraMap R K)).c₆ := by
  have hinj : Function.Injective (algebraMap R K) := IsFractionRing.injective R K
  have hlK : (l : K) ≠ 0 := by
    intro h
    have hR : ((l : ℕ) : R) = 0 := hinj (by rw [map_natCast, h, map_zero])
    exact (Nat.cast_ne_zero.2 hl.pos.ne' : ((l : ℕ) : R) ≠ 0) hR
  have h2K : (2 : K) ≠ 0 := by
    intro h
    have hR : (2 : R) = 0 := hinj (by rw [map_ofNat, h, map_zero])
    exact two_ne_zero hR
  have h4 := c4_velu_tateDF hl hζ q hq hql
  have h6 := c6_velu_tateDF hl hζ q hq hql
  rw [c6DFlhs] at h6
  set SV : R := ∑ i ∈ (range l).erase 0,
      veluV2DF l (tateCurveAt q hq)
        (tateXpairDF l (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        (tateYpairDF l (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq) with hSVdef
  set SW : R := ∑ i ∈ (range l).erase 0,
      ((l : R) ^ 2 * veluUDF l (tateCurveAt q hq)
            (tateXpairDF l (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            (tateYpairDF l (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        + 2 * (veluV2DF l (tateCurveAt q hq)
                (tateXpairDF l (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                (tateYpairDF l (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              * tateXpairDF l (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)) with hSWdef
  have h4K := congrArg (algebraMap R K) h4
  have h6K := congrArg (algebraMap R K) h6
  simp only [map_add, map_mul, map_pow, map_natCast, map_ofNat] at h4K h6K
  refine ⟨algebraMap R K SV / (l : K) ^ 6,
    algebraMap R K SW / (2 * (l : K) ^ 8), ?_, ?_⟩
  · rw [veluCurve_c₄, WeierstrassCurve.map_c₄, WeierstrassCurve.map_c₄]
    field_simp
    linear_combination h4K
  · rw [veluCurve_c₆, WeierstrassCurve.map_c₆, WeierstrassCurve.map_c₆,
      WeierstrassCurve.map_b₂, tateCurveAt_b₂, map_one]
    field_simp
    linear_combination 2 * h6K

/-- ★★★★★★★★★★★★★★**判別式の形**——`Δ(E_q/H) = l¹²·Δ(E_{q^l})`（第 1131）。

★`12·v_p(l)` のずれは「`u = 1/l` の変数変換」そのものであり、
**極小化でほどける**。☆したがって `Δ_min` は変わらない。 -/
theorem exists_vw_tate_mu_field_Delta {R : Type} [CommRing R] [IsDomain R] [CharZero R]
    {I : Ideal R} [IsAdicComplete I R]
    {K : Type} [Field K] [CharZero K] [Algebra R K] [IsFractionRing R K]
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (q : R) (hq : q ∈ I) (hql : q ^ l ∈ I) :
    ∃ v w : K,
      (veluCurve ((tateCurveAt q hq).map (algebraMap R K)) v w).Δ
        = (l : K) ^ 12 * ((tateCurveAt (q ^ l) hql).map (algebraMap R K)).Δ := by
  obtain ⟨v, w, h4, h6⟩ := exists_vw_tate_mu_field (K := K) hl hodd hζ q hq hql
  exact ⟨v, w, Delta_of_c4_c6 _ _ _ h4 h6⟩

/-! ## ★出典の紐付け(`.src`) -/

def exists_vw_tate_mu_field.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(p ∣ l でも Vélu の商は商体の上で E_{q^l} の l 倍スケールになる)",
    sectionId := "genell-lemma-3-5" }

def Delta_of_c4_c6.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(c₄・c₆ が l⁴・l⁶ 倍なら Δ は l¹² 倍。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_vw_tate_mu_field.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "c4_velu_tateDF(c₄ の分母なし版、第 1129、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.c4_velu_tateDF") 1,
    .citation "[ABC3]" "c6_velu_tateDF(c₆ の分母なし版、第 1130、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.c6_velu_tateDF") 1,
    .citation "[ABC3]" "veluCurve_c₄ / veluCurve_c₆(Vélu の不変量の変化、在庫)"
      (.inProject "ABC3" "ABC3.Found.GenEll.veluCurve_c₄") 1,
    .citation "[mathlib]" "WeierstrassCurve.c_relation(1728Δ = c₄³ − c₆²)"
      (.inMathlib "WeierstrassCurve.c_relation") 1,
    .implicitStep
      ("★★★★**2026-09-01（第 1131）の設計**——`exists_vw_tate_mu`（第 1003）は " ++
       "`hlu : IsUnit ((l : R))` と「`v`・`w` が `R` に取れる」を要求する。" ++
       "☆`p ∣ l` ではどちらも成り立たないが、**それは障害ではない**: " ++
       "商体で `l⁶`・`l⁸` で割れば `c₄ ↦ l⁴c₄`・`c₆ ↦ l⁶c₆` になり、" ++
       "`Δ ↦ l¹²Δ` すなわち `u = 1/l` の変数変換のぶんだけずれる。" ++
       "★`Δ_min` は極小モデルの判別式なのでそのずれは消える。" ++
       "☆残るのは `minDeltaExp_eq_mul_at_bad_prime`（第 1000）の `hvw` の口を " ++
       "「`R` の `v`・`w`」から「`K` の `v`・`w` ＋ `c₄`・`c₆` の関係」に" ++
       "付け替える配管である。") 10 ]

end ABC3.Found.GaloisRep
