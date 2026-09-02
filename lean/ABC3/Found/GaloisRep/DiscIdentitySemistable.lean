/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.DegInfBaseChange
import ABC3.Meta.Claim

/-!
# 第 1385 ブロック —— **判別式の恒等式から良い素点の半安定性が出る**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か——残る葉を**恒等式 1 本**に落とす

`Skeleton/GenEll/VeluSemistable.lean` の `semistableAt_veluQuotientFull` は
良い素点で `minDeltaExp p E′ = 0` を要求する。

★★★2026-09-02（第 1384）に**鍵になる恒等式**を見つけた:

    Δ(E)^l = Δ(E/C) · ( ∏_{P ∈ C∖{O}} (2 y_P + a₁ x_P + a₃) )^4

（`l = 3, 5, 7` について `tools/velu-disc-check.py` で厳密に確認。）

☆本ブロックはこの形の恒等式**だけ**を仮説に置いて、
良い素点での半安定性を**無条件に**導く。

★道は 4 行である:

1. `valAdd` は乗法的（`valAdd_mul`・`valAdd_pow`、在庫）なので
   `l·v(Δ(E)) = v(Δ(E′)) + 4·v(N)`
2. `E′` は `p` 上整なので `0 ≤ minDeltaExp p E′ ≤ v(Δ(E′))`
   （`minDeltaExp_nonneg`・`minDeltaExp_le_of_isIntegral`、在庫）
3. `N` も整なので `v(N) ≥ 0`（`valAdd_nonneg_iff`、在庫）
4. 良い素点では `v(Δ(E)) = 0` なので、和が `0` で両項が非負 ⟹ **両方 `0`**

★★★これで `minDeltaExp p E′ = 0`、すなわち `SemistableAt p E′` である。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField ABC3.Meta

open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-- ★★★★★★★★**整な曲線の判別式の付値は非負**（第 1385）。 -/
theorem valAdd_Delta_nonneg_of_isIntegral (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) (hΔ : W.Δ ≠ 0)
    (hint : WeierstrassCurve.IsIntegral (primeSubring p) W) :
    0 ≤ valAdd p (Units.mk0 W.Δ hΔ) := by
  have hone : ((1 : WeierstrassCurve.VariableChange L) • W) = W := one_smul _ _
  have hle := minDeltaExp_le_of_isIntegral p W hΔ (1 : WeierstrassCurve.VariableChange L)
    (by rwa [hone])
  have hEq : Units.mk0 (((1 : WeierstrassCurve.VariableChange L) • W).Δ)
      (variableChange_Delta_ne_zero W hΔ 1) = Units.mk0 W.Δ hΔ := by
    refine Units.ext ?_
    simp only [Units.val_mk0]
    rw [hone]
  rw [hEq] at hle
  exact le_trans (minDeltaExp_nonneg p W) hle

/-- ★★★★★★★★★★★★★★★★★★★★
**判別式の恒等式から良い素点の半安定性が出る**——★**無条件**（第 1385）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆仮説は「`Δ(E)^l = Δ(E′)·N^4` で `N` は `p` 上整」と
「`v_p(Δ(E)) = 0`（良い素点、極小モデル）」と「`E′` は `p` 上整」だけである。

★★★これで残る葉は**この恒等式 1 本**になる。 -/
theorem semistableAt_of_disc_pow_eq (p : HeightOneSpectrum (𝓞 L))
    (E E' : WeierstrassCurve L) (hΔ : E.Δ ≠ 0) (hΔ' : E'.Δ ≠ 0)
    (hint : WeierstrassCurve.IsIntegral (primeSubring p) E')
    (hgood : valAdd p (Units.mk0 E.Δ hΔ) = 0)
    {N : L} (hN : N ≠ 0) (hNint : N ∈ primeSubring p)
    {l : ℕ} (hid : E.Δ ^ l = E'.Δ * N ^ 4) :
    SemistableAt p E' := by
  -- ★段 1: `valAdd` は乗法的
  have hu : (Units.mk0 E.Δ hΔ) ^ l = (Units.mk0 E'.Δ hΔ') * (Units.mk0 N hN) ^ 4 := by
    refine Units.ext ?_
    simpa using hid
  have h1 := congrArg (valAdd p) hu
  rw [valAdd_pow, valAdd_mul, valAdd_pow, hgood, mul_zero] at h1
  -- ★段 2-3: 両項は非負
  have hNnn : 0 ≤ valAdd p (Units.mk0 N hN) :=
    (valAdd_nonneg_iff p _).2 ((mem_primeSubring_iff p N).1 hNint)
  have hD'nn : 0 ≤ valAdd p (Units.mk0 E'.Δ hΔ') :=
    valAdd_Delta_nonneg_of_isIntegral p E' hΔ' hint
  -- ★段 4: 和が `0` で両項が非負なら両方 `0`
  have hD'0 : valAdd p (Units.mk0 E'.Δ hΔ') = 0 := by omega
  -- ★`minDeltaExp p E′ = 0`
  have hone : ((1 : WeierstrassCurve.VariableChange L) • E') = E' := one_smul _ _
  have hle := minDeltaExp_le_of_isIntegral p E' hΔ' (1 : WeierstrassCurve.VariableChange L)
    (by rwa [hone])
  have hEq : Units.mk0 (((1 : WeierstrassCurve.VariableChange L) • E').Δ)
      (variableChange_Delta_ne_zero E' hΔ' 1) = Units.mk0 E'.Δ hΔ' := by
    refine Units.ext ?_
    simp only [Units.val_mk0]
    rw [hone]
  rw [hEq, hD'0] at hle
  exact Or.inl (le_antisymm hle (minDeltaExp_nonneg p E'))

/-! ## ★出典の紐付け(`.src`) -/

def valAdd_Delta_nonneg_of_isIntegral.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(整な曲線の判別式の付値は非負。★無条件)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_of_disc_pow_eq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(判別式の恒等式から良い素点の半安定性が出る。★無条件)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_of_disc_pow_eq.needs : List ProofObligation :=
  [ .citation "[ABC3]" "minDeltaExp_le_of_isIntegral(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.minDeltaExp_le_of_isIntegral") 1,
    .citation "[ABC3]" "valAdd_pow・valAdd_nonneg_iff(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.valAdd_pow") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1385）**——残る葉を**恒等式 1 本**に落とす段である。" ++
       "☆`Δ(E)^l = Δ(E/C)·(∏_{P∈C∖O}(2y_P+a₁x_P+a₃))^4` は " ++
       "`l = 3, 5, 7` について `tools/velu-disc-check.py` で厳密に確認した。" ++
       "★核でない ±閉集合では成り立たないので、部分群構造が本質である。") 17 ]

end ABC3.Found.GaloisRep
