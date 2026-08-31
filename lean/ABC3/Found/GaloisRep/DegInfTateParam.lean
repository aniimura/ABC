/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.DegInfLocal

/-!
# Galois (G6) 第 892 ブロック —— **★★★★★★★★`R` の側の母数で述べた `Δ_min` の関係**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★これは何か

`minDeltaExp_eq_mul_of_tateParam`（第 719）は仮説を **`Kˣ` の側**の母数
`tateParamK` で述べているが、葉 1 の連鎖（第 883・891）が与えるのは
**`R` の側**の `tateParamR` である。

★本ブロックはその隔たりを埋める——`tateParamK` は `tateParamR` の像なので
`Units.ext` と `map_pow` だけで済む。

| 定理 | 内容 |
|---|---|
| `tateParamK_pow` | ★`q_{E′} = q_E^l`（`R`）⇒ 同じ（`K`） |
| `minDeltaExp_eq_mul_of_tateParamR` | ★★★★★★★★**`R` の側で述べた `Δ_min` の関係** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField
open scoped Classical

variable {L : Type} [Field L] [NumberField L]
  {Lv : Type} [Field Lv] [Algebra L Lv]
  {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [Algebra R Lv] [IsFractionRing R Lv]

/-- ★**`R` の側の関係は `Kˣ` の側の関係を与える**。 -/
theorem tateParamK_pow [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (W W' : WeierstrassCurve Lv) [W.IsElliptic] [W.IsMinimal R]
    [W'.IsElliptic] [W'.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R) (h' : W'.HasSplitMultiplicativeReduction R)
    (l : ℕ) (hq : tateParamR W' h' = (tateParamR W h) ^ l) :
    tateParamK W' h' = (tateParamK W h) ^ l := by
  refine Units.ext ?_
  rw [Units.val_pow_eq_pow_val, tateParamK, tateParamK, Units.val_mk0, Units.val_mk0,
    ← map_pow, hq]

/-- ★★★★★★★★**`R` の側の母数で述べた `Δ_min` の関係**。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★これが葉 1（`tateParam_quot_velu`、第 891）の結論を
`Lemma 3.5` の入力 `hloc` に直接繋ぐ形である。 -/
theorem minDeltaExp_eq_mul_of_tateParamR [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (E E' : WeierstrassCurve L) (l : ℕ)
    [(E.baseChange Lv).IsElliptic] [(E.baseChange Lv).IsMinimal R]
    [(E'.baseChange Lv).IsElliptic] [(E'.baseChange Lv).IsMinimal R]
    (h : (E.baseChange Lv).HasSplitMultiplicativeReduction R)
    (h' : (E'.baseChange Lv).HasSplitMultiplicativeReduction R)
    (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = (HeightOneSpectrum.valuation L p) x)
    (C : VariableChange L) (hC : IsMinimal (primeSubring p) (C • E))
    (hc4ne : (C • E).c₄ ≠ 0) (hc4 : valAdd p (Units.mk0 ((C • E).c₄) hc4ne) = 0)
    (C' : VariableChange L) (hC' : IsMinimal (primeSubring p) (C' • E'))
    (hc4ne' : (C' • E').c₄ ≠ 0) (hc4' : valAdd p (Units.mk0 ((C' • E').c₄) hc4ne') = 0)
    (hq : tateParamR (E'.baseChange Lv) h' = (tateParamR (E.baseChange Lv) h) ^ l) :
    minDeltaExp p E' = l * minDeltaExp p E :=
  minDeltaExp_eq_mul_of_tateParam E E' l h h' p hp C hC hc4ne hc4 C' hC' hc4ne' hc4'
    (tateParamK_pow (R := R) _ _ h h' l hq)

/-! ## ★出典の紐付け(`.src`) -/

def tateParamK_pow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(R の側の q_{E′} = q_E^l から Kˣ の側へ。★無条件)",
    sectionId := "genell-lemma-3-2" }

def minDeltaExp_eq_mul_of_tateParamR.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(R の側の母数で述べた Δ_min の関係。★無条件)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GaloisRep
