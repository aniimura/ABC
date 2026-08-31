/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateParamJ
import ABC3.Found.GaloisRep.TateCurveNatural

/-!
# 第 944 ブロック —— **★★★★★★★★★★★★★★★★Tate 母数は体の拡大で変わらない**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★★★★★★★★★これは何か——**分岐指数の層を消す**

第 942 で測ったとおり、`Lemma 3.5` の組み立てには
**`Lv` の有限拡大 `L'`**が要る（`μ_l` や `l`-捉れ点が `Lv` にあるとは限らない）。

☆そこで第 943 の台帳には「(B) `Δ_min` が拡大でどう変わるか（分岐指数）」を
残る層として書いていた。★★**本ブロックはその層を不要にする**。

★見方を変える: 拡大の上で得た結論を `Δ_min` の段で降ろすのではなく、
**Tate 母数の段で降ろす**。Tate 母数は `R` の元であり、
`σ : R → R'` は単射だから、`σ(q_{E'}) = σ(q_E)^l` から
**その場で** `q_{E'} = q_E^l` が出る。
☆`Δ_min` はその後で `R` の中だけで計ればよいので、分岐指数は一切現れない。

★道は 2 段だけである:

| 段 | 中身 | 出どころ |
|---|---|---|
| 1 | `E_{σq} ⊗ K' = (E_q ⊗ K) ⊗ K'` | `tateCurveAt_map`（第 869）＋`map_map` |
| 2 | `j` が一致するので Tate 母数も一致 | `tateParamR_eq_of_j_tateCurveAt`（第 931） |

☆`j` の一致は `map_j`（体準同型と可換）と
`variableChange_j`（変数変換で不変）だけで出る。

| 定理 | 内容 |
|---|---|
| `tateCurveAt_map_comm` | ★Tate 曲線の底変換の可換図 |
| `j_tateModel_eq` | ★Tate モデルの `j` はもとの曲線の `j` |
| `tateParamR_map` | ★★★★★★★★★★★★★★★★**`q_{W ⊗ K'} = σ(q_W)`** |
| `tateParam_descend` | ★★★★★★★★**拡大の上の `q' = q^l` を `R` に降ろす** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsLocalRing

section ParamMap

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
  {K : Type} [Field K] [CharZero K] [Algebra R K] [IsFractionRing R K]
  {R' : Type} [CommRing R'] [IsDomain R'] [IsDiscreteValuationRing R']
  [IsAdicComplete (IsLocalRing.maximalIdeal R') R']
  {K' : Type} [Field K'] [CharZero K'] [Algebra R' K'] [IsFractionRing R' K']

/-! ## ★段 1——Tate 曲線の底変換の可換図 -/

/-- ★★**`E_{σq}` を `K'` に上げたものは、`E_q` を `K` に上げてから `K'` に送ったもの**。

☆四つの写像 `R → R' → K'` と `R → K → K'` が可換であればよい。 -/
theorem tateCurveAt_map_comm (σ : R →+* R') (φ : K →+* K')
    (hσI : ∀ x ∈ IsLocalRing.maximalIdeal R, σ x ∈ IsLocalRing.maximalIdeal R')
    (hcomm : ∀ x : R, algebraMap R' K' (σ x) = φ (algebraMap R K x))
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) :
    (tateCurveAt (σ q) (hσI q hq)).map (algebraMap R' K')
      = ((tateCurveAt q hq).map (algebraMap R K)).map φ := by
  rw [← tateCurveAt_map σ hσI q hq (hσI q hq), WeierstrassCurve.map_map,
    WeierstrassCurve.map_map]
  congr 1
  exact RingHom.ext hcomm

/-! ## ★段 2——Tate モデルの `j` -/

/-- ★★**Tate モデルを `K` に戻すと、変数変換を除いてもとの曲線**。

☆`tateParamR_spec` の `C • integralModel R W = tateCurveAt q hq` を
`K` に上げるだけである。 -/
theorem tateModel_baseChange (W : WeierstrassCurve K) [W.IsElliptic] [W.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R) {C : VariableChange R}
    {hq : tateParamR W h ∈ IsLocalRing.maximalIdeal R}
    (hCE : C • integralModel R W = tateCurveAt (tateParamR W h) hq) :
    (tateCurveAt (tateParamR W h) hq).map (algebraMap R K)
      = (C.map (algebraMap R K)) • W := by
  rw [← hCE, ← WeierstrassCurve.map_variableChange]
  congr 1
  exact WeierstrassCurve.baseChange_integralModel_eq R W

/-- ★★**Tate モデルの `j` はもとの曲線の `j`**。 -/
theorem j_tateModel_eq (W : WeierstrassCurve K) [W.IsElliptic] [W.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R) {C : VariableChange R}
    {hq : tateParamR W h ∈ IsLocalRing.maximalIdeal R}
    (hCE : C • integralModel R W = tateCurveAt (tateParamR W h) hq)
    [((tateCurveAt (tateParamR W h) hq).map (algebraMap R K)).IsElliptic] :
    ((tateCurveAt (tateParamR W h) hq).map (algebraMap R K)).j = W.j := by
  have hb := tateModel_baseChange (R := R) W h hCE
  rw [ABC3.Found.GenEll.j_congr_curve hb, WeierstrassCurve.variableChange_j]

/-! ## ★★★★★★★★★★★★★★★★Tate 母数は体の拡大で変わらない -/

/-- ★★★★★★★★★★★★★★★★**`q_{W ⊗ K'} = σ(q_W)`**——Tate 母数の自然性。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★★★**2026-09-01（第 944）**——これがあるので、`Lemma 3.5` の組み立てで
**分岐指数の層（第 943 の (B)）は不要になる**。
`μ_l` を含む拡大 `L'` の上で得た `q' = q^l` を、
`Δ_min` に直す前に **Tate 母数の段で** `R` に降ろせるからである。 -/
theorem tateParamR_map (σ : R →+* R') (φ : K →+* K')
    (hσI : ∀ x ∈ IsLocalRing.maximalIdeal R, σ x ∈ IsLocalRing.maximalIdeal R')
    (hcomm : ∀ x : R, algebraMap R' K' (σ x) = φ (algebraMap R K x))
    (W : WeierstrassCurve K) [W.IsElliptic] [W.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R)
    [(W.map φ).IsMinimal R']
    (h' : (W.map φ).HasSplitMultiplicativeReduction R') :
    tateParamR (W.map φ) h' = σ (tateParamR W h) := by
  obtain ⟨hq, C, hne, hCE⟩ := tateParamR_spec W h
  -- ★段 1: `E_{σq} ⊗ K'` は `(C • W) ⊗ K'` である
  have hEq : (tateCurveAt (σ (tateParamR W h)) (hσI _ hq)).map (algebraMap R' K')
      = ((C.map (algebraMap R K)) • W).map φ := by
    rw [tateCurveAt_map_comm σ φ hσI hcomm _ hq, tateModel_baseChange (R := R) W h hCE]
  haveI : ((tateCurveAt (σ (tateParamR W h)) (hσI _ hq)).map (algebraMap R' K')).IsElliptic := by
    rw [hEq]; infer_instance
  -- ★段 2: `j` が一致するので Tate 母数も一致
  refine tateParamR_eq_of_j_tateCurveAt (W.map φ) h' _ (hσI _ hq) ?_
  simp only [ABC3.Found.GenEll.j_congr_curve hEq, WeierstrassCurve.map_j,
    WeierstrassCurve.variableChange_j]

/-! ## ★★★★★★★★拡大の上の等式を `R` に降ろす -/

/-- ★★★★★★★★**`q_{W'} = q_W^l` を拡大から降ろす**。

☆`σ` が単射なら、拡大の上での等式はそのまま下の環の等式である。

★★★★**これが `Lemma 3.5` で実際に使う形**である——
`μ_l` を含む有限拡大 `L'` の上で `tateParam_quot_veluCurve_dvr`（第 927）を当て、
その結論を本定理で `R` に戻す。 -/
theorem tateParam_descend (σ : R →+* R') (φ : K →+* K')
    (hσ : Function.Injective σ)
    (hσI : ∀ x ∈ IsLocalRing.maximalIdeal R, σ x ∈ IsLocalRing.maximalIdeal R')
    (hcomm : ∀ x : R, algebraMap R' K' (σ x) = φ (algebraMap R K x))
    (W W' : WeierstrassCurve K) [W.IsElliptic] [W'.IsElliptic]
    [W.IsMinimal R] [W'.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R) (h' : W'.HasSplitMultiplicativeReduction R)
    [(W.map φ).IsMinimal R'] [(W'.map φ).IsMinimal R']
    (hL : (W.map φ).HasSplitMultiplicativeReduction R')
    (hL' : (W'.map φ).HasSplitMultiplicativeReduction R')
    (l : ℕ)
    (hquot : tateParamR (W'.map φ) hL' = (tateParamR (W.map φ) hL) ^ l) :
    tateParamR W' h' = (tateParamR W h) ^ l := by
  refine hσ ?_
  rw [map_pow]
  rw [← tateParamR_map σ φ hσI hcomm W h hL, ← tateParamR_map σ φ hσI hcomm W' h' hL']
  exact hquot

end ParamMap

/-! ## 出典 -/

def tateCurveAt_map_comm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 曲線の底変換の可換図)",
    sectionId := "genell-def-3-3" }

def tateModel_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate モデルを K に戻す)",
    sectionId := "genell-def-3-3" }

def j_tateModel_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate モデルの j はもとの曲線の j)",
    sectionId := "genell-def-3-3" }

def tateParamR_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.5(Tate 母数は体の拡大で変わらない)",
    sectionId := "genell-lemma-3-5" }

def tateParam_descend.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.5(拡大の上の q' = q^l を R に降ろす)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
