import ABC3.Found.GaloisRep.TateModel
import ABC3.Interface.GaloisRep.Reduction

/-!
# Galois (G6) 第 310 ブロック —— **★★★★★★★★★★★G6 達成 `TateCurveData`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★★到達点

> **`TateCurveData` は非空虚である**(`TateCurveData.nonvacuous`)——**G6 達成**

★★★これで Galois 表現論の 8 件のうち **5 件**(G1, G2, G3, G4, G6)が埋まった。

## ★★★★★★組み立ての 4 本

| 欄 | 中身 | 出どころ |
|---|---|---|
| `tateParam` | `q ∈ Kˣ` | 第 309 の `exists_tate_model` から選択 |
| `uniformization` | `W(K) ≅ Kˣ/q^ℤ` | 第 303(点への作用)+ 第 301(一意化) |
| `localHeight` | `v(q)` の自然数化 | 第 302 の橋 |
| `localHeight_eq`・`_pos` | 定義の同定と正値性 | `q ∈ 𝔪`、`q ≠ 0` |

## ★★★★★★★一意化の鎖

    W.Point  ≃+  ((C の像) • W).Point        第 303 `vcPointAddEquiv`
             =   (E_q の K への像).Point       第 309 `exists_tate_model` を `map` した
             ≃+  Additive (Kˣ/q^ℤ)            第 301 `tate_uniformization_dvr`

★★★★**mathlib の `map_variableChange`(`C.map φ • W.map φ = (C • W).map φ`)が
ちょうど継ぎ目になっている**。★整モデルの底変換が `W` に戻ること
(`baseChange_integralModel_eq`)も mathlib に在った。

## ★★界面の 2 つの訂正を経ている(記録)

本witness が通るのは、第 302(付値を正規化付値に固定)と第 304(`Δ ≠ 0` を要求)の
**2 つの界面訂正**を先に済ませてあるからである。
★★★訂正前の `TateCurveData` は**充足不能**だった——実装は仕様の検算でもある。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tateParamR`・`tateParamK` | ★★★Tate 母数の選択 |
| `uniformization_of_split` | ★★★★★★★★**`W(K) ≅ Kˣ/q^ℤ`** |
| `localHeightOf`・`localHeight_eq_of_split` | ★★★★局所高さ |
| `TateCurveData.nonvacuous` | ★★★★★★★★★★★**G6 達成** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve ABC3.Interface.GaloisRep IsDedekindDomain

section Param

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
  {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

/-- ★★★**Tate 母数**(`R` の元)——第 309 の存在から選ぶ。 -/
noncomputable def tateParamR (W : WeierstrassCurve K) [hell : W.IsElliptic] [W.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R) : R :=
  haveI : W.HasMultiplicativeReduction R := h.toHasMultiplicativeReduction
  (exists_tate_model W h hell.isUnit.ne_zero).choose

theorem tateParamR_spec (W : WeierstrassCurve K) [hell : W.IsElliptic] [W.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R) :
    ∃ (hq : tateParamR W h ∈ IsLocalRing.maximalIdeal R) (C : VariableChange R),
      tateParamR W h ≠ 0 ∧ C • integralModel R W = tateCurveAt (tateParamR W h) hq := by
  haveI : W.HasMultiplicativeReduction R := h.toHasMultiplicativeReduction
  exact (exists_tate_model W h hell.isUnit.ne_zero).choose_spec

theorem tateParamR_mem (W : WeierstrassCurve K) [W.IsElliptic] [W.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R) :
    tateParamR W h ∈ IsLocalRing.maximalIdeal R :=
  (tateParamR_spec W h).choose

theorem tateParamR_ne_zero (W : WeierstrassCurve K) [W.IsElliptic] [W.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R) : tateParamR W h ≠ 0 := by
  obtain ⟨hq, C, hne, _⟩ := tateParamR_spec W h
  exact hne

/-- ★★★**Tate 母数**(`Kˣ` の元)。 -/
noncomputable def tateParamK (W : WeierstrassCurve K) [W.IsElliptic] [W.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R) : Kˣ :=
  Units.mk0 (algebraMap R K (tateParamR W h))
    ((map_ne_zero_iff _ (IsFractionRing.injective R K)).2 (tateParamR_ne_zero W h))

/-! ## ★★★★局所高さ -/

/-- ★★★**局所高さ** `v_K(q_E)`。 -/
noncomputable def localHeightOf (W : WeierstrassCurve K) [W.IsElliptic] [W.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R) : ℕ :=
  (vAdd (tateDvrVal R K) (tateParamK W h)).toNat

theorem vAdd_tateParamK_pos (W : WeierstrassCurve K) [W.IsElliptic] [W.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R) :
    0 < vAdd (tateDvrVal R K) (tateParamK W h) :=
  tateDvrVal_pos_of_mem_max _ ⟨tateParamR W h, tateParamR_mem W h, rfl⟩

/-- ★★★★局所高さは正(原文が `∈ ℤ_{>0}` と書いているもの)。 -/
theorem localHeight_pos_of_split (W : WeierstrassCurve K) [W.IsElliptic] [W.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R) : 0 < localHeightOf W h := by
  have hpos := vAdd_tateParamK_pos W h
  rw [localHeightOf]
  omega

/-- ★★★★局所高さは Tate 母数の付値そのもの。 -/
theorem localHeight_eq_of_split (W : WeierstrassCurve K) [W.IsElliptic] [W.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R) :
    (IsDiscreteValuationRing.maximalIdeal R).valuation K ((tateParamK W h : Kˣ) : K)
      = (Multiplicative.ofAdd (-(localHeightOf W h : ℤ)) : Multiplicative ℤ) := by
  have hpos := vAdd_tateParamK_pos W h
  have hcast : ((localHeightOf W h : ℤ)) = vAdd (tateDvrVal R K) (tateParamK W h) := by
    rw [localHeightOf]
    omega
  rw [valuation_eq_ofAdd_neg_vAdd, hcast]

end Param

/-! ## ★★★★★★★★一意化の鎖 -/

section Uniformization

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
  {K : Type} [Field K] [DecidableEq K] [Algebra R K] [IsFractionRing R K]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★**`W(K) ≅ Kˣ/q^ℤ`**——第 303・第 309・第 301 を継ぐ。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem uniformization_of_split (W : WeierstrassCurve K) [W.IsElliptic] [W.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R) :
    Nonempty (W.toAffine.Point ≃+ Additive (Kˣ ⧸ Subgroup.zpowers (tateParamK W h))) := by
  obtain ⟨hq, C, hqne, hCE⟩ := tateParamR_spec W h
  have hWmap : (integralModel R W).map (algebraMap R K) = W :=
    WeierstrassCurve.baseChange_integralModel_eq R W
  have hkey : (C.map (algebraMap R K)) • W
      = (tateCurveAt (tateParamR W h) hq).map (algebraMap R K) := by
    conv_lhs => rw [← hWmap]
    rw [WeierstrassCurve.map_variableChange, hCE]
  have e1 := vcPointAddEquiv W (C.map (algebraMap R K))
  obtain ⟨e2⟩ := tate_uniformization_dvr (K := K) hq hqne
  exact ⟨e1.trans (hkey ▸ e2)⟩

end Uniformization

/-! ## ★★★★★★★★★★★G6 達成 -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**`TateCurveData` の実装**。 -/
noncomputable def tateCurveDataWitness : TateCurveData where
  tateParam := by
    intro R _ _ _ _ K _ _ _ W _ _ h
    exact tateParamK W h
  uniformization := by
    intro R _ _ _ _ K _ _ _ _ W _ _ h
    exact uniformization_of_split W h
  localHeight := by
    intro R _ _ _ _ K _ _ _ W _ _ h
    exact localHeightOf W h
  localHeight_eq := by
    intro R _ _ _ _ K _ _ _ W _ _ h
    exact localHeight_eq_of_split W h
  localHeight_pos := by
    intro R _ _ _ _ K _ _ _ W _ _ h
    exact localHeight_pos_of_split W h

/-- ★★★★★★★★★★★**`TateCurveData` は非空虚である**——G6 達成。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★★★これが Galois 表現論の 8 件のうち **G6** である。 -/
theorem TateCurveData.nonvacuous : Nonempty TateCurveData := ⟨tateCurveDataWitness⟩

/-! ## ★出典の紐付け(`.src`) -/

def uniformization_of_split.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化 E(K) ≅ Kˣ/q^ℤ、分裂乗法還元から)",
    sectionId := "genell-def-3-3" }

def TateCurveData.nonvacuous.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(TateCurveData の非空虚性——G6)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
