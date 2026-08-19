import ABC3.Found.Arakelov.ArcTorsorMetric
import ABC3.Found.Arakelov.PicWitness
import ABC3.Found.Arakelov.ArcSpaceInterface

/-!
# Arakelov (C3) 第 295 ブロック —— **★★★★★★★`HermitianMetricData` の非空虚 witness**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★これが Arakelov 理論の 9 件のうち **C3** である

| 欄 | 出所 |
|---|---|
| `Metric` | ★第 294 `TorsorMetric` |
| `metric_nonempty` | ★★第 285 `exists_contArcMetric`(コンパクト Hausdorff) |
| `logMetric` + 3 法則 | ★★★第 294(構成から) |
| `normSection` + 4 法則 | ★★★★第 283 の連続性 + 第 294 |
| `IsConjCompatible` | ★第 293 |
| `tensorMetric` + `logMetric_tensor` | ★★第 294(Green の和) |

## ★★可逆層は局所自明である(実測)

`IsInvertibleSheaf` の第 2 連言は **`IsLocallyTrivial` そのもの**であった
——`h.2` でそのまま通る(実測 2026-08-19)。
-/

namespace ABC3.Interface.Arakelov

open AlgebraicGeometry CategoryTheory ABC3.Found.Arakelov ABC3.Found.GenEll

/-- ★**可逆層は局所自明である**。 -/
theorem isLocallyTrivial_of_isInvertibleSheaf {X : Scheme.{0}} (F : X.Modules)
    (h : IsInvertibleSheaf F) : IsLocallyTrivial X F.val := h.2

/-- ★`Pic` の元の下にある層は局所自明である。 -/
theorem picSheaf_locallyTrivial (X : Scheme.{0}) (L : picardDataWitness.Pic X) :
    IsLocallyTrivial X (picardDataWitness.sheafOf X L).val :=
  isLocallyTrivial_of_isInvertibleSheaf _ (picardDataWitness.sheafOf_invertible X L)

/-- ★★★★★★★**`HermitianMetricData` の実装**。 -/
noncomputable def hermitianMetricDataWitness : HermitianMetricData where
  toPicardData := picardDataWitness
  toArcSpaceData := arcSpaceDataImpl
  Metric := fun X L => TorsorMetric X (picardDataWitness.sheafOf X L)
  metric_nonempty := fun X L hc ht =>
    nonempty_torsorMetric (picSheaf_locallyTrivial X L) (hasContMetrics_of_compact hc ht)
  logMetric := fun _ _ m p => m.green p
  logMetric_continuous := fun _ _ m => m.green_cont
  scale := fun _ _ c m => TorsorMetric.scale c m
  logMetric_scale := fun _ _ _ _ _ => rfl
  IsConjCompatible := fun _ _ m => m.IsConjCompatible
  isConjCompatible_iff := fun _ _ m => TorsorMetric.isConjCompatible_iff m
  normSection := fun _ _ m s p => m.norm s p
  normSection_nonneg := fun _ _ m s p => TorsorMetric.norm_nonneg m s p
  normSection_eq_zero_iff := fun _ _ m s p => TorsorMetric.norm_eq_zero_iff m s p
  normSection_scale := fun _ _ c m s p => TorsorMetric.norm_scale c m s p
  normSection_continuous := fun _ _ m s => TorsorMetric.norm_continuous m s
  tensorMetric := fun X _ _ mL mM =>
    TorsorMetric.tensor mL mM (picSheaf_locallyTrivial X _)
  logMetric_tensor := fun _ _ _ _ _ _ => rfl

/-- ★★★★★★★**`HermitianMetricData` は非空虚である**——C3 達成。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが Arakelov 理論の 9 件のうち **C3** である。 -/
theorem HermitianMetricData.nonvacuous : Nonempty HermitianMetricData :=
  ⟨hermitianMetricDataWitness⟩

/-! ## ★出典の紐付け(`.src`) -/

def hermitianMetricDataWitness.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——解析化と hermitian 計量)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Interface.Arakelov
