import ABC3.Found.Arakelov.ArcRatio

/-!
# Arakelov (C3) 第 290 ブロック —— **★★★★★Green 関数と平行移動則**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★`logMetric` は絶対量ではない——基準を明示する

Interface の `logMetric : Metric X L → Arc X → ℝ` を**絶対的な関数として定義することはできない**。
★Green 関数は基準計量 `m₀` に対する比の対数だからである。

★★本ブロックでは**基準を引数に出した** `logMetricOf m₀ m` を作る。
★★★witness を組むときは、`Nonempty` が Prop であることを使って
**型だけで決まる基準**(`Classical.choice`)を取ればよい。

## ★★平行移動則が退化を殺す

    logMetricOf m₀ (scaleMetric c m) p = logMetricOf m₀ m p + c

★`scale` をノルムの `exp (-c)` 倍と定めたので、`-log` は `+c` 平行移動になる。
★★★`logMetric ≡ 0` だと `0 = 0 + c` を全 `c` で要求して矛盾する——これが退化を殺す欄である。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `someNonzero` | ★各点で 0 でないベクトルを 1 本選ぶ |
| `gRatio` | ★★基準に対する比 |
| `gRatio_pos` | ★比は正 |
| `gRatio_indep` | ★★★選び方に依らない(第 289) |
| `logMetricOf` | ★★★★Green 関数 |
| `logMetricOf_scale` | ★★★★★**平行移動則** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite

variable {X : Scheme.{0}} {F : X.Modules} (hF : IsLocallyTrivial X F.val)

/-- ★**各複素点で 0 でないベクトルを 1 本選ぶ**。 -/
noncomputable def someNonzero (p : Spec (CommRingCat.of ℂ) ⟶ X) : ↥(arcFiber p F) :=
  Classical.choose (exists_ne_zero_arcFiber F hF p)

theorem someNonzero_ne_zero (p : Spec (CommRingCat.of ℂ) ⟶ X) : someNonzero hF p ≠ 0 :=
  Classical.choose_spec (exists_ne_zero_arcFiber F hF p)

/-- ★★**基準計量 `m'` に対する `m` の比**。 -/
noncomputable def gRatio (m m' : ContArcMetric X F) (p : Spec (CommRingCat.of ℂ) ⟶ X) : ℝ :=
  fiberRatio m m' p (someNonzero hF p)

/-- ★★★比はベクトルの選び方に依らない——第 289 そのもの。 -/
theorem gRatio_indep (m m' : ContArcMetric X F) (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (w : ↥(arcFiber p F)) (hw : w ≠ 0) :
    gRatio hF m m' p = fiberRatio m m' p w :=
  fiberRatio_indep hF m m' p _ w (someNonzero_ne_zero hF p) hw

/-- ★比は正である。 -/
theorem gRatio_pos (m m' : ContArcMetric X F) (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    0 < gRatio hF m m' p := by
  have hw := someNonzero_ne_zero hF p
  have h1 : 0 < m.nrm p (someNonzero hF p) :=
    lt_of_le_of_ne (m.nonneg _ _) (fun h => hw ((m.eq_zero_iff _ _).1 h.symm))
  have h2 : 0 < m'.nrm p (someNonzero hF p) :=
    lt_of_le_of_ne (m'.nonneg _ _) (fun h => hw ((m'.eq_zero_iff _ _).1 h.symm))
  exact div_pos h1 h2

/-- ★★★★**Green 関数**——基準 `m₀` に対する `m` の対数比。 -/
noncomputable def logMetricOf (m₀ m : ContArcMetric X F) (p : Spec (CommRingCat.of ℂ) ⟶ X) : ℝ :=
  - Real.log (gRatio hF m m₀ p)

/-- ★★★★★**平行移動則**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが Interface の `logMetric_scale` であり、**退化を殺す欄**である。 -/
theorem logMetricOf_scale (m₀ m : ContArcMetric X F) (c : ℝ)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    logMetricOf hF m₀ (scaleMetric c m) p = logMetricOf hF m₀ m p + c := by
  have hr : gRatio hF (scaleMetric c m) m₀ p = Real.exp (-c) * gRatio hF m m₀ p := by
    show (Real.exp (-c) * m.nrm p (someNonzero hF p)) / m₀.nrm p (someNonzero hF p)
      = Real.exp (-c) * (m.nrm p (someNonzero hF p) / m₀.nrm p (someNonzero hF p))
    rw [mul_div_assoc]
  show - Real.log (gRatio hF (scaleMetric c m) m₀ p) = - Real.log (gRatio hF m m₀ p) + c
  rw [hr, Real.log_mul (Real.exp_ne_zero _) (ne_of_gt (gRatio_pos hF m m₀ p)), Real.log_exp]
  ring

/-! ## ★出典の紐付け(`.src`) -/

def logMetricOf_scale.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C3——Green 関数と平行移動則)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
