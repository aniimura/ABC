import ABC3.Found.Arakelov.ArcMetricExists

/-!
# Arakelov (C3) 第 286 ブロック —— **定数倍と切断のノルム**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★Interface の欄をそのまま埋めにいく

| Interface の欄 | 本ブロック |
|---|---|
| `scale` | ★`scaleMetric` |
| `normSection` | ★`normSectionOf` |
| `normSection_nonneg` | ★`normSectionOf_nonneg` |
| `normSection_eq_zero_iff` | ★★`normSectionOf_eq_zero_iff` |
| `normSection_scale` | ★★`normSectionOf_scale` |
| `normSection_continuous` | ★★★`normSectionOf_continuous`(第 283 の `cont` そのもの) |

★★★`scale` を `exp (-c)` 倍と定めるのは Interface の `normSection_scale` の要求である
——`logMetric` を `+c` 平行移動させるための符号。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{0}} {F : X.Modules}

/-- ★**計量の定数倍**——ノルムを一様に `exp (-c)` 倍する。 -/
noncomputable def scaleMetric (c : ℝ) (m : ContArcMetric X F) : ContArcMetric X F where
  nrm := fun p w => Real.exp (-c) * m.nrm p w
  nonneg := fun p w => mul_nonneg (Real.exp_nonneg _) (m.nonneg p w)
  eq_zero_iff := fun p w => by
    rw [mul_eq_zero]
    exact ⟨fun h => h.elim (fun h => absurd h (Real.exp_ne_zero _))
        (fun h => (m.eq_zero_iff p w).1 h),
      fun h => Or.inr ((m.eq_zero_iff p w).2 h)⟩
  smul := fun p c' w => by rw [m.smul]; ring
  cont := fun s => by
    letI := arcTopology X
    exact continuous_const.mul (m.cont s)

@[simp] theorem scaleMetric_nrm (c : ℝ) (m : ContArcMetric X F)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (w : ↥(arcFiber p F)) :
    (scaleMetric c m).nrm p w = Real.exp (-c) * m.nrm p w := rfl

/-- ★**切断のノルム** `|s|` ——各点で切断を評価してノルムを取る。 -/
noncomputable def normSectionOf (m : ContArcMetric X F) (s : (F.val.obj (op ⊤) : Type))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) : ℝ :=
  m.nrm p (arcEval p F s)

theorem normSectionOf_nonneg (m : ContArcMetric X F) (s : (F.val.obj (op ⊤) : Type))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) : 0 ≤ normSectionOf m s p :=
  m.nonneg p _

/-- ★★**ノルムが 0 になるのは、切断がその点で消えるときだけ**
——Interface の `normSection_eq_zero_iff` そのもの。 -/
theorem normSectionOf_eq_zero_iff (m : ContArcMetric X F) (s : (F.val.obj (op ⊤) : Type))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    normSectionOf m s p = 0 ↔
      (((Scheme.Modules.pullbackPushforwardAdjunction p).unit.app F).val.app (op ⊤)).hom s = 0 :=
  m.eq_zero_iff p _

/-- ★★定数倍はノルムを一様に `exp (-c)` 倍する。 -/
theorem normSectionOf_scale (c : ℝ) (m : ContArcMetric X F) (s : (F.val.obj (op ⊤) : Type))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    normSectionOf (scaleMetric c m) s p = Real.exp (-c) * normSectionOf m s p := rfl

/-- ★★★**切断のノルムは連続である**——第 283 の `cont` そのもの。 -/
theorem normSectionOf_continuous (m : ContArcMetric X F) (s : (F.val.obj (op ⊤) : Type)) :
    @Continuous _ ℝ (arcTopology X) _ (normSectionOf m s) :=
  m.cont s

/-! ## ★出典の紐付け(`.src`) -/

def normSectionOf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C3——切断のノルムと計量の定数倍)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
