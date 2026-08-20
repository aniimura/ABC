import ABC3.Found.Arakelov.ArcPouPull

/-!
# Arakelov (C3) 第 292 ブロック —— **★★★★★★計量 = 基準計量 × exp(-Green)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★§9-335 の値段を測り直した——**設計を変えると 10 ブロック浮く**

`logMetric` の連続性を「任意の計量の比」として証明しようとすると、
局所枠に対する連続性が要り、**サイトの移送**(`V.toScheme` と `X` の開集合の圏の同値)まで
降りることになる。★測ったところ **10–15 ブロック**である。

★★★しかし数学的に正しい定式化は**別にある**:

    直線束の上の連続計量の全体は、C(X^arc, ℝ) 上の**捻れ集合(torsor)**である。

★したがって計量を「**基準計量 `base` と Green 関数 `green` の対**」として持てば、

| 欄 | 状態 |
|---|---|
| `logMetric_continuous` | ★★★**構成から自明**(`green` が連続) |
| `logMetric_scale` | ★★★**構成から自明**(`green + c`) |
| `logMetric_tensor` | ★★★**構成から自明**(`green` の和) |

となる。★★これは逃げではない——`|·| = base · exp(-green)` は Arakelov 理論の標準的な表示である。

## ★★退化しないこと

★`Check/Arakelov/MetricNondegenerate.lean` が落とした退化 witness は
`Metric X L := Green X`(`L` が現れない)であった。
★★本設計は `base : ContArcMetric X F` を**対の第 1 成分に持つ**ので、
`normSection` が `L` に依存する——退化 witness は通らない。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `GreenMetric` | ★★★基準計量と Green 関数の対 |
| `GreenMetric.norm` | ★切断のノルム |
| `norm_nonneg` / `_eq_zero_iff` / `_continuous` | ★★3 法則 |
| `GreenMetric.scale` | ★定数倍(Green を平行移動) |
| `logMetric_scale` / `norm_scale` | ★★平行移動則 |
| `nonempty_greenMetric` | ★★★★存在(第 285) |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{0}} {F : X.Modules}

/-- ★★★**計量 = 基準計量 × exp(-Green)**。 -/
structure GreenMetric (X : Scheme.{0}) (F : X.Modules) where
  /-- 基準となる連続計量。 -/
  base : ContArcMetric X F
  /-- Green 関数。 -/
  green : (Spec (CommRingCat.of ℂ) ⟶ X) → ℝ
  /-- ★Green 関数は連続である。 -/
  green_cont : @Continuous _ ℝ (arcTopology X) _ green

namespace GreenMetric

/-- ★**切断のノルム** `|s|`。 -/
noncomputable def norm (m : GreenMetric X F) (s : (F.val.obj (op ⊤) : Type))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) : ℝ :=
  Real.exp (-(m.green p)) * m.base.nrm p (arcEval p F s)

theorem norm_nonneg (m : GreenMetric X F) (s : (F.val.obj (op ⊤) : Type))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) : 0 ≤ m.norm s p :=
  mul_nonneg (Real.exp_nonneg _) (m.base.nonneg _ _)

/-- ★★ノルムが 0 になるのは切断がその点で消えるときだけ。 -/
theorem norm_eq_zero_iff (m : GreenMetric X F) (s : (F.val.obj (op ⊤) : Type))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    m.norm s p = 0 ↔
      (((Scheme.Modules.pullbackPushforwardAdjunction p).unit.app F).val.app (op ⊤)).hom s = 0 := by
  show Real.exp (-(m.green p)) * m.base.nrm p (arcEval p F s) = 0 ↔ _
  rw [mul_eq_zero]
  exact ⟨fun h => h.elim (fun h => absurd h (Real.exp_ne_zero _))
      (fun h => (m.base.eq_zero_iff p _).1 h),
    fun h => Or.inr ((m.base.eq_zero_iff p _).2 h)⟩

/-- ★★★ノルムは連続である——第 283 の `cont` と `green_cont` の積。 -/
theorem norm_continuous (m : GreenMetric X F) (s : (F.val.obj (op ⊤) : Type)) :
    @Continuous _ ℝ (arcTopology X) _ (m.norm s) := by
  letI := arcTopology X
  exact (Real.continuous_exp.comp m.green_cont.neg).mul (m.base.cont s)

/-- ★**定数倍**——Green 関数を `+c` 平行移動する。 -/
def scale (c : ℝ) (m : GreenMetric X F) : GreenMetric X F where
  base := m.base
  green := fun p => m.green p + c
  green_cont := by
    letI := arcTopology X
    exact m.green_cont.add continuous_const

@[simp] theorem scale_green (c : ℝ) (m : GreenMetric X F) (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    (scale c m).green p = m.green p + c := rfl

/-- ★★定数倍はノルムを一様に `exp (-c)` 倍する。 -/
theorem norm_scale (c : ℝ) (m : GreenMetric X F) (s : (F.val.obj (op ⊤) : Type))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    (scale c m).norm s p = Real.exp (-c) * m.norm s p := by
  show Real.exp (-(m.green p + c)) * m.base.nrm p (arcEval p F s)
    = Real.exp (-c) * (Real.exp (-(m.green p)) * m.base.nrm p (arcEval p F s))
  rw [show -(m.green p + c) = -c + -(m.green p) by ring, Real.exp_add, mul_assoc]

end GreenMetric

/-- ★★★★**存在**——第 285 の計量の存在から、Green 関数を 0 に取ればよい。 -/
theorem nonempty_greenMetric (hF : IsLocallyTrivial X F.val)
    (hc : @CompactSpace (Spec (CommRingCat.of ℂ) ⟶ X) (arcTopology X))
    (ht : @T2Space (Spec (CommRingCat.of ℂ) ⟶ X) (arcTopology X)) :
    Nonempty (GreenMetric X F) := by
  obtain ⟨m⟩ := exists_contArcMetric F hF hc ht
  letI := arcTopology X
  exact ⟨⟨m, fun _ => 0, continuous_const⟩⟩

/-! ## ★出典の紐付け(`.src`) -/

def GreenMetric.norm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C3——計量を基準計量と Green 関数の対で持つこと)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
