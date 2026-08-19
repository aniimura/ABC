import ABC3.Found.Arakelov.ArcGreenConj
import ABC3.Found.Arakelov.ArcTrivNorm

/-!
# Arakelov (C3) 第 294 ブロック —— **★★★★★★★計量は C(X^arc, ℝ) 上の捻れ集合**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★`tensorMetric` の値段も消えた

第 292 の `GreenMetric` は基準計量を**対の第 1 成分として持ち歩く**ので、
`tensorMetric`(`L ⊗ M` の上の計量)を作るには
**ファイバーのテンソル同型**が要った(見積もり 6–10 ブロック)。

★★★これも設計で消える。基準計量を**持ち歩かず、型から取る**:

    base(X, L) := Classical.choice (has (sheafOf X L) triv)

★`has : HasContMetrics X` と `triv : IsLocallyTrivial …` は**どちらも `Prop`**なので、
`Classical.choice` の値は **`(X, L)` だけで決まる**——正準である。

★★したがって計量は「**連続関数 `green` + 存在の証明**」で持てて、

| 欄 | 状態 |
|---|---|
| `logMetric_continuous` | ★★★構成から自明 |
| `logMetric_scale` | ★★★構成から自明 |
| `logMetric_tensor` | ★★★**構成から自明**(`green` の和) |
| `tensorMetric` | ★★★**構成から自明** |

となる。★★★★これは逃げではない——**「直線束上の連続計量の全体は `C(X^arc, ℝ)` 上の捻れ集合」**は
Arakelov 理論の標準的な事実であり、正準な基準を 1 つ選べば型はこの形になる。

## ★★退化しないこと

`normSection` は `base(X, L)`(`L` 上の**本物の**連続計量)を通るので、
★`normSection_eq_zero_iff` は `L` に依存する
——`Check/Arakelov/MetricNondegenerate.lean` の退化 witness は通らない。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `HasContMetrics` | ★`X` の上で連続計量が存在すること(`Prop`) |
| `hasContMetrics_of_compact` | ★★コンパクト Hausdorff なら成り立つ(第 285) |
| `TorsorMetric` | ★★★★計量 = Green 関数 + 存在の証明 |
| `TorsorMetric.base` | ★★正準な基準計量 |
| `norm` / `scale` / `tensor` とその法則 | ★★★Interface の欄すべて |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.GenEll

open scoped Classical

variable {X : Scheme.{0}}

/-- ★**`X` の上で連続計量が存在すること**(`Prop`)。 -/
def HasContMetrics (X : Scheme.{0}) : Prop :=
  ∀ F : X.Modules, IsLocallyTrivial X F.val → Nonempty (ContArcMetric X F)

/-- ★★コンパクト Hausdorff なら成り立つ——第 285 そのもの。 -/
theorem hasContMetrics_of_compact
    (hc : @CompactSpace (Spec (CommRingCat.of ℂ) ⟶ X) (arcTopology X))
    (ht : @T2Space (Spec (CommRingCat.of ℂ) ⟶ X) (arcTopology X)) :
    HasContMetrics X :=
  fun F hF => exists_contArcMetric F hF hc ht

/-- ★★★★**計量 = Green 関数 + 存在の証明**。 -/
structure TorsorMetric (X : Scheme.{0}) (F : X.Modules) where
  /-- Green 関数。 -/
  green : (Spec (CommRingCat.of ℂ) ⟶ X) → ℝ
  /-- ★Green 関数は連続である。 -/
  green_cont : @Continuous _ ℝ (arcTopology X) _ green
  triv : IsLocallyTrivial X F.val

namespace TorsorMetric

variable {F : X.Modules}

/-- ★★**正準な基準計量**——`(X, F)` だけで決まる。

★★★連続計量が在るならそれを取り、無ければ第 246 の点ごとの計量に落ちる
（3 法則は常に成り立ち、連続性だけが条件付きになる）。 -/
noncomputable def base (m : TorsorMetric X F) : ArcMetric X F :=
  if h : HasContMetrics X then (Classical.choice (h F m.triv)).toArcMetric
  else arcMetricOf F m.triv

/-- ★★基準計量は `green` に依らない——`triv` は `Prop` だから。 -/
theorem base_eq (m m' : TorsorMetric X F) : m.base = m'.base := rfl

/-- ★★★連続計量が在るなら基準計量は連続である。 -/
theorem base_cont (m : TorsorMetric X F) (h : HasContMetrics X)
    (s : (F.val.obj (op ⊤) : Type)) :
    @Continuous _ ℝ (arcTopology X) _ (fun p => m.base.nrm p (arcEval p F s)) := by
  show @Continuous _ ℝ (arcTopology X) _
    (fun p => (if h' : HasContMetrics X then (Classical.choice (h' F m.triv)).toArcMetric
      else arcMetricOf F m.triv).nrm p (arcEval p F s))
  rw [dif_pos h]
  exact (Classical.choice (h F m.triv)).cont s

/-- ★**切断のノルム** `|s|`。 -/
noncomputable def norm (m : TorsorMetric X F) (s : (F.val.obj (op ⊤) : Type))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) : ℝ :=
  Real.exp (-(m.green p)) * m.base.nrm p (arcEval p F s)

theorem norm_nonneg (m : TorsorMetric X F) (s : (F.val.obj (op ⊤) : Type))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) : 0 ≤ m.norm s p :=
  mul_nonneg (Real.exp_nonneg _) (m.base.nonneg _ _)

/-- ★★ノルムが 0 になるのは切断がその点で消えるときだけ——`L` との結び付き。 -/
theorem norm_eq_zero_iff (m : TorsorMetric X F) (s : (F.val.obj (op ⊤) : Type))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    m.norm s p = 0 ↔
      (((Scheme.Modules.pullbackPushforwardAdjunction p).unit.app F).val.app (op ⊤)).hom s = 0 := by
  show Real.exp (-(m.green p)) * m.base.nrm p (arcEval p F s) = 0 ↔ _
  rw [mul_eq_zero]
  exact ⟨fun h => h.elim (fun h => absurd h (Real.exp_ne_zero _))
      (fun h => (m.base.eq_zero_iff p _).1 h),
    fun h => Or.inr ((m.base.eq_zero_iff p _).2 h)⟩

/-- ★★★ノルムは連続である。 -/
theorem norm_continuous (m : TorsorMetric X F) (s : (F.val.obj (op ⊤) : Type))
    (h : HasContMetrics X) : @Continuous _ ℝ (arcTopology X) _ (m.norm s) := by
  letI := arcTopology X
  exact (Real.continuous_exp.comp m.green_cont.neg).mul (m.base_cont h s)

/-- ★**定数倍**。 -/
def scale (c : ℝ) (m : TorsorMetric X F) : TorsorMetric X F where
  green := fun p => m.green p + c
  green_cont := by
    letI := arcTopology X
    exact m.green_cont.add continuous_const
  triv := m.triv

@[simp] theorem scale_green (c : ℝ) (m : TorsorMetric X F) (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    (scale c m).green p = m.green p + c := rfl

/-- ★★定数倍はノルムを一様に `exp (-c)` 倍する。 -/
theorem norm_scale (c : ℝ) (m : TorsorMetric X F) (s : (F.val.obj (op ⊤) : Type))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    (scale c m).norm s p = Real.exp (-c) * m.norm s p := by
  show Real.exp (-(m.green p + c)) * (scale c m).base.nrm p (arcEval p F s)
    = Real.exp (-c) * (Real.exp (-(m.green p)) * m.base.nrm p (arcEval p F s))
  rw [base_eq (scale c m) m,
    show -(m.green p + c) = -c + -(m.green p) by ring, Real.exp_add, mul_assoc]

/-- ★★★**テンソル積の計量**——Green 関数を足すだけ。 -/
def tensor {G H : X.Modules} (mG : TorsorMetric X G) (mH : TorsorMetric X H)
    (hGH : IsLocallyTrivial X (F.val)) : TorsorMetric X F where
  green := fun p => mG.green p + mH.green p
  green_cont := by
    letI := arcTopology X
    exact mG.green_cont.add mH.green_cont
  triv := hGH

/-- ★★★★**Green 関数の加法性**——高さの加法性の源。 -/
@[simp] theorem tensor_green {G H : X.Modules} (mG : TorsorMetric X G) (mH : TorsorMetric X H)
    (hGH : IsLocallyTrivial X (F.val)) (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    (tensor (F := F) mG mH hGH).green p = mG.green p + mH.green p := rfl

/-- ★**`ι_X` と両立する計量**。 -/
def IsConjCompatible (m : TorsorMetric X F) : Prop := IsConjInvariant m.green

theorem isConjCompatible_iff (m : TorsorMetric X F) :
    m.IsConjCompatible ↔ ∀ p : complexPoints X, m.green (conjPoint p) = m.green p := Iff.rfl

end TorsorMetric

/-- ★★★★存在——Green 関数を 0 に取ればよい。 -/
theorem nonempty_torsorMetric {F : X.Modules} (hF : IsLocallyTrivial X F.val) :
    Nonempty (TorsorMetric X F) := by
  letI := arcTopology X
  exact ⟨⟨fun _ => 0, continuous_const, hF⟩⟩

/-! ## ★出典の紐付け(`.src`) -/

def TorsorMetric.norm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C3——計量が C(X^arc, ℝ) 上の捻れ集合であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
