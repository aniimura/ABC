import ABC3.Found.Arakelov.PicCartierWitness
import ABC3.Found.Arakelov.ArcSpaceInterface

/-!
# 負の対照 —— **旧 (C3) の 9 欄は「層に依らない計量」を通した**(`Check`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

## ★★★★★何を確かめるか

2026-08-19 以前の `HermitianMetricData` は 9 欄

    Metric / metric_nonempty / logMetric / logMetric_continuous
    scale / logMetric_scale / IsConjCompatible / isConjCompatible_iff
    tensorMetric / logMetric_tensor

だけを要求していた。★★本ファイルはそれが**弱すぎた**ことを機械で示す:

    Metric X L := { g : X^arc → ℝ // g は連続 }        (★`L` が現れない)

が 9 欄すべてを満たす。

## ★★★これは「嘘の witness」ではない——だから厄介だった

数学的にも「`L` 上の連続計量全体は `C(X^arc, ℝ)` 上の**捩れ集合**」であり、
各 `L` に基準計量を 1 つ選べば `Metric X L ≃ C(X^arc, ℝ)` は**正しい**。

★★★したがって `scale` / `tensorMetric` / `IsConjCompatible` の**形の欄**を
いくら足しても検出できない——どの欄も基準の選択と両立するからである。

## ★★★★★塞いだ手(2026-08-19)—— `normSection`

`Interface/Arakelov/ArcSpace.lean` に**切断のノルム** `|s|_L` と

    normSection X L m s p = 0  ↔  s を p で評価したものが 0

を足した。★右辺は `L` の切断を複素点で評価したものなので、
`Metric X L` が `L` を無視していると**満たせない**。

★★評価が代数的に書けること(解析化が要らないこと)は
`Found/Arakelov/ArcFiber.lean`(第 244 ブロック、§9-284)で確かめた。
-/

namespace ABC3.Check.Arakelov

open AlgebraicGeometry CategoryTheory ABC3.Found.Arakelov ABC3.Interface.Arakelov

/-- ★**旧 (C3) の 9 欄**(2026-08-19 以前の `HermitianMetricData`)。 -/
structure WeakMetricData where
  toPicardData : PicardData
  toArcSpaceData : ArcSpaceData
  Metric : (X : Scheme.{0}) → toPicardData.Pic X → Type
  metric_nonempty : ∀ (X : Scheme.{0}) (L : toPicardData.Pic X), Nonempty (Metric X L)
  logMetric : (X : Scheme.{0}) → (L : toPicardData.Pic X) → Metric X L →
    toArcSpaceData.Arc X → ℝ
  logMetric_continuous : ∀ (X : Scheme.{0}) (L : toPicardData.Pic X) (m : Metric X L),
    @Continuous (toArcSpaceData.Arc X) ℝ (toArcSpaceData.topology X) inferInstance
      (logMetric X L m)
  scale : (X : Scheme.{0}) → (L : toPicardData.Pic X) → ℝ → Metric X L → Metric X L
  logMetric_scale : ∀ (X : Scheme.{0}) (L : toPicardData.Pic X) (c : ℝ) (m : Metric X L)
    (p : toArcSpaceData.Arc X),
    logMetric X L (scale X L c m) p = logMetric X L m p + c
  IsConjCompatible : (X : Scheme.{0}) → (L : toPicardData.Pic X) → Metric X L → Prop
  isConjCompatible_iff : ∀ (X : Scheme.{0}) (L : toPicardData.Pic X) (m : Metric X L),
    IsConjCompatible X L m ↔
      ∀ p : toArcSpaceData.Arc X,
        logMetric X L m (toArcSpaceData.conj X p) = logMetric X L m p
  tensorMetric : (X : Scheme.{0}) → (L M : toPicardData.Pic X) → Metric X L → Metric X M →
    Metric X (@HMul.hMul _ _ _
      (@instHMul _ (toPicardData.group X).toDivInvMonoid.toMonoid.toMulOneClass.toMul) L M)
  logMetric_tensor : ∀ (X : Scheme.{0}) (L M : toPicardData.Pic X)
    (mL : Metric X L) (mM : Metric X M) (p : toArcSpaceData.Arc X),
    logMetric X _ (tensorMetric X L M mL mM) p
      = logMetric X L mL p + logMetric X M mM p

/-- ★**`L` に依らない「計量」**——連続な実数値関数(Green 関数)そのもの。 -/
abbrev Green (X : Scheme.{0}) : Type :=
  { g : arcSpaceDataImpl.Arc X → ℝ //
      @Continuous _ _ (arcSpaceDataImpl.topology X) _ g }

/-- ★定数 `0`。 -/
def Green.zero (X : Scheme.{0}) : Green X :=
  letI := arcSpaceDataImpl.topology X
  ⟨fun _ => 0, continuous_const⟩

/-- ★定数を足す。 -/
def Green.shift {X : Scheme.{0}} (m : Green X) (c : ℝ) : Green X :=
  letI := arcSpaceDataImpl.topology X
  ⟨fun p => m.val p + c, m.2.add continuous_const⟩

/-- ★和。 -/
def Green.plus {X : Scheme.{0}} (a b : Green X) : Green X :=
  letI := arcSpaceDataImpl.topology X
  ⟨fun p => a.val p + b.val p, a.2.add b.2⟩

/-- ★★★★★**退化 witness**——`Metric X L` が `L` を無視しても旧 9 欄は通る。 -/
noncomputable def degenerateMetricData : WeakMetricData where
  toPicardData := picardDataWitness
  toArcSpaceData := arcSpaceDataImpl
  Metric := fun X _ => Green X
  metric_nonempty := fun X _ => ⟨Green.zero X⟩
  logMetric := fun _ _ m => m.val
  logMetric_continuous := fun _ _ m => m.2
  scale := fun _ _ c m => Green.shift m c
  logMetric_scale := fun _ _ _ _ _ => rfl
  IsConjCompatible := fun X _ m =>
    ∀ p : arcSpaceDataImpl.Arc X, m.val (arcSpaceDataImpl.conj X p) = m.val p
  isConjCompatible_iff := fun _ _ _ => Iff.rfl
  tensorMetric := fun _ _ _ mL mM => Green.plus mL mM
  logMetric_tensor := fun _ _ _ _ _ _ => by rfl

/-- ★★★★★**旧 9 欄は計量と直線束の結び付きを要求していなかった**。

★この定理が通ること自体が、`normSection`(2026-08-19 追加)が要る理由の証明である。 -/
theorem weakMetricData_ignores_bundle :
    ∃ H : WeakMetricData,
      ∀ (X : Scheme.{0}) (L M : H.toPicardData.Pic X), H.Metric X L = H.Metric X M :=
  ⟨degenerateMetricData, fun _ _ _ => by rfl⟩

end ABC3.Check.Arakelov
