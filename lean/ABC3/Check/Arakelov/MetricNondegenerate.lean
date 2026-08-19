import ABC3.Found.Arakelov.PicCartierWitness
import ABC3.Found.Arakelov.ArcSpaceInterface

/-!
# 負の対照 —— **`HermitianMetricData` は「層に依らない計量」を通す**(`Check`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

## ★★★★★何を確かめるか

(C3) `HermitianMetricData` は `Metric : (X) → Pic X → Type` を posit する。
★★2026-08-17 に「計量を `L` ごとに 1 つ固定する」退化は塞いだが、
**`Metric X L` が `L` に依らない**退化は塞げていない。

★★★本ファイルはそれを機械で確かめる——

    Metric X L := { g : X^arc → ℝ // g は連続 }        (★`L` が現れない)
    logMetric  := g,  scale c g := g + c,  tensor := (+)

が **9 欄すべてを満たす**。したがって (C3) は現状、
**計量と直線束の結び付きを何も要求していない**。

## ★★★これは「嘘の witness」ではない——だから厄介である

数学的にも「`L` 上の連続計量全体は `C(X^arc, ℝ)` 上の**捩れ集合**(torsor)」であり、
計量は存在する。★したがって「各 `L` に基準計量を 1 つ選ぶ」ことで
`Metric X L ≃ C(X^arc, ℝ)` は**正しい**。

★★★欠けているのは「**基準の取り方が標準的でない**」ことを表す仕組みであり、
`scale` / `tensor` / `IsConjCompatible` の形の欄をいくら足しても検出できない
——どの欄も基準の選択と両立するからである。

★★★★★**結論**: 結び付きは (C3) の中では表せない。
`Metric` を `L` に縛るには **`L` の切断のノルム**が要り、
それは `PicardData.sheafOf`(B1、実装済)を経由して初めて書ける。
★あるいは (D2) の `degF` が**有限部分**(B2/B3)を見ることを要求する欄を足す。

★★したがって (C3) を「達成」と数える前に、この欄を足さねばならない。
-/

namespace ABC3.Check.Arakelov

open AlgebraicGeometry CategoryTheory ABC3.Found.Arakelov ABC3.Interface.Arakelov

/-- ★**`L` に依らない「計量」**——連続な実数値関数（Green 関数）そのもの。 -/
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

/-- ★★★★★**退化 witness**——`Metric X L` が `L` を無視しても 9 欄すべて通る。 -/
noncomputable def degenerateMetricData : HermitianMetricData where
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

/-- ★★★★★**(C3) は現状、計量と直線束の結び付きを要求していない**。

★この定理が通ること自体が、`HermitianMetricData` に欄が足りないことの証明である。 -/
theorem hermitianMetricData_ignores_bundle :
    ∃ H : HermitianMetricData,
      ∀ (X : Scheme.{0}) (L M : H.toPicardData.Pic X), H.Metric X L = H.Metric X M :=
  ⟨degenerateMetricData, fun _ _ _ => by rfl⟩

end ABC3.Check.Arakelov
