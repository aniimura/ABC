import ABC3.Found.Arakelov.ArcFiber
import Mathlib.Analysis.SpecialFunctions.Exp

/-!
# Arakelov (C3) 第 246 ブロック —— ★★★★★**ファイバーごとのノルムの族**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★正直な `Metric` の型

§9-286 で `Interface` に足した `normSection` を満たすには、
計量が **`L` のファイバーを知っている**必要がある。★その最小の形が

    ArcMetric X L := 各複素点 p でのファイバー `arcFiber p L` 上のノルムの族

である。★★`normSection X L m s p := ‖(arcEval p L s)‖` と置けば
`normSection_eq_zero_iff` は `eq_zero_iff` **そのもの**になる。

## ★摩擦 —— `(e : Type)` の型上書きは `∀` の中で潰れる

`arcEval` では `(arcFiber p L : Type)` と書けたが、構造体のフィールドの
`∀` 束縛の中では `nrm p : Type → ℝ` と読まれて落ちた。★`↥(arcFiber p L)` と書けば通る。
★★[[type-spelling-two-paths]]の新しい顔である。

| 定義・定理 | 内容 |
|---|---|
| `ArcMetric` | ★★★★ファイバーごとのノルムの族 |
| `ArcMetric.scale` | ★★定数倍(`exp (-c)` 倍) |
| `ArcMetric.normOf` | ★★★**切断のノルム** `|s|_L` |
| `normOf_eq_zero_iff` | ★★★★**退化封じの欄そのもの** |
| `normOf_scale` | ★`scale` との整合 |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite

variable {X : Scheme.{0}}

/-- ★★★★**複素点ごとのノルムの族**——正直な計量の型。 -/
structure ArcMetric (X : Scheme.{0}) (L : X.Modules) where
  /-- 各複素点のファイバー上のノルム。 -/
  nrm : ∀ (p : Spec (CommRingCat.of ℂ) ⟶ X), ↥(arcFiber p L) → ℝ
  /-- ★非負。 -/
  nonneg : ∀ (p : Spec (CommRingCat.of ℂ) ⟶ X) (v : ↥(arcFiber p L)), 0 ≤ nrm p v
  /-- ★★**0 になるのは 0 のときだけ**——これが `L` との結び付きを与える。 -/
  eq_zero_iff : ∀ (p : Spec (CommRingCat.of ℂ) ⟶ X) (v : ↥(arcFiber p L)),
    nrm p v = 0 ↔ v = 0
  /-- ★★スカラー倍で絶対値倍。 -/
  smul : ∀ (p : Spec (CommRingCat.of ℂ) ⟶ X) (c : ℂ) (v : ↥(arcFiber p L)),
    nrm p (c • v) = ‖c‖ * nrm p v

namespace ArcMetric

variable {L : X.Modules}

/-- ★★**計量を `exp (-c)` 倍する**。 -/
noncomputable def scale (c : ℝ) (m : ArcMetric X L) : ArcMetric X L where
  nrm p v := Real.exp (-c) * m.nrm p v
  nonneg p v := mul_nonneg (Real.exp_pos _).le (m.nonneg p v)
  eq_zero_iff p v := by
    rw [mul_eq_zero, m.eq_zero_iff]
    exact or_iff_right (Real.exp_ne_zero _)
  smul p c' v := by rw [m.smul]; ring

/-- ★★★**大域切断のノルム** `|s|_L`。 -/
noncomputable def normOf (m : ArcMetric X L) (s : (L.val.obj (op ⊤) : Type))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) : ℝ :=
  m.nrm p (arcEval p L s)

/-- ★非負。 -/
theorem normOf_nonneg (m : ArcMetric X L) (s : (L.val.obj (op ⊤) : Type))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) : 0 ≤ m.normOf s p :=
  m.nonneg p _

/-- ★★★★**切断が点で消えることとノルムが 0 になることは同値**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `Interface` の `normSection_eq_zero_iff`(§9-286 で追加)そのものであり、
「`Metric X L` が `L` を無視する」退化を殺す欄である。 -/
theorem normOf_eq_zero_iff (m : ArcMetric X L) (s : (L.val.obj (op ⊤) : Type))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    m.normOf s p = 0 ↔ arcEval p L s = 0 :=
  m.eq_zero_iff p _

/-- ★`scale` はノルムを一様に `exp (-c)` 倍する。 -/
theorem normOf_scale (c : ℝ) (m : ArcMetric X L) (s : (L.val.obj (op ⊤) : Type))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    (ArcMetric.scale c m).normOf s p = Real.exp (-c) * m.normOf s p :=
  rfl

end ArcMetric

/-! ## ★出典の紐付け(`.src`) -/

def ArcMetric.normOf_eq_zero_iff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——切断のノルムと退化封じ)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
