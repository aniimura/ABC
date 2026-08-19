import ABC3.Found.Arakelov.PicSquareLE
import ABC3.Found.Arakelov.PicAppIsoInv

/-!
# Arakelov (B2) 第 227 ブロック —— ★★★★★★**`hcompat` の核(座標を揃えた形)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★§9-256 の「座標系 `⊤` を捨てる」が完成した

    hB.fromSpec.app B ≫ (Spec.map appLE).appLE _ _ e₁
      = f.app B ≫ hA.fromSpec.appLE _ _ e₂

★**`⊤` が一度も出ない**——座標を `fromSpec ⁻¹ᵁ A` に取った。

## ★★機構は 3 本の合流

| 段 | 補題 |
|---|---|
| 可換正方形を任意座標で | ★第 226 |
| 左辺を分解 | ★`Scheme.Hom.comp_appLE`(mathlib) |
| 右辺を分解 | ★同上 |

★★★`comp_appLE` は `(f ≫ g).appLE U V e = g.app U ≫ f.appLE _ V e` を与えるので、
**両辺が同じ形に割れる**。あとは第 226 を当てるだけ。

## ★★★7 つの摩擦を越えて到達した

| # | 症状 | 逃げ道 | ブロック |
|---|---|---|---|
| 1 | `whnf` timeout | 抽象の側で書く | 219 |
| 2 | 始域が違う | 分解を挟む | 215 |
| 3 | `simp` の過剰正規化 | 元のレベルで書く | 217 |
| 4 | 依存位置 | grep し直す / `congr_app` | 222 / 225 |
| 5 | 転送の伝播 | 同型であることを使う | 223 |
| 6 | `rw` が証明部分まで照合 | 証明を引数に出す | 224 |
| 7 | 座標系が合わない | ★**`appLE` で依存しない座標を選ぶ** | ★226–227 |

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `square_at_pre` | ★第 226 を `fromSpec⁻¹A` で読む |
| `lhs_split` / `rhs_split` | ★`comp_appLE` による分解 |
| `hcompat_core` | ★★★★★★**`hcompat` の核** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-- ★`W := fromSpec ⁻¹ᵁ A` で読む。 -/
theorem square_at_pre {A : X.Opens} {B : Y.Opens}
    (hA : IsAffineOpen A) (hB : IsAffineOpen B) (i : A ≤ f ⁻¹ᵁ B)
    (e1 : (hA.fromSpec ⁻¹ᵁ A) ≤ (Spec.map (f.appLE B A i) ≫ hB.fromSpec) ⁻¹ᵁ B)
    (e2 : (hA.fromSpec ⁻¹ᵁ A) ≤ (hA.fromSpec ≫ f) ⁻¹ᵁ B) :
    (Spec.map (f.appLE B A i) ≫ hB.fromSpec).appLE B (hA.fromSpec ⁻¹ᵁ A) e1
      = (hA.fromSpec ≫ f).appLE B (hA.fromSpec ⁻¹ᵁ A) e2 :=
  square_appLE f hA hB i _ e1 e2

/-- ★左辺の分解。 -/
theorem lhs_split {A : X.Opens} {B : Y.Opens}
    (hA : IsAffineOpen A) (hB : IsAffineOpen B) (i : A ≤ f ⁻¹ᵁ B)
    (e1 : (hA.fromSpec ⁻¹ᵁ A) ≤ (Spec.map (f.appLE B A i) ≫ hB.fromSpec) ⁻¹ᵁ B) :
    (Spec.map (f.appLE B A i) ≫ hB.fromSpec).appLE B (hA.fromSpec ⁻¹ᵁ A) e1
      = hB.fromSpec.app B ≫ (Spec.map (f.appLE B A i)).appLE _ _ e1 :=
  Scheme.Hom.comp_appLE _ _ _ _ _




/-- ★右辺の分解。 -/
theorem rhs_split {A : X.Opens} {B : Y.Opens}
    (hA : IsAffineOpen A) (i : A ≤ f ⁻¹ᵁ B)
    (e2 : (hA.fromSpec ⁻¹ᵁ A) ≤ (hA.fromSpec ≫ f) ⁻¹ᵁ B) :
    (hA.fromSpec ≫ f).appLE B (hA.fromSpec ⁻¹ᵁ A) e2
      = f.app B ≫ hA.fromSpec.appLE _ _ e2 :=
  Scheme.Hom.comp_appLE _ _ _ _ _

/-- ★★★★★★hcompat の核(座標を揃えた形)。 -/
theorem hcompat_core {A : X.Opens} {B : Y.Opens}
    (hA : IsAffineOpen A) (hB : IsAffineOpen B) (i : A ≤ f ⁻¹ᵁ B)
    (e1 : (hA.fromSpec ⁻¹ᵁ A) ≤ (Spec.map (f.appLE B A i) ≫ hB.fromSpec) ⁻¹ᵁ B)
    (e2 : (hA.fromSpec ⁻¹ᵁ A) ≤ (hA.fromSpec ≫ f) ⁻¹ᵁ B) :
    hB.fromSpec.app B ≫ (Spec.map (f.appLE B A i)).appLE _ _ e1
      = f.app B ≫ hA.fromSpec.appLE _ _ e2 := by
  rw [← lhs_split f hA hB i e1, ← rhs_split f hA i e2]
  exact square_at_pre f hA hB i e1 e2

/-! ## ★出典の紐付け(`.src`) -/

def hcompat_core.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——hcompat の核、座標を揃えた形)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
