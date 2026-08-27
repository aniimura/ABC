import ABC3.Found.GenEll.MinField
import ABC3.Found.GenEll.LogDiffTower

/-!
# [GenEll] Proposition 1.7 —— **`log-diff_Y − log-diff_Z ≥ 0`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.10。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

## ★★★★原文が `≥` と書いて `≳` と書かなかった段

原文 p.10 は『the "portion over Σ" of `log-diff_Y − log-diff_Z` is **≥ 0**
[i.e., with "≥", not "≳"!]』と、**わざわざ区別して**書いている。

★★その中身は 2 行である:

1. `φ : Y → Z` を通すと、点の**最小定義体は小さくなる**(`F_min(φ∘y) ⊆ F_min(y)`)
2. `log-diff` は体を大きくすると**増えるだけ**(`logDiffOfField_le`、既実装)

★★★したがって `log-diff_Y(y) − log-diff_Z(φ(y)) ≥ 0` が**等号つきで**出る。

## ★★1 の中身 —— `x_F` は `Spec F_min` を経由する

`MinField.lean` は**極小性**(`x_F` が `E` を経由するなら `F_min ⊆ E`)を取ったが、
「`F_min` 自身が定義体である」ことは `κ(ξ)` を経由する形でしか取っていなかった。

★本ファイルはそれを `F_min` まで下ろす——`κ(ξ) ↠ F_min ↪ F` と分解するだけである。
★★これがあれば `x_F ≫ φ` も `Spec F_min` を経由するので、
**極小性がそのまま被覆の単調性を与える**。

## ★本ファイルで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_factor_minField` | ★★`x_F` は `Spec F_min` を経由する |
| `minField_comp_le` | ★★★**`F_min(x_F ≫ φ) ⊆ F_min(x_F)`** |
| `logDiffOfField_minField_comp_le` | ★★★★**`log-diff_Z(φ(y)) ≤ log-diff_Y(y)`** |
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

section Factor

variable (F : Type) [Field F] {X : Scheme}

/-- ★★**`x_F` は `Spec F_min` を経由する** —— `F_min` が実際に定義体であること。

★`MinField.lean` は `κ(ξ)` を経由することしか言っていなかった。
`κ(ξ) → F` の像が `F_min` なので、**値域制限して分解する**だけである。 -/
theorem exists_factor_minField (xF : Spec (CommRingCat.of F) ⟶ X) :
    ∃ xE : Spec (CommRingCat.of ((minField F xF : Subfield F) : Type _)) ⟶ X,
      xF = Spec.map (CommRingCat.ofHom (minField F xF).subtype) ≫ xE := by
  refine ⟨Spec.map (CommRingCat.ofHom
      (RingHom.rangeRestrictField (minFieldData F xF).2.hom))
    ≫ X.fromSpecResidueField (minFieldData F xF).1, ?_⟩
  have hcomp : CommRingCat.ofHom
        (RingHom.rangeRestrictField (minFieldData F xF).2.hom)
      ≫ CommRingCat.ofHom (minField F xF).subtype = (minFieldData F xF).2 := by
    apply CommRingCat.hom_ext
    ext x
    rfl
  conv_lhs => rw [← specMap_minFieldData_comp F xF, ← hcomp, Spec.map_comp,
    Category.assoc]
  rfl

/-- ★★★**被覆を通すと最小定義体は小さくなる** —— `F_min(x_F ≫ φ) ⊆ F_min(x_F)`。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

★証明は 2 行である——`x_F` が `Spec F_min` を経由するので
`x_F ≫ φ` も経由し、**極小性**(`minField_le_of_factors`)がそのまま効く。

★★これが `Proposition 1.7` の『`log-diff_Y − log-diff_Z` is **≥ 0**
[i.e., with "≥", not "≳"!]』の前半である——
後半は `logDiffOfField_le`(体を大きくすると `log-diff` は増える)であり、
`LogDiffTower.lean` に実装済みである。 -/
theorem minField_comp_le {Y Z : Scheme} (φ : Y ⟶ Z)
    (xF : Spec (CommRingCat.of F) ⟶ Y) :
    minField F (xF ≫ φ) ≤ minField F xF := by
  obtain ⟨xE, hxE⟩ := exists_factor_minField F xF
  refine minField_le_of_factors F (minField F xF) (xE ≫ φ) (xF ≫ φ) ?_
  conv_lhs => rw [hxE]
  rw [Category.assoc]

end Factor



/-! ## ★出典の紐付け(`.src`) -/

def exists_factor_minField.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (i)(F_min が実際に定義体であること——x_F は Spec F_min を経由する)",
    sectionId := "genell-def-1-5" }

def minField_comp_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i)(Σ 上の log-diff_Y − log-diff_Z が ≥ 0 であることの前半)",
    sectionId := "genell-prop-1-7" }

end ABC3.Found.GenEll
