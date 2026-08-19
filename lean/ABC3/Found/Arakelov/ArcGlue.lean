import ABC3.Found.Arakelov.ArcLocalMetric

/-!
# Arakelov (C3) 第 262 ブロック —— **貼り合わせたノルムの定義と代数的法則**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★1 の分割で貼る —— 代数の段

局所ノルム `genNorm` は開集合 `U i` の上でしか定義されない。★**0 で延長**して

    extNorm i p w := if p⁻¹(U i) = ⊤ then genNorm … else 0
    gluedNorm p w := ∑ᶠ i, ρ i p * extNorm i p w

と置く(`ρ` は 1 の分割)。

★★代数的法則(非負・`‖c·v‖ = |c|‖v‖`)は**各項ごとに成り立つ**ので、
`finsum` の線形性(mathlib `mul_finsum`)でそのまま上がる。

★★★残るのは
- `gluedNorm p w = 0 → w = 0`(★どこかで `ρ i p > 0` かつ `p ∈ U i`)
- 連続性(★mathlib `PartitionOfUnity.IsSubordinate.continuous_finsum_smul`)

## ★摩擦 —— `dite` には `split` が効かない

`if h : … then … else …`(依存 `if`)は `split` で場合分けできない。
★`by_cases` + `simp only [dif_pos h]` / `[dif_neg h]` を使う。
★★`Decidable` が要るので `open scoped Classical` を先に置く。

| 定義・定理 | 内容 |
|---|---|
| `extNorm` | ★0 で延長した局所ノルム |
| `gluedNorm` | ★★1 の分割で貼ったノルム |
| `extNorm_nonneg` / `_smul` | ★各項の法則 |
| `gluedNorm_nonneg` / `_smul` | ★★★和の法則(`finsum_nonneg` / `mul_finsum`) |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite
open scoped Classical

variable {X : Scheme.{0}} {ι : Type}

/-- ★局所ノルムを 0 で延長する。 -/
noncomputable def extNorm (F : X.Modules) (hF : IsLocallyTrivial X F.val)
    (U : ι → X.Opens) (g : ∀ i, (F.val.obj (op (U i)) : Type)) (i : ι)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (w : ↥(arcFiber p F)) : ℝ :=
  if h : p ⁻¹ᵁ (U i) = ⊤ then genNorm F hF (U i) (g i) p h w else 0

/-- ★★貼り合わせたノルム。 -/
noncomputable def gluedNorm (F : X.Modules) (hF : IsLocallyTrivial X F.val)
    (U : ι → X.Opens) (g : ∀ i, (F.val.obj (op (U i)) : Type))
    (ρ : ι → (Spec (CommRingCat.of ℂ) ⟶ X) → ℝ)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (w : ↥(arcFiber p F)) : ℝ :=
  ∑ᶠ i, ρ i p * extNorm F hF U g i p w

theorem extNorm_nonneg (F : X.Modules) (hF : IsLocallyTrivial X F.val)
    (U : ι → X.Opens) (g : ∀ i, (F.val.obj (op (U i)) : Type)) (i : ι)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (w : ↥(arcFiber p F)) :
    0 ≤ extNorm F hF U g i p w := by
  by_cases h : p ⁻¹ᵁ (U i) = ⊤
  · simp only [extNorm, dif_pos h]
    exact genNorm_nonneg F hF _ _ _ _ _
  · simp only [extNorm, dif_neg h]
    exact le_rfl

theorem extNorm_smul (F : X.Modules) (hF : IsLocallyTrivial X F.val)
    (U : ι → X.Opens) (g : ∀ i, (F.val.obj (op (U i)) : Type)) (i : ι)
    (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (c : (CommRingCat.of ℂ : Type)) (w : ↥(arcFiber p F)) :
    extNorm F hF U g i p (c • w) = ‖c‖ * extNorm F hF U g i p w := by
  by_cases h : p ⁻¹ᵁ (U i) = ⊤
  · simp only [extNorm, dif_pos h]
    exact genNorm_smul F hF _ _ _ _ _ _
  · simp only [extNorm, dif_neg h]
    rw [mul_zero]


theorem gluedNorm_nonneg (F : X.Modules) (hF : IsLocallyTrivial X F.val)
    (U : ι → X.Opens) (g : ∀ i, (F.val.obj (op (U i)) : Type))
    (ρ : ι → (Spec (CommRingCat.of ℂ) ⟶ X) → ℝ) (hρ : ∀ i p, 0 ≤ ρ i p)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (w : ↥(arcFiber p F)) :
    0 ≤ gluedNorm F hF U g ρ p w :=
  finsum_nonneg (fun i => mul_nonneg (hρ i p) (extNorm_nonneg F hF U g i p w))

theorem gluedNorm_smul (F : X.Modules) (hF : IsLocallyTrivial X F.val)
    (U : ι → X.Opens) (g : ∀ i, (F.val.obj (op (U i)) : Type))
    (ρ : ι → (Spec (CommRingCat.of ℂ) ⟶ X) → ℝ)
    (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (c : (CommRingCat.of ℂ : Type)) (w : ↥(arcFiber p F)) :
    gluedNorm F hF U g ρ p (c • w) = ‖c‖ * gluedNorm F hF U g ρ p w := by
  show (∑ᶠ i, ρ i p * extNorm F hF U g i p (c • w)) = ‖c‖ * ∑ᶠ i, ρ i p * extNorm F hF U g i p w
  rw [mul_finsum]
  refine finsum_congr (fun i => ?_)
  rw [extNorm_smul]
  ring


/-! ## ★出典の紐付け(`.src`) -/

def gluedNorm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——1 の分割で貼ったノルム)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
