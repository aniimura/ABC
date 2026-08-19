import ABC3.Found.Arakelov.ArcSemilinear

/-!
# Arakelov (C3) 第 259 ブロック —— **開集合を経由する点の分解**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★連続性を `V^arc` へ移すための足場

第 258 で `genNorm (arcEvalOnTop (c·g)) = ‖evalOn c‖` が出た。
★残るのは **`p ↦ evalOn p V h c` の連続性**である。

★★ところが `h : p ⁻¹ᵁ V = ⊤` は `p` に依存するので、
定義域が**部分型**になってしまい、そのままでは連続性を述べにくい。

★★★逃げ道: **`V.toScheme` の点でパラメータ付ける**。

    q : Spec ℂ ⟶ V.toScheme   ↦   q ≫ V.ι : Spec ℂ ⟶ X

★このとき `(q ≫ V.ι) ⁻¹ᵁ V = ⊤` は**自動**である(`comp_preimage_eq_top`)。
★★逆に `p ⁻¹ᵁ V = ⊤` なら `p` は `V.ι` を経由して分解する(`liftToOpenOfTop`)。

★★★★これで連続性は `arcTopology V.toScheme` の上の話になり、
第 252(大域の正則関数は連続)を `V.toScheme` に適用すればよい。

| 定理 | 内容 |
|---|---|
| `range_subset` | ★`p⁻¹V = ⊤` なら像は `V` に入る |
| `liftToOpenOfTop` / `liftToOpen_fac` | ★★分解(mathlib `IsOpenImmersion.lift`) |
| `comp_preimage_eq_top` | ★★経由する点は自動的に条件を満たす |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite

variable {X : Scheme.{0}} (V : X.Opens) (p : Spec (CommRingCat.of ℂ) ⟶ X)

/-- ★`p⁻¹V = ⊤` なら `p` の像は `V` に入る。 -/
theorem range_subset (h : p ⁻¹ᵁ V = ⊤) : Set.range p.base ⊆ Set.range V.ι.base := by
  rw [V.range_ι]
  rintro _ ⟨z, rfl⟩
  have : z ∈ p ⁻¹ᵁ V := by rw [h]; trivial
  exact this

/-- ★★`p` は `V` を経由して分解する。 -/
noncomputable def liftToOpenOfTop (h : p ⁻¹ᵁ V = ⊤) : Spec (CommRingCat.of ℂ) ⟶ V.toScheme :=
  IsOpenImmersion.lift V.ι p (range_subset V p h)

theorem liftToOpen_fac (h : p ⁻¹ᵁ V = ⊤) : liftToOpenOfTop V p h ≫ V.ι = p :=
  IsOpenImmersion.lift_fac _ _ _


/-- ★★`V` を経由する点は自動的に `p⁻¹V = ⊤` を満たす。 -/
theorem comp_preimage_eq_top (q : Spec (CommRingCat.of ℂ) ⟶ V.toScheme) :
    (q ≫ V.ι) ⁻¹ᵁ V = ⊤ := by
  ext z
  simp only [TopologicalSpace.Opens.coe_top, Set.mem_univ, iff_true]
  have hr : (q ≫ V.ι).base z ∈ Set.range V.ι.base := ⟨q.base z, rfl⟩
  rw [V.range_ι] at hr
  exact hr


/-! ## ★出典の紐付け(`.src`) -/

def liftToOpenOfTop.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——開集合を経由する点の分解)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
