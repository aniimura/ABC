import ABC3.Found.Arakelov.PicIdealPre

/-!
# Arakelov (B2) 第 150 ブロック —— **イデアル層の切断は局所的**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★層であることの中身

`idealSections` は `𝒪_X` の部分前層なので、貼り合わせは `𝒪_X` から来る。
★残るのは「貼った切断が条件を満たす」——すなわち**局所性**である:

    U を Wᵢ が覆い、各 s|_{Wᵢ} ∈ idealSections D (Wᵢ) ⟹ s ∈ idealSections D U

## ★★筋——**茎で判定する**

アフィン開 `A ≤ U` を取り、`mem_ideal_iff`(mathlib、`IsAffineOpen`)で
**茎ごと**に判定する。★各点 `x ∈ A` に対し

| 段 | 内容 |
|---|---|
| 1 | `x ∈ Wᵢ` なる `i` を取る |
| 2 | `A` の基本開集合 `B := D(f) ∋ x` で `B ≤ Wᵢ ⊓ A` を取る(`exists_basicOpen_le`) |
| 3 | `s|_B ∈ D.ideal B`(`B ≤ Wᵢ` はアフィン開だから仮定から出る) |
| 4 | `map_ideal` で `D.ideal B = (D.ideal A).map (制限)` |
| 5 | 茎は `B` を経由する(`germ_res`)ので `Ideal.map` が合う |

★★★これで `idealSections` が層の条件を満たす。

## ★★逃げ道——`germ_res` は `erw`

`X.presheaf.germ (↑(X.affineBasicOpen f)) x hxf` の `↑(affineBasicOpen f)` と
`X.basicOpen f` は `rfl` だが**簡約透明度では合わない**ので `rw` が当たらない。
★`erw` で通る([[exact-term-over-rw]] の系統)。
-/

universe u v

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{u}} (D : X.IdealSheafData)

theorem idealSections_of_local {U : X.Opens} {ι : Type v} (W : ι → X.Opens)
    (hle : ∀ i, W i ≤ U) (hcover : U ≤ ⨆ i, W i) (s : (Γ(X, U) : Type u))
    (h : ∀ i, (X.presheaf.map (homOfLE (hle i)).op).hom s ∈ idealSections D (W i)) :
    s ∈ idealSections D U := by
  intro A hAU
  rw [(A.2).mem_ideal_iff]
  intro x hxA
  obtain ⟨i, hxi⟩ := Opens.mem_iSup.1 (hcover (hAU hxA))
  obtain ⟨f, hfle, hxf⟩ :=
    (A.2).exists_basicOpen_le (V := W i ⊓ A.1) ⟨x, ⟨hxi, hxA⟩⟩ hxA
  have hBA : (X.affineBasicOpen f : X.Opens) ≤ A.1 := X.basicOpen_le f
  have hBW : (X.affineBasicOpen f : X.Opens) ≤ W i := le_trans hfle inf_le_left
  have hmem : (X.presheaf.map (homOfLE hBW).op).hom
      ((X.presheaf.map (homOfLE (hle i)).op).hom s) ∈ D.ideal (X.affineBasicOpen f) :=
    h i (X.affineBasicOpen f) hBW
  rw [← CommRingCat.comp_apply, ← Functor.map_comp, ← op_comp] at hmem
  have key : Ideal.map (X.presheaf.germ (A : X.Opens) x hxA).hom (D.ideal A)
      = Ideal.map (X.presheaf.germ (X.affineBasicOpen f : X.Opens) x hxf).hom
        (D.ideal (X.affineBasicOpen f)) := by
    rw [← D.map_ideal hBA, Ideal.map_map, ← CommRingCat.hom_comp]
    erw [X.presheaf.germ_res]
  have heq : (X.presheaf.germ (A : X.Opens) x hxA).hom
        ((X.presheaf.map (homOfLE hAU).op).hom s)
      = (X.presheaf.germ (X.affineBasicOpen f : X.Opens) x hxf).hom
        ((X.presheaf.map (homOfLE (le_trans hBW (hle i))).op).hom s) := by
    rw [← CommRingCat.comp_apply, ← CommRingCat.comp_apply]
    congr 1
    erw [X.presheaf.germ_res, X.presheaf.germ_res]
  rw [key, heq]
  exact Ideal.mem_map_of_mem _ hmem


/-! ## ★出典の紐付け(`.src`) -/

def idealSections_of_local.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——イデアル層の切断は局所的)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
