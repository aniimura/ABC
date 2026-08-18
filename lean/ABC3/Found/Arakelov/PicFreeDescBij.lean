import ABC3.Found.Arakelov.PicFreeUnique

/-!
# Arakelov (B1) 第 106 ブロック —— **一点添字の `freeDesc` の全単射性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★局所自明性の全単射判定が閉じた

    freeDesc φ : free(ι) ⟶ M   (ι は一点)

は、**生成元の像による乗法**と同じ全単射性を持つ:

    freeDesc φ ∘ (階数 1 の同型).symm  =  (c ↦ c • φ default)

★★これで第 103(生成元の乗法は全単射)が直に効く。

## ★★詰まりの回避(記録)

★`ModuleCat.freeDesc` は `φ : ι ⟶ ↑M`(**圏の射**)を要求する。
`φ : ι → ↑M`(**関数**)で書くと**型が合わない**——`Type` の圏では同じものだが、
Lean は `⟶` の形を要求する(実測)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `freeDesc_uniqueLinearEquiv_symm` | ★★合成は「生成元の倍」 |
| `bijective_freeDesc_of_unique` | ★★★★**一点添字の `freeDesc` の全単射性** |
-/

universe u

namespace ABC3.Found.Arakelov

open CategoryTheory

variable {A : Type u} [CommRing A] {ι : Type u} [Unique ι] {M : ModuleCat.{u} A}

/-- ★★**合成は「生成元の倍」である**。 -/
theorem freeDesc_uniqueLinearEquiv_symm (φ : ι ⟶ (M : Type u)) (c : A) :
    (ModuleCat.freeDesc (X := ι) (M := M) φ)
        ((Finsupp.uniqueLinearEquiv A A (default : ι)).symm c)
      = c • φ default := by
  rw [uniqueLinearEquiv_symm_apply]
  rw [show (c • (Finsupp.single (default : ι) (1 : A)))
      = c • (ModuleCat.freeMk (default : ι) : (ModuleCat.free A).obj ι) from rfl]
  rw [map_smul, ModuleCat.freeDesc_apply]

/-- ★★★★**一点添字の `freeDesc` は、生成元の像による乗法と同じ全単射性を持つ**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで第 103 ブロック(生成元の乗法は全単射)が直に効く。 -/
theorem bijective_freeDesc_of_unique (φ : ι ⟶ (M : Type u))
    (h : Function.Bijective (fun c : A => c • φ default)) :
    Function.Bijective (ModuleCat.freeDesc (X := ι) (M := M) φ) := by
  have he : (fun c : A => c • φ default)
      = (fun x => (ModuleCat.freeDesc (X := ι) (M := M) φ) x)
        ∘ (Finsupp.uniqueLinearEquiv A A (default : ι)).symm := by
    funext c
    exact (freeDesc_uniqueLinearEquiv_symm φ c).symm
  rw [he] at h
  have := (Finsupp.uniqueLinearEquiv A A (default : ι)).symm.bijective
  exact (Function.Bijective.of_comp_iff _ this).1 h

/-! ## ★出典の紐付け(`.src`) -/

def bijective_freeDesc_of_unique.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——一点添字の freeDesc の全単射性)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
