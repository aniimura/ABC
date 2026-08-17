import Mathlib.AlgebraicGeometry.Modules.Sheaf
import Mathlib.Algebra.Category.ModuleCat.Presheaf.Monoidal
import Mathlib.Algebra.Category.ModuleCat.Sheaf.LocallyFree
import ABC3.Meta.Claim

/-!
# Arakelov (B1) 第 1 ブロック —— **スキーム上の前層加群のテンソル積**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★なぜこれが要るか

B1(`PicardData`)は `Pic X` に **`CommGroup`** を要求する。その乗法は
**可逆層のテンソル積**である。★2026-08-17 の実測で、

| 何 | mathlib |
|---|---|
| **前層**加群の対称モノイダル構造 | ★**在る**(`ModuleCat/Presheaf/Monoidal.lean`) |
| **層**加群のモノイダル構造 | ★★★**無い**(`SheafOfModules` に `⊗` は 0 件) |

したがって B1 の底は「前層のテンソル積を層に持ち上げる」ことである。

## ★★本ファイルが渡す橋

★★★**インスタンスが自動では見つからない**——`X.PresheafOfModules` は
`PresheafOfModules X.ringCatSheaf.obj` であり、mathlib のモノイダル構造は
`PresheafOfModules (R ⋙ forget₂ CommRingCat RingCat)` の形に付いている。
★両者は **`rfl` で等しい**が、インスタンス探索は `X.ringCatSheaf.obj` を
`⋙ forget₂` の形まで展開しない(2026-08-17 実測)。

★本ファイルは `inferInstanceAs` で橋を渡し、`X.PresheafOfModules` の上で
`⊗` が使えるようにする。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory

variable (X : Scheme.{u})

/-! ## ★★モノイダル構造の橋渡し -/

/-- ★**スキーム上の前層加群の圏はモノイダル圏である**。

★機構は mathlib の `PresheafOfModules.Monoidal.monoidalCategory` を
`X.presheaf ⋙ forget₂ CommRingCat RingCat` に当てるだけ。
★★ただし**インスタンス探索は自動では届かない**ので、明示的に渡す。 -/
noncomputable instance presheafOfModulesMonoidalStruct :
    MonoidalCategoryStruct X.PresheafOfModules :=
  inferInstanceAs (MonoidalCategoryStruct
    (PresheafOfModules.{u} (X.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u})))

noncomputable instance presheafOfModulesMonoidal :
    MonoidalCategory X.PresheafOfModules :=
  inferInstanceAs (MonoidalCategory
    (PresheafOfModules.{u} (X.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u})))

/-- ★★**しかも対称モノイダル圏である**(テンソル積は可換)。

★★★これが `Pic X` の **`CommGroup`**(可換群)を出すもとになる。 -/
noncomputable instance presheafOfModulesSymmetric :
    SymmetricCategory X.PresheafOfModules :=
  inferInstanceAs (SymmetricCategory
    (PresheafOfModules.{u} (X.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u})))

/-! ## ★切断の計算 -/

variable {X}

/-- ★**テンソル積の切断は切断のテンソル積である**(前層だから各開集合ごと)。

★★★これが「層であること」を確かめる足場になる——
右辺は各 `U` ごとの加群のテンソル積なので、局所自由なら計算できる。 -/
theorem presheafTensor_obj (M N : X.PresheafOfModules) (U : X.Opensᵒᵖ) :
    (M ⊗ N).obj U
      = MonoidalCategory.tensorObj (C := ModuleCat (X.presheaf.obj U)) (M.obj U) (N.obj U) :=
  PresheafOfModules.Monoidal.tensorObj_obj (R := X.presheaf) M N U

/-! ## ★★層の側へ渡す道具 -/

/-- ★**層加群からその台の前層加群を取る**(記法の短縮)。 -/
noncomputable abbrev val (M : X.Modules) : X.PresheafOfModules := M.val

/-- ★★**2 つの層加群の前層テンソル積**。

★★★**これが層であるかどうかが B1 の底である。**
一般には層ではないが、**局所自由有限階数なら層である**
(局所的に `𝒪^n` の有限直和で、層であることは局所的性質だから)。 -/
noncomputable abbrev presheafTensor (M N : X.Modules) : X.PresheafOfModules :=
  M.val ⊗ N.val

/-- ★前層テンソル積の切断も、切断のテンソル積である。 -/
theorem presheafTensor_obj_eq (M N : X.Modules) (U : X.Opensᵒᵖ) :
    (presheafTensor M N).obj U
      = MonoidalCategory.tensorObj (C := ModuleCat (X.presheaf.obj U))
          (M.val.obj U) (N.val.obj U) :=
  presheafTensor_obj M.val N.val U

/-! ## ★出典の紐付け(`.src`) -/

def presheafOfModulesMonoidal.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——スキーム上の前層加群のテンソル積)",
    sectionId := "genell-def-1-1-i" }

def presheafTensor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——2 つの層加群の前層テンソル積)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
