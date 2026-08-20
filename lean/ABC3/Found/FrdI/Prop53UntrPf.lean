/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55PfUnit

/-!
# `⟨A,n⟩` の End は `𝒞^pf` の End の共役(鎖 `prop55` の `p53-untrpf` の (a))

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.103。

原文 (FrdI p.103):
> the Frobenioid C the model Frobenioid [cf. Theorem 5.2, (ii)] associated to

## ★★`compRoot` は `compPf` の共役である

`Prop55PfUnit.lean` の `otimes_pfRoot_eq_bot` は `⟨A,1⟩` についてしか言っていない。
★すべての `⟨A,n⟩` へ広げるには、**`compRoot` の定義を開く**のが要点である。

`compRoot_eq_lift`(`Def31Pf.lean`)を `c = 1`、`PA = PB = PE = n·n`、
`ef = eg = er = n` で使うと、`X = Y = Z = ⟨A,n⟩` のとき

  `compRoot f g = Θ.hom (compPf (Θ.inv f) (Θ.inv g))`,
  `Θ := rtRootIso A A (n·n = n·n) (n·n = n·n)`

★★つまり **`Θ.inv` は単系の同型 `End ⟨A,n⟩ ≅ Hom^pf(A^{(n·n)}, A^{(n·n)})`** である。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `rootSelfIso` | ★`Θ`(根 `n·n` から根 `n` へ) |
| `rootSelfIso_inv_compRoot` | ★★**`Θ.inv` は合成を `compPf` へ移す** |
| `rootSelfIso_inv_id` | `Θ.inv` は単位を単位へ |
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u v2 u2 w

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

/-- ★**`Θ`** —— 根 `n·n` の `Hom^pf` から根 `n` の `Hom^pf` への同型。 -/
noncomputable def rootSelfIso (A : C) (n : ℕ+) :
    HomPf P F (rtObj P F A (n * n)) (rtObj P F A (n * n))
      ≅ HomPf P F (rtObj P F A n) (rtObj P F A n) :=
  rtRootIso P F A A (show n * n = n * n from rfl) (show n * n = n * n from rfl)

/-- ★★★**`Θ.inv` は `compRoot` を `compPf` へ移す**。 -/
theorem rootSelfIso_inv_compRoot (A : C) (n : ℕ+)
    (f g : HomRoot P F (⟨A, n⟩ : PfRootObj P F) ⟨A, n⟩) :
    (rootSelfIso (F := F) A n).inv (compRoot P F f g)
      = compPf P F ((rootSelfIso (F := F) A n).inv f) ((rootSelfIso (F := F) A n).inv g) := by
  have hlift := compRoot_eq_lift (P := P) (F := F) (X := (⟨A, n⟩ : PfRootObj P F))
    (Y := ⟨A, n⟩) (Z := ⟨A, n⟩) f g
    (c := 1) (PA := n * n) (PB := n * n) (PE := n * n)
    (hcA := (one_mul _).symm) (hcB := (one_mul _).symm) (hcE := (one_mul _).symm)
    (ef := n) (eg := n) (er := n)
    (hfA := rfl) (hfB := rfl) (hgA := rfl) (hgE := rfl) (hrA := rfl) (hrE := rfl)
  rw [hlift]
  exact Iso.hom_inv_id_apply _ _

/-- ★**`Θ.inv` は単位を単位へ**。 -/
theorem rootSelfIso_inv_id (A : C) (n : ℕ+) :
    (rootSelfIso (F := F) A n).inv (idRoot P F (⟨A, n⟩ : PfRootObj P F))
      = toHomPf (F := F) (𝟙 (rtObj P F A (n * n))) :=
  rtRootIso_inv_id P F A (show n * n = n * n from rfl)

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Proposition 5.3` の「`(𝒞^un-tr)^pf` は unit-trivial」をすべての
`⟨A,n⟩` へ広げる段(第 1 段)。 -/
def rootSelfIso_inv_compRoot.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 103,
    item := "Proposition 5.3 — ⟨A,n⟩ の End は 𝒞^pf の End の共役",
    sectionId := "frdi-prop-5-3" }

end ABC3.Found.FrdI
