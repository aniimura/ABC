import ABC3.Found.Arakelov.PicTrivialNoSheaf
import ABC3.Found.Arakelov.PicDualPre
import ABC3.Found.Arakelov.PicResTrans

/-!
# Arakelov (B2) 第 171 ブロック —— ★★★★★★★**双対も局所自明**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★逆層の局所自明性が取れた

    F|_V ≅ 𝟙_  ⟹  (dualPresheaf F)|_V ≅ 𝟙_

★★切断として **`e.hom` 自身**を取るのが要点である
——`e.hom : F|_V ⟶ 𝟙_` は**双対の切断そのもの**だからである。

## ★★全単射性の筋

    c ↦ c • (e.hom|_W) = (e.hom|_W) ≫ unitMul W c

は 2 つの全単射の合成である:

| 段 | 全単射 |
|---|---|
| `c ↦ unitMul W c` | ★第 166(`unitEndEquiv`) |
| `ψ ↦ (e.hom|_W) ≫ ψ` | ★`e.hom|_W` が同型だから(`Iso.homCongr`) |

★★★層の仮定は要らない(★第 170)——すべての `W` で全単射だからである。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} (F : X.PresheafOfModules)

theorem bijective_dual_smul (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V)) (W : Over V) :
    Function.Bijective
      (fun c : ((((Over.forget V).op ⋙ X.presheaf)
          ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj (op W) : Type u) =>
        c • restrictSec V ((restrictPresheafFunctor X V).obj (dualPresheaf F)) e.hom W) := by
  have hrs : restrictSec V ((restrictPresheafFunctor X V).obj (dualPresheaf F)) e.hom W
      = (restrictOnFunctor (leOfHom W.hom)).map e.hom := rfl
  have hfun : (fun c : ((((Over.forget V).op ⋙ X.presheaf)
        ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj (op W) : Type u) =>
      c • restrictSec V ((restrictPresheafFunctor X V).obj (dualPresheaf F)) e.hom W)
      = (fun ψ : (𝟙_ (PresheafModulesOn X W.left) ⟶ 𝟙_ (PresheafModulesOn X W.left)) =>
          (restrictSec V ((restrictPresheafFunctor X V).obj (dualPresheaf F)) e.hom W) ≫ ψ)
        ∘ (unitMul W.left) := by
    funext c
    exact dual_smul_eq F W.left c _
  rw [hfun]
  refine Function.Bijective.comp ?_ ?_
  · haveI hiso : IsIso (restrictSec V
        ((restrictPresheafFunctor X V).obj (dualPresheaf F)) e.hom W) := by
      rw [hrs]
      exact ((restrictOnFunctor (leOfHom W.hom)).mapIso e).isIso_hom
    have hc : (fun ψ : (𝟙_ (PresheafModulesOn X W.left) ⟶ 𝟙_ (PresheafModulesOn X W.left)) =>
        (restrictSec V ((restrictPresheafFunctor X V).obj (dualPresheaf F)) e.hom W) ≫ ψ)
        = (Iso.homCongr (asIso (restrictSec V
            ((restrictPresheafFunctor X V).obj (dualPresheaf F)) e.hom W)) (Iso.refl _)).symm := by
      funext ψ
      simp [Iso.homCongr]
    rw [hc]
    exact Equiv.bijective _
  · have hu : unitMul W.left = fun c => (unitEndEquiv W.left).symm c := by
      funext c
      exact (unitEndEquiv_symm_apply W.left c).symm
    rw [hu]
    exact (unitEndEquiv W.left).symm.bijective

/-- ★★★★双対も局所自明。 -/
noncomputable def dualTrivialIso (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V)) :
    𝟙_ (PresheafModulesOn X V) ≅ (restrictPresheafFunctor X V).obj (dualPresheaf F) :=
  trivialIsoOfSection' V ((restrictPresheafFunctor X V).obj (dualPresheaf F)) e.hom
    (bijective_dual_smul F V e)


/-! ## ★出典の紐付け(`.src`) -/

def dualTrivialIso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——双対も局所自明)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
