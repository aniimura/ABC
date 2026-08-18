import ABC3.Found.Arakelov.PicDualHom

/-!
# Arakelov (B2) 第 179 ブロック —— ★★★★★★★**双対の前層は層である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★B2 の山を越えた

§9-200 で採った **A の道**が完走した。

    Hom(F|_U, 𝒪|_U) は U について層である

mathlib には**加群の層の内部 Hom が無い**(`Sheaf/Monoidal.lean` は存在しない)ので、
第 163–179 の 17 ブロックで**双対を一から作った**ことになる。

## ★★証明の骨

| 歩 | 内容 | ブロック |
|---|---|---|
| 1 | 局所値が貼り合わせ可能 | 第 176 |
| 2 | 貼り合わせた値が加法的・線型・自然 | 第 177 |
| 3 | 束ねて**射**にする(`homMk`) | 第 178 |
| 4 | ★**貼り合わせ性**と**一意性** | 本ブロック |

★★★どちらも「`Z.left ⊓ Uᵢ` へ制限して局所値に落とす」だけである
——制限は `dual_app_res`(切断の自然性)、一致は `dual_eq_at`(仮定)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `isGluing_dualGlueHom` | ★★★貼り合わせた射は制限すると `φ i` に戻る |
| `dualGlue_uniq` | ★★★そのような射は一意 |
| `isSheaf_dualPresheaf` | ★★★★★★★**双対の前層は層である** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} {F : X.PresheafOfModules} {ι : Type u} {U : ι → X.Opens}
  (φ : ∀ i, ((dualPresheaf F).obj (op (U i)) : Type u))
  (hφ : TopCat.Presheaf.IsCompatible (dualPresheaf F).presheaf U φ)

/-- ★★★**貼り合わせた射は制限すると `φ i` に戻る**。 -/
theorem isGluing_dualGlueHom :
    TopCat.Presheaf.IsGluing (dualPresheaf F).presheaf U φ
      (dualGlueHom φ hφ (le_refl (⨆ i, U i))) := by
  intro i
  refine dualHom_ext _ _ (fun Z x => ?_)
  have key : ∀ (h : Z.left ≤ ⨆ i, U i),
      dualGlueVal φ hφ h x = ((φ i).app (op Z)).hom x := by
    intro h
    refine dualGlue_ext h _ _ (fun j => ?_)
    have hij : Z.left ⊓ U j ≤ U i := inf_le_left.trans (leOfHom Z.hom)
    refine Eq.trans (dualGlueVal_res φ hφ h x j) ?_
    refine Eq.trans (dualLocVal_apply φ Z.left x j) ?_
    refine Eq.trans (dual_eq_at φ hφ j i inf_le_right hij _) ?_
    exact (dual_app_res (φ i)
      (Over.homMk (homOfLE (inf_le_left : Z.left ⊓ U j ≤ Z.left)) :
        Over.mk (homOfLE hij) ⟶ Z) x).symm
  exact key ((leOfHom Z.hom).trans (le_iSup U i))

/-- ★★★**貼り合わせた射は一意である**。 -/
theorem dualGlue_uniq (t : ((dualPresheaf F).obj (op (⨆ i, U i)) : Type u))
    (ht : TopCat.Presheaf.IsGluing (dualPresheaf F).presheaf U φ t) :
    t = dualGlueHom φ hφ (le_refl (⨆ i, U i)) := by
  refine dualHom_ext _ _ (fun Z x => ?_)
  rw [dualGlueHom_app]
  refine dualGlue_ext ((leOfHom Z.hom).trans (le_refl (⨆ i, U i))) _ _ (fun i => ?_)
  rw [dualGlueVal_res, dualLocVal_apply]
  have hZi : Z.left ⊓ U i ≤ ⨆ j, U j := (inf_le_left.trans (leOfHom Z.hom)).trans le_rfl
  refine Eq.trans (dual_app_res t
    (Over.homMk (homOfLE (inf_le_left : Z.left ⊓ U i ≤ Z.left)) :
      Over.mk (homOfLE hZi) ⟶ Z) x) ?_
  exact congrArg
    (fun (ψ : ((dualPresheaf F).obj (op (U i)) : Type u)) =>
      ((ψ.app (op (infOver U Z.left i))).hom
        ((F.map (homOfLE (inf_le_left : Z.left ⊓ U i ≤ Z.left)).op).hom x))) (ht i)

/-- ★★★★★★★**双対の前層は層である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★★★mathlib には加群の層の**内部 Hom が無い**ので、第 163–179 の 17 ブロックで
双対を一から作った。これが B2 の律速であった。 -/
theorem isSheaf_dualPresheaf (F : X.PresheafOfModules) :
    Presheaf.IsSheaf (Opens.grothendieckTopology X) (dualPresheaf F).presheaf := by
  refine (TopCat.Presheaf.isSheaf_iff_isSheafUniqueGluing _).2 ?_
  intro ι U sf hsf
  exact ⟨dualGlueHom sf hsf (le_refl (⨆ i, U i)), isGluing_dualGlueHom sf hsf,
    fun t ht => dualGlue_uniq sf hsf t ht⟩

/-! ## ★出典の紐付け(`.src`) -/

def isSheaf_dualPresheaf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——双対の前層が層であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
