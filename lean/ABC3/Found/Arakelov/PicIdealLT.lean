import ABC3.Found.Arakelov.PicMulCover
import ABC3.Found.Arakelov.PicIdealSheaf
import ABC3.Found.Arakelov.PicPointwise
import ABC3.Found.Arakelov.PicPresieveIso
import ABC3.Found.Arakelov.PicTrivialSheaf

/-!
# Arakelov (B2) 第 162 ブロック —— ★★★★★★★★★**Cartier なら局所自明**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★★B2 の律速が抜けた

    (∀ アフィン開 A, D.ideal A が可逆) ⟹ IsLocallyTrivial X (idealSheaf D).val

★これで `ofDivisor` が書ける——`idealSheaf D` が `InvSheaf` の台になる。

## ★★組み上げ(第 148–162、15 ブロック)

| # | 内容 |
|---|---|
| 148–151 | `IdealSheafData → X.Modules`(mathlib に無い接続) |
| 152–153 | 基本開集合での切断 = イデアルの局所化 |
| 154–155 | 点と素イデアルの対応、アフィン基本開集合の被覆 |
| 156–158 | 切断と局所化加群の同型、全単射性 |
| 159–161 | 可逆性、可換図式、`D(h·g)` 形の被覆 |
| 162 | ★★組み立て |

★★★`tilde M` のとき(第 93–132、**40 ブロック**)に比べて **15 ブロック**で済んだ:
第 130(`bijective_smul_liftGen`)を `R`・`M` について一般に書いておいたのと、
mathlib の `Submodule.toLocalized'` を**測ってから**使ったためである。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} (D : X.IdealSheafData)

set_option maxHeartbeats 2000000 in
theorem isLocallyTrivial_idealSheaf
    (hcart : ∀ A : X.affineOpens, Module.Invertible (Γ(X, A.1) : Type u) (D.ideal A)) :
    IsLocallyTrivial X (idealSheaf D).val := by
  refine isLocallyTrivial_of_pointwise _ (fun x => ?_)
  obtain ⟨U, hU, hxU, -⟩ :=
    Opens.isBasis_iff_nbhd.1 X.isBasis_affineOpens (U := ⊤) (x := x) trivial
  set A : X.affineOpens := ⟨U, hU⟩ with hA
  haveI := hcart A
  obtain ⟨g, hxg, ⟨e⟩⟩ := exists_basicOpen_free D A ⟨x, hxU⟩
  refine ⟨X.basicOpen g, hxg, ⟨?_⟩⟩
  set s0 : ((((restrictPresheafFunctor X (X.basicOpen g)).obj
      (idealSheaf D).val).obj (op (Over.mk (𝟙 (X.basicOpen g))))) : Type u) :=
    idealAwayEquiv D A g (e.symm 1) with hs0
  refine (trivialIsoOfSectionPresieve (X := X) (X.basicOpen g)
    ((restrictPresheafFunctor X (X.basicOpen g)).obj (idealSheaf D).val)
    (isSheaf_restrictModules _ (idealSheaf D)) (isSheaf_unitOn _)
    (s := s0) (h := ?_)).symm
  intro W
  refine ⟨Presieve.functorPullback (Over.forget (X.basicOpen g))
      (fun (Z : X.Opens) (_ : Z ⟶ W.left) =>
        ∃ hh : (Γ(X, A.1) : Type u), Z = X.basicOpen (hh * g)),
    overAffineBasicMulPresieve_mem A g _ le_rfl W, ?_⟩
  rintro Z i ⟨hh, hZ⟩
  obtain ⟨Zl, Zr, Zhom⟩ := Z
  change Zl = X.basicOpen (hh * g) at hZ
  subst hZ
  have hkey2 : @restrictSec X (X.basicOpen g)
      ((restrictPresheafFunctor X (X.basicOpen g)).obj (idealSheaf D).val) s0
      { left := X.basicOpen (hh * g), right := Zr, hom := Zhom }
      (overTermInst _ _)
      = idealAwayEquiv D A (hh * g)
        (liftAwayMap (Γ(X, A.1) : Type u) g hh (D.ideal A) (e.symm 1)) := by
    rw [hs0]
    exact DFunLike.congr_fun (idealAwayEquiv_res D A g hh) (e.symm 1)
  rw [hkey2]
  exact bijective_smul_idealGen D A g hh e


/-! ## ★出典の紐付け(`.src`) -/

def isLocallyTrivial_idealSheaf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——Cartier なら局所自明)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
