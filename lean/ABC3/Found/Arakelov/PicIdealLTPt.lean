import ABC3.Found.Arakelov.PicMulCover

/-!
# Arakelov (B2) 第 196 ブロック —— **点ごとの仮定で足りる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★第 162 の仮定は**強すぎた**

第 162 は `∀ A : X.affineOpens, Module.Invertible Γ(X,A) (D.ideal A)` を仮定するが、
★証明が使うのは**各点につきアフィン開 1 つ**だけである。

    ∀ x, ∃ A アフィン開, x ∈ A ∧ D.ideal A が可逆   ⟹   IsLocallyTrivial

## ★★★これが「可逆性の局所性」への入口になる

`(g_i)` が単位イデアルを生成し各 `D(g_i)` で可逆なら、点ごとの仮定は満たされる。
★★そこから `IsLocallyTrivial`(本ブロック)→ `InvSheaf`(第 182)→ 大域切断が可逆、
と辿れば**局所 ⟹ 大域**が出る筋である。
★★★これは mathlib の TODO(`PicardGroup.lean`「可逆加群の他の特徴づけ」)に当たる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `isLocallyTrivial_idealSheaf_of_pointwise` | ★★★**点ごとの仮定で局所自明** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} (D : X.IdealSheafData)

set_option maxHeartbeats 2000000 in
/-- ★★★**点ごとの仮定で局所自明である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★第 162 の仮定を弱めたもの——証明が使うのは各点につきアフィン開 1 つだけである。 -/
theorem isLocallyTrivial_idealSheaf_of_pointwise
    (hcart : ∀ x : X, ∃ A : X.affineOpens, ∃ _ : x ∈ A.1,
      Module.Invertible (Γ(X, A.1) : Type u) (D.ideal A)) :
    IsLocallyTrivial X (idealSheaf D).val := by
  refine isLocallyTrivial_of_pointwise _ (fun x => ?_)
  obtain ⟨A, hxA, hinv⟩ := hcart x
  haveI := hinv
  obtain ⟨g, hxg, ⟨e⟩⟩ := exists_basicOpen_free D A ⟨x, hxA⟩
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

def isLocallyTrivial_idealSheaf_of_pointwise.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——点ごとの仮定で局所自明)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
