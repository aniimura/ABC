import ABC3.Found.Arakelov.PicGammaInv

/-!
# Arakelov (B1) 第 145 ブロック —— ★★★★★★★★★★**`Pic(Spec R) ≃* Pic R`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★★★`PicardData` の**最後の欄**

`Interface` の `PicardData` は 14 欄あり、13 欄は既に埋まっていた。
★残る 1 欄が `equivPicRing`——「アフィンでは mathlib の `Pic R` と一致する」である。

★★これは**退化を殺す**条件である:`Pic := PUnit` では
`Pic (Spec ℤ[√-5]) ≃* Pic ℤ[√-5]`(位数 2)が成り立たない。

## ★★両向き

| 向き | 中身 |
|---|---|
| `toPicRing` | `L ↦ Pic.mk (Γ L.carrier)`(★第 144 で可逆性) |
| `ofPicRing` | `M ↦ mk (tilde M)`(★第 133) |
| 左逆 | `tilde (Γ Lc) ≅ Lc`(★第 143) |
| 右逆 | `Γ (tilde M) ≅ M`(mathlib の `tildeGammaIso`) |
| 乗法性 | `Γ(A ⊗ B) ≅ Γ A ⊗ Γ B`(★第 143 + 第 91) |

## ★★★★ここまでの依存

    第 91(tilde はテンソルを保つ)      ← mathlib に無い
    第 132(tilde M は局所自明)         ← 40 ブロック
    第 143(局所自明 ⟹ F ≅ (Γ F)~)     ← mathlib の TODO
      ⟹ 第 145
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace
open scoped TensorProduct

variable (R : CommRingCat.{u})


/-- ★★局所自明な層のテンソル積の大域切断。 -/
noncomputable def gammaTensorIsoGen (A B : (Spec R).Modules)
    (hA : IsLocallyTrivial (Spec R) A.val) (hB : IsLocallyTrivial (Spec R) B.val) :
    AlgebraicGeometry.moduleSpecΓFunctor.obj (tensorModules A B)
      ≅ ModuleCat.of (R : Type u)
        ((AlgebraicGeometry.moduleSpecΓFunctor.obj A : Type u) ⊗[(R : Type u)]
         (AlgebraicGeometry.moduleSpecΓFunctor.obj B : Type u)) :=
  AlgebraicGeometry.moduleSpecΓFunctor.mapIso
      (tensorModulesIso (tildeGammaIsoOfTrivial R A hA).symm (tildeGammaIsoOfTrivial R B hB).symm
        ≪≫ tildeTensorIso R _ _)
    ≪≫ ((tildeGammaIso R).app _).symm

/-- ★★可逆層の類から可逆加群の類へ。 -/
noncomputable def toPicRing : PicSheaf (Spec R) → CommRing.Pic (R : Type u) :=
  Quotient.lift
    (fun L : InvSheaf (Spec R) =>
      letI := invertible_gammaCarrier R L
      CommRing.Pic.mk (R : Type u)
        (AlgebraicGeometry.moduleSpecΓFunctor.obj L.carrier : Type u))
    (by
      rintro L M ⟨e⟩
      letI := invertible_gammaCarrier R L
      letI := invertible_gammaCarrier R M
      exact CommRing.Pic.mk_eq_mk_iff.2
        ⟨(AlgebraicGeometry.moduleSpecΓFunctor.mapIso e).toLinearEquiv⟩)

theorem toPicRing_one : toPicRing R 1 = 1 := by
  letI := invertible_gammaCarrier R (InvSheaf.one (Spec R))
  refine CommRing.Pic.mk_eq_one_iff.2 ⟨?_⟩
  exact (AlgebraicGeometry.moduleSpecΓFunctor.mapIso (tildeUnit R).symm
    ≪≫ ((tildeGammaIso R).app (ModuleCat.of (R : Type u) (R : Type u))).symm).toLinearEquiv

theorem toPicRing_mul (x y : PicSheaf (Spec R)) :
    toPicRing R (x * y) = toPicRing R x * toPicRing R y := by
  induction x using Quotient.inductionOn with
  | h L =>
    induction y using Quotient.inductionOn with
    | h M =>
      letI := invertible_gammaCarrier R L
      letI := invertible_gammaCarrier R M
      letI := invertible_gammaCarrier R (InvSheaf.mul L M)
      show CommRing.Pic.mk (R : Type u)
          (AlgebraicGeometry.moduleSpecΓFunctor.obj (InvSheaf.mul L M).carrier : Type u)
        = CommRing.Pic.mk (R : Type u)
            (AlgebraicGeometry.moduleSpecΓFunctor.obj L.carrier : Type u)
          * CommRing.Pic.mk (R : Type u)
            (AlgebraicGeometry.moduleSpecΓFunctor.obj M.carrier : Type u)
      rw [← CommRing.Pic.mk_tensor]
      exact CommRing.Pic.mk_eq_mk_iff.2
        ⟨(gammaTensorIsoGen R L.carrier M.carrier L.trivial M.trivial).toLinearEquiv⟩

/-- ★★可逆加群の類から可逆層の類へ。 -/
noncomputable def ofPicRing (M : CommRing.Pic (R : Type u)) : PicSheaf (Spec R) :=
  PicSheaf.mk (invSheafOfModule R (ModuleCat.of (R : Type u) (M : Type u)))

theorem toPicRing_ofPicRing (M : CommRing.Pic (R : Type u)) :
    toPicRing R (ofPicRing R M) = M := by
  letI := invertible_gammaCarrier R (invSheafOfModule R (ModuleCat.of (R : Type u) (M : Type u)))
  refine Eq.trans ?_ CommRing.Pic.mk_eq_self
  exact CommRing.Pic.mk_eq_mk_iff.2
    ⟨((tildeGammaIso R).app (ModuleCat.of (R : Type u) (M : Type u))).symm.toLinearEquiv⟩

theorem ofPicRing_toPicRing (x : PicSheaf (Spec R)) :
    ofPicRing R (toPicRing R x) = x := by
  induction x using Quotient.inductionOn with
  | h L =>
    letI := invertible_gammaCarrier R L
    refine PicSheaf.mk_eq_mk ?_
    exact (tildeFunctor R).mapIso
        (LinearEquiv.toModuleIso
          (CommRing.Pic.mk.linearEquiv (R : Type u)
            (AlgebraicGeometry.moduleSpecΓFunctor.obj L.carrier : Type u)))
      ≪≫ tildeGammaIsoOfTrivial R L.carrier L.trivial

/-- ★★★★★★★**`Pic(Spec R) ≃* Pic R`**。 -/
noncomputable def equivPicRingSheaf : PicSheaf (Spec R) ≃* CommRing.Pic (R : Type u) where
  toFun := toPicRing R
  invFun := ofPicRing R
  left_inv := ofPicRing_toPicRing R
  right_inv := toPicRing_ofPicRing R
  map_mul' := toPicRing_mul R


/-! ## ★出典の紐付け(`.src`) -/

def equivPicRingSheaf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——Pic(Spec R) ≃* Pic R)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
