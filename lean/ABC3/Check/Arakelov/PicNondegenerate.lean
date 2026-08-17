import ABC3.Interface.Arakelov.LineBundle

/-!
# 負の対照 —— **`PicardData` の穴が実在したことの証明**(`Check`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

## ★★★★何を確かめるか

2026-08-17 に `PicardData` の穴を見つけた。★当初の要求

    Pic / group / pullback / pullback_mul / pullback_id / pullback_comp / equivPicRing

は、**大域切断の Picard 群**

    Pic X := CommRing.Pic Γ(X, ⊤)

で**すべて満たされる**。しかしこれは★★**非アフィンな `X` では数学的に誤り**である
——`Pic(ℙ¹) = ℤ` だが `Γ(ℙ¹, 𝒪) = k` なので `CommRing.Pic k = 0`。

★★★**本ファイルはそれを機械で確かめる**——下の 4 定理が通ること自体が、
「アフィンでの一致だけでは `Pic` は決まらない」ことの証明である。

★塞いだ手は `sheafOf` / `sheafOf_injective` / `sheafOf_surjective`
(`Interface/Arakelov/LineBundle.lean`)。
-/

universe u

namespace ABC3.Check.Arakelov

open AlgebraicGeometry CategoryTheory

/-! ## ★★誤った候補 —— 大域切断の Picard 群 -/

/-- ★**大域切断の Picard 群**(誤った `Pic` の候補)。 -/
def gammaPic (X : Scheme.{0}) : Type :=
  CommRing.Pic (X.presheaf.obj (Opposite.op ⊤))

noncomputable instance (X : Scheme.{0}) : CommGroup (gammaPic X) :=
  inferInstanceAs (CommGroup (CommRing.Pic (X.presheaf.obj (Opposite.op ⊤))))

/-- ★引き戻し(`Γ` は反変なので向きは正しい)。 -/
noncomputable def gammaPicPullback {X Y : Scheme.{0}} (f : X ⟶ Y) :
    gammaPic Y →* gammaPic X :=
  CommRing.Pic.mapRingHom f.appTop.hom

/-! ## ★★★誤った候補が旧要求をすべて満たすこと -/

/-- ★**関手性(恒等射)**——旧 `pullback_id` を満たす。 -/
theorem gammaPicPullback_id (X : Scheme.{0}) (L : gammaPic X) :
    gammaPicPullback (𝟙 X) L = L := by
  show CommRing.Pic.mapRingHom (Scheme.Hom.appTop (𝟙 X)).hom L = L
  have h : (Scheme.Hom.appTop (𝟙 X)).hom = RingHom.id _ := by
    simp
  rw [h, CommRing.Pic.mapRingHom_id_apply]

/-- ★**関手性(合成)**——旧 `pullback_comp` を満たす。 -/
theorem gammaPicPullback_comp {X Y Z : Scheme.{0}} (f : X ⟶ Y) (g : Y ⟶ Z)
    (L : gammaPic Z) :
    gammaPicPullback (f ≫ g) L = gammaPicPullback f (gammaPicPullback g L) := by
  show CommRing.Pic.mapRingHom (Scheme.Hom.appTop (f ≫ g)).hom L
      = CommRing.Pic.mapRingHom (Scheme.Hom.appTop f).hom
          (CommRing.Pic.mapRingHom (Scheme.Hom.appTop g).hom L)
  rw [CommRing.Pic.mapRingHom_mapRingHom]
  congr 1

/-- ★**群準同型**——旧 `pullback_mul` を満たす(`mapRingHom` が `MonoidHom` だから)。 -/
theorem gammaPicPullback_mul {X Y : Scheme.{0}} (f : X ⟶ Y) (L M : gammaPic Y) :
    gammaPicPullback f (L * M) = gammaPicPullback f L * gammaPicPullback f M :=
  map_mul _ _ _

/-! ## ★★★★アフィンでは正しい値を与えてしまう -/

/-- ★★★★**`Spec R` では mathlib の `Pic R` と一致する**——旧 `equivPicRing` を満たす。

★★★**これが穴の核心である。**アフィンで一致するだけでは
非アフィンでの値は何も決まらない。 -/
noncomputable def gammaPicEquivSpec (R : CommRingCat.{0}) :
    gammaPic (Spec R) ≃* CommRing.Pic R where
  toFun := CommRing.Pic.mapRingHom (Scheme.ΓSpecIso R).hom.hom
  invFun := CommRing.Pic.mapRingHom (Scheme.ΓSpecIso R).inv.hom
  left_inv L := by
    show CommRing.Pic.mapRingHom _ (CommRing.Pic.mapRingHom _ L) = L
    rw [CommRing.Pic.mapRingHom_mapRingHom]
    have h : ((Scheme.ΓSpecIso R).inv.hom).comp ((Scheme.ΓSpecIso R).hom.hom)
        = RingHom.id _ := by
      ext x
      exact congrArg (fun m => CommRingCat.Hom.hom m x) (Scheme.ΓSpecIso R).hom_inv_id
    rw [h, CommRing.Pic.mapRingHom_id_apply]
  right_inv L := by
    show CommRing.Pic.mapRingHom _ (CommRing.Pic.mapRingHom _ L) = L
    rw [CommRing.Pic.mapRingHom_mapRingHom]
    have h : ((Scheme.ΓSpecIso R).hom.hom).comp ((Scheme.ΓSpecIso R).inv.hom)
        = RingHom.id _ := by
      ext x
      exact congrArg (fun m => CommRingCat.Hom.hom m x) (Scheme.ΓSpecIso R).inv_hom_id
    rw [h, CommRing.Pic.mapRingHom_id_apply]
  map_mul' _ _ := map_mul _ _ _

end ABC3.Check.Arakelov
