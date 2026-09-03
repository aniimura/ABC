/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.LocalizationBridge
import ABC3.Found.GenEll.ComapAffine
import ABC3.Found.GenEll.PullbackNatural

/-!
# PullbackLocalization —— `[GenEll] Definition 1.2` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open CategoryTheory AlgebraicGeometry

/-! ## ★★★★★★任意の環の上の引き戻しイデアル -/

/-- ★★★★**任意の可換環 `B` について、`Spec B` の点に沿った因子の引き戻し**。

★`pullbackIdeal F D xF` はこの `B = 𝓞_F` の場合であり、**定義が同じ**である。 -/
noncomputable def pullbackIdealOf (B : CommRingCat.{0}) {X : Scheme.{0}}
    (D : X.IdealSheafData) (x : Spec B ⟶ X) : Ideal B :=
  (Scheme.IdealSheafData.equivOfIsAffine (D.comap x)).comap (Scheme.ΓSpecIso B).inv.hom

/-- ★既存の `pullbackIdeal` はその特別な場合。 -/
theorem pullbackIdealOf_eq_pullbackIdeal (F : Type) [Field F] [NumberField F]
    {X : Scheme.{0}} (D : X.IdealSheafData) (xF : specRingOfIntegers F ⟶ X) :
    pullbackIdealOf (CommRingCat.of (NumberField.RingOfIntegers F)) D xF
      = pullbackIdeal F D xF := rfl

/-- ★★★★★★★★**任意の環準同型に沿って拡大イデアルになる**。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★`pullbackIdeal_specMap`（数体しか扱えない）の一般化である。
★★証明は `pullbackIdeal_comp` の骨格そのまま——
`ideal_comap_eq_map_of_isAffine` が数体を仮定していないので通る。 -/
theorem pullbackIdealOf_specMap {B B' : CommRingCat.{0}} {X : Scheme.{0}}
    (D : X.IdealSheafData) (x : Spec B ⟶ X) (φ : B ⟶ B') :
    pullbackIdealOf B' D (Spec.map φ ≫ x) = (pullbackIdealOf B D x).map φ.hom := by
  have keyB : ∀ J : Ideal Γ(Spec B, ⊤),
      Ideal.comap (Scheme.ΓSpecIso B).inv.hom J = J.map (Scheme.ΓSpecIso B).hom.hom :=
    fun J => Ideal.comap_symm (Scheme.ΓSpecIso B).commRingCatIsoToRingEquiv
  have keyB' : ∀ J : Ideal Γ(Spec B', ⊤),
      Ideal.comap (Scheme.ΓSpecIso B').inv.hom J = J.map (Scheme.ΓSpecIso B').hom.hom :=
    fun J => Ideal.comap_symm (Scheme.ΓSpecIso B').commRingCatIsoToRingEquiv
  show Ideal.comap _ (Scheme.IdealSheafData.equivOfIsAffine (D.comap (Spec.map φ ≫ x))) = _
  rw [Scheme.IdealSheafData.comap_comp]
  show Ideal.comap _ (((D.comap x).comap (Spec.map φ)).ideal ⟨⊤, isAffineOpen_top _⟩) = _
  rw [ideal_comap_eq_map_of_isAffine, keyB']
  show _ = Ideal.map _ (Ideal.comap _ (Scheme.IdealSheafData.equivOfIsAffine (D.comap x)))
  rw [keyB]
  simp only [Ideal.map_map, Scheme.IdealSheafData.equivOfIsAffine_apply]
  refine Eq.trans (Ideal.map_map _ _) ?_
  have hmor : (Spec.map φ).app ⊤ ≫ (Scheme.ΓSpecIso B').hom
      = (Scheme.ΓSpecIso B).hom ≫ φ := Scheme.ΓSpecIso_naturality φ
  have hring := congrArg CommRingCat.Hom.hom hmor
  rw [CommRingCat.hom_comp, CommRingCat.hom_comp] at hring
  exact congrArg (fun ψ => Ideal.map ψ ((D.comap x).ideal ⟨⊤, isAffineOpen_top _⟩)) hring

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` はまだ置かない。** 残っているのは
「`ConductorDescent.lean` の `X_n` 上の因子の一致から、本定理の仮定 `h` を作る」段
——`Spec 𝓞_F[1/N] ⟶ X_n` を作る配線（`𝓞_F[1/N]` は `ℤ[1/n!]`-代数である）である。 -/

def pullbackIdealOf_specMap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2, (i)(引き戻しの底変換——任意の可換環の上で)",
    sectionId := "genell-def-1-2-i" }

end ABC3.Found.GenEll
