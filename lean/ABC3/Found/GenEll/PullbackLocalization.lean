/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.LocalizationBridge
import ABC3.Found.GenEll.ComapAffine
import ABC3.Found.GenEll.PullbackNatural

/-!
# [GenEll] Remark 1.5.1 —— **引き戻しイデアルを任意の環の上で**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★★★★★★なぜ一般化するのか

`pullbackIdeal` は `Spec 𝓞_F` の点についてしか定義されていない。
★`LocalizationBridge.lean` が要求するのは **`𝓞_F[1/N]` 上のイデアルの一致**であり、
`𝓞_F[1/N]` は数体の整数環ではないので、既存の `pullbackIdeal_specMap`
（`K` が数体であることを要求する）は当たらない。

★★そこで**任意の可換環の上**で引き戻しイデアルを定義し直す（`pullbackIdealOf`）。
`pullbackIdeal F D xF` はその特別な場合であり、**`rfl` で一致する**。

## ★★★★機構

`pullbackIdeal_comp`（`PullbackBase.lean`）の証明がそのまま一般化する。
★鍵は `ideal_comap_eq_map_of_isAffine`（`ComapAffine.lean`）が**既に一般**で、
数体を仮定していないこと。
★★`ΓSpecIso` の往復は `Scheme.ΓSpecIso_naturality` 1 本で吸収される。
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

/-! ## ★★★★★★★★★★`LocalizationBridge` が要求する形 -/

variable (F : Type) [Field F] [NumberField F]

/-- ★★★★★★★★**`𝓞_F[1/N]` へ拡大したイデアルは、`Spec 𝓞_F[1/N]` の点での引き戻しである**。

★★これが `LocalizationBridge.lean` の `conductorADiv_fin_eq_of_map_eq` が要求する
仮定 `h` を作るための配線である。 -/
theorem pullbackIdeal_map_algebraMap (A : Type) [CommRing A]
    [Algebra (NumberField.RingOfIntegers F) A]
    {X : Scheme.{0}} (D : X.IdealSheafData) (xF : specRingOfIntegers F ⟶ X) :
    (pullbackIdeal F D xF).map (algebraMap (NumberField.RingOfIntegers F) A)
      = pullbackIdealOf (CommRingCat.of A) D
          (Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) A)) ≫ xF) := by
  rw [← pullbackIdealOf_eq_pullbackIdeal F D xF]
  have h := (pullbackIdealOf_specMap D xF
    (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) A))).symm
  simpa using h

/-- ★★★★★★★★★★**[GenEll] Remark 1.5.1 —— `hagree` を作る形**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

`Spec 𝓞_F[1/N]` の水準で 2 つの因子の引き戻しが一致するなら、
**`N` を割らない素点 `v` で導手の係数が一致する**。

★★★これが `Skeleton/GenEll/Section1.lean` の `remark_1_5_1` が受けている
仮定 `hagree` そのものである。 -/
theorem conductorADiv_fin_eq_of_localized
    (M : Submonoid (NumberField.RingOfIntegers F)) (A : Type) [CommRing A]
    [Algebra (NumberField.RingOfIntegers F) A] [IsLocalization M A]
    (v : FinitePlace F) (hM : ∀ m ∈ M, m ∉ v.asIdeal)
    {X X' : Scheme.{0}} (D : X.IdealSheafData) (D' : X'.IdealSheafData)
    (xF : specRingOfIntegers F ⟶ X) (xF' : specRingOfIntegers F ⟶ X')
    (hI : pullbackIdeal F D xF ≠ 0) (hJ : pullbackIdeal F D' xF' ≠ 0)
    (h : pullbackIdealOf (CommRingCat.of A) D
          (Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) A)) ≫ xF)
        = pullbackIdealOf (CommRingCat.of A) D'
          (Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) A)) ≫ xF')) :
    (conductorADiv F D xF).fin v = (conductorADiv F D' xF').fin v := by
  refine conductorADiv_fin_eq_of_map_eq F M A v hM D D' xF xF' hI hJ ?_
  rw [pullbackIdeal_map_algebraMap, pullbackIdeal_map_algebraMap, h]

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` はまだ置かない。** 残っているのは
「`ConductorDescent.lean` の `X_n` 上の因子の一致から、本定理の仮定 `h` を作る」段
——`Spec 𝓞_F[1/N] ⟶ X_n` を作る配線（`𝓞_F[1/N]` は `ℤ[1/n!]`-代数である）である。 -/

def pullbackIdealOf_specMap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2, (i)(引き戻しの底変換——任意の可換環の上で)",
    sectionId := "genell-def-1-2-i" }

def conductorADiv_fin_eq_of_localized.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(hagree を作る形——X_n の点への配線は含まない)",
    sectionId := "genell-rem-1-5-1" }

def conductorADiv_fin_eq_of_localized.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "conductorADiv_fin_eq_of_map_eq(局所化の橋)"
      (.inProject "ABC3" "ABC3.Found.GenEll.conductorADiv_fin_eq_of_map_eq") 9,
    .citation "[ABC3]" "ideal_comap_eq_map_of_isAffine(アフィンでの comap は map)"
      (.inProject "ABC3" "ABC3.Found.GenEll.ideal_comap_eq_map_of_isAffine") 5,
    .citation "[mathlib]" "Scheme.ΓSpecIso_naturality(Spec を取って Γ を取ると戻る)"
      (.inMathlib "AlgebraicGeometry.Scheme.ΓSpecIso_naturality") 5,
    .implicitStep
      ("★既存の pullbackIdeal_specMap は K が数体であることを要求するので、" ++
       "𝓞_F[1/N](数体の整数環ではない)には当たらない。" ++
       "★★そこで任意の可換環の上で定義し直した(pullbackIdealOf)。" ++
       "pullbackIdeal はその特別な場合であり rfl で一致する") 5,
    .implicitStep
      ("★★★残る配線: ConductorDescent.lean の X_n 上の因子の一致から本定理の仮定 h を作る。" ++
       "★機構は Spec 𝓞_F[1/N] ⟶ X_n を作ること" ++
       "(𝓞_F[1/N] は ℤ[1/n!]-代数なので Spec 𝓞_F[1/N] ⟶ Spec ℤ[1/n!] があり、" ++
       "X への射と合わせて引き戻しの普遍性で作れる)") 9 ]

end ABC3.Found.GenEll
