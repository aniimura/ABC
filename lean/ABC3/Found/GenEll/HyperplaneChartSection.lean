/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.KerHyperplaneChart
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★超平面イデアルを**スキームの切断で**読む —— 段 C2c-1d（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★これは何か —— 代数の側とスキームの側を繋ぐ

`§9-863`（段 C2c-1 の (a)(b)(c)）が取ったのは**代数の側**である:

    `ker (Away.map (hyperplaneHom) (x_i)) = (x₀/x_i)`   （`A⁰_{x_i}` の中）

★しかし消費側（`pullbackIdealOf`）が読むのは**スキームの側**

    `hyperplaneIdeal.ideal U = ker ((hyperplaneι).app U)`   （`§9-857`）

であり、`U = D₊(x_i)` の切断環 `Γ(Proj 𝒜, D₊(x_i))` は `A⁰_{x_i}` と
**同型ではあるが同じ型ではない**。★★本ファイルはその橋を架ける。

## ★★★機構 —— mathlib の可換な四角を `app` に落とす

★mathlib の `Proj.awayι_comp_map`:

    `awayι ℬ (g s) ≫ Proj.map g = Spec.map (Away.map g s) ≫ awayι 𝒜 s`

★★これに `.app (D₊(x_i))` を当てる。左の合成の第 2 因子（`awayι ℬ` の `app`）は
**同型**なので核に効かず、右は `chartA.app` と `(Spec.map φ).app` に分かれる。

★★★配管の要点は 3 つ:

| 障害 | 直し方 |
|---|---|
| `f = g` から `f.app U = g.app U` が**直接は書けない**（終域の型が `f` に依る） | ★核 `RingHom.ker (f.app U).hom` は `Ideal Γ(Y,U)` で**`f` に依らない**ので、そこで `subst` する（`ker_app_congr`） |
| `chartA ⁻¹ᵁ (chartA ''ᵁ ⊤) = ⊤` が型に現れる | ★`V` を変数にして `subst`（`ker_specMap_app_eqToHom`） |
| `CommRingCat.hom_comp` の `rw` が通らない | ★`set_option backward.isDefEq.respectTransparency false`（mathlib 自身が `IdealSheaf` で使っている） |
-/

namespace ABC3.Found.GenEll

open MvPolynomial AlgebraicGeometry CategoryTheory HomogeneousLocalization

attribute [local instance] MvPolynomial.gradedAlgebra

variable (N : ℕ) (R : Type) [CommRing R] (i : Fin (N+1))

/-! ## ★チャートと超平面の射 -/

/-- ★**`ℙᴺ` の標準チャート** `Spec A⁰_{x_i} ⟶ Proj 𝒜`。 -/
noncomputable abbrev chartA : Spec (.of <| HomogeneousLocalization.Away
      (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R) (MvPolynomial.X i))
      ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R) :=
  Proj.awayι _ (MvPolynomial.X i) (MvPolynomial.isHomogeneous_X R i) one_pos

/-- ★**超平面の側のチャート** `Spec B⁰_{g(x_i)} ⟶ Proj ℬ`。 -/
noncomputable abbrev chartB : Spec (.of <| HomogeneousLocalization.Away
      (MvPolynomial.homogeneousSubmodule (Fin N) R) (hyperplaneHom N R (MvPolynomial.X i)))
      ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin N) R) :=
  Proj.awayι _ (hyperplaneHom N R (MvPolynomial.X i))
    ((hyperplaneHom N R).map_mem (MvPolynomial.isHomogeneous_X R i)) one_pos

/-- ★**チャートの上での超平面の引き戻し**（環の射）。 -/
noncomputable abbrev awayPhi : CommRingCat.of (HomogeneousLocalization.Away
      (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R) (MvPolynomial.X i))
    ⟶ CommRingCat.of (HomogeneousLocalization.Away
      (MvPolynomial.homogeneousSubmodule (Fin N) R) (hyperplaneHom N R (MvPolynomial.X i))) :=
  CommRingCat.ofHom (Away.map (hyperplaneHom N R) (MvPolynomial.X i))

/-! ## ★★配管の 2 本 -/

/-- ★**単射を後置しても核は変わらない**。 -/
theorem ker_comp_of_injective {A B C : Type} [Ring A] [Ring B] [Ring C]
    (a : A →+* B) (b : B →+* C) (hb : Function.Injective b) :
    RingHom.ker (b.comp a) = RingHom.ker a := by
  rw [← RingHom.comap_ker, (RingHom.injective_iff_ker_eq_bot b).1 hb]
  rfl

/-- ★★**等しい射の `app` の核は等しい**。

★`f.app U : Γ(Y,U) ⟶ Γ(X, f⁻¹ᵁU)` は**終域が `f` に依る**ので `f = g` から
`f.app U = g.app U` は直接書けない。★★しかし `RingHom.ker (f.app U).hom` の型は
`Ideal Γ(Y,U)` で**`f` に依らない**ので、そこでなら `subst` できる。 -/
theorem ker_app_congr {X Y : Scheme} {f g : X ⟶ Y} (h : f = g) (U : Y.Opens) :
    RingHom.ker (Scheme.Hom.app f U).hom = RingHom.ker (Scheme.Hom.app g U).hom := by
  subst h; rfl

set_option backward.isDefEq.respectTransparency false in
/-- ★★★**`Spec.map φ` の切断の核は `ker φ` の引き戻しである**。

★`V = ⊤` を**変数のまま**受け取って `subst` するのが要点である
——`Γ(Spec B, V)` は型なので `rw` では動かせない。 -/
theorem ker_specMap_app_eqToHom {B B' : CommRingCat.{0}} (φ : B ⟶ B')
    (V : (Spec B).Opens) (h : (⊤ : (Spec B).Opens) = V) :
    RingHom.ker ((Spec.map φ).app V).hom
      = Ideal.comap (((Spec B).presheaf.map (eqToHom h).op ≫ (Scheme.ΓSpecIso B).hom).hom)
          (RingHom.ker φ.hom) := by
  subst h
  have hinj : Function.Injective (Scheme.ΓSpecIso B').hom.hom :=
    (ConcreteCategory.bijective_of_isIso (Scheme.ΓSpecIso B').hom).1
  have hnat := Scheme.ΓSpecIso_naturality φ
  calc RingHom.ker ((Spec.map φ).app ⊤).hom
      = RingHom.ker ((Scheme.ΓSpecIso B').hom.hom.comp ((Spec.map φ).app ⊤).hom) :=
        (ker_comp_of_injective _ _ hinj).symm
    _ = RingHom.ker (((Spec.map φ).app ⊤ ≫ (Scheme.ΓSpecIso B').hom).hom) := by
        rw [CommRingCat.hom_comp]
    _ = RingHom.ker (((Scheme.ΓSpecIso B).hom ≫ φ).hom) := by rw [hnat]
    _ = RingHom.ker (φ.hom.comp (Scheme.ΓSpecIso B).hom.hom) := by rw [CommRingCat.hom_comp]
    _ = Ideal.comap (Scheme.ΓSpecIso B).hom.hom (RingHom.ker φ.hom) :=
        (RingHom.comap_ker _ _).symm
    _ = Ideal.comap ((((Spec B).presheaf.map (eqToHom rfl).op) ≫
          (Scheme.ΓSpecIso B).hom).hom) (RingHom.ker φ.hom) := by
        simp

/-! ## ★★★★★★★★★★段 C2c-1d の本体 -/

set_option maxHeartbeats 2000000 in
set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★**チャートの切断環で読んだ超平面イデアル** —— 段 C2c-1d。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

    `ker ((hyperplaneι).app (D₊(x_i))) = (x₀/x_i) の引き戻し`

★`§9-863` の代数側（`A⁰_{x_i}` の中の `ker (Away.map …)`）を、
`Γ(Proj 𝒜, D₊(x_i))` の言葉に翻訳したものである。
★★同型 `Γ(Proj 𝒜, D₊(x_i)) ≅ A⁰_{x_i}` は
`chartA.app ≫ presheaf.map (eqToHom) ≫ ΓSpecIso` である。 -/
theorem ker_app_hyperplaneι :
    RingHom.ker ((hyperplaneι N R).app ((chartA N R i) ''ᵁ ⊤)).hom
      = Ideal.comap (((chartA N R i).app ((chartA N R i) ''ᵁ ⊤)
          ≫ (Spec (CommRingCat.of (HomogeneousLocalization.Away
              (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R)
              (MvPolynomial.X i)))).presheaf.map
                (eqToHom (Scheme.Hom.preimage_image_eq (chartA N R i) ⊤).symm).op
          ≫ (Scheme.ΓSpecIso _).hom).hom) (Ideal.span {projCoord N R i 0}) := by
  have hV : hyperplaneι N R ⁻¹ᵁ ((chartA N R i) ''ᵁ ⊤) = (chartB N R i) ''ᵁ ⊤ := by
    rw [Scheme.Hom.image_top_eq_opensRange, Proj.opensRange_awayι, hyperplaneι,
      Proj.map_preimage_basicOpen, Scheme.Hom.image_top_eq_opensRange, Proj.opensRange_awayι]
  haveI : IsIso (Scheme.Hom.app (chartB N R i)
      (hyperplaneι N R ⁻¹ᵁ ((chartA N R i) ''ᵁ ⊤))) := by
    rw [hV]; infer_instance
  have hinj : Function.Injective (Scheme.Hom.app (chartB N R i)
      (hyperplaneι N R ⁻¹ᵁ ((chartA N R i) ''ᵁ ⊤))).hom :=
    (ConcreteCategory.bijective_of_isIso (Scheme.Hom.app (chartB N R i)
      (hyperplaneι N R ⁻¹ᵁ ((chartA N R i) ''ᵁ ⊤)))).1
  have hsq : chartB N R i ≫ hyperplaneι N R = Spec.map (awayPhi N R i) ≫ chartA N R i :=
    Proj.awayι_comp_map (hyperplaneHom N R) (hyperplane_irrelevant_le N R) one_pos
      (MvPolynomial.X i) (MvPolynomial.isHomogeneous_X R i)
  have hker : RingHom.ker ((hyperplaneι N R).app ((chartA N R i) ''ᵁ ⊤)).hom
      = Ideal.comap ((chartA N R i).app ((chartA N R i) ''ᵁ ⊤)).hom
          (RingHom.ker ((Spec.map (awayPhi N R i)).app
              ((chartA N R i) ⁻¹ᵁ ((chartA N R i) ''ᵁ ⊤))).hom) := by
    calc RingHom.ker ((hyperplaneι N R).app ((chartA N R i) ''ᵁ ⊤)).hom
        = RingHom.ker ((Scheme.Hom.app (chartB N R i)
            (hyperplaneι N R ⁻¹ᵁ ((chartA N R i) ''ᵁ ⊤))).hom.comp
              ((hyperplaneι N R).app ((chartA N R i) ''ᵁ ⊤)).hom) :=
          (ker_comp_of_injective _ _ hinj).symm
      _ = RingHom.ker ((chartB N R i ≫ hyperplaneι N R).app ((chartA N R i) ''ᵁ ⊤)).hom := by
          rw [Scheme.Hom.comp_app, CommRingCat.hom_comp]
      _ = RingHom.ker ((Spec.map (awayPhi N R i) ≫ chartA N R i).app
            ((chartA N R i) ''ᵁ ⊤)).hom := ker_app_congr hsq _
      _ = RingHom.ker ((Scheme.Hom.app (Spec.map (awayPhi N R i))
            ((chartA N R i) ⁻¹ᵁ ((chartA N R i) ''ᵁ ⊤))).hom.comp
              ((chartA N R i).app ((chartA N R i) ''ᵁ ⊤)).hom) := by
          rw [Scheme.Hom.comp_app, CommRingCat.hom_comp]
      _ = Ideal.comap ((chartA N R i).app ((chartA N R i) ''ᵁ ⊤)).hom
            (RingHom.ker ((Spec.map (awayPhi N R i)).app
                ((chartA N R i) ⁻¹ᵁ ((chartA N R i) ''ᵁ ⊤))).hom) :=
          (RingHom.comap_ker _ _).symm
  rw [hker, ker_specMap_app_eqToHom (awayPhi N R i) _
    (Scheme.Hom.preimage_image_eq (chartA N R i) ⊤).symm, Ideal.comap_comap,
    ← CommRingCat.hom_comp, CommRingCat.hom_ofHom, ker_awayMap_hyperplane]

/-! ## ★出典の紐付け(`.src`) -/

def ker_app_congr.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(等しい射の app の核は等しい)",
    sectionId := "genell-prop-1-4" }

def ker_specMap_app_eqToHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(Spec.map φ の切断の核は ker φ の引き戻し)",
    sectionId := "genell-prop-1-4" }

def ker_app_hyperplaneι.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(段 C2c-1d——チャートの切断環で読んだ超平面イデアル)",
    sectionId := "genell-prop-1-4" }

def ker_app_hyperplaneι.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "ker_awayMap_hyperplane(代数の側、段 C2c-1、§9-863)"
      (.inProject "ABC3" "ABC3.Found.GenEll.ker_awayMap_hyperplane") 2,
    .citation "[mathlib]" "Proj.awayι_comp_map(チャートの上で Proj.map は Away.map)"
      (.inMathlib "AlgebraicGeometry.Proj.awayι_comp_map") 2,
    .citation "[mathlib]" "Scheme.ΓSpecIso_naturality"
      (.inMathlib "AlgebraicGeometry.Scheme.ΓSpecIso_naturality") 2,
    .implicitStep
      ("★配管 1: f.app U : Γ(Y,U) ⟶ Γ(X, f⁻¹ᵁU) は**終域が f に依る**ので " ++
       "f = g から f.app U = g.app U は直接書けない。" ++
       "★★しかし RingHom.ker (f.app U).hom の型は Ideal Γ(Y,U) で**f に依らない**ので、" ++
       "そこでなら subst できる(ker_app_congr)") 3,
    .implicitStep
      ("★配管 2: chartA ⁻¹ᵁ (chartA ''ᵁ ⊤) = ⊤ は**型に現れる**ので rw では動かせない。" ++
       "V を変数のまま受け取って subst する(ker_specMap_app_eqToHom)") 3,
    .implicitStep
      ("★配管 3: CommRingCat.hom_comp の rw が「型が instances 透明度で正しくない」で" ++
       "通らない。set_option backward.isDefEq.respectTransparency false を置く" ++
       "(mathlib 自身が IdealSheaf/Basic.lean で同じことをしている)") 3,
    .implicitStep
      ("★★残るのは pullbackIdealOf(点に沿った引き戻し)への接続である: " ++
       "IdealSheafData.ideal_comap_of_isOpenImmersion で D.comap chartA の ⊤ での値が" ++
       "本補題の左辺の comap になり、pullbackIdealOf_specMap で点まで延びる") 4 ]

end ABC3.Found.GenEll
