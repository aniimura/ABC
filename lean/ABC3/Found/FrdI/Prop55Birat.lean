/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop32Frob
import ABC3.Found.FrdI.Prop44
import ABC3.Found.FrdI.Prop55UntrIdx

/-!
# [FrdI] Proposition 5.5, (ii) の `birat` の側 —— 添字圏の出発点

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.105。

原文 (FrdI p.105):
> tween the respective sets of morphisms between the images of two given objects of C

## ★★測って分かった構造

`(𝒞^pf)^birat` と `(𝒞^birat)^pf` はどちらも**二重の余極限**である:

* 左 = `colim_{Z ∈ IdxBirat(𝒞^pf)(A)} Hom^pf(Z, B)`
* 右 = `colim_{W ∈ IdxPf(𝒞^birat)(A,B)} Hom^birat(W)`

★内側の添字圏が外側の対象に依るので**単純な `colim_{I×J}` にはならない**。

## ★★★本ファイルの一歩 —— 左の添字圏が簡単になること

★`𝒞^pf` では**すべての射が co-angular**(`pfRoot_isCoAngular`、`Proposition 1.4, (i)`)
なので、`IdxBirat(𝒞^pf)(A)` の定義にある「co-angular pre-step」は
**「pre-step」だけ**になる。
★さらに `Proposition 3.2, (ii)` の判定(`isPreStep_mk_iff`)により、
それは `𝒞` の pre-step で代表される。

★★**対象は両側とも `𝒞` の対象である**(`PfRootObj` は根つきだが `BiratCat _ _ = C`)ので、
残るのは射の集合の全単射だけである。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

/-- ★★**`𝒞^pf` では「co-angular pre-step」＝「pre-step」**
(すべての射が co-angular だから)。

★★これで `(𝒞^pf)^birat` の添字圏(`IdxBirat`)は
**`𝒞^pf` の pre-step のスライスそのもの**になる。 -/
theorem coaPreProp_pfRoot_iff (hfi : IsOfFrobeniusIsotropicType P) {X Y : PfRootObj P F}
    (f : X ⟶ Y) : coaPreProp (pfRootPre P F) f ↔ IsPreStep (pfRootPre P F) f :=
  ⟨fun h => h.2, fun h => ⟨pfRoot_isCoAngular hfi f, h⟩⟩

/-- ★★★**`𝒞^pf` の co-angular pre-step は `𝒞` の pre-step で代表される**。 -/
theorem coaPreProp_pfRoot_mk_iff (hfi : IsOfFrobeniusIsotropicType P) {X Y : PfRootObj P F}
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    coaPreProp (pfRootPre P F) (show HomRoot P F X Y from HomPf.mk Z φ) ↔ IsPreStep P φ :=
  (coaPreProp_pfRoot_iff hfi _).trans (isPreStep_mk_iff Z φ)

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Proposition 5.5, (ii)` の `birat` の側の添字圏。 -/
def coaPreProp_pfRoot_mk_iff.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — (𝒞^pf)^birat の添字圏は 𝒞^pf の pre-step のスライス",
    sectionId := "frdi-prop-5-5" }

variable (P) (G : Frobenioid P)

/-! ## ★★右の側 —— `𝒞^birat` の Frobenius 型を代表元で判定する

★★**`biratUp` を置くのが実務上の要点**である ——
`BiratCat P G` は定義上 `C` そのものだが、`show BiratCat P G from A` と書くと
`rw` と instance 合成が「`instances` 透明度で型が合わない」と言って落ちる
(`PfCat` のときと同じ事情、2026-08-21 に踏んだ)。★**別名の `def` を置くこと。** -/

/-- ★`𝒞` の対象を `𝒞^birat` の対象として見る(`biratDown` の逆向き)。 -/
def biratUp (A : C) : BiratCat P G := A

@[simp] theorem biratDown_biratUp (A : C) : biratDown P G (biratUp P G A) = A := rfl

@[simp] theorem biratUp_biratDown (A : BiratCat P G) : biratUp P G (biratDown P G A) = A := rfl

/-- ★★★**`𝒞^birat` の Frobenius 型は代表元で判定できる**。

★`𝒞^birat` ではすべての射が isometric(`birat_isIsometric`)なので、
判定は **co-angular ＋ base-isomorphism** に落ちる。 -/
theorem birat_isFrobeniusType_repr {X Y : BiratCat P G} (f : X ⟶ Y) {A' : C}
    (φ : A' ⟶ biratDown P G Y)
    (aa : biratUp P G A' ⟶ X) (a' : X ⟶ biratUp P G A')
    (h1 : aa ≫ a' = 𝟙 (biratUp P G A')) (h2 : a' ≫ aa = 𝟙 X)
    (heq : f = a' ≫ (toBiratCat P G).map φ) :
    IsFrobeniusType (biratPre P G) f ↔ (IsCoAngular P φ ∧ IsBaseIsomorphism P φ) := by
  haveI hiso : IsIso a' := ⟨aa, h2, h1⟩
  haveI hba : IsIso ((biratPre P G).Base a') := isIso_Base_of_isIso a'
  have hcomp : (biratPre P G).Base f
      = (biratPre P G).Base a' ≫ (biratPre P G).Base ((toBiratCat P G).map φ) := by
    rw [heq]
    exact (biratPre P G).Base_comp _ _
  constructor
  · rintro ⟨⟨hc, _⟩, hb⟩
    have hcφ : IsCoAngular P φ := (birat_isCoAngular_repr P G f φ aa a' h1 h2 heq).mp hc
    refine ⟨hcφ, ?_⟩
    haveI hc2 : IsIso ((biratPre P G).Base a' ≫ (biratPre P G).Base ((toBiratCat P G).map φ)) :=
      hcomp ▸ hb
    have hbb : IsIso ((biratPre P G).Base ((toBiratCat P G).map φ)) :=
      IsIso.of_isIso_comp_left ((biratPre P G).Base a')
        ((biratPre P G).Base ((toBiratCat P G).map φ))
    exact ((birat_isFrobeniusType_iff P G φ).mp
      ⟨⟨(birat_isCoAngular_iff P G φ).mpr hcφ, birat_isIsometric _⟩, hbb⟩).2
  · rintro ⟨hc, hb⟩
    have hft := (birat_isFrobeniusType_iff P G φ).mpr ⟨hc, hb⟩
    refine ⟨⟨(birat_isCoAngular_repr P G f φ aa a' h1 h2 heq).mpr hc, birat_isIsometric _⟩, ?_⟩
    show IsIso ((biratPre P G).Base f)
    rw [hcomp]
    exact @IsIso.comp_isIso _ _ _ _ _ _ _ hba hft.2

/-- ★★★locator —— `Proposition 5.5, (ii)` の `birat` の側の Frobenius 型判定。 -/
def birat_isFrobeniusType_repr.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — 𝒞^birat の Frobenius 型は代表元で判定できる",
    sectionId := "frdi-prop-5-5" }

/-! ## ★★★添字圏のあいだの関手 -/

variable (F) (F' : FrobenioidCore (biratPre P G))

/-- ★★**`𝒞^{bi-Fr}` から `(𝒞^birat)^{bi-Fr}` への関手**。

★`𝒞` の Frobenius 型射は `𝒞^birat` でも Frobenius 型である
(`birat_isFrobeniusType_iff`、co-angular と base-isomorphism がそのまま渡る)。
★次数が合うのは `biratDeg_toHomBirat`。 -/
noncomputable def biFrToBirat : BiFr P F ⥤ BiFr (biratPre P G) F' where
  obj Z := ⟨((toBiratCat P G).obj Z.obj.1, (toBiratCat P G).obj Z.obj.2)⟩
  map {Z W} u :=
    ⟨((toBiratCat P G).map u.hom.1, (toBiratCat P G).map u.hom.2),
      (birat_isFrobeniusType_iff P G _).mpr ⟨u.property.1.1.1, u.property.1.2⟩,
      (birat_isFrobeniusType_iff P G _).mpr ⟨u.property.2.1.1.1, u.property.2.1.2⟩,
      (biratDeg_toHomBirat (P := P) (G := G) u.hom.1).trans
        (u.property.2.2.trans (biratDeg_toHomBirat (P := P) (G := G) u.hom.2).symm)⟩
  map_id Z := by
    apply InducedWideCategory.Hom.ext
    exact Prod.ext ((toBiratCat P G).map_id _) ((toBiratCat P G).map_id _)
  map_comp u v := by
    apply InducedWideCategory.Hom.ext
    exact Prod.ext ((toBiratCat P G).map_comp _ _) ((toBiratCat P G).map_comp _ _)

/-- ★★★locator —— `Proposition 5.5, (ii)` の添字圏のあいだの関手。 -/
def biFrToBirat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — 𝒞^{bi-Fr} から (𝒞^birat)^{bi-Fr} への関手",
    sectionId := "frdi-prop-5-5" }

/-! ## ★★★★添字圏のあいだの**終尾関手** -/

/-- ★★**添字圏のあいだの関手** —— コスライスへ持ち上げただけ。 -/
noncomputable def idxToBirat (A B : C) :
    IdxPf P F A B ⥤ IdxPf (biratPre P G) F' (biratUp P G A) (biratUp P G B) :=
  Under.post (biFrToBirat P F G F')

/-- ★★★★**本質的全射**。

★★**pre-step 分解は要らなかった** —— `𝒞^birat` 側の Frobenius 型射
`g : A → A″` に対し、**同じ次数の Frobenius 型射を `𝒞` の中で選び直す**
(`F.frobDegSurj`)。★2 本は `𝒞^birat` の中で同じ次数の Frobenius 型射なので、
**`Definition 1.3, (ii)` の本質的一意性**(`F'.frobDegUniq`)で同型が 1 本出る。
★これが `Under` 圏の同型を与える。 -/
instance idxToBirat_essSurj (A B : C) : (idxToBirat P F G F' A B).EssSurj where
  mem_essImage d := by
    have hg₁ : IsFrobeniusType (biratPre P G) d.hom.hom.1 := d.hom.property.1
    have hg₂ : IsFrobeniusType (biratPre P G) d.hom.hom.2 := d.hom.property.2.1
    have hdg : (biratPre P G).degFr d.hom.hom.1 = (biratPre P G).degFr d.hom.hom.2 :=
      d.hom.property.2.2
    obtain ⟨A₁, α, hα, hαd⟩ := F.frobDegSurj A ((biratPre P G).degFr d.hom.hom.1)
    obtain ⟨B₁, β, hβ, hβd⟩ := F.frobDegSurj B ((biratPre P G).degFr d.hom.hom.1)
    have hαb : IsFrobeniusType (biratPre P G) ((toBiratCat P G).map α) :=
      (birat_isFrobeniusType_iff P G α).mpr ⟨hα.1.1, hα.2⟩
    have hβb : IsFrobeniusType (biratPre P G) ((toBiratCat P G).map β) :=
      (birat_isFrobeniusType_iff P G β).mpr ⟨hβ.1.1, hβ.2⟩
    have hαdb : (biratPre P G).degFr ((toBiratCat P G).map α)
        = (biratPre P G).degFr d.hom.hom.1 :=
      (biratDeg_toHomBirat (P := P) (G := G) α).trans hαd
    have hβdb : (biratPre P G).degFr ((toBiratCat P G).map β)
        = (biratPre P G).degFr d.hom.hom.2 :=
      ((biratDeg_toHomBirat (P := P) (G := G) β).trans hβd).trans hdg
    obtain ⟨θ₁, hθ₁iso, hθ₁⟩ := F'.frobDegUniq (biratUp P G A) (biratUp P G A₁)
      d.right.obj.1 ((toBiratCat P G).map α) d.hom.hom.1 hαb hg₁ hαdb
    obtain ⟨θ₂, hθ₂iso, hθ₂⟩ := F'.frobDegUniq (biratUp P G B) (biratUp P G B₁)
      d.right.obj.2 ((toBiratCat P G).map β) d.hom.hom.2 hβb hg₂ hβdb
    haveI := hθ₁iso
    haveI := hθ₂iso
    have hd1 : (biratPre P G).degFr θ₁ = 1 := degFr_of_isIso _ θ₁
    have hd2 : (biratPre P G).degFr θ₂ = 1 := degFr_of_isIso _ θ₂
    have hdi1 : (biratPre P G).degFr (inv θ₁) = 1 := degFr_of_isIso _ (inv θ₁)
    have hdi2 : (biratPre P G).degFr (inv θ₂) = 1 := degFr_of_isIso _ (inv θ₂)
    let W : BiFr P F := ⟨(A₁, B₁)⟩
    let u : biFrObj P F A B ⟶ W := ⟨(α, β), hα, hβ, hαd.trans hβd.symm⟩
    let W' : BiFr (biratPre P G) F' := ⟨(biratUp P G A₁, biratUp P G B₁)⟩
    let e : W' ≅ d.right :=
      { hom := ⟨(θ₁, θ₂), isFrobeniusType_of_isIso _ θ₁, isFrobeniusType_of_isIso _ θ₂,
          hd1.trans hd2.symm⟩
        inv := ⟨(inv θ₁, inv θ₂), isFrobeniusType_of_isIso _ _, isFrobeniusType_of_isIso _ _,
          hdi1.trans hdi2.symm⟩
        hom_inv_id := by
          apply InducedWideCategory.Hom.ext
          exact Prod.ext (IsIso.hom_inv_id θ₁) (IsIso.hom_inv_id θ₂)
        inv_hom_id := by
          apply InducedWideCategory.Hom.ext
          exact Prod.ext (IsIso.inv_hom_id θ₁) (IsIso.inv_hom_id θ₂) }
    refine ⟨Under.mk u, ⟨Under.isoMk e ?_⟩⟩
    apply InducedWideCategory.Hom.ext
    exact Prod.ext hθ₁ hθ₂

/-- ★★★★★**終尾(final)** —— 一般補題(`final_of_essSurj_of_thin`)を当てるだけ。

★細さ(`idx_hom_ext`)は `𝒞^birat` が totally epimorphic であること、
すなわち **pre-Frobenioid の定義そのもの**から来る。 -/
instance idxToBirat_final (A B : C) : (idxToBirat P F G F' A B).Final :=
  final_of_essSurj_of_thin _ (fun f g => idx_hom_ext f g)

/-- ★★★locator —— `Proposition 5.5, (ii)` の添字圏の終尾性。 -/
def idxToBirat_final.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — 添字圏のあいだの終尾関手(birat の側)",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
