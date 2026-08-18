/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop44Equiv

/-!
# [FrdI] `Proposition 4.4, (ii)` —— `𝒞^birat` は **group-like 型**

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.83。

原文 (FrdI p.83):
> (ii) The functor Cbirat →F0D of (i) determines a structure of Frobenioid

★原文は「**of group-like type**」と言う。★★その中身は
**`0_𝒟` の値がすべて 1 元単系である**ことに尽きる ——
`Div` が情報を持たないので、`Φ^birat` の値は群(実は自明群)になる。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w v2 u2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (G : Frobenioid P)

/-- ★★★**[FrdI] Proposition 4.4, (ii)** —— `𝒞^birat` は **group-like 型**。

★`biratPre` の単系は `0_𝒟`(値がすべて 1 元)なので、
すべての元が `0` に等しく、したがって加法的可逆である。 -/
theorem birat_isOfGroupLikeType : IsOfGroupLikeType (biratPre P G) := by
  intro A
  refine (isGroupLike_iff _).mpr (fun a => ?_)
  rw [Subsingleton.elim a 0]
  exact isAddUnit_zero

/-- ★locator —— `Proposition 4.4, (ii)` の「group-like 型」。 -/
def birat_isOfGroupLikeType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (ii) — group-like 型",
    sectionId := "frdi-prop-4-4" }


/-! ## ★辞書の pull-back の条

原文 (FrdI p.83):
> of a given Frobenius degree; isometry; pre-step; base-isomorphism) of

★★原文の辞書は「pull-back 射 ⇔ **co-angular linear** 射」と言う。
★`𝒞` では pull-back は **LB-invertible ＋ linear**＝
**co-angular ＋ 等長 ＋ linear** だが、
★★★`𝒞^birat` では**等長が自動**(`birat_isIsometric`)なので
その 1 条が落ちる——これが「等長を忘れる」という
birationalization の働きそのものである。 -/

/-- ★★★★**[FrdI] Proposition 4.4, (iv) の pull-back の条** —— 辞書。

★`⟹` は `birat_pullBackLB` と co-angular / linear の辞書から。
★`⟸` は `𝒞^birat` で等長が自動なので LB-invertible になり、
`Proposition 1.4, (ii)` を `𝒞^birat` で使う。

★★**逆向きに `𝒞^birat` の `FrobenioidCore` が要る**ので、
分類 (B)(birat-Frobenius-normalized 型を仮定)の道になる
(`birat_frobenioidCore_of_frobNormalized`)。 -/
theorem birat_isPullBack_iff (Fc : FrobenioidCore (biratPre P G)) {A B : C} (φ : A ⟶ B) :
    IsPullBack (biratPre P G) ((toBiratCat P G).map φ) ↔
      (IsCoAngular P φ ∧ IsLinear P φ) := by
  constructor
  · intro hpb
    obtain ⟨hlb, hlin⟩ := birat_pullBackLB P G _ hpb
    exact ⟨(birat_isCoAngular_iff P G φ).mp hlb.1,
      (birat_isLinear_iff (P := P) (G := G) φ).mp hlin⟩
  · rintro ⟨hc, hlin⟩
    refine prop_1_4_ii_mpr (biratPre P G) Fc _ ?_ ?_
    · exact ⟨(birat_isCoAngular_iff P G φ).mpr hc, birat_isIsometric _⟩
    · exact (birat_isLinear_iff (P := P) (G := G) φ).mpr hlin

/-- ★locator —— `Proposition 4.4, (iv)` の pull-back の条。 -/
def birat_isPullBack_iff.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (iv) — pull-back の条",
    sectionId := "frdi-prop-4-4" }


/-! ## ★辞書の base-identity 自己射の条

原文 (FrdI p.83):
> a pair (α : A → A; φ : A → A), where α is a co-angular pre-step in the in-

★★代表元で計算すれば
`biratBase (mk Z φ) = inv (Base α) ≫ Base φ` なので、
これが `𝟙` になることと `Base α = Base φ`(＝base-equivalent)は同じことである。

★★★**罠の記録**: 目標式の `inv` は `sliceBaseOf` が抱える `IsIso` の
**実例**を使うので、`haveI` で同じ命題を入れても
`rw` / `IsIso.inv_comp_eq` が照合しない。
★**`asIso` で包んで `Iso` の API に乗せ換える**と実例引数が消えて通る。 -/

/-- ★★**抽象化した補題** —— `inv f ≫ g = 𝟙 ⇔ f = g`。

★★★実例を**引数で固定**するために切り出している。
目標式の `inv` は `sliceBaseOf` が抱える実例を使うので、
`rw` や `IsIso.inv_comp_eq` は照合しないが、
**`@` で実例を渡して `exact` すれば defeq で通る**。 -/
theorem inv_comp_eq_id_iff {K : Type u} [Category.{v} K] {X Y : K}
    (f : X ⟶ Y) [IsIso f] (g : X ⟶ Y) : inv f ≫ g = 𝟙 Y ↔ f = g := by
  constructor
  · intro h
    calc f = f ≫ 𝟙 Y := (Category.comp_id f).symm
      _ = f ≫ (inv f ≫ g) := by rw [h]
      _ = (f ≫ inv f) ≫ g := (Category.assoc _ _ _).symm
      _ = g := by rw [IsIso.hom_inv_id, Category.id_comp]
  · rintro rfl
    exact IsIso.inv_hom_id f

variable {P G} in
/-- ★★★★**[FrdI] Proposition 4.4, (iv) の base-identity の条**(代表元の形)。 -/
theorem birat_isBaseIdentity_mk {A : C} (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ A) :
    IsBaseIdentity (biratPre P G) (HomBirat.mk Z φ)
      ↔ BaseEquivalent P Z.unop.hom.hom φ := by
  have hid : biratBase (𝟙 (show BiratCat P G from A))
      = 𝟙 (P.toElem.obj A).base := by
    show biratBase (toHomBirat (P := P) (G := G) (𝟙 A)) = 𝟙 _
    rw [biratBase_toHomBirat, P.Base_id]
  constructor
  · intro h
    have h1 : sliceBaseOf (P := P) Z.unop.hom.hom Z.unop.hom.property.2.2 φ
        = 𝟙 (P.toElem.obj A).base := by
      rw [← biratBase_mk]
      exact h.trans hid
    exact (@inv_comp_eq_id_iff _ _ _ _ (P.Base Z.unop.hom.hom)
      Z.unop.hom.property.2.2 (P.Base φ)).mp h1
  · intro h
    show biratBase (HomBirat.mk Z φ) = biratBase (𝟙 (show BiratCat P G from A))
    refine Eq.trans ?_ hid.symm
    rw [biratBase_mk]
    exact (@inv_comp_eq_id_iff _ _ _ _ (P.Base Z.unop.hom.hom)
      Z.unop.hom.property.2.2 (P.Base φ)).mpr h

variable {P G} in
/-- ★★★★**[FrdI] Proposition 4.4, (iv) の base-identity の条**(存在の形)。

★原文の「arises from a pair (α; φ) … base-equivalent」そのもの。 -/
theorem birat_isBaseIdentity_iff {A : C}
    (f : (show BiratCat P G from A) ⟶ (show BiratCat P G from A)) :
    IsBaseIdentity (biratPre P G) f ↔
      ∃ (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ A),
        f = HomBirat.mk Z φ ∧ BaseEquivalent P Z.unop.hom.hom φ := by
  constructor
  · intro h
    obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep f
    exact ⟨Z, φ, rfl, (birat_isBaseIdentity_mk Z φ).mp h⟩
  · rintro ⟨Z, φ, rfl, hbe⟩
    exact (birat_isBaseIdentity_mk Z φ).mpr hbe

/-- ★locator —— `Proposition 4.4, (iv)` の base-identity 自己射の条。 -/
def birat_isBaseIdentity_iff.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (iv) — base-identity 自己射の条",
    sectionId := "frdi-prop-4-4" }

/-- ★★★★**[FrdI] Proposition 4.4, (iv)** の locator。

★辞書の 10 項がすべて揃った:

| 項 | 宣言 |
|---|---|
| co-angular | `birat_isCoAngular_iff` |
| 同型 ⇔ co-angular pre-step | `birat_isIso_iff` |
| Frobenius 型 ⇔ co-angular 底同型 | `birat_isFrobeniusType_iff` |
| pull-back ⇔ co-angular linear | `birat_isPullBack_iff` |
| 与えられた Frobenius 次数 | `birat_degFr_iff` |
| isometry ⇔ 任意の射 | `birat_isIsometric` |
| pre-step | `birat_isPreStep_iff` |
| 底同型 | `birat_isBaseIsomorphism_iff` |
| base-identity 自己射 | `birat_isBaseIdentity_iff` |
| isotropic 対象 | `birat_isIsotropic_iff` | -/
def prop_4_4_iv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83, item := "Proposition 4.4, (iv)",
    sectionId := "frdi-prop-4-4" }

/-! ## ★birat-Frobenius-normalized なら Frobenius-normalized

原文 (FrdI p.86):
> normalized object of C is Frobenius-normalized — cf. Proposition 4.4, (ii), (iv).]

★★原文はこれを `Proposition 4.4, (ii), (iv)` に送る。
★実際に使うのは (ii) の**忠実性**と (iv) の**辞書**だけである。

★★これが要るのは、`𝒪^▷(A)^gp` が意味を持つには
`𝒪^▷(A)` の**可換性**が要り(`otri_mul_comm`)、
それが Frobenius-normalized から出るからである。 -/

variable {P G} in
/-- ★★★★**birationally Frobenius-normalized な対象は Frobenius-normalized**。

★手: 両辺を `𝒞^birat` へ送り、そこで仮定を使ってから
**忠実性**で `𝒞` へ戻す。★`α ^ n` と関手の交換は `mapEnd` の `map_pow`。 -/
theorem isFrobeniusNormalized_of_birat {A : C}
    (h : IsFrobeniusNormalized (biratPre P G) ((toBiratCat P G).obj A)) :
    IsFrobeniusNormalized P A := by
  haveI : (toBiratCat P G).Faithful := toBiratCat_faithful
  intro φ hφ α hα
  set T := toBiratCat P G with hT
  set e := CategoryTheory.Functor.mapEnd A T with he
  have hb : ∀ δ : End A, IsBaseIdentity P δ →
      IsBaseIdentity (biratPre P G) (e δ) := by
    intro δ hδ
    show biratBase (toHomBirat (P := P) (G := G) (δ : A ⟶ A))
      = biratBase (toHomBirat (P := P) (G := G) (𝟙 A))
    rw [biratBase_toHomBirat, biratBase_toHomBirat]
    exact hδ
  have hd : ∀ δ : End A, (biratPre P G).degFr (e δ) = P.degFr δ := by
    intro δ
    show biratDeg (toHomBirat (P := P) (G := G) (δ : A ⟶ A)) = P.degFr δ
    rw [biratDeg_toHomBirat]
  have hα2 : e α ∈ OTri (biratPre P G) (T.obj A) :=
    ⟨hb α hα.1, by
      show (biratPre P G).degFr (e α) = 1
      rw [hd]
      exact hα.2⟩
  have hres := h (e φ) (hb φ hφ) (e α) hα2
  rw [hd] at hres
  refine T.map_injective ?_
  have hp : (e α ^ ((P.degFr φ : ℕ+) : ℕ) : End (T.obj A))
      = e ((α ^ ((P.degFr φ : ℕ+) : ℕ) : End A)) := (map_pow e α _).symm
  rw [hp] at hres
  exact (T.map_comp (φ : A ⟶ A) _).trans (hres.trans (T.map_comp _ _).symm)

/-- ★★**その帰結** —— birat-Frobenius-normalized なら `𝒪^▷(A)` は可換。 -/
theorem otri_mul_comm_of_birat {A : C}
    (h : IsFrobeniusNormalized (biratPre P G) ((toBiratCat P G).obj A))
    (x y : OTri P A) : x * y = y * x :=
  otri_mul_comm P (isFrobeniusNormalized_of_birat h) x y

/-- ★locator —— `Definition 4.5, (i)` の括弧書き。 -/
def isFrobeniusNormalized_of_birat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 86,
    item := "Definition 4.5, (i) — birat-Frobenius-normalized なら Frobenius-normalized",
    sectionId := "frdi-def-4-5" }

/-! ## ★`𝒪^▷(A)` から `𝒪^×(A^birat)` へ

原文 (FrdI p.83):
> determines an injection of groups O (A)gp →O×(Abirat). We shall refer to the

★★★**測定(2026-08-18)**: `Definition 1.3, (iii)(b)`(`coAngularOfPreStep`)は
「**co-angular pre-step が 1 本でもあれば、その 2 対象間の射はすべて co-angular**」
と言う。★`𝟙 A` がそれなので、**`A` の自己射はすべて co-angular** である。
★★これで `𝒪^▷(A)` の元は co-angular pre-step になり、
`𝒞^birat` で**同型**になる(`birat_isIso_iff`)。
★自己射の co-angular 性は既に `Prop110.lean` の `endo_isCoAngular` にある。 -/

variable {P G} in
/-- ★★`𝒪^▷(A)` の元は `𝒞^birat` で同型になる。 -/
theorem otri_isIso_birat {A : C} (α : OTri P A) :
    IsIso ((toBiratCat P G).map ((α : End A) : A ⟶ A)) :=
  birat_isIso_of_coaPre _ (endo_isCoAngular P G.core _)
    ⟨α.2.2, by
      show IsIso (P.Base ((α : End A) : A ⟶ A))
      rw [show P.Base ((α : End A) : A ⟶ A) = P.Base (𝟙 A) from α.2.1, P.Base_id]
      infer_instance⟩

variable {P G} in
/-- ★★★**`𝒪^▷(A) → End(A^birat)` は単射**(忠実性から)。

★これが原文の「an injection of groups」の**単系の層**である。 -/
theorem otri_toBirat_injective {A : C} :
    Function.Injective (fun α : OTri P A =>
      (toBiratCat P G).map ((α : End A) : A ⟶ A)) := by
  haveI : (toBiratCat P G).Faithful := toBiratCat_faithful
  intro x y h
  exact Subtype.ext ((toBiratCat P G).map_injective h)

/-- ★locator —— `Proposition 4.4, (ii)` の単射の条(単系の層)。 -/
def otri_toBirat_injective.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (ii) — 𝒪^▷(A) の像への単射性",
    sectionId := "frdi-prop-4-4" }

end ABC3.Found.FrdI
