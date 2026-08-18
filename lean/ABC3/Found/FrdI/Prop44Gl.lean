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

variable {P G} in
/-- ★`sliceBaseOf` を `Iso` の形で書き直す(実例引数を消す)。 -/
theorem sliceBaseOf_asIso {A B A' : C} (a : A' ⟶ A) (ha : IsIso (P.Base a)) (φ : A' ⟶ B) :
    sliceBaseOf (P := P) a ha φ = (@asIso _ _ _ _ (P.Base a) ha).inv ≫ P.Base φ := rfl

variable {P G} in
/-- ★★★★**[FrdI] Proposition 4.4, (iv) の base-identity の条**(代表元の形)。 -/
theorem birat_isBaseIdentity_mk {A : C} (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ A) :
    IsBaseIdentity (biratPre P G) (HomBirat.mk Z φ)
      ↔ BaseEquivalent P Z.unop.hom.hom φ := by
  have hid : biratBase (𝟙 (A : BiratCat P G)) = 𝟙 (P.toElem.obj A).base := by
    show biratBase (toHomBirat (P := P) (G := G) (𝟙 A)) = 𝟙 _
    rw [biratBase_toHomBirat, P.Base_id]
  show biratBase (HomBirat.mk Z φ) = biratBase (𝟙 (A : BiratCat P G)) ↔ _
  rw [biratBase_mk, hid, sliceBaseOf_asIso, Iso.inv_comp_eq, Category.comp_id]
  exact eq_comm

variable {P G} in
/-- ★★★★**[FrdI] Proposition 4.4, (iv) の base-identity の条**(存在の形)。

★原文の「arises from a pair (α; φ) … base-equivalent」そのもの。 -/
theorem birat_isBaseIdentity_iff {A : C} (f : (A : BiratCat P G) ⟶ (A : BiratCat P G)) :
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

end ABC3.Found.FrdI
