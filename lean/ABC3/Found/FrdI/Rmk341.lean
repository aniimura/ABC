import ABC3.Found.FrdI.Thm34

/-!
# [FrdI] Remark 3.4.1 —— 次数を保てば底同型も保つ

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.69。

原文 (FrdI p.69):
> Remark 3.4.1. With regard to assumption (b) of Theorem 3.4, (iii), (iv), (v), we

原文 (FrdI p.69):
> observe the following: Suppose, in the situation of Theorem 3.4, that C1, C2 are of

★この `Remark` は `Theorem 3.4, (iii), (iv), (v)` の**仮定 (b) を外すための道具**である。
仮定 (b)「group-like 型なら `Ψ` とその擬逆が底同型を保つ」を、
**「次数を保つ」というより検査しやすい条件に置き換える**。

## ★★原文の証明の 3 段

原文 (FrdI p.69):
> 3.4, (i), we may assume, without loss of generality, that C1, C2 are of isotropic

| 段 | 内容 | 我々の実装 |
|---|---|---|
| 1 | `Theorem 3.4, (i)` で **isotropic 型に帰着** | `isotropificationCommute` ＋ `isotropification_baseIso_iff` |
| 2 | 次数を保つ ⟹ linear を保つ ⟹ **Frobenius 型を保つ** | `isFrobeniusType_map` |
| 3 | group-like ＋ isotropic では **底同型 = Frobenius 型**(同型を除く) | `isBaseIsomorphism_map_of_isotropic` |

原文 (FrdI p.69):
> preserve linear morphisms, hence morphisms of Frobenius type [cf. Proposition 1.7,

★★段 2 の中身は **`Proposition 1.7, (iii)` の圏論的特徴づけ**
(`Frobenius 型 ⟺ linear に minimal-coadjoint`)である。
★圏同値は分解を反射する(`exists_factor_of_map_factor`)ので、
**linear を保てばそのまま minimal-coadjoint 性が移る**。

## ★★測ったこと —— 原文の仮定「およびある擬逆」は**自動**である

原文は「`Ψ` **とその擬逆のどれか**が次数を保つなら」と 2 つ仮定するが、
★**`Ψ` が保てば擬逆も保つ**(`degFr_map_inv`)。同型は linear だからである。
★したがって仮定は 1 つでよい。**原文が弱い主張を述べているのではなく、
読者に確認の手間を掛けさせない書き方**と見るのが妥当である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v₁ u₁ v₂ u₂

variable {C : Type u₁} [Category.{v₁} C] {D : Type u₂} [Category.{v₂} D]
  (Ψ : C ⥤ D) [Ψ.IsEquivalence]

section Rmk341

universe v₃ u₃ w₃ v₄ u₄ w₄

variable {E₁ : Type u₃} [Category.{v₃} E₁] {Φ₁ : MonoidOn.{v₃, u₃, w₃} E₁}
  (P₁ : PreFrobenioid C Φ₁)
  {E₂ : Type u₄} [Category.{v₄} E₂] {Φ₂ : MonoidOn.{v₄, u₄, w₄} E₂}
  (P₂ : PreFrobenioid D Φ₂)

/-- ★**次数を保つ**という条件。原文の「preserve Frobenius degrees」。 -/
def PreservesDegFr : Prop := ∀ {A B : C} (φ : A ⟶ B), P₂.degFr (Ψ.map φ) = P₁.degFr φ

/-- ★★**次数を保てば linear 射を(両向きに)保つ**。★`IsLinear` は `degFr = 1` だから
**書き換え 1 回**である。 -/
theorem isLinear_map_iff (hdeg : PreservesDegFr Ψ P₁ P₂) {A B : C} (φ : A ⟶ B) :
    IsLinear P₂ (Ψ.map φ) ↔ IsLinear P₁ φ := by
  show P₂.degFr (Ψ.map φ) = 1 ↔ P₁.degFr φ = 1
  rw [hdeg φ]

/-- ★★★**次数を保てば Frobenius 型射を保つ**。

原文 (FrdI p.69):
> preserve linear morphisms, hence morphisms of Frobenius type [cf. Proposition 1.7,

★★鍵は `Proposition 1.7, (iii)` の**圏論的特徴づけ**
`IsFrobeniusType φ ↔ IsMinimalCoadjoint (linProp) φ` である。
★右辺は「linear 射」という語彙だけで書けているので、
**linear の保存(`isLinear_map_iff`)と分解の反射(`exists_factor_of_map_factor`)**
だけで移る。★**isotropy も group-like 性も要らない。** -/
theorem isFrobeniusType_map (F₁ : FrobenioidCore P₁) (F₂ : FrobenioidCore P₂)
    (hdeg : PreservesDegFr Ψ P₁ P₂) {A B : C} {φ : A ⟶ B}
    (h : IsFrobeniusType P₁ φ) : IsFrobeniusType P₂ (Ψ.map φ) := by
  rw [prop_1_7_iii_frobType_minimalCoadjoint P₂ F₂]
  rw [prop_1_7_iii_frobType_minimalCoadjoint P₁ F₁] at h
  intro X α β hfac hβ
  obtain ⟨X₀, α₀, β₀, e, hfac₀, hα, hβe⟩ := exists_factor_of_map_factor Ψ φ α β hfac
  -- ★`β = e.inv ≫ Ψ.map β₀` であり、`e.inv` は同型ゆえ linear。よって `Ψ.map β₀` も linear。
  have hmapβ₀ : IsLinear P₂ (Ψ.map β₀) := by
    have hcomp : P₂.degFr β = P₂.degFr (Ψ.map β₀) * P₂.degFr e.inv := by
      rw [hβe, P₂.degFr_comp]
    have hone : P₂.degFr e.inv = 1 := isLinear_of_isIso P₂ e.inv
    rw [hone, mul_one] at hcomp
    show P₂.degFr (Ψ.map β₀) = 1
    rw [← hcomp]
    exact hβ
  haveI : IsIso β₀ := h X₀ α₀ β₀ hfac₀ ((isLinear_map_iff Ψ P₁ P₂ hdeg β₀).mp hmapβ₀)
  rw [hβe]
  infer_instance

/-- ★★★**域が isotropic で、`𝒞₁` が group-like 型なら、`Ψ` は底同型を保つ**。

原文 (FrdI p.69):
> and base-isomorphisms [i.e., morphisms of Frobenius type, since C1, C2 are of

★★**筋**: `Proposition 1.7, (ii)` で底同型を `(Frobenius 型) ≫ (pre-step)` に分ける。
★group-like ＋ isotropic では **pre-step はすべて同型**
(`isIso_of_preStep_of_isGroupLikeObj`)なので、残るのは Frobenius 型の因子だけ。
★**isotropy は `isotropicClosed` で域から中間対象へ伝わる**——
原文が言う「isotropic 型に帰着してよい」の、対象 1 つぶんの中身である。 -/
theorem isBaseIsomorphism_map_of_isotropic (F₁ : FrobenioidCore P₁) (F₂ : FrobenioidCore P₂)
    (hg₁ : IsOfGroupLikeType P₁) (hdeg : PreservesDegFr Ψ P₁ P₂)
    {X Y : C} (hX : IsIsotropic P₁ X) {g : X ⟶ Y}
    (hg : IsBaseIsomorphism P₁ g) : IsBaseIsomorphism P₂ (Ψ.map g) := by
  obtain ⟨Z, γ, β, rfl, hγ, hβ⟩ := (prop_1_7_ii_baseIso_factor P₁ F₁ g).mp hg
  have hZ : IsIsotropic P₁ Z := F₁.isotropicClosed γ hX
  haveI : IsIso β :=
    isIso_of_preStep_of_isGroupLikeObj P₁ F₁ (fun _ f => F₁.isotropicClosed f hZ)
      (hg₁ Z) β hβ
  rw [Ψ.map_comp]
  exact isBaseIsomorphism_comp P₂ (isFrobeniusType_map Ψ P₁ P₂ F₁ F₂ hdeg hγ).2
    (isBaseIsomorphism_of_isIso P₂ (Ψ.map β))

/-- ★★★★★**[FrdI] Remark 3.4.1** —— `𝒞₁`・`𝒞₂` が group-like かつ
quasi-isotropic 型で、`Ψ` が Frobenius 次数を保つなら、**`Ψ` は底同型を保つ**。

原文 (FrdI p.69):
> 3.4, (i), we may assume, without loss of generality, that C1, C2 are of isotropic

★★**「isotropic 型に帰着してよい」の中身は 3 手**である:
1. `isotropification_baseIso_iff` —— 底同型性は isotropification で**両向きに移る**
2. `isotropificationCommute`(`Theorem 3.4, (i)` の 1-可換性)——
   `Ψ` の像の isotropification と、isotropification の `Ψ` 像が**同型で結ばれる**
3. その同型は `Base` を同型に送る(同型は底同型)ので、**両端で消える**

★★**3 が要点である**——原文の「without loss of generality」は、
形式化すると**この 1-可換図式を `Base` で潰す操作**になる。 -/
theorem rmk_3_4_1 (F₁ : FrobenioidCore P₁) (F₂ : FrobenioidCore P₂)
    (h₁ : IsOfQuasiIsotropicType C P₁) (h₂ : IsOfQuasiIsotropicType D P₂)
    (hg₁ : IsOfGroupLikeType P₁) (hdeg : PreservesDegFr Ψ P₁ P₂)
    {A B : C} {φ : A ⟶ B} (hφ : IsBaseIsomorphism P₁ φ) :
    IsBaseIsomorphism P₂ (Ψ.map φ) := by
  -- ★手 1: `𝒞₁` 側で isotropification に持ち上げる
  have h1 : IsBaseIsomorphism P₁ (istrMap P₁ F₁ φ) :=
    (isotropification_baseIso_iff P₁ F₁ φ).mpr hφ
  have h2 : IsBaseIsomorphism P₂ (Ψ.map (istrMap P₁ F₁ φ)) :=
    isBaseIsomorphism_map_of_isotropic Ψ P₁ P₂ F₁ F₂ hg₁ hdeg
      (hullMap_spec P₁ F₁ A).2.2.1 h1
  -- ★手 2: 1-可換性の自然性で `𝒟` 側の isotropification に移す
  have hnat : Ψ.map (istrMap P₁ F₁ φ)
        ≫ ((isotropificationCommute Ψ P₁ P₂ F₁ F₂ h₁ h₂).hom.app B).hom
      = ((isotropificationCommute Ψ P₁ P₂ F₁ F₂ h₁ h₂).hom.app A).hom
        ≫ istrMap P₂ F₂ (Ψ.map φ) :=
    congrArg (fun m => InducedCategory.Hom.hom m)
      ((isotropificationCommute Ψ P₁ P₂ F₁ F₂ h₁ h₂).hom.naturality φ)
  -- ★手 3: 両端の同型は `Base` で潰れる
  haveI : IsIso ((isotropificationCommute Ψ P₁ P₂ F₁ F₂ h₁ h₂).hom.app A) := inferInstance
  haveI : IsIso ((isotropificationCommute Ψ P₁ P₂ F₁ F₂ h₁ h₂).hom.app B) := inferInstance
  haveI : IsIso ((isotropificationCommute Ψ P₁ P₂ F₁ F₂ h₁ h₂).hom.app A).hom := inferInstance
  haveI : IsIso ((isotropificationCommute Ψ P₁ P₂ F₁ F₂ h₁ h₂).hom.app B).hom := inferInstance
  haveI hbA : IsIso (P₂.Base ((isotropificationCommute Ψ P₁ P₂ F₁ F₂ h₁ h₂).hom.app A).hom) :=
    isBaseIsomorphism_of_isIso P₂ _
  haveI hbB : IsIso (P₂.Base ((isotropificationCommute Ψ P₁ P₂ F₁ F₂ h₁ h₂).hom.app B).hom) :=
    isBaseIsomorphism_of_isIso P₂ _
  haveI hb2 : IsIso (P₂.Base (Ψ.map (istrMap P₁ F₁ φ))) := h2
  have hbase := congrArg P₂.Base hnat
  erw [P₂.Base_comp, P₂.Base_comp] at hbase
  have hcomp : IsIso (P₂.Base (Ψ.map (istrMap P₁ F₁ φ))
      ≫ P₂.Base ((isotropificationCommute Ψ P₁ P₂ F₁ F₂ h₁ h₂).hom.app B).hom) :=
    IsIso.comp_isIso' hb2 hbB
  haveI hright : IsIso (P₂.Base ((isotropificationCommute Ψ P₁ P₂ F₁ F₂ h₁ h₂).hom.app A).hom
      ≫ P₂.Base (istrMap P₂ F₂ (Ψ.map φ))) := hbase ▸ hcomp
  have hfin : IsIso (P₂.Base (istrMap P₂ F₂ (Ψ.map φ))) :=
    IsIso.of_isIso_comp_left
      (P₂.Base ((isotropificationCommute Ψ P₁ P₂ F₁ F₂ h₁ h₂).hom.app A).hom)
      (P₂.Base (istrMap P₂ F₂ (Ψ.map φ)))
  exact (isotropification_baseIso_iff P₂ F₂ (Ψ.map φ)).mp hfin

/-- ★★**原文の仮定「およびある擬逆が次数を保つ」は自動である**。

★`Ψ` が次数を保てば、擬逆 `Ψ⁻¹` も保つ。★★**同型が linear であること**だけが効く
——`Ψ (Ψ⁻¹ ψ)` は `ψ` と同型を挟んで等しく、同型の次数は 1 だからである。 -/
theorem degFr_map_inv (hdeg : PreservesDegFr Ψ P₁ P₂) {X Y : D} (ψ : X ⟶ Y) :
    P₁.degFr (Ψ.inv.map ψ) = P₂.degFr ψ := by
  have hnat := (Ψ.asEquivalence.counitIso).hom.naturality ψ
  -- hnat : (Ψ.inv ⋙ Ψ).map ψ ≫ ε.hom.app Y = ε.hom.app X ≫ ψ
  have h := congrArg P₂.degFr hnat
  rw [P₂.degFr_comp, P₂.degFr_comp] at h
  have e1 : P₂.degFr ((Ψ.asEquivalence.counitIso).hom.app X) = 1 :=
    isLinear_of_isIso P₂ _
  have e2 : P₂.degFr ((Ψ.asEquivalence.counitIso).hom.app Y) = 1 :=
    isLinear_of_isIso P₂ _
  rw [e1, e2, one_mul, mul_one] at h
  rw [← hdeg (Ψ.inv.map ψ)]
  exact h

/-- ★★★**擬逆の側も同じ** —— 原文の「**they** also preserve base-isomorphisms」の
「they」は `Ψ` **と擬逆の両方**を指す。

★`Ψ.inv` も圏同値であり、`degFr_map_inv` により次数を保つ。
★★要るのは `𝒞₂` が group-like 型であること —— `Ψ` の側で `𝒞₁` が要ったのと**対称**である。
★原文が「`C1, C2` are of group-like type」と**両方**を仮定する理由がこれである。 -/
theorem rmk_3_4_1_inv (F₁ : FrobenioidCore P₁) (F₂ : FrobenioidCore P₂)
    (h₁ : IsOfQuasiIsotropicType C P₁) (h₂ : IsOfQuasiIsotropicType D P₂)
    (hg₂ : IsOfGroupLikeType P₂) (hdeg : PreservesDegFr Ψ P₁ P₂)
    {X Y : D} {ψ : X ⟶ Y} (hψ : IsBaseIsomorphism P₂ ψ) :
    IsBaseIsomorphism P₁ (Ψ.inv.map ψ) :=
  rmk_3_4_1 Ψ.inv P₂ P₁ F₂ F₁ h₂ h₁ hg₂
    (fun φ => degFr_map_inv Ψ P₁ P₂ hdeg φ) hψ

/-- ★**`Remark 3.4.1` 全体**(`Ψ` 側 `rmk_3_4_1` ＋ 擬逆側 `rmk_3_4_1_inv`
＋ 擬逆の次数保存が自動であること `degFr_map_inv`)。 -/
def rmk_3_4_1.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 69, item := "Remark 3.4.1",
    sectionId := "frdi-remark-3-4-1" }

end Rmk341

end ABC3.Found.FrdI
