/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Rmk451

/-!
# 「型」の圏をまたぐ移送

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]。

★★`Definition 1.2, (iv)` の性質のうち、**`End A`・`𝒪^▷(A)`・`Base`・`degFr` の
4 つだけで書かれているもの**は、その 4 つが対応すれば
**圏も単系も違ってよい**——そのまま移送できる。

★在庫の `isFrobeniusCompact_transport`(`Rmk451.lean`)が同じ形の先例である
(そちらは `End`・`𝒪^×`・自己同型の 3 つ)。本ファイルは
**Frobenius-normalized 版**を置く。

## ★なぜ要るか

`𝒞^un-tr` / `𝒞^rlf` が **standard 型**であること(`Proposition 5.5, (iii)` の最後の文)を
言うには `IsOfStandardType` の 6 条を移す必要があり、そのうち
`frobNormalized` は model Frobenioid 側では無条件に成り立つ
(`ModelData.model_frobNormalizedType`)。★`Theorem 5.2, (iv)` の圏同値で
それを `𝒞^un-tr` へ引き戻すのがこの補題の用途である。

## ★★`End` の積の向き

mathlib の `End X` は **`x * y = y ≫ x`** である。したがって
単系の同型 `e : End A₁ ≃* End A₂` は
**`e (f ≫ g) = e f ≫ e g`**(向きはそのまま)を満たす —— 両側で反転するから。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2 u3 v3 w3 u4 v4

/-! ## ★0. 充満忠実な関手が与える `End` の単系同型 -/

section EndEquiv

variable {C : Type u2} [Category.{v2} C] {D : Type u4} [Category.{v4} D]

/-- ★★**充満忠実な関手は `End` の単系の同型を与える**。

★`End` の積は `x * y = y ≫ x` なので、`map_mul'` はそのまま `Functor.map_comp` である。 -/
noncomputable def endMulEquiv (Ψ : C ⥤ D) [Ψ.Full] [Ψ.Faithful] (A : C) :
    End A ≃* End (Ψ.obj A) where
  toFun f := Ψ.map f
  invFun g := (Functor.FullyFaithful.ofFullyFaithful Ψ).preimage g
  left_inv f := (Functor.FullyFaithful.ofFullyFaithful Ψ).map_injective (by simp)
  right_inv g := by simp [Functor.FullyFaithful.ofFullyFaithful]
  map_mul' x y := by
    show Ψ.map ((y : A ⟶ A) ≫ (x : A ⟶ A)) = _
    show _ = (Ψ.map (y : A ⟶ A)) ≫ (Ψ.map (x : A ⟶ A))
    exact Ψ.map_comp _ _

end EndEquiv

/-! ## ★1. `Frobenius-normalized` の移送 -/

/-- ★★★**`Frobenius-normalized` は「`End`・`𝒪^▷`・`Base`・`degFr`」の対応で移送できる**。

★`isFrobeniusCompact_transport`(`Rmk451.lean`)と同じ形である。
★★仮定は 3 本だけ:
* `hbi` —— base-identity 性が対応する
* `hdeg` —— Frobenius 次数が対応する
* `hot` —— `𝒪^▷` が対応する -/
theorem isFrobeniusNormalized_transport
    {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
    {Φ₁ : MonoidOn.{v, u, w} D₁} (P₁ : PreFrobenioid C₁ Φ₁)
    {D₂ : Type u3} [Category.{v3} D₂] {C₂ : Type u4} [Category.{v4} C₂]
    {Φ₂ : MonoidOn.{v3, u3, w3} D₂} (P₂ : PreFrobenioid C₂ Φ₂)
    {A₁ : C₁} {A₂ : C₂}
    (e : End A₁ ≃* End A₂)
    (hbi : ∀ φ : End A₁, IsBaseIdentity P₁ φ ↔ IsBaseIdentity P₂ (e φ))
    (hdeg : ∀ φ : End A₁, P₂.degFr ((e φ : End A₂) : A₂ ⟶ A₂)
      = P₁.degFr ((φ : End A₁) : A₁ ⟶ A₁))
    (hot : ∀ α : End A₁, α ∈ OTri P₁ A₁ ↔ e α ∈ OTri P₂ A₂)
    (h : IsFrobeniusNormalized P₁ A₁) : IsFrobeniusNormalized P₂ A₂ := by
  intro φ hφ α hα
  have hφ₁ : IsBaseIdentity P₁ (e.symm φ) := by
    rw [hbi (e.symm φ), MulEquiv.apply_symm_apply]; exact hφ
  have hα₁ : e.symm α ∈ OTri P₁ A₁ := by
    rw [hot (e.symm α), MulEquiv.apply_symm_apply]; exact hα
  have key := h (e.symm φ) hφ₁ (e.symm α) hα₁
  have hd : P₁.degFr ((e.symm φ : End A₁) : A₁ ⟶ A₁) = P₂.degFr ((φ : End A₂) : A₂ ⟶ A₂) := by
    have := hdeg (e.symm φ)
    rw [MulEquiv.apply_symm_apply] at this
    exact this.symm
  have hmul : ∀ x y : End A₁,
      e ((x : A₁ ⟶ A₁) ≫ (y : A₁ ⟶ A₁)) = ((e x : A₂ ⟶ A₂) ≫ (e y : A₂ ⟶ A₂)) := by
    intro x y
    show e (y * x) = _
    rw [map_mul]
    rfl
  have h2 := congrArg e key
  rw [hmul, hmul, map_pow, MulEquiv.apply_symm_apply, MulEquiv.apply_symm_apply, hd] at h2
  exact h2

/-! ## ★2. 関手の形に束ねる

★★実際に使うときは「圏同値」が与えられるので、
**充満忠実で `Base` と `degFr` を保つ関手**という形にしておくと当てやすい。
★`𝒪^▷` / `𝒪^×` の対応はその 2 つから**自動で出る**(定義がその 2 つで書かれているから)。 -/

section MapForm

variable {C₁ : Type u2} [Category.{v2} C₁] {C₂ : Type u4} [Category.{v4} C₂]
  {D₁ : Type u} [Category.{v} D₁] {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u3} [Category.{v3} D₂] {Φ₂ : MonoidOn.{v3, u3, w3} D₂}
  {P₂ : PreFrobenioid C₂ Φ₂}

/-- ★★★★**充満忠実で「底恒等性」と「Frobenius 次数」を保つ関手は
Frobenius-normalized 性を移す**。 -/
theorem isFrobeniusNormalized_map (Ψ : C₁ ⥤ C₂) [Ψ.Full] [Ψ.Faithful]
    (hdeg : ∀ {X Y : C₁} (φ : X ⟶ Y), P₂.degFr (Ψ.map φ) = P₁.degFr φ)
    (hbase : ∀ {X : C₁} (φ : End X),
      IsBaseIdentity P₁ φ ↔ IsBaseIdentity P₂ (Ψ.map ((φ : X ⟶ X))))
    (A : C₁) (h : IsFrobeniusNormalized P₁ A) :
    IsFrobeniusNormalized P₂ (Ψ.obj A) := by
  refine isFrobeniusNormalized_transport P₁ P₂ (endMulEquiv Ψ A) (fun φ => hbase φ)
    (fun φ => hdeg ((φ : A ⟶ A))) (fun α => ?_) h
  show (IsBaseIdentity P₁ α ∧ P₁.degFr ((α : A ⟶ A)) = 1)
    ↔ (IsBaseIdentity P₂ (Ψ.map ((α : A ⟶ A))) ∧ P₂.degFr (Ψ.map ((α : A ⟶ A))) = 1)
  rw [hdeg ((α : A ⟶ A))]
  exact and_congr (hbase α) Iff.rfl

/-- ★★★★**同じ 2 条件で Frobenius-compact 性も移る**。

★★`𝒪^×` の対応は `𝒪^▷` ＋ `IsUnit` であり、`IsUnit` は**単系の同型で移る**。
★自己同型の対応は充満忠実性(`Functor.FullyFaithful.isoEquiv`)、
共役との両立は `Functor.map_comp` 2 回である。 -/
theorem isFrobeniusCompact_map (Ψ : C₁ ⥤ C₂) [Ψ.Full] [Ψ.Faithful]
    (hdeg : ∀ {X Y : C₁} (φ : X ⟶ Y), P₂.degFr (Ψ.map φ) = P₁.degFr φ)
    (hbase : ∀ {X : C₁} (φ : End X),
      IsBaseIdentity P₁ φ ↔ IsBaseIdentity P₂ (Ψ.map ((φ : X ⟶ X))))
    (A : C₁) (h : IsFrobeniusCompact P₁ A) :
    IsFrobeniusCompact P₂ (Ψ.obj A) := by
  refine isFrobeniusCompact_transport P₁ P₂ (endMulEquiv Ψ A)
    ((Functor.FullyFaithful.ofFullyFaithful Ψ).isoEquiv) (fun u => ?_) (fun θ u => ?_) h
  · show ((IsBaseIdentity P₁ u ∧ P₁.degFr ((u : A ⟶ A)) = 1) ∧ IsUnit u)
      ↔ ((IsBaseIdentity P₂ (Ψ.map ((u : A ⟶ A)))
        ∧ P₂.degFr (Ψ.map ((u : A ⟶ A))) = 1) ∧ IsUnit (endMulEquiv Ψ A u))
    rw [hdeg ((u : A ⟶ A))]
    refine and_congr (and_congr (hbase u) Iff.rfl) ?_
    exact ⟨fun hu => hu.map (endMulEquiv Ψ A).toMonoidHom,
      fun hu => by simpa using hu.map (endMulEquiv Ψ A).symm.toMonoidHom⟩
  · show Ψ.map ((endConj θ u : End A) : A ⟶ A)
      = ((endConj ((Functor.FullyFaithful.ofFullyFaithful Ψ).isoEquiv θ)
          (Ψ.map ((u : A ⟶ A))) : End (Ψ.obj A)) : Ψ.obj A ⟶ Ψ.obj A)
    show Ψ.map (θ.inv ≫ ((u : A ⟶ A)) ≫ θ.hom)
      = (Ψ.mapIso θ).inv ≫ Ψ.map ((u : A ⟶ A)) ≫ (Ψ.mapIso θ).hom
    rw [Ψ.map_comp, Ψ.map_comp]
    rfl

end MapForm

/-! ## ★3. 仮定を「`𝔽_Φ` への関手との 1-可換性」1 本にまとめる

★★★実務では圏同値が **`𝒞 → 𝔽_Φ` と 1-可換**な形で与えられる
(`Theorem 5.2, (iv)` の `modelTypeEquiv_comp_toElem` がそれ)。
★その自然同型 1 本から **`degFr` の保存も「底恒等性」の保存も出る** ——
`Base` と `degFr` は `𝔽_Φ` の射の成分そのものだからである。

★鍵は「`𝔽_Φ` の同型は Frobenius 次数 1」(`elem_deg_of_iso`)である。 -/

section ToElemIso

variable {C₁ : Type u2} [Category.{v2} C₁] {C₂ : Type u4} [Category.{v4} C₂]
  {D₀ : Type u} [Category.{v} D₀] {Φ₀ : MonoidOn.{v, u, w} D₀}
  {P₁ : PreFrobenioid C₁ Φ₀} {P₂ : PreFrobenioid C₂ Φ₀}

/-- ★★**`𝔽_Φ` の同型は Frobenius 次数 1**。 -/
theorem elem_deg_of_iso {A B : ElemFrobCat Φ₀} (e : A ≅ B) : (e.hom : A ⟶ B).deg = 1 := by
  have h1 : ((e.hom ≫ e.inv : A ⟶ A)).deg = ((𝟙 A : A ⟶ A)).deg :=
    congrArg (fun t : A ⟶ A => t.deg) e.hom_inv_id
  rw [ElemFrobCat.comp_deg] at h1
  have h2 : (e.inv : B ⟶ A).deg * (e.hom : A ⟶ B).deg = 1 := h1
  have h3 : (((e.inv : B ⟶ A).deg : ℕ)) * (((e.hom : A ⟶ B).deg : ℕ)) = 1 := by
    exact_mod_cast congrArg (fun n : ℕ+ => (n : ℕ)) h2
  exact PNat.coe_injective (Nat.eq_one_of_mul_eq_one_left h3)

/-- ★★★**`𝔽_Φ` への関手と 1-可換なら Frobenius 次数を保つ**。 -/
theorem degFr_map_of_toElem_iso (Ψ : C₁ ⥤ C₂) (η : Ψ ⋙ P₂.toElem ≅ P₁.toElem)
    {X Y : C₁} (f : X ⟶ Y) : P₂.degFr (Ψ.map f) = P₁.degFr f := by
  have hnat := congrArg (fun t : (Ψ ⋙ P₂.toElem).obj X ⟶ P₁.toElem.obj Y => t.deg)
    (η.hom.naturality f)
  simp only [Functor.comp_map, ElemFrobCat.comp_deg] at hnat
  have hY : (η.hom.app Y).deg = 1 := elem_deg_of_iso (η.app Y)
  have hX : (η.hom.app X).deg = 1 := elem_deg_of_iso (η.app X)
  rw [hY, hX, one_mul, mul_one] at hnat
  exact hnat

/-- ★★★**`𝔽_Φ` への関手と 1-可換なら「底恒等性」を(両向きに)保つ**。

★`η` の成分の `base` は同型なので、`Base` の四角形が両向きに割れる。 -/
theorem isBaseIdentity_map_iff_of_toElem_iso (Ψ : C₁ ⥤ C₂)
    (η : Ψ ⋙ P₂.toElem ≅ P₁.toElem) {X : C₁} (φ : End X) :
    IsBaseIdentity P₁ φ ↔ IsBaseIdentity P₂ (Ψ.map ((φ : X ⟶ X))) := by
  have hnat := congrArg (fun t : (Ψ ⋙ P₂.toElem).obj X ⟶ P₁.toElem.obj X => t.base)
    (η.hom.naturality ((φ : X ⟶ X)))
  simp only [Functor.comp_map, ElemFrobCat.comp_base] at hnat
  have hb : P₂.Base (Ψ.map ((φ : X ⟶ X))) ≫ (η.hom.app X).base
      = (η.hom.app X).base ≫ P₁.Base ((φ : X ⟶ X)) := hnat
  haveI : IsIso ((η.hom.app X).base) := by
    refine ⟨(η.inv.app X).base, ?_, ?_⟩
    · have h := congrArg (fun t : (Ψ ⋙ P₂.toElem).obj X ⟶ (Ψ ⋙ P₂.toElem).obj X => t.base)
        (η.hom_inv_id_app X)
      rw [ElemFrobCat.comp_base] at h
      exact h
    · have h := congrArg (fun t : P₁.toElem.obj X ⟶ P₁.toElem.obj X => t.base)
        (η.inv_hom_id_app X)
      rw [ElemFrobCat.comp_base] at h
      exact h
  constructor
  · intro h
    have hP₁ : P₁.Base ((φ : X ⟶ X)) = 𝟙 ((P₁.toElem.obj X).base) :=
      (show P₁.Base ((φ : X ⟶ X)) = P₁.Base (𝟙 X) from h).trans (P₁.Base_id X)
    refine (cancel_mono ((η.hom.app X).base)).mp ?_
    refine hb.trans ?_
    refine (congrArg (fun t : (P₁.toElem.obj X).base ⟶ (P₁.toElem.obj X).base =>
      (η.hom.app X).base ≫ t) hP₁).trans ?_
    exact (Category.comp_id _).trans
      ((Category.id_comp _).symm.trans
        (congrArg (fun t : (P₂.toElem.obj (Ψ.obj X)).base ⟶ (P₂.toElem.obj (Ψ.obj X)).base =>
          t ≫ (η.hom.app X).base) (P₂.Base_id (Ψ.obj X)).symm))
  · intro h
    have hP₂ : P₂.Base (Ψ.map ((φ : X ⟶ X))) = 𝟙 ((P₂.toElem.obj (Ψ.obj X)).base) :=
      (show P₂.Base (Ψ.map ((φ : X ⟶ X))) = P₂.Base (𝟙 (Ψ.obj X)) from h).trans
        (P₂.Base_id (Ψ.obj X))
    refine (cancel_epi ((η.hom.app X).base)).mp ?_
    refine hb.symm.trans ?_
    refine (congrArg (fun t : (P₂.toElem.obj (Ψ.obj X)).base ⟶ (P₂.toElem.obj (Ψ.obj X)).base =>
      t ≫ (η.hom.app X).base) hP₂).trans ?_
    exact (Category.id_comp _).trans
      ((Category.comp_id _).symm.trans
        (congrArg (fun t : (P₁.toElem.obj X).base ⟶ (P₁.toElem.obj X).base =>
          (η.hom.app X).base ≫ t) (P₁.Base_id X).symm))

/-- ★★★★★**1-可換な充満忠実関手は Frobenius-normalized 性を移す**(使いやすい形)。 -/
theorem isFrobeniusNormalized_map_of_toElem_iso (Ψ : C₁ ⥤ C₂) [Ψ.Full] [Ψ.Faithful]
    (η : Ψ ⋙ P₂.toElem ≅ P₁.toElem) (A : C₁) (h : IsFrobeniusNormalized P₁ A) :
    IsFrobeniusNormalized P₂ (Ψ.obj A) :=
  isFrobeniusNormalized_map Ψ (fun f => degFr_map_of_toElem_iso Ψ η f)
    (fun φ => isBaseIdentity_map_iff_of_toElem_iso Ψ η φ) A h

/-- ★★★★★**1-可換な充満忠実関手は Frobenius-compact 性を移す**(使いやすい形)。 -/
theorem isFrobeniusCompact_map_of_toElem_iso (Ψ : C₁ ⥤ C₂) [Ψ.Full] [Ψ.Faithful]
    (η : Ψ ⋙ P₂.toElem ≅ P₁.toElem) (A : C₁) (h : IsFrobeniusCompact P₁ A) :
    IsFrobeniusCompact P₂ (Ψ.obj A) :=
  isFrobeniusCompact_map Ψ (fun f => degFr_map_of_toElem_iso Ψ η f)
    (fun φ => isBaseIdentity_map_iff_of_toElem_iso Ψ η φ) A h

end ToElemIso

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Proposition 5.5, (iii)` の standard 型の移送に要る一本
(★**条つき**: 4 つの対応を仮定として受け取っている)。 -/
def isFrobeniusNormalized_transport.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — Frobenius-normalized 型の圏をまたぐ移送",
    sectionId := "frdi-prop-5-5" }

/-- ★★★★locator —— `Proposition 5.5, (iii)` の「standard 型が `𝒞^un-tr`・`𝒞^rlf` へ移る」に
要る移送を、**充満忠実で `Base` と `degFr` を保つ関手**の形に束ねたもの。 -/
def isFrobeniusCompact_map.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — Frobenius-compact / -normalized 性の関手による移送",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
