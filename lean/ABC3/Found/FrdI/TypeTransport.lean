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
