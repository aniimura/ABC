/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.MorphismTypes

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

/-- ★★★locator —— `Proposition 5.5, (iii)` の standard 型の移送に要る一本
(★**条つき**: 4 つの対応を仮定として受け取っている)。 -/
def isFrobeniusNormalized_transport.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — Frobenius-normalized 型の圏をまたぐ移送",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
