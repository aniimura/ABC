/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor57Unit
import ABC3.Found.FrdI.Rmk341

/-!
# [FrdI] Corollary 5.7, (iii)(iv) —— base-Frobenius 対は圏同値で移る

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.108。

原文 (FrdI p.108):
> quasi-base-Frobenius pair of a Frobenius-trivial object A2 ∈Ob(C2).

## ★★対象 1 つの上の対

`Proposition 5.6` の `(σ, φ)` は**対象 1 つの上のデータ**である。そこで
`BaseSection` を経由せずに `BaseFrobPair P A` として切り出す:

* `σ : Aut_𝒟(A_𝒟) →* End_𝒞(A)`(底を取ると恒等に戻る)
* `φ : ℕ≥1 →* End_𝒞(A)`(base-identity・Frobenius 型・次数 `n`)

★★`Aut` と `End` はどちらも `f * g = g ≫ f` の規約なので、
`σ` は**そのまま単系準同型**として書ける。

## ★★底の四角形は「素の型」で持つ(逸脱の記録、2026-08-20)

在庫の `otimes_map_of_baseSquare` などは四角形を
`P₁.proj ⋙ ΨBase ≅ Ψ ⋙ P₂.proj` という**合成関手の自然同型**で受け取る。
★これを使うと `(P₁.proj ⋙ ΨBase).obj A` と `ΨBase.obj ((P₁.toElem.obj A).base)` が
**構文上別物**になり、`rw` と `simp` が軒並み外れる(2026-08-20 に実測)。

★そこで**素の型で書いた `BaseSquare`** を置き、合成関手版からは
`BaseSquare.ofNatIso` で作る(中身は同じで、型の見え方だけが違う)。

## ★輸送の中身

| 段 | 使うもの |
|---|---|
| `σ` の底の付け替え | `bs.app A` による共役 ＋ `ΨBase` の充満忠実性(`preimageIso`) |
| `σ` の値の移送 | `Functor.mapEnd A Ψ` |
| 底に戻ること | 四角形の自然性 ＋ `bs.app A` が同型(epi) |
| `φ` の base-identity | 四角形の自然性 |
| `φ` の Frobenius 型 | `isFrobeniusType_map`(`Remark 3.4.1`) |
| `φ` の次数 | `PreservesDegFr` |

★★★**(iv) の「quasi- が外れる」のは `phi_deg` の行**である ——
`Ψ` が次数を保つなら `φ₂ n` の次数もちょうど `n` になるので、
`ℕ≥1` の自己同型で補正する必要がない。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section Pair

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D}

/-- ★★★**[FrdI] Proposition 5.6 の「base-Frobenius pair of `A`」**。

★`BaseSection` を経由せずに、対象 1 つの上のデータとして切り出したもの。 -/
structure BaseFrobPair (P : PreFrobenioid C Φ) (A : C) where
  /-- **`σ`** —— `Aut_𝒟(A_𝒟) → Aut_𝒞(A) ⊆ End_𝒞(A)`。 -/
  sigma : Aut ((P.toElem.obj A).base) →* End A
  /-- **`σ` は底の全射の切断** —— 底を取ると元に戻る。 -/
  sigma_base : ∀ α : Aut ((P.toElem.obj A).base),
    P.Base ((sigma α : End A) : A ⟶ A) = α.hom
  /-- **`φ`** —— `ℕ≥1 → End_𝒞(A)`。 -/
  phi : ℕ+ →* End A
  /-- `φ` の像は base-identity。 -/
  phi_baseId : ∀ n : ℕ+, IsBaseIdentity P ((phi n : End A) : A ⟶ A)
  /-- `φ` の像は Frobenius 型。 -/
  phi_frobType : ∀ n : ℕ+, IsFrobeniusType P ((phi n : End A) : A ⟶ A)
  /-- `φ n` の Frobenius 次数はちょうど `n`(★**(iv) の「quasi- が外れる」条**)。 -/
  phi_deg : ∀ n : ℕ+, P.degFr ((phi n : End A) : A ⟶ A) = n

end Pair

section Transport

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

/-- ★★**底の 1-可換図式を「素の型」で持ったもの**。 -/
structure BaseSquare (Ψ : C₁ ⥤ C₂) (ΨBase : D₁ ⥤ D₂)
    (P₁ : PreFrobenioid C₁ Φ₁) (P₂ : PreFrobenioid C₂ Φ₂) where
  /-- 各対象での同型。 -/
  app : ∀ A : C₁, ΨBase.obj ((P₁.toElem.obj A).base) ≅ (P₂.toElem.obj (Ψ.obj A)).base
  /-- 自然性。 -/
  naturality : ∀ {A B : C₁} (f : A ⟶ B),
    ΨBase.map (P₁.Base f) ≫ (app B).hom = (app A).hom ≫ P₂.Base (Ψ.map f)

/-- ★合成関手の自然同型から `BaseSquare` を作る(中身は同じ)。 -/
def BaseSquare.ofNatIso (Ψ : C₁ ⥤ C₂) (ΨBase : D₁ ⥤ D₂)
    (sq : P₁.proj ⋙ ΨBase ≅ Ψ ⋙ P₂.proj) : BaseSquare Ψ ΨBase P₁ P₂ where
  app A := sq.app A
  naturality f := sq.hom.naturality f

variable {Ψ : C₁ ⥤ C₂} {ΨBase : D₁ ⥤ D₂}

/-- ★★底の四角形から base-identity が移る。 -/
theorem BaseSquare.isBaseIdentity_map (bs : BaseSquare Ψ ΨBase P₁ P₂)
    {A : C₁} (φ : End A) (h : IsBaseIdentity P₁ φ) :
    IsBaseIdentity P₂ (Ψ.map ((φ : End A) : A ⟶ A)) := by
  have hnat := bs.naturality ((φ : End A) : A ⟶ A)
  have hb : P₁.Base ((φ : End A) : A ⟶ A) = 𝟙 _ := by
    rw [show P₁.Base ((φ : End A) : A ⟶ A) = P₁.Base (𝟙 A) from h, P₁.Base_id]
  rw [hb, ΨBase.map_id, Category.id_comp] at hnat
  show P₂.Base (Ψ.map ((φ : End A) : A ⟶ A)) = P₂.Base (𝟙 (Ψ.obj A))
  rw [P₂.Base_id]
  refine (cancel_epi ((bs.app A).hom)).mp ?_
  rw [Category.comp_id]
  exact hnat.symm

/-! ## ★1. 底の自己同型を降ろす -/

/-- ★★底の自己同型を `𝒟₂` から `𝒟₁` へ降ろす —— 四角形で共役を取り、`ΨBase` の
充満忠実性で原像を取る。 -/
noncomputable def autDown [ΨBase.Full] [ΨBase.Faithful]
    (bs : BaseSquare Ψ ΨBase P₁ P₂) (A : C₁) :
    Aut ((P₂.toElem.obj (Ψ.obj A)).base) →* Aut ((P₁.toElem.obj A).base) where
  toFun α := ΨBase.preimageIso (bs.app A ≪≫ α ≪≫ (bs.app A).symm)
  map_one' := by
    refine Iso.ext (ΨBase.map_injective ?_)
    show ΨBase.map (ΨBase.preimageIso (bs.app A ≪≫ (1 : Aut _) ≪≫ (bs.app A).symm)).hom
      = ΨBase.map (𝟙 _)
    rw [Functor.preimageIso_hom, ΨBase.map_preimage, ΨBase.map_id]
    show (bs.app A).hom ≫ 𝟙 _ ≫ (bs.app A).inv = 𝟙 _
    rw [Category.id_comp, Iso.hom_inv_id]
  map_mul' α β := by
    refine Iso.ext (ΨBase.map_injective ?_)
    show ΨBase.map (ΨBase.preimageIso (bs.app A ≪≫ (α * β) ≪≫ (bs.app A).symm)).hom
      = ΨBase.map ((ΨBase.preimageIso (bs.app A ≪≫ β ≪≫ (bs.app A).symm)).hom
          ≫ (ΨBase.preimageIso (bs.app A ≪≫ α ≪≫ (bs.app A).symm)).hom)
    rw [Functor.preimageIso_hom, ΨBase.map_preimage, ΨBase.map_comp,
      Functor.preimageIso_hom, Functor.preimageIso_hom,
      ΨBase.map_preimage, ΨBase.map_preimage]
    show (bs.app A).hom ≫ (β.hom ≫ α.hom) ≫ (bs.app A).inv
      = ((bs.app A).hom ≫ β.hom ≫ (bs.app A).inv)
        ≫ ((bs.app A).hom ≫ α.hom ≫ (bs.app A).inv)
    simp

@[simp] theorem autDown_hom [ΨBase.Full] [ΨBase.Faithful]
    (bs : BaseSquare Ψ ΨBase P₁ P₂) (A : C₁)
    (α : Aut ((P₂.toElem.obj (Ψ.obj A)).base)) :
    ΨBase.map ((autDown bs A) α).hom
      = (bs.app A).hom ≫ α.hom ≫ (bs.app A).inv := by
  show ΨBase.map (ΨBase.preimageIso (bs.app A ≪≫ α ≪≫ (bs.app A).symm)).hom = _
  rw [Functor.preimageIso_hom, ΨBase.map_preimage]
  rfl

/-! ## ★2. 対の輸送 -/

/-- ★★★★★**[FrdI] Corollary 5.7, (iii)(iv)** —— base-Frobenius 対は圏同値で移る。

★★★**(iv) の「quasi- が外れる」のは `phi_deg` の行**である ——
`Ψ` が次数を保つなら `φ₂ n` の次数もちょうど `n` になる。

原文 (FrdI p.108):
> quasi-base-Frobenius pair of a Frobenius-trivial object A2 ∈Ob(C2). -/
noncomputable def BaseFrobPair.map {A : C₁} (S : BaseFrobPair P₁ A)
    [Ψ.IsEquivalence] [ΨBase.Full] [ΨBase.Faithful]
    (bs : BaseSquare Ψ ΨBase P₁ P₂)
    (F₁ : FrobenioidCore P₁) (F₂ : FrobenioidCore P₂)
    (hdeg : PreservesDegFr Ψ P₁ P₂) :
    BaseFrobPair P₂ (Ψ.obj A) where
  sigma := (Functor.mapEnd A Ψ).comp (S.sigma.comp (autDown bs A))
  sigma_base α := by
    set f : End A := S.sigma ((autDown bs A) α) with hf
    have hnat := bs.naturality ((f : End A) : A ⟶ A)
    have hb : P₁.Base ((f : End A) : A ⟶ A) = ((autDown bs A) α).hom := S.sigma_base _
    rw [hb, autDown_hom] at hnat
    have hL : ((bs.app A).hom ≫ α.hom ≫ (bs.app A).inv) ≫ (bs.app A).hom
        = (bs.app A).hom ≫ α.hom := by simp
    rw [hL] at hnat
    show P₂.Base (Ψ.map ((f : End A) : A ⟶ A)) = α.hom
    exact ((cancel_epi ((bs.app A).hom)).mp hnat).symm
  phi := (Functor.mapEnd A Ψ).comp S.phi
  phi_baseId n := bs.isBaseIdentity_map _ (S.phi_baseId n)
  phi_frobType n := isFrobeniusType_map Ψ P₁ P₂ F₁ F₂ hdeg (S.phi_frobType n)
  phi_deg n := by
    show P₂.degFr (Ψ.map ((S.phi n : End A) : A ⟶ A)) = n
    rw [hdeg]
    exact S.phi_deg n

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Corollary 5.7, (iii)(iv)`。 -/
def BaseFrobPair.map.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 108,
    item := "Corollary 5.7, (iii)(iv) — base-Frobenius 対は圏同値で移る",
    sectionId := "frdi-cor-5-7" }

end Transport

end ABC3.Found.FrdI
