/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm52

/-!
# [FrdI] model Frobenioid どうしを結ぶ関手 —— `ModelData` の射

原文 (FrdI p.104):
> [where the vertical arrows are the natural functors of Proposition 5.3; the horizontal

★`Proposition 5.3` の末尾に置かれた 1-可換図式

```
𝒞      ⟶ 𝒞^istr    ⟶ 𝒞^pf
↓                      ↓
𝒞^un-tr ⟶ (𝒞^un-tr)^pf ⟶ 𝒞^rlf
```

の**縦・横の矢印はすべて「因子の単系と有理関数の単系の自然変換」から来る**。
★★そこで、`ModelData` の射(2 つの単系の自然変換で `Div_B` と両立するもの)を定め、
それが model Frobenioid のあいだの関手を誘導することを一度だけ示しておく。

## ★これが効くところ

* `Proposition 5.3` の 1-可換図式の縦の矢印(`Prop53Rlf.lean` の `untrToCrlf`)
* `Corollary 5.4` の `Ψ^rlf`(実化どうしの比較)
* `Proposition 5.5, (iv)` の各種の比較
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w

/-- ★`gpMap` は `toGpHom` と可換。 -/
theorem gpMap_toGpHom {M N : Type w} [AddCommMonoid M] [AddCommMonoid N] (f : M →+ N) (m : M) :
    gpMap _ f (toGpHom M m) = toGpHom N (f m) := gpMap_toGp M f m

variable {D : Type u} [Category.{v} D]

/-- ★★**`ModelData` の射** —— 因子の単系と有理関数の単系の、`Div_B` と両立する自然変換。 -/
structure ModelDataHom (M M' : ModelData.{v, u, w} D) where
  /-- `Φ → Φ'` -/
  phiHom : ∀ d : D, M.phi.val d →+ M'.phi.val d
  /-- その自然性 -/
  phiNat : ∀ {A B : D} (f : B ⟶ A) (x : M.phi.val A),
    phiHom B (M.phi.map f x) = M'.phi.map f (phiHom A x)
  /-- `B → B'` -/
  bmonHom : ∀ d : D, M.bmon.val d →+ M'.bmon.val d
  /-- その自然性 -/
  bmonNat : ∀ {A B : D} (f : B ⟶ A) (x : M.bmon.val A),
    bmonHom B (M.bmon.map f x) = M'.bmon.map f (bmonHom A x)
  /-- `Div_B` との両立 -/
  divCompat : ∀ (d : D) (u : M.bmon.val d),
    M'.divB d (bmonHom d u) = gpMap _ (phiHom d) (M.divB d u)

namespace ModelDataHom

variable {M M' : ModelData.{v, u, w} D} (F : ModelDataHom M M')

/-- ★`Φ^gp` の側の自然性。 -/
theorem gpPhiNat {A B : D} (f : B ⟶ A) (x : Gp (M.phi.val A)) :
    gpMap _ (F.phiHom B) (M.phi.gpMapOn f x) = M'.phi.gpMapOn f (gpMap _ (F.phiHom A) x) := by
  have hc : (F.phiHom B).comp (M.phi.map f) = (M'.phi.map f).comp (F.phiHom A) := by
    ext y; exact F.phiNat f y
  show gpMap _ _ (gpMap _ _ x) = gpMap _ _ (gpMap _ _ x)
  rw [← AddMonoidHom.comp_apply, ← gpMap_comp, ← AddMonoidHom.comp_apply, ← gpMap_comp, hc]

/-- ★対象の対応 —— 底はそのまま、類は `Φ → Φ'` で送る。 -/
noncomputable def obj (A : ModelData.Obj M) : ModelData.Obj M' :=
  ⟨A.base, gpMap _ (F.phiHom A.base) A.cls⟩

@[simp] theorem obj_base (A : ModelData.Obj M) : (F.obj A).base = A.base := rfl

@[simp] theorem obj_cls (A : ModelData.Obj M) :
    (F.obj A).cls = gpMap _ (F.phiHom A.base) A.cls := rfl

/-- ★射の対応 —— 4 成分のうち `Base` と `deg_Fr` はそのまま。 -/
noncomputable def hom {A B : ModelData.Obj M} (φ : A ⟶ B) : F.obj A ⟶ F.obj B where
  base := φ.base
  div := F.phiHom _ φ.div
  deg := φ.deg
  u := F.bmonHom _ φ.u
  cond := by
    have h := congrArg (gpMap _ (F.phiHom A.base)) φ.cond
    rw [map_add, map_add, map_nsmul, gpMap_toGpHom, F.gpPhiNat, ← F.divCompat] at h
    exact h

@[simp] theorem hom_base {A B : ModelData.Obj M} (φ : A ⟶ B) : (F.hom φ).base = φ.base := rfl

@[simp] theorem hom_deg {A B : ModelData.Obj M} (φ : A ⟶ B) : (F.hom φ).deg = φ.deg := rfl

@[simp] theorem hom_div {A B : ModelData.Obj M} (φ : A ⟶ B) :
    (F.hom φ).div = F.phiHom _ φ.div := rfl

@[simp] theorem hom_u {A B : ModelData.Obj M} (φ : A ⟶ B) :
    (F.hom φ).u = F.bmonHom _ φ.u := rfl

/-- ★★★**`ModelData` の射が誘導する関手**。

★`Base` と `deg_Fr` をそのまま保つので、`𝔽_Φ` への関手とも(底の側を `Φ → Φ'` で
送るぶんを除いて)両立する。 -/
noncomputable def functor : ModelData.Obj M ⥤ ModelData.Obj M' where
  obj := F.obj
  map := F.hom
  map_id A := by
    refine ModelData.Hom.ext ?_ ?_ ?_ ?_
    · rfl
    · show F.phiHom _ (0 : M.phi.val A.base) = 0
      exact map_zero _
    · rfl
    · show F.bmonHom _ (0 : M.bmon.val A.base) = 0
      exact map_zero _
  map_comp {A B E} φ ψ := by
    refine ModelData.Hom.ext ?_ ?_ ?_ ?_
    · rfl
    · show F.phiHom _ (M.phi.map φ.base ψ.div + (ψ.deg : ℕ) • φ.div)
        = M'.phi.map φ.base (F.phiHom _ ψ.div) + (ψ.deg : ℕ) • F.phiHom _ φ.div
      rw [map_add, map_nsmul, F.phiNat]
    · rfl
    · show F.bmonHom _ (M.bmon.map φ.base ψ.u + (ψ.deg : ℕ) • φ.u)
        = M'.bmon.map φ.base (F.bmonHom _ ψ.u) + (ψ.deg : ℕ) • F.bmonHom _ φ.u
      rw [map_add, map_nsmul, F.bmonNat]

@[simp] theorem functor_obj (A : ModelData.Obj M) : F.functor.obj A = F.obj A := rfl

@[simp] theorem functor_map {A B : ModelData.Obj M} (φ : A ⟶ B) :
    F.functor.map φ = F.hom φ := rfl

end ModelDataHom

/-! ## ★2. 底の圏が違う場合 —— `Corollary 5.4` の `Ψ^rlf` の形

原文 (FrdI p.104):
> there exists a 1-unique functor Ψrlf : Crlf

★`Corollary 5.4` は **2 つの Frobenioid の実化のあいだ**の関手を作る。
底の圏が `𝒟₁`, `𝒟₂` と違うので、`ModelDataHom` を**底の関手 `Ψ_𝒟` の上で**一般化する。
★`Ψ_𝒟 = 𝟭` の場合が §1 の `ModelDataHom` である。 -/

variable {D₂ : Type u} [Category.{v} D₂]

/-- ★★**底の関手 `Ψ_𝒟 : 𝒟₁ ⥤ 𝒟₂` の上の `ModelData` の射**。 -/
structure ModelDataHomOver (ΨB : D ⥤ D₂) (M : ModelData.{v, u, w} D)
    (M₂ : ModelData.{v, u, w} D₂) where
  /-- `Φ₁ → Ψ_𝒟^* Φ₂` -/
  phiHom : ∀ d : D, M.phi.val d →+ M₂.phi.val (ΨB.obj d)
  /-- その自然性 -/
  phiNat : ∀ {A B : D} (f : B ⟶ A) (x : M.phi.val A),
    phiHom B (M.phi.map f x) = M₂.phi.map (ΨB.map f) (phiHom A x)
  /-- `B₁ → Ψ_𝒟^* B₂` -/
  bmonHom : ∀ d : D, M.bmon.val d →+ M₂.bmon.val (ΨB.obj d)
  /-- その自然性 -/
  bmonNat : ∀ {A B : D} (f : B ⟶ A) (x : M.bmon.val A),
    bmonHom B (M.bmon.map f x) = M₂.bmon.map (ΨB.map f) (bmonHom A x)
  /-- `Div_B` との両立 -/
  divCompat : ∀ (d : D) (u : M.bmon.val d),
    M₂.divB (ΨB.obj d) (bmonHom d u) = gpMap _ (phiHom d) (M.divB d u)

namespace ModelDataHomOver

variable {ΨB : D ⥤ D₂} {M : ModelData.{v, u, w} D} {M₂ : ModelData.{v, u, w} D₂}
  (F : ModelDataHomOver ΨB M M₂)

theorem gpPhiNat {A B : D} (f : B ⟶ A) (x : Gp (M.phi.val A)) :
    gpMap _ (F.phiHom B) (M.phi.gpMapOn f x)
      = M₂.phi.gpMapOn (ΨB.map f) (gpMap _ (F.phiHom A) x) := by
  have hc : (F.phiHom B).comp (M.phi.map f) = (M₂.phi.map (ΨB.map f)).comp (F.phiHom A) := by
    ext y; exact F.phiNat f y
  show gpMap _ _ (gpMap _ _ x) = gpMap _ _ (gpMap _ _ x)
  rw [← AddMonoidHom.comp_apply, ← gpMap_comp, ← AddMonoidHom.comp_apply, ← gpMap_comp, hc]

/-- ★対象の対応。 -/
noncomputable def obj (A : ModelData.Obj M) : ModelData.Obj M₂ :=
  ⟨ΨB.obj A.base, gpMap _ (F.phiHom A.base) A.cls⟩

@[simp] theorem obj_base (A : ModelData.Obj M) : (F.obj A).base = ΨB.obj A.base := rfl

/-- ★射の対応。 -/
noncomputable def hom {A B : ModelData.Obj M} (φ : A ⟶ B) : F.obj A ⟶ F.obj B where
  base := ΨB.map φ.base
  div := F.phiHom _ φ.div
  deg := φ.deg
  u := F.bmonHom _ φ.u
  cond := by
    have h := congrArg (gpMap _ (F.phiHom A.base)) φ.cond
    rw [map_add, map_add, map_nsmul, gpMap_toGpHom, F.gpPhiNat, ← F.divCompat] at h
    exact h

@[simp] theorem hom_base {A B : ModelData.Obj M} (φ : A ⟶ B) :
    (F.hom φ).base = ΨB.map φ.base := rfl

@[simp] theorem hom_deg {A B : ModelData.Obj M} (φ : A ⟶ B) : (F.hom φ).deg = φ.deg := rfl

/-- ★★★**底の関手の上の `ModelData` の射が誘導する関手**。

★`deg_Fr` はそのまま、`Base` は `Ψ_𝒟` で送る。 -/
noncomputable def functor : ModelData.Obj M ⥤ ModelData.Obj M₂ where
  obj := F.obj
  map := F.hom
  map_id A := by
    refine ModelData.Hom.ext ?_ ?_ ?_ ?_
    · exact ΨB.map_id A.base
    · show F.phiHom _ (0 : M.phi.val A.base) = 0
      exact map_zero _
    · rfl
    · show F.bmonHom _ (0 : M.bmon.val A.base) = 0
      exact map_zero _
  map_comp {A B E} φ ψ := by
    refine ModelData.Hom.ext ?_ ?_ ?_ ?_
    · exact ΨB.map_comp _ _
    · show F.phiHom _ (M.phi.map φ.base ψ.div + (ψ.deg : ℕ) • φ.div)
        = M₂.phi.map (ΨB.map φ.base) (F.phiHom _ ψ.div) + (ψ.deg : ℕ) • F.phiHom _ φ.div
      rw [map_add, map_nsmul, F.phiNat]
    · rfl
    · show F.bmonHom _ (M.bmon.map φ.base ψ.u + (ψ.deg : ℕ) • φ.u)
        = M₂.bmon.map (ΨB.map φ.base) (F.bmonHom _ ψ.u) + (ψ.deg : ℕ) • F.bmonHom _ φ.u
      rw [map_add, map_nsmul, F.bmonNat]

end ModelDataHomOver

end ABC3.Found.FrdI
