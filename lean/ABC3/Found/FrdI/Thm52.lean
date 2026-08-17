/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.ElementaryFrobenioid

/-!
# [FrdI] Theorem 5.2 (Model Frobenioids) —— (i) 圏の構成

原文 (FrdI p.100):
> (i) A well-defined category C may be constructed in the following fashion.
> The objects of C are pairs of the form

★**入力**は `𝒟` 上のモノイド `Φ`(divisorial)・`B`(group-like)と
準同型 `Div_B : B → Φ^gp` の 3 つである。

★**対象**は `(A_𝒟, α)`(`A_𝒟 ∈ Ob(𝒟)`、`α ∈ Φ(A_𝒟)^gp`)。
★**射** `φ : (A_𝒟, α) → (B_𝒟, β)` は 4 つ組
`(deg_Fr(φ), Base(φ), Div(φ), u_φ)` で、関係式

  `deg_Fr(φ) · α + Div(φ) = Base(φ)^*(β) + Div_B(u_φ)`

を `Φ(A_𝒟)^gp` で満たすもの。

★★**合成は `𝔽_Φ` と同じ形**(`Div` の合成則がそのまま)で、
`u` 成分がそれに並走する。したがって
**忘却関手 `𝒞 → 𝔽_Φ` が自動的に得られる**——これが原文の最後の一文である。

原文 (FrdI p.100):
> [cf. Remark 1.1.1]. Moreover, the Frobenius degree, projection to D, and zero
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w

variable {D : Type u} [Category.{v} D]

/-! ## ★0. 下準備 —— `M → M^gp` を加法準同型として、`Φ^gp` の関手性 -/

/-- ★`M → M^gp` を **加法準同型**として見たもの。 -/
def toGpHom (M : Type w) [AddCommMonoid M] : M →+ Gp M where
  toFun := toGp M
  map_zero' := by
    have h : toGp M 0 + toGp M 0 = toGp M 0 + 0 := by
      rw [add_zero]
      show AddLocalization.mk (0 : M) 0 + AddLocalization.mk (0 : M) 0
        = AddLocalization.mk (0 : M) 0
      rw [AddLocalization.mk_add]
      simp
    exact add_left_cancel h
  map_add' a b := by
    show AddLocalization.mk (a + b) 0 = AddLocalization.mk a 0 + AddLocalization.mk b 0
    rw [AddLocalization.mk_add]
    simp

theorem toGpHom_apply {M : Type w} [AddCommMonoid M] (m : M) :
    toGpHom M m = toGp M m := rfl

namespace MonoidOn

/-- ★**`Φ^gp` の誘導射** —— `α : B ⟶ A` に対する `Gp (Φ(A)) → Gp (Φ(B))`。 -/
noncomputable def gpMapOn (Φ : MonoidOn.{v, u, w} D) {A B : D} (α : B ⟶ A) :
    Gp (Φ.val A) →+ Gp (Φ.val B) := gpMap _ (Φ.map α)

@[simp] theorem gpMapOn_toGpHom (Φ : MonoidOn.{v, u, w} D) {A B : D} (α : B ⟶ A)
    (m : Φ.val A) : Φ.gpMapOn α (toGpHom _ m) = toGpHom _ (Φ.map α m) :=
  gpMap_toGp (Φ.val A) (Φ.map α) m

@[simp] theorem gpMapOn_id (Φ : MonoidOn.{v, u, w} D) (A : D) (x : Gp (Φ.val A)) :
    Φ.gpMapOn (𝟙 A) x = x := by
  have h : Φ.map (𝟙 A) = AddMonoidHom.id (Φ.val A) := by
    ext y; exact Φ.map_id A y
  show gpMap _ (Φ.map (𝟙 A)) x = x
  rw [h, gpMap_id]
  rfl

theorem gpMapOn_comp (Φ : MonoidOn.{v, u, w} D) {A B E : D} (α : B ⟶ A) (β : E ⟶ B)
    (x : Gp (Φ.val A)) : Φ.gpMapOn (β ≫ α) x = Φ.gpMapOn β (Φ.gpMapOn α x) := by
  have h : Φ.map (β ≫ α) = (Φ.map β).comp (Φ.map α) := by
    ext y; exact Φ.map_comp α β y
  show gpMap _ (Φ.map (β ≫ α)) x = _
  rw [h, gpMap_comp]
  rfl

end MonoidOn

/-! ## ★1. 入力データ

原文 (FrdI p.100):
> monoid on D; DivB : B →
-/

/-- ★★★**[FrdI] Theorem 5.2 の入力** —— `Φ`(零因子のモノイド)・
`B`(有理関数のモノイド)・`Div_B : B → Φ^gp`。

★`Div_B` は **`𝒟` 上のモノイドの準同型**なので、`𝒟` の射に沿う自然性を課す。 -/
structure ModelData (D : Type u) [Category.{v} D] where
  /-- 零因子のモノイド `Φ` -/
  phi : MonoidOn.{v, u, w} D
  /-- 有理関数のモノイド `B` -/
  bmon : MonoidOn.{v, u, w} D
  /-- `Div_B : B → Φ^gp` -/
  divB : ∀ A : D, bmon.val A →+ Gp (phi.val A)
  /-- `Div_B` の自然性 -/
  divB_nat : ∀ {A B : D} (f : A ⟶ B) (x : bmon.val B),
    divB A (bmon.map f x) = phi.gpMapOn f (divB B x)

namespace ModelData

variable (M : ModelData.{v, u, w} D)

/-! ## ★2. 対象と射

原文 (FrdI p.100):
> of data as follows: (a) an element degFr
-/

/-- ★★**model Frobenioid の対象** —— `(A_𝒟, α)`。 -/
structure Obj where
  /-- `Base(A) ∈ Ob(𝒟)` -/
  base : D
  /-- `α ∈ Φ(A_𝒟)^gp` -/
  cls : Gp (M.phi.val base)

/-- ★★**model Frobenioid の射**。

原文 (FrdI p.100):
> refer to as the projection to D to

原文 (FrdI p.100):
> shall refer to as the zero divisor of
-/
@[ext]
structure Hom (A B : Obj M) where
  /-- `Base(φ) : A_𝒟 ⟶ B_𝒟` -/
  base : A.base ⟶ B.base
  /-- `Div(φ) ∈ Φ(A_𝒟)` -/
  div : M.phi.val A.base
  /-- `deg_Fr(φ) ∈ ℕ≥1` -/
  deg : ℕ+
  /-- `u_φ ∈ B(A_𝒟)` -/
  u : M.bmon.val A.base
  /-- `deg_Fr(φ) · α + Div(φ) = Base(φ)^*(β) + Div_B(u_φ)` -/
  cond : (deg : ℕ) • A.cls + toGpHom _ div
    = M.phi.gpMapOn base B.cls + M.divB _ u

variable {M}

/-- ★恒等射。 -/
def Hom.id (A : Obj M) : Hom M A A where
  base := 𝟙 A.base
  div := 0
  deg := 1
  u := 0
  cond := by simp

/-- ★合成 —— `Div` は `𝔽_Φ` と同じ形、`u` はそれに並走する。 -/
def Hom.comp {A B E : Obj M} (φ : Hom M A B) (ψ : Hom M B E) : Hom M A E where
  base := φ.base ≫ ψ.base
  div := M.phi.map φ.base ψ.div + (ψ.deg : ℕ) • φ.div
  deg := ψ.deg * φ.deg
  u := M.bmon.map φ.base ψ.u + (ψ.deg : ℕ) • φ.u
  cond := by
    have hφ := φ.cond
    have hψ := ψ.cond
    have hnat := M.divB_nat φ.base ψ.u
    have h1 : M.phi.gpMapOn (φ.base ≫ ψ.base) E.cls
        = M.phi.gpMapOn φ.base (M.phi.gpMapOn ψ.base E.cls) :=
      M.phi.gpMapOn_comp _ _ _
    -- ★右辺を `φ.base^*` でくくり出し、`ψ`・`φ` の関係式を順に入れる
    refine Eq.symm ?_
    calc M.phi.gpMapOn (φ.base ≫ ψ.base) E.cls
            + M.divB _ (M.bmon.map φ.base ψ.u + (ψ.deg : ℕ) • φ.u)
        = M.phi.gpMapOn φ.base (M.phi.gpMapOn ψ.base E.cls + M.divB _ ψ.u)
            + (ψ.deg : ℕ) • M.divB _ φ.u := by
          rw [h1, map_add, map_add, hnat, map_nsmul]
          abel
      _ = M.phi.gpMapOn φ.base ((ψ.deg : ℕ) • B.cls + toGpHom _ ψ.div)
            + (ψ.deg : ℕ) • M.divB _ φ.u := by rw [hψ]
      _ = (ψ.deg : ℕ) • (M.phi.gpMapOn φ.base B.cls + M.divB _ φ.u)
            + toGpHom _ (M.phi.map φ.base ψ.div) := by
          rw [map_add, map_nsmul, MonoidOn.gpMapOn_toGpHom, smul_add]
          abel
      _ = (ψ.deg : ℕ) • ((φ.deg : ℕ) • A.cls + toGpHom _ φ.div)
            + toGpHom _ (M.phi.map φ.base ψ.div) := by rw [hφ]
      _ = ((ψ.deg * φ.deg : ℕ+) : ℕ) • A.cls
            + toGpHom _ (M.phi.map φ.base ψ.div + (ψ.deg : ℕ) • φ.div) := by
          rw [map_add, map_nsmul, smul_add, smul_smul]
          push_cast
          abel

instance : Category (Obj M) where
  Hom A B := Hom M A B
  id := Hom.id
  comp := Hom.comp
  id_comp φ := by
    refine Hom.ext ?_ ?_ ?_ ?_
    · exact Category.id_comp _
    · show M.phi.map (𝟙 _) φ.div + (φ.deg : ℕ) • (0 : M.phi.val _) = φ.div
      simp
    · show φ.deg * 1 = φ.deg
      simp
    · show M.bmon.map (𝟙 _) φ.u + (φ.deg : ℕ) • (0 : M.bmon.val _) = φ.u
      simp
  comp_id φ := by
    refine Hom.ext ?_ ?_ ?_ ?_
    · exact Category.comp_id _
    · show M.phi.map φ.base (0 : M.phi.val _) + ((1 : ℕ+) : ℕ) • φ.div = φ.div
      simp
    · show (1 : ℕ+) * φ.deg = φ.deg
      simp
    · show M.bmon.map φ.base (0 : M.bmon.val _) + ((1 : ℕ+) : ℕ) • φ.u = φ.u
      simp
  assoc φ ψ χ := by
    refine Hom.ext ?_ ?_ ?_ ?_
    · exact Category.assoc _ _ _
    · show M.phi.map (φ.base ≫ ψ.base) χ.div
          + (χ.deg : ℕ) • (M.phi.map φ.base ψ.div + (ψ.deg : ℕ) • φ.div)
        = M.phi.map φ.base (M.phi.map ψ.base χ.div + (χ.deg : ℕ) • ψ.div)
          + ((χ.deg * ψ.deg : ℕ+) : ℕ) • φ.div
      rw [MonoidOn.map_comp, map_add, map_nsmul, smul_add, ← add_assoc, ← mul_smul]
      norm_cast
    · show χ.deg * (ψ.deg * φ.deg) = χ.deg * ψ.deg * φ.deg
      rw [mul_assoc]
    · show M.bmon.map (φ.base ≫ ψ.base) χ.u
          + (χ.deg : ℕ) • (M.bmon.map φ.base ψ.u + (ψ.deg : ℕ) • φ.u)
        = M.bmon.map φ.base (M.bmon.map ψ.base χ.u + (χ.deg : ℕ) • ψ.u)
          + ((χ.deg * ψ.deg : ℕ+) : ℕ) • φ.u
      rw [MonoidOn.map_comp, map_add, map_nsmul, smul_add, ← add_assoc, ← mul_smul]
      norm_cast

@[simp] theorem comp_base {A B E : Obj M} (φ : A ⟶ B) (ψ : B ⟶ E) :
    (φ ≫ ψ).base = φ.base ≫ ψ.base := rfl
@[simp] theorem comp_deg {A B E : Obj M} (φ : A ⟶ B) (ψ : B ⟶ E) :
    (φ ≫ ψ).deg = ψ.deg * φ.deg := rfl
@[simp] theorem comp_div {A B E : Obj M} (φ : A ⟶ B) (ψ : B ⟶ E) :
    (φ ≫ ψ).div = M.phi.map φ.base ψ.div + (ψ.deg : ℕ) • φ.div := rfl
@[simp] theorem comp_u {A B E : Obj M} (φ : A ⟶ B) (ψ : B ⟶ E) :
    (φ ≫ ψ).u = M.bmon.map φ.base ψ.u + (ψ.deg : ℕ) • φ.u := rfl
@[simp] theorem id_base (A : Obj M) : Hom.base (𝟙 A) = 𝟙 A.base := rfl
@[simp] theorem id_deg (A : Obj M) : Hom.deg (𝟙 A) = 1 := rfl
@[simp] theorem id_div (A : Obj M) : Hom.div (𝟙 A) = 0 := rfl
@[simp] theorem id_u (A : Obj M) : Hom.u (𝟙 A) = 0 := rfl

/-! ## ★3. `𝒞 → 𝔽_Φ`

原文 (FrdI p.100):
> divisor determine a functor
-/

/-- ★★★**忘却関手 `𝒞 → 𝔽_Φ`** —— `α` と `u` を落とす。

★★合成則が `𝔽_Φ` と**同じ形**なので、`map_id`・`map_comp` はどちらも `rfl`。 -/
def toElem : Obj M ⥤ ElemFrobCat M.phi where
  obj A := ⟨A.base⟩
  map φ := ⟨φ.base, φ.div, φ.deg⟩
  map_id _ := rfl
  map_comp _ _ := rfl

@[simp] theorem toElem_obj_base (A : Obj M) : ((toElem (M := M)).obj A).base = A.base := rfl
@[simp] theorem toElem_map_base {A B : Obj M} (φ : A ⟶ B) :
    ((toElem (M := M)).map φ).base = φ.base := rfl
@[simp] theorem toElem_map_div {A B : Obj M} (φ : A ⟶ B) :
    ((toElem (M := M)).map φ).div = φ.div := rfl
@[simp] theorem toElem_map_deg {A B : Obj M} (φ : A ⟶ B) :
    ((toElem (M := M)).map φ).deg = φ.deg := rfl

end ModelData

/-- ★★★**[FrdI] Theorem 5.2, (i)** —— model Frobenioid の圏が構成でき、
Frobenius 次数・`𝒟` への射影・零因子が関手 `𝒞 → 𝔽_Φ` を定める。 -/
def ModelData.Obj.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 100, item := "Theorem 5.2, (i) — 圏の構成と 𝔽_Φ への関手",
    sectionId := "frdi-thm-5-2" }

end ABC3.Found.FrdI
