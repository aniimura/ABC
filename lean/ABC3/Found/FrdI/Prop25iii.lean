import ABC3.Found.FrdI.Prop25
import ABC3.Found.FrdI.Def24
import ABC3.Found.FrdI.PlBkShuffle

/-!
# [FrdI] Proposition 2.5, (iii) へ向けた 4 重分解

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.49。

原文 (FrdI p.49):
> Finally, we consider assertion (iii). By applying the factorizations of Definition

原文 (FrdI p.49):
> alences of categories of Definition 1.3, (iii), (d)], we conclude that every morphism

## ★この段の内容

`Proposition 2.5, (iii)`(unit-linear Frobenius 函手 `Ψ : 𝒞 ≃ 𝒞(d)`)の
証明が最初に用意する **4 重分解**

  `φ = δ ≫ γ ≫ β ≫ α`
  (`δ` Frobenius 型、`γ` 等長 pre-step、`β` **base-identity な pre-step 自己射**、
   `α` pull-back 射)

をここで取る。★**原文が「`Definition 1.3, (iv), (a)`; `(v), (c)` と
assertion (i) の全単射(および `Definition 1.3, (iii), (d)` の圏同値)を当てる」と
一言で述べる箇所**である。

★★**要点は `β` を「自己射」に、しかも「base-identity」にすること。**
`Definition 1.3` の分解が与えるのは co-angular pre-step `β₀ : X ⟶ Y` までで、
- **metrically trivial 型**が `Y ≅ X` を与えて**自己射**にし、
- **Aut-ample 型**がずれた底の自己同型を打ち消して **base-identity** にする。

★`δ`・`γ`・`α` は**どれも等長射**である(Frobenius 型は LB-invertible、
pull-back も LB-invertible)。★**これが (a)「`Ψ` は等長射の上で恒等」から
`Ψ` が一意に決まる理由**でもある —— 非等長な部分は `β` に集まっている。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) {A : C}

include P in
/-- ★★★**4 重分解** —— metrically trivial ＋ Aut-ample 型のもとで、
どの射も `δ ≫ γ ≫ β ≫ α` と分解する。

| 因子 | 型 |
|---|---|
| `δ` | Frobenius 型 |
| `γ` | 等長 pre-step |
| `β` | **base-identity な pre-step 自己射**(＝ `𝒪^▷(X)` の元) |
| `α` | pull-back 射 | -/
theorem quadFactor (Fc : FrobenioidCore P)
    (hmt : ∀ X : C, IsMetricallyTrivial P X) (haa : IsOfAutAmpleType P)
    {A B : C} (φ : A ⟶ B) :
    ∃ (X Y : C) (δ : A ⟶ X) (γ : X ⟶ Y) (β : Y ⟶ Y) (α : Y ⟶ B),
      φ = δ ≫ γ ≫ β ≫ α ∧ IsFrobeniusType P δ ∧
        IsIsometric P γ ∧ IsPreStep P γ ∧
        IsBaseIdentity P β ∧ IsPreStep P β ∧ IsCoAngular P β ∧ IsPullBack P α := by
  -- 手 1: `Definition 1.3, (iv), (a)`
  obtain ⟨X, Z, δ, p, α₀, hfac, hδ, hp, hα₀⟩ := Fc.arbFactor φ
  -- 手 2: `Definition 1.3, (v), (c)` —— pre-step を「等長 ≫ co-angular」に
  obtain ⟨Y, γ, β₀, hpfac, hγi, hγs, hβ₀c, hβ₀s⟩ := Fc.preStepFactor' p hp
  -- 手 3: metrically trivial 型で `β₀` の終域を `Y` に戻す
  obtain ⟨e⟩ := hmt Y Z β₀ hβ₀c hβ₀s
  -- 手 4: Aut-ample 型で底を `𝟙` にする
  have hβ₁s : IsPreStep P (β₀ ≫ e.hom) :=
    IsPreStep.comp P hβ₀s (isPreStep_of_isIso P e.hom)
  haveI hb₁ : IsIso (P.Base (β₀ ≫ e.hom)) := hβ₁s.2
  obtain ⟨u, hui, hub⟩ := haa Y (P.Base (β₀ ≫ e.hom)) hb₁
  haveI := hui
  obtain ⟨w, hw1, hw2⟩ : ∃ w : Y ⟶ Y, (u : Y ⟶ Y) ≫ w = 𝟙 Y ∧ w ≫ (u : Y ⟶ Y) = 𝟙 Y :=
    ⟨inv (u : Y ⟶ Y), IsIso.hom_inv_id _, IsIso.inv_hom_id _⟩
  haveI hwi : IsIso w := ⟨(u : Y ⟶ Y), hw2, hw1⟩
  refine ⟨X, Y, δ, γ, (β₀ ≫ e.hom) ≫ w, ((u : Y ⟶ Y) ≫ e.inv) ≫ α₀,
    ?_, hδ, hγi, hγs, ?_, ?_, ?_, ?_⟩
  · -- 分解が元の射に戻る
    have hkey : ((β₀ ≫ e.hom) ≫ w) ≫ (((u : Y ⟶ Y)) ≫ e.inv) ≫ α₀ = β₀ ≫ α₀ := by
      simp only [Category.assoc]
      rw [← Category.assoc w, hw2, Category.id_comp, e.hom_inv_id_assoc]
    rw [hfac, hpfac, hkey]
    simp
  · -- `β` は base-identity
    show P.Base ((β₀ ≫ e.hom) ≫ w) = P.Base (𝟙 Y)
    rw [P.Base_comp, P.Base_id, ← hub, ← P.Base_comp, hw1, P.Base_id]
  · -- `β` は pre-step
    exact IsPreStep.comp P hβ₁s (isPreStep_of_isIso P w)
  · -- ★`β` は co-angular —— `β₀` が co-angular で、残りは同型
    exact Fc.coAngularComp _ w
      (Fc.coAngularComp _ e.hom hβ₀c (isCoAngular_of_isIso P e.hom))
      (isCoAngular_of_isIso P w)
  · -- `α` は pull-back
    exact IsPullBack.comp P
      (IsPullBack.comp P (isPullBack_of_isIso P (u : Y ⟶ Y)) (isPullBack_of_isIso P e.inv))
      hα₀

include P in
/-- ★**4 重分解の `δ`・`γ`・`α` はどれも等長射**。

★**これが原文の (a)「`Ψ` は対象と等長射の上で恒等」から
`Ψ` が一意に決まる理由**である —— 非等長な部分は `β` にすべて集まる。 -/
theorem quadFactor_isometric (Fc : FrobenioidCore P) {A X Y B : C}
    {δ : A ⟶ X} {γ : X ⟶ Y} {α : Y ⟶ B}
    (hδ : IsFrobeniusType P δ) (hγ : IsIsometric P γ) (hα : IsPullBack P α) :
    IsIsometric P δ ∧ IsIsometric P γ ∧ IsIsometric P α :=
  ⟨hδ.1.2, hγ, (Fc.pullBackLB α hα).1.2⟩

include P in
/-- ★**`β` は `𝒪^▷(Y)` の元**である。 -/
theorem quadFactor_mem_otri {Y : C} {β : Y ⟶ Y}
    (hβb : IsBaseIdentity P β) (hβs : IsPreStep P β) : β ∈ OTri P Y :=
  ⟨hβb, hβs.1⟩

include P in
/-- ★★**分解から `Div` が読める** —— `Div φ = Φ.map (Base (δ ≫ γ)) (Div β)`。

★`δ`・`γ`・`α` が等長で `γ`・`β`・`α` が linear なので、
**`Div` はすべて `β` から来る**。 -/
theorem quadFactor_div (Fc : FrobenioidCore P) {A X Y B : C}
    {δ : A ⟶ X} {γ : X ⟶ Y} {β : Y ⟶ Y} {α : Y ⟶ B}
    (hδ : IsFrobeniusType P δ) (hγi : IsIsometric P γ) (hγs : IsPreStep P γ)
    (hβb : IsBaseIdentity P β) (hβs : IsPreStep P β) (hα : IsPullBack P α) :
    P.Div (δ ≫ γ ≫ β ≫ α) = Φ.map (P.Base (δ ≫ γ)) (P.Div β) := by
  have hαd : P.Div α = 0 := (Fc.pullBackLB α hα).1.2
  have hαdeg : P.degFr α = 1 := (Fc.pullBackLB α hα).2
  have hβbb : P.Base β = 𝟙 _ := by
    have h : P.Base β = P.Base (𝟙 _) := hβb
    rwa [P.Base_id] at h
  have h1 : P.Div (β ≫ α) = P.Div β := by
    rw [P.Div_comp, hαd, hβbb, MonoidOn.map_id, hαdeg]
    simp
  have h2 : P.Div (γ ≫ β ≫ α) = Φ.map (P.Base γ) (P.Div β) := by
    rw [P.Div_comp, h1, hγi]
    simp
  rw [P.Div_comp, h2, show P.Div δ = 0 from hδ.1.2, P.Base_comp]
  simp

include P in
/-- ★★★**底が同型な射の 4 重分解では `α` が同型になる**。

★★段 6 で `plBk_shuffle` の出す `ρ̃`(`Base ρ̃` が同型)を分解するときに、
**pull-back 部分が消えて `δ₂ ≫ γ₂ ≫ β₂` になる**ことを保証する。 -/
theorem baseIso_quadFactor (Fc : FrobenioidCore P)
    (hmt : ∀ X : C, IsMetricallyTrivial P X) (haa : IsOfAutAmpleType P)
    {A B : C} (φ : A ⟶ B) (hb : IsIso (P.Base φ)) :
    ∃ (X Y : C) (δ : A ⟶ X) (γ : X ⟶ Y) (β : Y ⟶ Y) (α : Y ⟶ B),
      φ = δ ≫ γ ≫ β ≫ α ∧ IsFrobeniusType P δ ∧
        IsIsometric P γ ∧ IsPreStep P γ ∧
        IsBaseIdentity P β ∧ IsPreStep P β ∧ IsCoAngular P β ∧ IsIso α := by
  obtain ⟨X, Y, δ, γ, β, α, hfac, hδ, hγi, hγs, hβb, hβs, hβc, hα⟩ :=
    quadFactor P Fc hmt haa φ
  haveI := hb
  haveI : IsIso (P.Base δ) := hδ.2
  haveI : IsIso (P.Base γ) := hγs.2
  haveI : IsIso (P.Base β) := hβs.2
  haveI hbα : IsIso (P.Base α) := by
    have h0 : P.Base φ = P.Base δ ≫ P.Base γ ≫ P.Base β ≫ P.Base α := by
      rw [hfac, P.Base_comp, P.Base_comp, P.Base_comp]
    have h1 : IsIso (P.Base δ ≫ P.Base γ ≫ P.Base β ≫ P.Base α) := h0 ▸ hb
    have h2 : IsIso (P.Base γ ≫ P.Base β ≫ P.Base α) :=
      IsIso.of_isIso_comp_left (P.Base δ) _
    have h3 : IsIso (P.Base β ≫ P.Base α) := IsIso.of_isIso_comp_left (P.Base γ) _
    exact IsIso.of_isIso_comp_left (P.Base β) _
  exact ⟨X, Y, δ, γ, β, α, hfac, hδ, hγi, hγs, hβb, hβs, hβc,
    isIso_of_isPullBack_of_baseIso P hα hbα⟩

/-! ## ★★`𝒪^▷(A)` は Frobenius-normalized 型なら**可換**

★`Definition 1.2, (iv)` の `Frobenius-normalized` は
`φ ≫ α^{degFr φ} = α ≫ φ`(`φ` base-identity、`α ∈ 𝒪^▷`)である。
★★**`φ` も `𝒪^▷` の元に取ると `degFr φ = 1` なので `φ ≫ α = α ≫ φ`** ——
すなわち `𝒪^▷(A)` は可換モノイドになる。

★**これが `Ψ(β) = β₀ · β₁^d` が乗法的である根拠**である。
-/

include P in
/-- ★★**Frobenius-normalized 型なら `𝒪^▷(A)` は可換**。 -/
theorem otri_mul_comm (hfn : IsFrobeniusNormalized P A) (x y : OTri P A) :
    x * y = y * x := by
  have hx1 : P.degFr ((x : End A) : A ⟶ A) = 1 := x.2.2
  have h := hfn (x : End A) x.2.1 (y : End A) y.2
  rw [hx1] at h
  refine Subtype.ext ?_
  show ((y : End A) : A ⟶ A) ≫ ((x : End A) : A ⟶ A)
    = ((x : End A) : A ⟶ A) ≫ ((y : End A) : A ⟶ A)
  simpa using h.symm

/-! ## ★`Div` と冪 -/

include P in
/-- ★**`𝒪^▷` の元の冪の `Div`** —— `Div (t^k) = k • Div t`。

★`t` は base-identity かつ linear なので、合成則が単純な足し算に潰れる。 -/
theorem otri_div_pow (t : OTri P A) (k : ℕ) :
    P.Div (((t ^ k : OTri P A) : End A) : A ⟶ A) = k • P.Div ((t : End A) : A ⟶ A) := by
  have htb : P.Base ((t : End A) : A ⟶ A) = 𝟙 _ := by
    have h : P.Base ((t : End A) : A ⟶ A) = P.Base (𝟙 A) := t.2.1
    rwa [P.Base_id] at h
  induction k with
  | zero => simpa using P.Div_id A
  | succ n ih =>
    have hstep : (((t ^ (n + 1) : OTri P A) : End A) : A ⟶ A)
        = (((t ^ n : OTri P A) : End A) : A ⟶ A) ≫ ((t : End A) : A ⟶ A) := by
      rw [pow_succ']
      rfl
    have htnb : P.Base (((t ^ n : OTri P A) : End A) : A ⟶ A) = 𝟙 _ := by
      have h : P.Base (((t ^ n : OTri P A) : End A) : A ⟶ A) = P.Base (𝟙 A) := (t ^ n).2.1
      rwa [P.Base_id] at h
    rw [hstep, P.Div_comp, htnb, MonoidOn.map_id, ih,
      show P.degFr ((t : End A) : A ⟶ A) = 1 from t.2.2, succ_nsmul]
    simp [add_comm]

include P in
/-- ★**単元を掛けても `Div` は変わらない**。 -/
theorem otri_div_unit_mul (u : OTimes P A) (x : OTri P A) :
    P.Div ((((⟨(u : End A), OTimes_le_OTri P A u.2⟩ : OTri P A) * x : OTri P A) : End A) : A ⟶ A)
      = P.Div ((x : End A) : A ⟶ A) := by
  haveI : IsIso (((u : End A)) : A ⟶ A) := (CategoryTheory.isUnit_iff_isIso _).mp u.2.2
  have hxb : P.Base ((x : End A) : A ⟶ A) = 𝟙 _ := by
    have h : P.Base ((x : End A) : A ⟶ A) = P.Base (𝟙 A) := x.2.1
    rwa [P.Base_id] at h
  show P.Div (((x : End A) : A ⟶ A) ≫ (((u : End A)) : A ⟶ A)) = _
  rw [P.Div_comp, hxb, MonoidOn.map_id,
    show P.Div (((u : End A)) : A ⟶ A) = 0 from isIsometric_of_isIso P _,
    show P.degFr (((u : End A)) : A ⟶ A) = 1 from u.2.1.2]
  simp

/-! ## ★★`Ψ` の `𝒪^▷(A)` 上の定義 —— `β = β₀ · β₁ ↦ β₀ · β₁^d`

原文 (FrdI p.49):
> we obtain a factorization

★**characteristic splitting `τ` が `β` を「単元部分」と「`τ` 部分」に分け、
`τ` 部分だけを `d` 乗する。**
-/

section Psi

variable (F : FrobenioidCore P) {τ : ∀ X : C, Submonoid (End X)}
  (hτ : IsCharacteristicSplitting P F τ) (hA : IsIsotropic P A)
  (hfn : IsFrobeniusNormalized P A)

/-- ★`𝒪^×(A)` の元を `𝒪^▷(A)` の元と見る。 -/
def uOf (u : OTimes P A) : OTri P A := ⟨(u : End A), OTimes_le_OTri P A u.2⟩

/-- ★`τ(A)` の元を `𝒪^▷(A)` の元と見る。 -/
def tOf (t : τ A) : OTri P A := ⟨(t : End A), hτ.le_otri A hA t.2⟩

/-- ★★**分裂の逆** —— `Definition 2.3, (a)` の全単射から。 -/
noncomputable def splitEquiv : (OTimes P A × τ A) ≃ OTri P A :=
  Equiv.ofBijective _ (charSplitting_bijective P F hτ hA)

theorem splitEquiv_apply (p : OTimes P A × τ A) :
    splitEquiv P F hτ hA p = uOf P p.1 * tOf P F hτ hA p.2 := rfl

/-- ★★★**`Ψ` の `𝒪^▷(A)` 上の定義** —— `β₀ · β₁ ↦ β₀ · β₁^d`。 -/
noncomputable def psiOTri (d : ℕ+) (β : OTri P A) : OTri P A :=
  uOf P ((splitEquiv P F hτ hA).symm β).1
    * (tOf P F hτ hA ((splitEquiv P F hτ hA).symm β).2) ^ ((d : ℕ+) : ℕ)

include hfn in
/-- ★★**`Div Ψ(β) = d • Div β`** —— 単元部分は `Div` に効かず、
`τ` 部分の `d` 乗が `Div` を `d` 倍する。

★**これが原文の (b)「`d` の Frobenius 函手と 1-compatible」の中身**である。 -/
theorem psiOTri_div (d : ℕ+) (β : OTri P A) :
    P.Div (((psiOTri P F hτ hA d β : OTri P A) : End A) : A ⟶ A)
      = ((d : ℕ+) : ℕ) • P.Div ((β : End A) : A ⟶ A) := by
  set p := (splitEquiv P F hτ hA).symm β with hp
  have hβ : β = uOf P p.1 * tOf P F hτ hA p.2 := by
    rw [hp, ← splitEquiv_apply P F hτ hA]
    exact ((splitEquiv P F hτ hA).apply_symm_apply β).symm
  have h1 : P.Div ((β : End A) : A ⟶ A) = P.Div (((tOf P F hτ hA p.2 : OTri P A) : End A)) := by
    rw [hβ]
    exact otri_div_unit_mul P p.1 (tOf P F hτ hA p.2)
  have h2 : P.Div (((psiOTri P F hτ hA d β : OTri P A) : End A))
      = P.Div ((((tOf P F hτ hA p.2) ^ ((d : ℕ+) : ℕ) : OTri P A) : End A)) :=
    otri_div_unit_mul P p.1 _
  rw [h2, otri_div_pow P (tOf P F hτ hA p.2) ((d : ℕ+) : ℕ), h1]

/-! ### ★★分裂はモノイド同型 —— そこから `Ψ` の乗法性が出る -/

include hfn in
/-- ★**可換性から出る並べ替え** —— `a * b * (c * e) = a * c * (b * e)`。 -/
theorem otri_mul_mul_comm (a b c e : OTri P A) : a * b * (c * e) = a * c * (b * e) := by
  have hcomm : ∀ x y : OTri P A, x * y = y * x := otri_mul_comm P hfn
  rw [mul_assoc a, ← mul_assoc b, hcomm b c, mul_assoc c, ← mul_assoc a]

include hfn in
/-- ★★**分裂はモノイド同型**。

★`(u₁u₂)(t₁t₂) = (u₁t₁)(u₂t₂)` は**可換性そのもの**である。 -/
noncomputable def splitMulEquiv : (OTimes P A × τ A) ≃* OTri P A where
  toEquiv := splitEquiv P F hτ hA
  map_mul' := fun p q => by
    show uOf P p.1 * uOf P q.1 * (tOf P F hτ hA p.2 * tOf P F hτ hA q.2)
      = uOf P p.1 * tOf P F hτ hA p.2 * (uOf P q.1 * tOf P F hτ hA q.2)
    exact otri_mul_mul_comm P hfn _ _ _ _

include hfn in
/-- ★★**`Ψ` は乗法的** —— 分裂がモノイド同型で、`ᵒ^▷` が可換だから。 -/
theorem psiOTri_mul (d : ℕ+) (x y : OTri P A) :
    psiOTri P F hτ hA d (x * y)
      = psiOTri P F hτ hA d x * psiOTri P F hτ hA d y := by
  have hcomm : ∀ a b : OTri P A, a * b = b * a := otri_mul_comm P hfn
  have hsymm : (splitEquiv P F hτ hA).symm (x * y)
      = (splitEquiv P F hτ hA).symm x * (splitEquiv P F hτ hA).symm y :=
    (splitMulEquiv P F hτ hA hfn).symm.map_mul x y
  show uOf P (((splitEquiv P F hτ hA).symm (x * y)).1)
      * (tOf P F hτ hA (((splitEquiv P F hτ hA).symm (x * y)).2)) ^ ((d : ℕ+) : ℕ) = _
  rw [hsymm]
  show uOf P (((splitEquiv P F hτ hA).symm x).1) * uOf P (((splitEquiv P F hτ hA).symm y).1)
      * (tOf P F hτ hA (((splitEquiv P F hτ hA).symm x).2)
          * tOf P F hτ hA (((splitEquiv P F hτ hA).symm y).2)) ^ ((d : ℕ+) : ℕ) = _
  rw [Commute.mul_pow (hcomm _ _), otri_mul_mul_comm P hfn]
  rfl

include hfn in
/-- ★**`Ψ` は単位元を保つ**。 -/
theorem psiOTri_one (d : ℕ+) : psiOTri P F hτ hA d 1 = 1 := by
  have h1 : (splitEquiv P F hτ hA).symm 1 = 1 :=
    (splitMulEquiv P F hτ hA hfn).symm.map_one
  show uOf P (((splitEquiv P F hτ hA).symm 1).1)
      * (tOf P F hτ hA (((splitEquiv P F hτ hA).symm 1).2)) ^ ((d : ℕ+) : ℕ) = 1
  rw [h1]
  show (1 : OTri P A) * (1 : OTri P A) ^ ((d : ℕ+) : ℕ) = 1
  rw [one_pow, one_mul]

include hfn in
/-- ★★★**`Ψ` は `ᵒ^▷(A)` のモノイド自己準同型**。 -/
noncomputable def psiOTriHom (d : ℕ+) : OTri P A →* OTri P A where
  toFun := psiOTri P F hτ hA d
  map_one' := psiOTri_one P F hτ hA hfn d
  map_mul' := psiOTri_mul P F hτ hA hfn d

include hfn in
/-- ★★**原文が「elementary computation」と呼ぶ等式** ——
`Ψ(β^{d'}) = Ψ(β₀^{d'} · β₁^{d'}) = β₀^{d'} · (β₁^{d'})^d = (β₀ · β₁^d)^{d'} = Ψ(β)^{d'}`。

★★これは `psiOTriHom` が**モノイド準同型**であることに尽きる ——
そしてその根拠は `𝒪^▷(A)` の可換性、すなわち **Frobenius-normalized 型**である。
★**合成との両立(関手性)を `𝒞^istr` 上で出すときに使う。** -/
theorem psiOTri_pow (d : ℕ+) (β : OTri P A) (k : ℕ) :
    psiOTri P F hτ hA d (β ^ k) = (psiOTri P F hτ hA d β) ^ k :=
  map_pow (psiOTriHom P F hτ hA hfn d) β k

include hfn in
/-- ★★**`𝒪^▷` を Frobenius 型射の右へ通す** —— `Definition 1.2, (iv)` の
Frobenius-normalized の定義そのもの。

  `β ≫ ζ = ζ ≫ β^{degFr ζ}`

★★**`ζ` が base-identity である必要がある**。`𝒞^istr` では
`Proposition 2.5, (ii)`(すべての対象が Frobenius-trivial)により
**Frobenius 型射を base-identity 自己射に取れる**ので、この形が使える。
★これが原文が関手性の議論を `𝒞^istr` に限る理由である。 -/
theorem otri_comp_frobBaseId {ζ : End A} (hζb : IsBaseIdentity P ζ) (β : OTri P A) :
    ((β : End A) : A ⟶ A) ≫ ((ζ : End A) : A ⟶ A)
      = ((ζ : End A) : A ⟶ A)
        ≫ (((β ^ (P.degFr ((ζ : End A) : A ⟶ A) : ℕ) : OTri P A) : End A) : A ⟶ A) :=
  (hfn ζ hζb (β : End A) β.2).symm

/-! ### ★★`Ψ` の自然性 —— well-defined 性の核

★4 重分解には**選び方の自由**がある(原文の
`(α∘ε, ε⁻¹∘β∘ζ, ζ⁻¹∘γ∘θ, θ⁻¹∘δ)`)。★**その自由度は 2 種類**である:

1. **単元倍** —— `β` が `𝒪^×` の元だけずれる
2. **同型による共役** —— `β` が別の対象へ移る

★**`Ψ` はどちらにも可換である**ことをここで示す。
-/

/-- ★★**単元倍との可換性** —— `Ψ(u · β) = u · Ψ(β)`。

★分裂の単元成分に `u` が乗るだけで、`τ` 成分は変わらない。 -/
theorem psiOTri_unit_mul (d : ℕ+) (u : OTimes P A) (β : OTri P A) :
    psiOTri P F hτ hA d (uOf P u * β) = uOf P u * psiOTri P F hτ hA d β := by
  set p := (splitEquiv P F hτ hA).symm β with hp
  have hβ : uOf P p.1 * tOf P F hτ hA p.2 = β := by
    rw [hp, ← splitEquiv_apply P F hτ hA]
    exact (splitEquiv P F hτ hA).apply_symm_apply β
  have hsplit : (splitEquiv P F hτ hA).symm (uOf P u * β) = (u * p.1, p.2) := by
    refine (splitEquiv P F hτ hA).symm_apply_eq.mpr ?_
    rw [splitEquiv_apply]
    show uOf P u * β = uOf P (u * p.1) * tOf P F hτ hA p.2
    rw [show uOf P (u * p.1) = uOf P u * uOf P p.1 from rfl, mul_assoc, hβ]
  show uOf P (((splitEquiv P F hτ hA).symm (uOf P u * β)).1)
      * (tOf P F hτ hA (((splitEquiv P F hτ hA).symm (uOf P u * β)).2)) ^ ((d : ℕ+) : ℕ) = _
  rw [hsplit]
  show uOf P (u * p.1) * (tOf P F hτ hA p.2) ^ ((d : ℕ+) : ℕ) = _
  rw [show uOf P (u * p.1) = uOf P u * uOf P p.1 from rfl, mul_assoc]
  rfl

/-- ★★**等長な元では `Ψ` は恒等** —— 原文 (a)「`Ψ` は等長射の上で恒等」の中身。

★`Div β = 0` なら `Definition 2.3, (a)` の全単射性から **`τ` 成分が `1`** になる。 -/
theorem psiOTri_of_div_zero (d : ℕ+) (β : OTri P A)
    (h : P.Div (((β : End A)) : A ⟶ A) = 0) : psiOTri P F hτ hA d β = β := by
  set p := (splitEquiv P F hτ hA).symm β with hp
  have hβ : β = uOf P p.1 * tOf P F hτ hA p.2 := by
    rw [hp, ← splitEquiv_apply P F hτ hA]
    exact ((splitEquiv P F hτ hA).apply_symm_apply β).symm
  have hdt : P.Div (((tOf P F hτ hA p.2 : OTri P A) : End A) : A ⟶ A) = 0 := by
    rw [← h, hβ]
    exact (otri_div_unit_mul P p.1 (tOf P F hτ hA p.2)).symm
  have hdt' : P.Div (((p.2 : End A)) : A ⟶ A) = 0 := hdt
  have hone : p.2 = 1 := by
    obtain ⟨t₀, -, huniq⟩ := hτ.charBij A hA (1 : OTri P A)
    have e1 : P.Div (((p.2 : End A)) : A ⟶ A)
        = P.Div ((((1 : OTri P A) : End A)) : A ⟶ A) := by
      rw [hdt']; exact (P.Div_id A).symm
    have e2 : P.Div ((((1 : τ A) : End A)) : A ⟶ A)
        = P.Div ((((1 : OTri P A) : End A)) : A ⟶ A) := rfl
    rw [huniq p.2 e1, huniq 1 e2]
  show uOf P p.1 * (tOf P F hτ hA p.2) ^ ((d : ℕ+) : ℕ) = β
  conv_rhs => rw [hβ]
  rw [hone]
  have h2 : tOf P F hτ hA (1 : τ A) = 1 := rfl
  rw [h2, one_pow, mul_one]

end Psi

/-! ## ★★同型による共役との可換性

★`otriLin`(`Proposition 2.2, (ii)` の射の対応)は、同型 `j` に沿えば
`𝒪^▷` の間のモノイド同型を与える。★**`τ` は部分関手なのでそれで保たれ**
(`IsCharacteristicSplitting.map_mem`)、★**`𝒪^×` も保たれる**
(`otriLin_otimes_mem`)。★したがって分裂が可換になり、`Ψ` も可換になる。
-/

section Conj

variable (F : FrobenioidCore P) {τ : ∀ X : C, Submonoid (End X)}
  (hτ : IsCharacteristicSplitting P F τ)

include hτ in
/-- ★★**`otriLin` は分裂と可換** —— `𝒪^×` も `τ` も保たれるから。 -/
theorem splitEquiv_otriLin {Y Y' : C} (hY : IsIsotropic P Y) (hY' : IsIsotropic P Y')
    {j : Y ⟶ Y'} (hjl : IsLinear P j) (β : OTri P Y') :
    (splitEquiv P F hτ hY).symm (otriLin P F hY hjl β)
      = (⟨((otriLin P F hY hjl (uOf P ((splitEquiv P F hτ hY').symm β).1) : OTri P Y) : End Y),
          otriLin_otimes_mem P F hY hjl _ ((splitEquiv P F hτ hY').symm β).1.2⟩,
         ⟨((otriLin P F hY hjl (tOf P F hτ hY' ((splitEquiv P F hτ hY').symm β).2) :
              OTri P Y) : End Y),
          hτ.map_mem hY hjl _ ((splitEquiv P F hτ hY').symm β).2.2⟩) := by
  refine (splitEquiv P F hτ hY).symm_apply_eq.mpr ?_
  rw [splitEquiv_apply]
  show otriLin P F hY hjl β
    = otriLin P F hY hjl (uOf P _) * otriLin P F hY hjl (tOf P F hτ hY' _)
  rw [← map_mul]
  congr 1
  show β = uOf P ((splitEquiv P F hτ hY').symm β).1 * tOf P F hτ hY' _
  rw [← splitEquiv_apply P F hτ hY']
  exact ((splitEquiv P F hτ hY').apply_symm_apply β).symm

include hτ in
/-- ★★★**`Ψ` は `otriLin` と可換** —— `Ψ(otriLin j β) = otriLin j (Ψ β)`。

★**これが 4 重分解の「対象の取り替え」に対する不変性**である。 -/
theorem psiOTri_otriLin {Y Y' : C} (hY : IsIsotropic P Y) (hY' : IsIsotropic P Y')
    {j : Y ⟶ Y'} (hjl : IsLinear P j) (d : ℕ+) (β : OTri P Y') :
    psiOTri P F hτ hY d (otriLin P F hY hjl β)
      = otriLin P F hY hjl (psiOTri P F hτ hY' d β) := by
  show uOf P (((splitEquiv P F hτ hY).symm (otriLin P F hY hjl β)).1)
      * (tOf P F hτ hY (((splitEquiv P F hτ hY).symm (otriLin P F hY hjl β)).2))
        ^ ((d : ℕ+) : ℕ) = _
  rw [splitEquiv_otriLin P F hτ hY hY' hjl β]
  show (otriLin P F hY hjl (uOf P _) : OTri P Y)
      * (otriLin P F hY hjl (tOf P F hτ hY' _) : OTri P Y) ^ ((d : ℕ+) : ℕ) = _
  rw [← map_pow, ← map_mul]
  rfl

end Conj

/-! ## ★★★`Ψ` の well-defined 性

★4 重分解は一意ではない。★**しかし曖昧性はちょうど 2 種類しかない**——
`arbFactorUniq` と `preStepFactorUniq'` を続けて使うと、
2 つの分解 `δᵢ ≫ γᵢ ≫ βᵢ ≫ αᵢ` の間に同型 `eX : X₁ ≅ X₂`、`eY : Y₁ ≅ Y₂`、
`g : Y₂ ≅ Y₁` が取れて

- `δ₂ = δ₁ ≫ eX.hom`
- `γ₂ ≫ g.hom = eX.inv ≫ γ₁`
- `β₂ = g.hom ≫ β₁ ≫ eY.hom`
- `α₂ = eY.inv ≫ α₁`

となる。★★ここで **`k := g.hom ≫ eY.hom` は `𝒪^×(Y₂)` の元**である
(`β₁`・`β₂` が base-identity だから底が `𝟙` になり、同型だから単元)。
★したがって `β₂ = k · (g による β₁ の共役)` で、**上の自然性 2 本がそのまま効き**

  `Ψ(β₂) = g.hom ≫ Ψ(β₁) ≫ eY.hom`

となる。★★★**あとは代入すると同型がすべて打ち消える。**
-/

section WellDef

variable (F : FrobenioidCore P) {τ : ∀ X : C, Submonoid (End X)}
  (hτ : IsCharacteristicSplitting P F τ) (hiso : IsOfIsotropicType P)

/-- ★**「`ψ` は `φ` の `Ψ` 値である」** —— 4 重分解を 1 つ選んで計算した結果。 -/
def IsPsiValue (d : ℕ+) {A B : C} (φ ψ : A ⟶ B) : Prop :=
  ∃ (X Y : C) (δ : A ⟶ X) (γ : X ⟶ Y) (β : Y ⟶ Y) (α : Y ⟶ B) (hβ : β ∈ OTri P Y),
    φ = δ ≫ γ ≫ β ≫ α ∧ IsFrobeniusType P δ ∧ IsIsometric P γ ∧ IsPreStep P γ ∧
      IsPreStep P β ∧ IsCoAngular P β ∧ IsPullBack P α ∧
      ψ = δ ≫ γ ≫ (((psiOTri P F hτ (hiso Y) d ⟨β, hβ⟩ : OTri P Y) : End Y) : Y ⟶ Y) ≫ α

include hτ hiso in
/-- ★**存在** —— 4 重分解が取れるから。 -/
theorem isPsiValue_exists (hmt : ∀ X : C, IsMetricallyTrivial P X)
    (haa : IsOfAutAmpleType P) (d : ℕ+) {A B : C} (φ : A ⟶ B) :
    ∃ ψ : A ⟶ B, IsPsiValue P F hτ hiso d φ ψ := by
  obtain ⟨X, Y, δ, γ, β, α, hfac, hδ, hγi, hγs, hβb, hβs, hβc, hα⟩ :=
    quadFactor P F hmt haa φ
  exact ⟨_, X, Y, δ, γ, β, α, ⟨hβb, hβs.1⟩, hfac, hδ, hγi, hγs, hβs, hβc, hα, rfl⟩

include hτ hiso in
/-- ★★★**well-defined 性** —— 4 重分解の取り方に依らない。

★**曖昧性が「同型による共役」と「単元倍」に尽きる**ことを、
`arbFactorUniq` ＋ `preStepFactorUniq'` で取り出して使う。 -/
theorem isPsiValue_unique (d : ℕ+) {A B : C} {φ ψ₁ ψ₂ : A ⟶ B}
    (h₁ : IsPsiValue P F hτ hiso d φ ψ₁) (h₂ : IsPsiValue P F hτ hiso d φ ψ₂) :
    ψ₁ = ψ₂ := by
  obtain ⟨X₁, Y₁, δ₁, γ₁, β₁, α₁, hm₁, hf₁, hδ₁, hγi₁, hγs₁, hβs₁, hβc₁, hα₁, he₁⟩ := h₁
  obtain ⟨X₂, Y₂, δ₂, γ₂, β₂, α₂, hm₂, hf₂, hδ₂, hγi₂, hγs₂, hβs₂, hβc₂, hα₂, he₂⟩ := h₂
  -- ★手 1: `arbFactorUniq` —— pre-step 部分をまとめて比較する
  have hps₁ : IsPreStep P (γ₁ ≫ β₁) := IsPreStep.comp P hγs₁ hβs₁
  have hps₂ : IsPreStep P (γ₂ ≫ β₂) := IsPreStep.comp P hγs₂ hβs₂
  have hcomp : δ₁ ≫ (γ₁ ≫ β₁) ≫ α₁ = δ₂ ≫ (γ₂ ≫ β₂) ≫ α₂ := by
    simp only [Category.assoc]
    rw [← hf₁, ← hf₂]
  obtain ⟨eY, eX, hαe, hβe, hδe⟩ :=
    F.arbFactorUniq X₁ Y₁ X₂ Y₂ δ₁ (γ₁ ≫ β₁) α₁ δ₂ (γ₂ ≫ β₂) α₂ hcomp
      hδ₁ hps₁ hα₁ hδ₂ hps₂ hα₂
  -- ★手 2: `preStepFactorUniq'` —— 中間対象 `Y` を比較する
  have hsplit : γ₂ ≫ β₂ = (eX.inv ≫ γ₁) ≫ (β₁ ≫ eY.hom) := by
    rw [hβe]; simp only [Category.assoc]
  obtain ⟨g, hga, hgb⟩ :=
    F.preStepFactorUniq' Y₂ Y₁ γ₂ β₂ (eX.inv ≫ γ₁) (β₁ ≫ eY.hom) hsplit
      hγi₂ hγs₂ hβc₂ hβs₂
      (by simpa using IsIsometric.comp P (isIsometric_of_isIso P eX.inv) hγi₁)
      (IsPreStep.comp P (isPreStep_of_isIso P eX.inv) hγs₁)
      (F.coAngularComp _ eY.hom hβc₁ (isCoAngular_of_isIso P eY.hom))
      (IsPreStep.comp P hβs₁ (isPreStep_of_isIso P eY.hom))
  -- `hga : β₁ ≫ eY.hom = g.inv ≫ β₂`、`hgb : eX.inv ≫ γ₁ = γ₂ ≫ g.hom`
  have hβ₂ : β₂ = g.hom ≫ β₁ ≫ eY.hom := by
    rw [hga, ← Category.assoc, g.hom_inv_id, Category.id_comp]
  -- ★手 3: `k := g.hom ≫ eY.hom` は `𝒪^×(Y₂)` の元
  have hb₁ : P.Base β₁ = P.Base (𝟙 Y₁) := hm₁.1
  have hb₂ : P.Base β₂ = P.Base (𝟙 Y₂) := hm₂.1
  have hkb : P.Base (g.hom ≫ eY.hom) = P.Base (𝟙 Y₂) := by
    rw [hβ₂, P.Base_comp, P.Base_comp, hb₁, P.Base_id, Category.id_comp] at hb₂
    rw [P.Base_comp]; exact hb₂
  have hku : IsUnit (M := End Y₂) (g.hom ≫ eY.hom) := by
    refine isUnit_iff_exists.mpr ⟨(eY.inv ≫ g.inv : End Y₂), ?_, ?_⟩
    · show (eY.inv ≫ g.inv) ≫ (g.hom ≫ eY.hom) = 𝟙 Y₂
      simp
    · show (g.hom ≫ eY.hom) ≫ (eY.inv ≫ g.inv) = 𝟙 Y₂
      simp
  let k : OTimes P Y₂ :=
    ⟨(g.hom ≫ eY.hom : End Y₂),
      ⟨hkb, isLinear_of_isIso P (g.hom ≫ eY.hom)⟩, hku⟩
  -- ★手 4: `β₂ = k · (g による β₁ の共役)`
  have hgl : IsLinear P g.hom := isLinear_of_isIso P g.hom
  have hconj : ((otriLin P F (hiso Y₂) hgl ⟨β₁, hm₁⟩ : OTri P Y₂) : End Y₂)
      = g.hom ≫ β₁ ≫ g.inv := by
    have hs : g.hom ≫ β₁
        = ((otriLin P F (hiso Y₂) hgl ⟨β₁, hm₁⟩ : OTri P Y₂) : End Y₂) ≫ g.hom :=
      otriLin_spec P F (hiso Y₂) hgl (⟨β₁, hm₁⟩ : OTri P Y₁)
    rw [← Category.assoc, hs]
    simp
  have hsplit₂ : (⟨β₂, hm₂⟩ : OTri P Y₂)
      = uOf P k * otriLin P F (hiso Y₂) hgl ⟨β₁, hm₁⟩ := by
    apply Subtype.ext
    show β₂ = ((otriLin P F (hiso Y₂) hgl ⟨β₁, hm₁⟩ : OTri P Y₂) : End Y₂)
      ≫ (g.hom ≫ eY.hom)
    rw [hconj, hβ₂]
    simp
  -- ★手 5: 自然性 2 本で `Ψ(β₂) = g.hom ≫ Ψ(β₁) ≫ eY.hom`
  have hpsi₂ : ((psiOTri P F hτ (hiso Y₂) d ⟨β₂, hm₂⟩ : OTri P Y₂) : End Y₂)
      = g.hom ≫ ((psiOTri P F hτ (hiso Y₁) d ⟨β₁, hm₁⟩ : OTri P Y₁) : End Y₁)
        ≫ eY.hom := by
    have hs : g.hom ≫ ((psiOTri P F hτ (hiso Y₁) d ⟨β₁, hm₁⟩ : OTri P Y₁) : End Y₁)
        = ((otriLin P F (hiso Y₂) hgl (psiOTri P F hτ (hiso Y₁) d ⟨β₁, hm₁⟩) :
            OTri P Y₂) : End Y₂) ≫ g.hom :=
      otriLin_spec P F (hiso Y₂) hgl _
    rw [hsplit₂, psiOTri_unit_mul P F hτ (hiso Y₂) d k,
      psiOTri_otriLin P F hτ (hiso Y₂) (hiso Y₁) hgl d ⟨β₁, hm₁⟩]
    show ((otriLin P F (hiso Y₂) hgl (psiOTri P F hτ (hiso Y₁) d ⟨β₁, hm₁⟩) :
        OTri P Y₂) : End Y₂) ≫ (g.hom ≫ eY.hom) = _
    rw [← Category.assoc, ← hs]
    simp
  -- ★手 6: 代入すると同型がすべて打ち消える
  have hγ₂ : γ₂ = (eX.inv ≫ γ₁) ≫ g.inv := by
    rw [hgb]; simp
  rw [he₁, he₂, hδe, hαe, hpsi₂, hγ₂]
  simp

end WellDef

/-! ## ★★`Ψ` を射の上の写像として取り出す -/

section PsiHom

variable (F : FrobenioidCore P) {τ : ∀ X : C, Submonoid (End X)}
  (hτ : IsCharacteristicSplitting P F τ) (hiso : IsOfIsotropicType P)
  (hmt : ∀ X : C, IsMetricallyTrivial P X) (haa : IsOfAutAmpleType P)
  (hfn : IsOfFrobeniusNormalizedType P)

include hτ hiso hmt haa in
/-- ★★★**`Ψ` は射の上でちょうど一つの値をとる**。

★これで `Ψ : Hom(A,B) → Hom(A,B)` が**関数として定まった**。 -/
theorem psiHom_existsUnique (d : ℕ+) {A B : C} (φ : A ⟶ B) :
    ∃! ψ : A ⟶ B, IsPsiValue P F hτ hiso d φ ψ := by
  obtain ⟨ψ, hψ⟩ := isPsiValue_exists P F hτ hiso hmt haa d φ
  exact ⟨ψ, hψ, fun ψ' hψ' => isPsiValue_unique P F hτ hiso d hψ' hψ⟩

/-- ★**`Ψ` の射の上の値**。 -/
noncomputable def psiMap (d : ℕ+) {A B : C} (φ : A ⟶ B) : A ⟶ B :=
  (psiHom_existsUnique P F hτ hiso hmt haa d φ).choose

include hτ hiso hmt haa in
theorem psiMap_spec (d : ℕ+) {A B : C} (φ : A ⟶ B) :
    IsPsiValue P F hτ hiso d φ (psiMap P F hτ hiso hmt haa d φ) :=
  (psiHom_existsUnique P F hτ hiso hmt haa d φ).choose_spec.1

include hτ hiso hmt haa in
/-- ★**分解から `Ψ` の値が読める**。 -/
theorem psiMap_eq (d : ℕ+) {A B : C} {φ ψ : A ⟶ B}
    (h : IsPsiValue P F hτ hiso d φ ψ) : psiMap P F hτ hiso hmt haa d φ = ψ :=
  isPsiValue_unique P F hτ hiso d (psiMap_spec P F hτ hiso hmt haa d φ) h

/-! ### ★原文 (a)(b) —— `Ψ` の基本性質 -/

include hτ hiso hmt haa in
/-- ★★**`Ψ` は `Base` を変えない**。 -/
theorem psiMap_base (d : ℕ+) {A B : C} (φ : A ⟶ B) :
    P.Base (psiMap P F hτ hiso hmt haa d φ) = P.Base φ := by
  obtain ⟨X, Y, δ, γ, β, α, hm, hf, hδ, hγi, hγs, hβs, hβc, hα, he⟩ :=
    psiMap_spec P F hτ hiso hmt haa d φ
  have hb : P.Base (((psiOTri P F hτ (hiso Y) d ⟨β, hm⟩ : OTri P Y) : End Y) : Y ⟶ Y)
      = P.Base β := by
    rw [show P.Base β = P.Base (𝟙 Y) from hm.1]
    exact (psiOTri P F hτ (hiso Y) d ⟨β, hm⟩).2.1
  rw [he, hf, P.Base_comp, P.Base_comp, P.Base_comp, P.Base_comp, P.Base_comp, P.Base_comp, hb]

include hτ hiso hmt haa in
/-- ★★**`Ψ` は Frobenius 次数を変えない**。 -/
theorem psiMap_degFr (d : ℕ+) {A B : C} (φ : A ⟶ B) :
    P.degFr (psiMap P F hτ hiso hmt haa d φ) = P.degFr φ := by
  obtain ⟨X, Y, δ, γ, β, α, hm, hf, hδ, hγi, hγs, hβs, hβc, hα, he⟩ :=
    psiMap_spec P F hτ hiso hmt haa d φ
  have hb : P.degFr (((psiOTri P F hτ (hiso Y) d ⟨β, hm⟩ : OTri P Y) : End Y) : Y ⟶ Y)
      = P.degFr β := by
    rw [show P.degFr β = 1 from hm.2]
    exact (psiOTri P F hτ (hiso Y) d ⟨β, hm⟩).2.2
  rw [he, hf]
  simp only [P.degFr_comp, hb]

include hτ hiso hmt haa hfn in
/-- ★★★**`Ψ` は `Div` を `d` 倍する** —— 原文 (b)「`d` の Frobenius 函手と 1-compatible」。 -/
theorem psiMap_div (d : ℕ+) {A B : C} (φ : A ⟶ B) :
    P.Div (psiMap P F hτ hiso hmt haa d φ) = ((d : ℕ+) : ℕ) • P.Div φ := by
  obtain ⟨X, Y, δ, γ, β, α, hm, hf, hδ, hγi, hγs, hβs, hβc, hα, he⟩ :=
    psiMap_spec P F hτ hiso hmt haa d φ
  have hbi : IsPreStep P (((psiOTri P F hτ (hiso Y) d ⟨β, hm⟩ : OTri P Y) : End Y) : Y ⟶ Y) := by
    refine ⟨(psiOTri P F hτ (hiso Y) d ⟨β, hm⟩).2.2, ?_⟩
    show IsIso (P.Base _)
    rw [show P.Base (((psiOTri P F hτ (hiso Y) d ⟨β, hm⟩ : OTri P Y) : End Y) : Y ⟶ Y)
      = P.Base (𝟙 Y) from (psiOTri P F hτ (hiso Y) d ⟨β, hm⟩).2.1, P.Base_id]
    infer_instance
  have hψ := quadFactor_div P F hδ hγi hγs
    (psiOTri P F hτ (hiso Y) d ⟨β, hm⟩).2.1 hbi hα
  rw [he, hψ, psiOTri_div P F hτ (hiso Y) (hfn Y) d ⟨β, hm⟩, map_nsmul, hf,
    quadFactor_div P F hδ hγi hγs hm.1 hβs hα]

include hτ hiso hmt haa in
/-- ★★★**等長射の上で `Ψ` は恒等** —— 原文 (a)。

★★`Φ.map` は **`Definition 1.1, (a)` の characteristic 単射性から常に単射**なので、
`Div φ = 0` から直ちに `Div β = 0` が出る。 -/
theorem psiMap_of_isometric (d : ℕ+) {A B : C} (φ : A ⟶ B) (hφ : IsIsometric P φ) :
    psiMap P F hτ hiso hmt haa d φ = φ := by
  obtain ⟨X, Y, δ, γ, β, α, hm, hf, hδ, hγi, hγs, hβs, hβc, hα, he⟩ :=
    psiMap_spec P F hτ hiso hmt haa d φ
  have hd0 : P.Div β = 0 := by
    refine Φ.map_injective (P.Base (δ ≫ γ)) ?_
    rw [map_zero, ← quadFactor_div P F hδ hγi hγs hm.1 hβs hα, ← hf]
    exact hφ
  rw [he, psiOTri_of_div_zero P F hτ (hiso Y) d ⟨β, hm⟩ hd0, ← hf]

include hτ hiso hmt haa in
/-- ★**`Ψ(𝟙) = 𝟙`** —— 恒等射は等長だから。 -/
theorem psiMap_id (d : ℕ+) (A : C) :
    psiMap P F hτ hiso hmt haa d (𝟙 A) = 𝟙 A :=
  psiMap_of_isometric P F hτ hiso hmt haa d (𝟙 A) (P.Div_id A)

/-! ### ★合成則のうち、分解を組み替えずに出る 2 つ

★4 重分解 `δ ≫ γ ≫ β ≫ α` の**両端**——左の Frobenius 型と右の pull-back ——は
合成しても形が変わらない。★したがってこの 2 方向の合成則は**ただちに出る**。
残るのは真ん中(等長 pre-step と base-identity 自己射)を通す場合である。
-/

include hτ hiso hmt haa in
/-- ★★**右から pull-back を掛ける場合** —— `Ψ(φ ≫ α') = Ψ(φ) ≫ α'`。

★pull-back どうしの合成は pull-back なので、分解の右端が伸びるだけ。 -/
theorem psiMap_comp_pullBack (d : ℕ+) {A B E : C} (φ : A ⟶ B) {α' : B ⟶ E}
    (hα' : IsPullBack P α') :
    psiMap P F hτ hiso hmt haa d (φ ≫ α')
      = psiMap P F hτ hiso hmt haa d φ ≫ α' := by
  obtain ⟨X, Y, δ, γ, β, α, hm, hf, hδ, hγi, hγs, hβs, hβc, hα, he⟩ :=
    psiMap_spec P F hτ hiso hmt haa d φ
  refine psiMap_eq P F hτ hiso hmt haa d ?_
  refine ⟨X, Y, δ, γ, β, α ≫ α', hm, ?_, hδ, hγi, hγs, hβs, hβc,
    IsPullBack.comp P hα hα', ?_⟩
  · rw [hf]; simp
  · rw [he]; simp

include hτ hiso hmt haa in
/-- ★★**左から Frobenius 型射を掛ける場合** —— `Ψ(δ' ≫ φ) = δ' ≫ Ψ(φ)`。

★Frobenius 型射どうしの合成は Frobenius 型なので、分解の左端が伸びるだけ。 -/
theorem psiMap_frob_comp (d : ℕ+) {A A' B : C} {δ' : A' ⟶ A}
    (hδ' : IsFrobeniusType P δ') (φ : A ⟶ B) :
    psiMap P F hτ hiso hmt haa d (δ' ≫ φ)
      = δ' ≫ psiMap P F hτ hiso hmt haa d φ := by
  obtain ⟨X, Y, δ, γ, β, α, hm, hf, hδ, hγi, hγs, hβs, hβc, hα, he⟩ :=
    psiMap_spec P F hτ hiso hmt haa d φ
  refine psiMap_eq P F hτ hiso hmt haa d ?_
  refine ⟨X, Y, δ' ≫ δ, γ, β, α, hm, ?_, IsFrobeniusType.comp P F hδ' hδ,
    hγi, hγs, hβs, hβc, hα, ?_⟩
  · rw [hf]; simp
  · rw [he]; simp

include hτ hiso hmt haa hfn in
/-- ★★★**右から `𝒪^▷` の元を掛ける場合** —— `Ψ(φ ≫ β') = Ψ(φ) ≫ Ψ(β')`。

★★**pull-back `α` は(isotropic 型のもとで)co-angular かつ linear なので
`otriLin` が使える** —— `α ≫ β' = otriLin(β') ≫ α`。したがって `β'` は
分解の中の `β` に**掛け算として吸収される**。

★★★あとは
- `psiOTri` がモノイド準同型であること(`psiOTri_mul`)と
- `Ψ` が `otriLin` と可換であること(`psiOTri_otriLin`)

の 2 本で閉じる。★**先に用意した自然性がそのまま効いた。** -/
theorem psiMap_comp_otri (d : ℕ+) {A B : C} (φ : A ⟶ B) (β' : OTri P B) :
    psiMap P F hτ hiso hmt haa d (φ ≫ (((β' : End B)) : B ⟶ B))
      = psiMap P F hτ hiso hmt haa d φ
        ≫ (((psiOTri P F hτ (hiso B) d β' : OTri P B) : End B) : B ⟶ B) := by
  obtain ⟨X, Y, δ, γ, β, α, hm, hf, hδ, hγi, hγs, hβs, hβc, hα, he⟩ :=
    psiMap_spec P F hτ hiso hmt haa d φ
  have hαlin : IsLinear P α := (F.pullBackLB α hα).2
  have hspec : α ≫ (((β' : End B)) : B ⟶ B)
      = (((otriLin P F (hiso Y) hαlin β' : OTri P Y) : End Y) : Y ⟶ Y) ≫ α :=
    otriLin_spec P F (hiso Y) hαlin β'
  have hpsispec : α ≫ (((psiOTri P F hτ (hiso B) d β' : OTri P B) : End B) : B ⟶ B)
      = (((otriLin P F (hiso Y) hαlin (psiOTri P F hτ (hiso B) d β') :
          OTri P Y) : End Y) : Y ⟶ Y) ≫ α :=
    otriLin_spec P F (hiso Y) hαlin _
  set s : OTri P Y := otriLin P F (hiso Y) hαlin β' with hsdef
  set b : OTri P Y := ⟨β, hm⟩ with hbdef
  have hnew : (((s * b : OTri P Y) : End Y) : Y ⟶ Y) = β ≫ (((s : End Y)) : Y ⟶ Y) := rfl
  refine psiMap_eq P F hτ hiso hmt haa d ?_
  refine ⟨X, Y, δ, γ, (((s * b : OTri P Y) : End Y) : Y ⟶ Y), α, (s * b).2, ?_, hδ,
    hγi, hγs, ?_, ?_, hα, ?_⟩
  · -- 分解が `φ ≫ β'` に戻る
    rw [hf, hnew]
    simp only [Category.assoc]
    rw [← hspec]
  · -- pre-step
    refine ⟨(s * b).2.2, ?_⟩
    show IsIso (P.Base _)
    rw [show P.Base ((((s * b : OTri P Y) : End Y) : Y ⟶ Y)) = P.Base (𝟙 Y) from (s * b).2.1,
      P.Base_id]
    infer_instance
  · -- co-angular —— isotropic 型なので `Proposition 1.4, (i)`
    rw [hnew]
    exact F.coAngularComp _ _ hβc (prop_1_4_i P _ (fun X _ => hiso X))
  · -- `Ψ` の値
    rw [he, psiOTri_mul P F hτ (hiso Y) (hfn Y) d s b, hsdef,
      psiOTri_otriLin P F hτ (hiso Y) (hiso B) hαlin d β']
    show _ ≫ _ = δ ≫ γ ≫ ((((psiOTri P F hτ (hiso Y) d b : OTri P Y) : End Y) : Y ⟶ Y)
      ≫ (((otriLin P F (hiso Y) hαlin (psiOTri P F hτ (hiso B) d β') :
        OTri P Y) : End Y) : Y ⟶ Y)) ≫ α
    simp only [Category.assoc]
    rw [← hpsispec]

end PsiHom

/-! ## ★★★残り —— 関手性と圏同値(段取り)

原文 (FrdI p.50):
> for d′ ∈N≥1 — it follows that the assignment φ →Ψ(φ) is compatible with com-

★**ここまでで揃ったもの**: `psiMap` が**関数として**定まり(well-defined)、
`Base`・`degFr` を保ち、`Div` を `d` 倍し、等長射の上で恒等である。
★残るのは**合成との両立**と、そこから出る `𝒞 ≃ 𝒞(d)` である。

### ★段取り(原文の筋を我々の道具に翻訳したもの)

| 段 | 内容 | 使う道具 | 状態 |
|---|---|---|---|
| 1 | pull-back を右へ寄せる | `plBk_shuffle` | ✅ |
| 2 | 底が同型な pull-back は同型 | `isIso_of_isPullBack_of_baseIso` | ✅ |
| 3 | `𝒪^▷` を Frobenius 型射の右へ通す | `otri_comp_frobBaseId` | ✅ |
| 4 | `𝒪^▷` を等長 pre-step の右へ通す | `otriFwd` / `otriLin` | ✅(既存) |
| 5 | `Ψ(β^k) = Ψ(β)^k` | `psiOTri_pow` | ✅ |
| 6 | 合成の 4 重分解を組み立てる | 上の 1–4 | ★**未** |
| 7 | `Ψ(φ ≫ φ') = Ψφ ≫ Ψφ'`(`𝒞^istr` 上) | 6 ＋ 5 ＋ 自然性 2 本 | ★**未** |
| 8 | isotropic hull が単射 ⟹ 一般の `𝒞` へ | `Definition 1.3, (v), (a)` | ★**未** |
| 9 | 関手 `𝒞 ⥤ 𝒞(d)`(`Cd` は `Definition 2.4` にある) | 7,8 | ★**未** |
| 10 | 本質的全射・忠実・充満 | `prop_2_5_i_bijective` | ★**未** |

### ★段 6 の具体的な形(導出済み、実装待ち)

`φ = δ≫γ≫β≫α`、`φ' = δ'≫γ'≫β'≫α'` とすると

  `φ≫φ' = δ≫γ≫β≫(α ≫ δ'≫γ'≫β')≫α'`

で、`ρ := δ'≫γ'≫β'` に段 1 を当てると `α≫ρ = ρ̃≫α̃`(`α̃` は pull-back、
`Base ρ̃` は同型)。★`α̃≫α'` は pull-back なので

  `φ≫φ' = δ≫γ≫(β≫ρ̃)≫(α̃≫α')`

となり、あとは `ρ̃` を `δ₂≫γ₂≫β₂` に分解して(段 2 でその pull-back 部分が
同型になる)、`β` を段 3・段 4 で `δ₂`・`γ₂` の右へ通せばよい。

★★**段 3 が `𝒞^istr` を要求する**(Frobenius 型射を base-identity に取るため)
——これが原文が `𝒞^istr` で議論し、段 8 で一般へ移す理由である。

### ★★我々の実装が原文より強い点(記録)

原文は特性分裂を「**`A` が isotropic でなくても使える**」と明言し、根拠に
`Definition 2.3, (b)`(我々の `hullMem`)を挙げている。★我々の
`charSplitting_bijective` は `IsIsotropic` を要求するので、`psiMap` は現在
`IsOfIsotropicType` を仮定していて**原文より強い**。`hullMem` で
非 isotropic へ延ばせば外せる。★段 8 と同時に解消するのが自然である。
-/

end ABC3.Found.FrdI
