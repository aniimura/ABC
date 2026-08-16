import ABC3.Found.FrdI.Prop25
import ABC3.Found.FrdI.Def24

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
        IsBaseIdentity P β ∧ IsPreStep P β ∧ IsPullBack P α := by
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
    ?_, hδ, hγi, hγs, ?_, ?_, ?_⟩
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

/-! ### ★★`Ψ` の自然性 —— well-defined 性の核

★4 重分解には**選び方の自由**がある(原文の
`(α∘ε, ε⁻¹∘β∘ζ, ζ⁻¹∘γ∘θ, θ⁻¹∘δ)`)。★**その自由度は 2 種類**である:

1. **単元倍** —— `β` が `𝒪^×` の元だけずれる
2. **同型による共役** —— `β` が別の対象へ移る

★**`Ψ` はどちらにも可換である**ことをここで示す。
-/

include hfn in
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

end ABC3.Found.FrdI
