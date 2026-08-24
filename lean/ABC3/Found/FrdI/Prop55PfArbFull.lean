/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55PfArbStd

/-!
# [FrdI] Proposition 5.5, (iv) の単系の同定 —— `𝒞^birat` 側の仮定をすべて外す

原文 (FrdI p.105):
> Finally, assertion (iv) is immediate from the definitions [cf. also assertions (i), (ii);

★★`Prop55PfRatFn.lean` / `Prop55PfArbStd.lean` の単系の同定
`(Φ^pf)^birat = ℚ·Φ^birat` は、`𝒞^birat` 側の性質を**仮引数で受けていた**:

| 仮引数 | 内容 |
|---|---|
| `hisoB` | `𝒞^birat` は isotropic |
| `hfnBirat` | `𝒞^birat` は Frobenius-normalized |
| `hAB` / `ζ` / `hdeg` / `hprop` | `A^birat` は Frobenius-trivial |
| `hfnB'` | 根 `(A^birat)^{(1)}` も Frobenius-normalized |
| `hfnPfBirat` | `(𝒞^pf)^birat` も Frobenius-normalized |

★本ファイルは **`𝒞` が isotropic ＋ unit-trivial 型**であれば
これらが**すべて在庫から出る**ことを示す。

## ★鍵は 2 本

1. **`𝒞^birat` の零因子の単系は自明**(`trivialOn`)なので、
   そこでは**すべての射が isometric** である。★したがって
   「co-angular pre-step」は「isometric pre-step」であり、
   isotropic 性がそれを同型にする —— つまり `𝒞^birat` は
   **metrically trivial 型**であり、`Proposition 2.5, (ii)`
   (`isFrobeniusTrivial_of_isotropic`)で**すべての対象が Frobenius-trivial** になる。
2. **unit-trivial 型なら `𝒞^birat` は Frobenius-normalized**
   (在庫 `birat_frobNormalized_of_unitTrivial`)。
   ★`𝒞^pf` も unit-trivial(`pfRoot_isOfUnitTrivialType_of_arb`)なので
   `(𝒞^pf)^birat` にも同じことが当たる。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section BiratTypes

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-- ★★**`𝒞^birat` ではすべての射が isometric**(零因子の単系が `PUnit` だから)。 -/
theorem birat_isIsometric_all (G : Frobenioid P) {X Y : BiratCat P G} (f : X ⟶ Y) :
    IsIsometric (biratPre P G) f :=
  Subsingleton.elim (α := PUnit.{w + 1}) _ _

/-- ★★★**`𝒞^birat` は metrically trivial 型** ——
co-angular pre-step は自動的に isometric pre-step であり、isotropic 性で同型になる。 -/
theorem birat_isMetricallyTrivial (G : Frobenioid P)
    (hisoB : ∀ X : BiratCat P G, IsIsotropic (biratPre P G) X) (A : BiratCat P G) :
    IsMetricallyTrivial (biratPre P G) A := by
  intro Dd φ _ hs
  haveI : IsIso φ := hisoB A Dd φ (birat_isIsometric_all G φ) hs
  exact ⟨(asIso φ).symm⟩

/-- ★★★★**`𝒞^birat` のすべての対象は Frobenius-trivial**。

★`Proposition 2.5, (ii)`(`isFrobeniusTrivial_of_isotropic`)を
`𝒞^birat` に当てるだけである —— 仮定 `hmt` は上で出た。 -/
theorem birat_isFrobTrivial_of_isotropic (G : Frobenioid P) (F' : FrobenioidCore (biratPre P G))
    (hisoB : ∀ X : BiratCat P G, IsIsotropic (biratPre P G) X) (A : BiratCat P G) :
    IsFrobeniusTrivial (biratPre P G) A :=
  isFrobeniusTrivial_of_isotropic (biratPre P G) F'
    (birat_isMetricallyTrivial G hisoB) A (hisoB A)

end BiratTypes

/-! ## ★2. `Proposition 5.5, (iv)` の単系の同定 —— 条なしの版 -/

section Full

variable {D : Type u} [Category.{v} D] [IsConnected D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★**[FrdI] Proposition 5.5, (iv) の単系の同定**(★`𝒞^birat` 側の仮定なし)——
`𝒞` が Frobenius-isotropic ＋ isotropic ＋ Frobenius-normalized ＋ unit-trivial 型なら

  `(Φ^pf)^birat = ℚ·Φ^birat`

が `𝒟` 上の部分群の等号として成り立つ。

★`𝒞^birat` 側の 8 つの仮引数はすべて `birat_isFrobTrivial_of_isotropic` と
`birat_frobNormalized_of_unitTrivial` から出る。 -/
theorem phiBiratOn_pf_eq_qPhiBiratOn_full (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hfnC : ∀ X : C, IsFrobeniusNormalized P X)
    (hut : ∀ X : C, IsUnitTrivial P X)
    (G : Frobenioid P) (Gpf : Frobenioid (pfRootPre P F))
    (d : D) :
    phiBiratOn Gpf d = qPhiBiratOn P G d := by
  -- `𝒞^birat` の型
  have hisoB : ∀ X : BiratCat P G, IsIsotropic (biratPre P G) X :=
    birat_isOfIsotropicType hiso
  have hfnBirat : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X :=
    birat_frobNormalized_of_unitTrivial G hiso hut
  set F' : FrobenioidCore (biratPre P G) := biratCoreOf G hfnBirat with hF'
  have hAB : ∀ A : C, IsFrobeniusTrivial (biratPre P G) (biratUp P G A) :=
    fun A => birat_isFrobTrivial_of_isotropic G F' hisoB _
  -- `(𝒞^pf)^birat` の Frobenius-normalized 性(`𝒞^pf` は unit-trivial)
  have hfnPfBirat : ∀ X : BiratCat (pfRootPre P F) Gpf,
      IsFrobeniusNormalized (biratPre (pfRootPre P F) Gpf) X :=
    birat_frobNormalized_of_unitTrivial Gpf (pfRoot_isOfIsotropicType (F := F) hfi)
      (pfRoot_isOfUnitTrivialType_of_arb hfi hiso hfnC hut Gpf)
  exact phiBiratOn_pf_eq_qPhiBiratOn_all hfi hiso Gpf F' hisoB hfnBirat hAB
    (fun A => hfnBirat _) (fun A => (hAB A).choose) (fun A n => (hAB A).choose_spec.1 n)
    (fun A n => (hAB A).choose_spec.2 n) hfnPfBirat d

/-- ★★★★★★★**[FrdI] Proposition 5.5, (iv)**(★条なし)——
`𝒞^pf` の有理関数の単系は各 `d ∈ Ob(𝒟)` で **`ℚ·Φ^birat(d)`** である。 -/
theorem pfRoot_ratFnData_bmon_val_full (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hfnC : ∀ X : C, IsFrobeniusNormalized P X)
    (hut : ∀ X : C, IsUnitTrivial P X)
    (G : Frobenioid P) (Gpf : Frobenioid (pfRootPre P F)) (hfsmD : IsOfFSMType D)
    (d : D) :
    (pfRoot_ratFnData_of_arb (F := F) hfi hiso hfnC hut Gpf hfsmD).bmon.val d
      = ↥(qPhiBiratOn P G d) := by
  show ↥(phiBiratOn Gpf d) = ↥(qPhiBiratOn P G d)
  rw [phiBiratOn_pf_eq_qPhiBiratOn_full hfi hiso hfnC hut G Gpf d]
  rfl

set_option maxHeartbeats 800000 in
/-- ★★★★★★★**`Proposition 5.3` の図式の右の縦の矢印**(★条なし)`𝒞^pf ⟶ 𝒞^rlf`。 -/
noncomputable def pfRootToRlfFunctor_full (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hfnC : ∀ X : C, IsFrobeniusNormalized P X)
    (hut : ∀ X : C, IsUnitTrivial P X)
    (G : Frobenioid P) (Gpf : Frobenioid (pfRootPre P F)) (hfsmD : IsOfFSMType D)
    (hcharInj' : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := NNReal) (Φ.map α)))
    (hint' : ∀ A : D, IsIntegralMonoid (RlfT (Φ.val A))) :
    ModelData.Obj (pfRoot_ratFnData_of_arb (F := F) hfi hiso hfnC hut Gpf hfsmD).model
      ⥤ Crlf G hiso (birat_frobNormalized_of_unitTrivial G hiso hut) hcharInj' hint' hfsmD :=
  pfRootToRlfFunctor_of_arb hfi hiso hfnC hut Gpf hfsmD
    (birat_frobNormalized_of_unitTrivial G hiso hut) hcharInj' hint'
    (fun d => phiBiratOn_pf_eq_qPhiBiratOn_full hfi hiso hfnC hut G Gpf d)

end Full

/-! ### ★出典の紐付け -/

/-- ★★★★★★★locator —— `Proposition 5.5, (iv)` の単系の同定(条なし)。 -/
def pfRoot_ratFnData_bmon_val_full.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iv) — 𝒞^pf の有理関数の単系は ℚ·Φ^birat(条なし)",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
