/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55Pf
import ABC3.Found.FrdI.Cor411Otimes

/-!
# [FrdI] `𝒞^pf` は unit-trivial 性を保つ

原文 (FrdI p.103):
> tion monoid Φbirat (respectively, Q · Φbirat = Φbirat ⊗Z Q = (Φbirat)pf). In particular,

★★`Proposition 5.3` は「`𝒞^un-tr` と `(𝒞^un-tr)^pf` は model 型」と言う。
その根拠は `Theorem 5.1, (iv)`(unit-trivial 型はつねに model 型)なので、
要るのは **`𝒞^pf` が unit-trivial 性を保つこと**である。

★★本ファイルはそれを **`Proposition 5.5, (i)` の同型から**出す:

* `isUnit_of_pf_isAddUnit` —— `M^pf` の元が可逆なら分子が可逆
  (`m/a + n/b = 0` から `(m^b·n^a)^k = 1`、そこから `m` が可逆)
* ★`otimes_pfRoot_eq_bot` —— `𝒪^▷(A)^pf ≅ 𝒪^▷(A^pf)`(`otriPfEquiv`)を
  **単元に制限する**だけ。分子は `𝒪^×(A) = ⊥` から `1` になる。

★これは `Proposition 5.5, (ii)` の `(𝒞^pf)^un-tr ≃ (𝒞^un-tr)^pf` の材料でもある。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

/-! ## ★1. 単系の一般論 -/

/-- ★★**`M^pf` の元が可逆なら、分子も可逆**(乗法の言葉で)。

★`m/a + n/b = 0` は `∃k, (m^b·n^a)^k = 1` を意味する。`k ≥ 1` なので
`m^b·n^a` は可逆、可換だから `m^b` も可逆、同じ理屈で `m` も可逆。 -/
theorem isUnit_of_pf_isAddUnit {N : Type w} [CommMonoid N] {n : N} {a : ℕ+}
    (h : IsAddUnit (@Pf.mk (Additive N) _ (Additive.ofMul n) a)) : IsUnit n := by
  obtain ⟨y, hy⟩ := isAddUnit_iff_exists_neg.mp h
  induction y using Pf.inductionOn with
  | _ m b =>
    have h0 : @Pf.mk (Additive N) _ ((b : ℕ) • Additive.ofMul n + (a : ℕ) • m) (a * b)
        = (0 : Pf (Additive N)) := hy
    obtain ⟨k, hk⟩ := Quotient.exact h0
    have hk' : ((k : ℕ) * 1) • ((b : ℕ) • Additive.ofMul n + (a : ℕ) • m)
        = ((k : ℕ) * ((a * b : ℕ+) : ℕ)) • (0 : Additive N) := hk
    rw [mul_one, smul_zero] at hk'
    have hz : (Additive.toMul ((b : ℕ) • Additive.ofMul n + (a : ℕ) • m)) ^ ((k : ℕ)) = 1 := hk'
    have hu : IsUnit ((Additive.toMul ((b : ℕ) • Additive.ofMul n + (a : ℕ) • m))) := by
      refine IsUnit.of_mul_eq_one
        ((Additive.toMul ((b : ℕ) • Additive.ofMul n + (a : ℕ) • m)) ^ ((k : ℕ) - 1)) ?_
      rw [← pow_succ', show ((k : ℕ) - 1 + 1) = (k : ℕ) from Nat.sub_add_cancel k.2]
      exact hz
    have hnb : IsUnit (n ^ ((b : ℕ))) :=
      isUnit_of_mul_isUnit_left (y := (Additive.toMul ((a : ℕ) • m))) hu
    refine isUnit_of_mul_isUnit_left (y := n ^ ((b : ℕ) - 1)) ?_
    rw [← pow_succ', show ((b : ℕ) - 1 + 1) = (b : ℕ) from Nat.sub_add_cancel b.2]
    exact hnb

/-! ## ★2. `𝒞^pf` は unit-trivial 性を保つ -/

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

/-- ★★★★★**完全化は unit-trivial 性を保つ**(`A` が Frobenius-trivial かつ
Frobenius-normalized のとき)。

★`Proposition 5.5, (i)` の同型 `𝒪^▷(A)^pf ≅ 𝒪^▷(A^pf)` を単元に制限するだけ。
★★これで `Theorem 5.1, (iv)` が `(𝒞^un-tr)^pf` に当たるようになる —— それが
`Proposition 5.3` の「`(𝒞^un-tr)^pf` は model 型」の根拠である。 -/
theorem otimes_pfRoot_eq_bot (hiso : ∀ X : C, IsIsotropic P X) {A : C}
    (hA : IsFrobeniusTrivial P A)
    (hfn : IsFrobeniusNormalized P A) (hfn' : IsFrobeniusNormalized P (rtObj P F A 1))
    (ζ : ℕ+ →* End A) (hdeg : ∀ n : ℕ+, P.degFr ((ζ n : End A) : A ⟶ A) = n)
    (hprop : ∀ n : ℕ+, IsBaseIdentity P (ζ n) ∧ IsFrobeniusType P ((ζ n : End A) : A ⟶ A))
    (hut : IsUnitTrivial P A) :
    IsUnitTrivial (pfRootPre P F) (⟨A, 1⟩ : PfRootObj P F) := by
  refine le_antisymm (fun f hf => ?_) bot_le
  rw [Submonoid.mem_bot]
  set E := otriPfEquiv (F := F) hiso hA hfn hfn' ζ hdeg hprop with hE
  set g : OTri (pfRootPre P F) (⟨A, 1⟩ : PfRootObj P F) := ⟨f, hf.1⟩ with hg
  have hgu : IsUnit g := isUnit_otri_of_otimes _ hf
  obtain ⟨x, hx⟩ := otriPfMap_surjective (F := F) hiso hA hfn hfn' ζ hdeg hprop g
  have hEx : E x = Additive.ofMul g := hx
  have hxu : IsAddUnit x := by
    obtain ⟨w, hw⟩ := hgu.exists_right_inv
    refine isAddUnit_iff_exists_neg.mpr ⟨E.symm (Additive.ofMul w), ?_⟩
    refine E.injective ?_
    rw [map_add, hEx, E.apply_symm_apply, map_zero]
    show Additive.ofMul (g * w) = 0
    rw [hw]
    rfl
  obtain ⟨⟨α, a⟩, hrep⟩ := Quotient.exists_rep x
  have hxm : x = @Pf.mk (Additive (OTri P A)) (otriAddCommMonoid hfn) α a := hrep.symm
  have hαu : IsUnit (Additive.toMul α) :=
    @isUnit_of_pf_isAddUnit (OTri P A) (otriCommMonoid hfn) (Additive.toMul α) a (hxm ▸ hxu)
  have hmem : ((Additive.toMul α : OTri P A) : End A) ∈ OTimes P A :=
    ⟨(Additive.toMul α : OTri P A).2, IsUnit.map (OTri P A).subtype hαu⟩
  have hα1 : (Additive.toMul α : OTri P A) = 1 :=
    Subtype.ext (Submonoid.mem_bot.mp (hut ▸ hmem))
  have hzero : x = 0 := by
    rw [hxm]
    refine @Pf.sound (Additive (OTri P A)) (otriAddCommMonoid hfn) α 0 a 1 1 ?_
    rw [show α = 0 from hα1]
    simp
  have hg1 : g = 1 := by
    have h0 : E x = E 0 := by rw [hzero]
    rw [hEx, map_zero] at h0
    exact congrArg Additive.toMul h0
  exact congrArg Subtype.val hg1

/-! ### ★出典の紐付け -/

/-- ★locator —— `Proposition 5.3` が使う「完全化は unit-trivial 性を保つ」。 -/
def otimes_pfRoot_eq_bot.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 103,
    item := "Proposition 5.3 — (𝒞^un-tr)^pf も unit-trivial(ゆえに model 型)",
    sectionId := "frdi-prop-5-3" }

end ABC3.Found.FrdI
