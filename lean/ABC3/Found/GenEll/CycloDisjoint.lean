/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.NumberTheory.NumberField.Cyclotomic.Basic
import Mathlib.NumberTheory.NumberField.Discriminant.Different
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★`L` と `ℚ(ζ_{l^n})` の線型無関連（`Found`）

`Skeleton/GenEll/GaloisLocal.lean` の **葉 4**（円分指標の全射性）の核。

## ★原文が畳んでいるもの

`GenEll` p.20 の証明は `det ρ = χ_cyc` が全射であることを
「`ℚ(ζ_{l^∞})/ℚ` は `l` で完全分岐するので `L/ℚ` と線型無関連」で済ませている。

## ★★★★★★★★2026-08-31 の道（第 780）

**「完全分岐」を経由しない**。判別式だけで済む:

* `disc ℚ(ζ_{l^k}) = ±l^m` —— `IsCyclotomicExtension.Rat.discr_prime_pow`
* `l ∤ disc L`（`l` が `L` で不分岐）
* ゆえに `IsCoprime (disc ℚ(ζ)) (disc L)`
* `NumberField.linearDisjoint_of_isGalois_isCoprime_discr` が線型無関連を出す

★★どれも mathlib にある。
-/

open NumberField Polynomial

namespace ABC3.Found.GenEll

open ABC3.Meta

variable {M : Type*} [Field M] [NumberField M]

/-- ★単元は何とでも互いに素。 -/
theorem isCoprime_of_isUnit_left {R : Type*} [CommSemiring R] {u b : R} (hu : IsUnit u) :
    IsCoprime u b := by
  obtain ⟨v, hv⟩ := hu.exists_left_inv
  exact ⟨v, 0, by simp [hv]⟩

/-- ★★★★★★★★**`disc ℚ(ζ_{l^k})` は `l` 以外の素因数を持たない**。

`IsCyclotomicExtension.Rat.discr_prime_pow` は
`disc = (-1)^(φ(l^k)/2) · l^(l^(k-1)·((l-1)k-1))` を与える。
★`(-1)^t` は単元なので、`l ∤ d` なら `disc` と `d` は互いに素。 -/
theorem discr_cyclo_isCoprime (l k : ℕ) [Fact l.Prime] (A : Type*) [Field A] [CharZero A]
    [IsCyclotomicExtension {l ^ k} ℚ A] (d : ℤ) (hd : ¬ (l : ℤ) ∣ d) :
    haveI : NumberField A := IsCyclotomicExtension.numberField {l ^ k} ℚ A
    IsCoprime (NumberField.discr A) d := by
  haveI : NumberField A := IsCyclotomicExtension.numberField {l ^ k} ℚ A
  rw [IsCyclotomicExtension.Rat.discr_prime_pow l k A]
  have hsign : IsUnit ((-1 : ℤ) ^ ((l ^ k).totient / 2)) :=
    (isUnit_one.neg).pow _
  refine IsCoprime.mul_left (isCoprime_of_isUnit_left hsign) (IsCoprime.pow_left ?_)
  have hlp : Prime (l : ℤ) := by
    rw [Int.prime_iff_natAbs_prime]
    simpa using (Fact.out : l.Prime)
  exact (hlp.coprime_iff_not_dvd).2 hd

/-- ★★★★★★★★★★★★**`ℚ(ζ_{l^k})` と `L` は `ℚ` 上線型無関連**。

仮説 `hL : ¬ (l : ℤ) ∣ discr B` は「**`l` は `B` で不分岐**」である
（`NumberField.not_dvd_discr_iff_forall_mem`）。 -/
theorem cyclo_linearDisjoint (A B : IntermediateField ℚ M) (l k : ℕ) [Fact l.Prime]
    [IsCyclotomicExtension {l ^ k} ℚ A]
    (hL : ¬ (l : ℤ) ∣ NumberField.discr B) :
    A.LinearDisjoint B := by
  haveI : IsGalois ℚ A := IsCyclotomicExtension.isGalois {l ^ k} ℚ A
  exact NumberField.linearDisjoint_of_isGalois_isCoprime_discr (L := M) A B
    (discr_cyclo_isCoprime l k A _ hL)

/-- ★★★★★★★★★★★★★★★★**`[L(ζ_{l^k}) : L] = φ(l^k)`**。

★★これが葉 4 の核心である。原文が「完全分岐だから線型無関連」で畳んだ一行が、
判別式だけで（完全分岐を経由せずに）出た。 -/
theorem finrank_adjoin_cyclo (A B : IntermediateField ℚ M) (l k : ℕ) [Fact l.Prime]
    [NeZero (l ^ k)] [IsCyclotomicExtension {l ^ k} ℚ A]
    (hL : ¬ (l : ℤ) ∣ NumberField.discr B) :
    Module.finrank B (IntermediateField.adjoin B (A : Set M)) = (l ^ k).totient := by
  have hLD : A.LinearDisjoint B := cyclo_linearDisjoint A B l k hL
  have hrank : Module.rank B (IntermediateField.adjoin B (A : Set M)) = Module.rank ℚ A :=
    hLD.adjoin_rank_eq_rank_left_of_isAlgebraic (Or.inl inferInstance)
  have := congrArg Cardinal.toNat hrank
  rw [← Module.finrank, ← Module.finrank] at this
  rw [this]
  exact IsCyclotomicExtension.Rat.finrank (l ^ k) A

def isCoprime_of_isUnit_left.src : Source :=
  { paper := "GenEll", pdfPage := 20,
    item := "Theorem 3.8(配管——単元は何とでも互いに素)",
    sectionId := "genell-thm-3-8" }

def discr_cyclo_isCoprime.src : Source :=
  { paper := "GenEll", pdfPage := 20,
    item := "Theorem 3.8(disc ℚ(ζ_{l^k}) は l 以外の素因数を持たない)",
    sectionId := "genell-thm-3-8" }

def finrank_adjoin_cyclo.src : Source :=
  { paper := "GenEll", pdfPage := 20,
    item := "Theorem 3.8([L(ζ_{l^k}) : L] = φ(l^k))",
    sectionId := "genell-thm-3-8" }

def cyclo_linearDisjoint.src : Source :=
  { paper := "GenEll", pdfPage := 20,
    item := "Theorem 3.8(ℚ(ζ_{l^k}) と L は ℚ 上線型無関連)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GenEll
