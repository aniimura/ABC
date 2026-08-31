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

/-- ★★★★★★★★★★★★★★**`[L(ζ) : L] = φ(l^k)`**（原始根の形）。 -/
theorem finrank_adjoin_zeta (B : IntermediateField ℚ M) (l k : ℕ) [Fact l.Prime]
    [NeZero (l ^ k)] {ζ : M} (hζ : IsPrimitiveRoot ζ (l ^ k))
    (hL : ¬ (l : ℤ) ∣ NumberField.discr B) :
    Module.finrank B (IntermediateField.adjoin B ({ζ} : Set M)) = (l ^ k).totient := by
  set A : IntermediateField ℚ M := IntermediateField.adjoin ℚ {ζ} with hA
  haveI : IsCyclotomicExtension {l ^ k} ℚ A := hζ.intermediateField_adjoin_isCyclotomicExtension ℚ
  have hset : IntermediateField.adjoin B (A : Set M)
      = IntermediateField.adjoin B ({ζ} : Set M) := by
    apply le_antisymm
    · rw [IntermediateField.adjoin_le_iff]
      intro x hx
      have : A ≤ (IntermediateField.adjoin B ({ζ} : Set M)).restrictScalars ℚ := by
        rw [hA, IntermediateField.adjoin_le_iff]
        rintro y rfl
        exact IntermediateField.subset_adjoin B _ rfl
      exact this hx
    · rw [IntermediateField.adjoin_le_iff]
      rintro y rfl
      exact IntermediateField.subset_adjoin B _ (IntermediateField.subset_adjoin ℚ _ rfl)
  rw [← hset]
  exact finrank_adjoin_cyclo A B l k hL

/-- ★★★★★★★★★★★★★★★★★★**`cyclotomic (l^k)` は `L` 上既約**。

★★これが葉 4 の目標である。`IsCyclotomicExtension.autEquivPow` の仮説がこれ。 -/
theorem cyclotomic_irreducible_of_not_dvd_discr (B : IntermediateField ℚ M) (l k : ℕ)
    [Fact l.Prime] [NeZero (l ^ k)] {ζ : M} (hζ : IsPrimitiveRoot ζ (l ^ k))
    (hL : ¬ (l : ℤ) ∣ NumberField.discr B) :
    Irreducible (Polynomial.cyclotomic (l ^ k) B) := by
  have hint : IsIntegral B ζ := IsIntegral.tower_top (Algebra.IsIntegral.isIntegral (R := ℚ) ζ)
  have hroot : Polynomial.aeval ζ (Polynomial.cyclotomic (l ^ k) (↥B)) = 0 := by
    rw [Polynomial.aeval_def, ← Polynomial.eval_map, Polynomial.map_cyclotomic]
    exact hζ.isRoot_cyclotomic (NeZero.pos _)
  have hdvd : minpoly (↥B) ζ ∣ Polynomial.cyclotomic (l ^ k) (↥B) :=
    minpoly.dvd (↥B) ζ hroot
  have hdeg : (minpoly (↥B) ζ).natDegree = (l ^ k).totient := by
    rw [← IntermediateField.adjoin.finrank hint]
    exact finrank_adjoin_zeta B l k hζ hL
  have hcyc : (Polynomial.cyclotomic (l ^ k) (↥B)).natDegree = (l ^ k).totient :=
    Polynomial.natDegree_cyclotomic _ _
  have heq : minpoly (↥B) ζ = Polynomial.cyclotomic (l ^ k) (↥B) := by
    refine Polynomial.eq_of_monic_of_associated (minpoly.monic hint)
      (Polynomial.cyclotomic.monic _ _) ?_
    exact (Polynomial.associated_of_dvd_of_natDegree_le hdvd
      (Polynomial.cyclotomic_ne_zero _ _) (by rw [hdeg, hcyc]))
  rw [← heq]
  exact minpoly.irreducible hint

/-- ★★★★★★★★★★★★★★★★★★★★**`ζ ↦ ζ^c` を実現する `K`-自己同型が存在する**。

`cyclotomic n K` が既約なら、`ζ` と `ζ^c`（`c` は単元）は同じ最小多項式を持つので、
冪基底の間の同型がそのまま体の自己同型になる。 -/
theorem exists_algEquiv_pow {K N : Type*} [Field K] [Field N] [Algebra K N] (n : ℕ) [NeZero n]
    [NeZero ((n : ℕ) : K)]
    [IsCyclotomicExtension {n} K N] (h : Irreducible (Polynomial.cyclotomic n K))
    {ζ : N} (hζ : IsPrimitiveRoot ζ n) (c : (ZMod n)ˣ) :
    ∃ σ : N ≃ₐ[K] N, σ ζ = ζ ^ ((c : ZMod n)).val := by
  have hcop : ((c : ZMod n)).val.Coprime n := ZMod.val_coe_unit_coprime c
  have hμ : IsPrimitiveRoot (ζ ^ ((c : ZMod n)).val) n := hζ.pow_of_coprime _ hcop
  have hmin : minpoly K (hζ.powerBasis K).gen = minpoly K (hμ.powerBasis K).gen := by
    rw [IsPrimitiveRoot.powerBasis_gen, IsPrimitiveRoot.powerBasis_gen]
    exact (hζ.minpoly_eq_cyclotomic_of_irreducible h).symm.trans
      (hμ.minpoly_eq_cyclotomic_of_irreducible h)
  refine ⟨(hζ.powerBasis K).equivOfMinpoly (hμ.powerBasis K) hmin, ?_⟩
  have := (hζ.powerBasis K).equivOfMinpoly_gen (hμ.powerBasis K) hmin
  simpa [IsPrimitiveRoot.powerBasis_gen] using this

/-- ★★★★★★★★★★★★★★★★★★★★★★**`l` が `L` で不分岐なら `ζ ↦ ζ^c` が `L` 上実現できる**。

★★★これが葉 4 の到達点である。`L(ζ)` の自己同型として得られる。 -/
theorem exists_algEquiv_pow_adjoin (B : IntermediateField ℚ M) (l k : ℕ) [Fact l.Prime]
    [NeZero (l ^ k)] {ζ : M} (hζ : IsPrimitiveRoot ζ (l ^ k))
    (hL : ¬ (l : ℤ) ∣ NumberField.discr B) (c : (ZMod (l ^ k))ˣ) :
    ∃ σ : (IntermediateField.adjoin (↥B) ({ζ} : Set M))
            ≃ₐ[↥B] (IntermediateField.adjoin (↥B) ({ζ} : Set M)),
      ((σ ⟨ζ, IntermediateField.subset_adjoin _ _ rfl⟩ : _) : M)
        = ζ ^ ((c : ZMod (l ^ k))).val := by
  set N : IntermediateField (↥B) M := IntermediateField.adjoin (↥B) ({ζ} : Set M) with hN
  have hmem : ζ ∈ N := IntermediateField.subset_adjoin _ _ rfl
  haveI : IsCyclotomicExtension {l ^ k} (↥B) N :=
    hζ.intermediateField_adjoin_isCyclotomicExtension (↥B)
  have hζN : IsPrimitiveRoot (⟨ζ, hmem⟩ : N) (l ^ k) := by
    refine IsPrimitiveRoot.of_map_of_injective (f := N.val.toMonoidHom) ?_ N.val.injective
    exact hζ
  haveI : NeZero (((l ^ k : ℕ) : ℕ) : ↥B) := ⟨by
    simpa using (Nat.cast_ne_zero (R := ↥B)).2 (NeZero.ne (l ^ k))⟩
  obtain ⟨σ, hσ⟩ := exists_algEquiv_pow (K := ↥B) (N := N) (l ^ k)
    (cyclotomic_irreducible_of_not_dvd_discr B l k hζ hL) hζN c
  exact ⟨σ, by simpa using congrArg (fun x : N => (x : M)) hσ⟩

section BaseFree

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**周囲体を要らない形**——
`L` が数体で `l ∤ disc L` なら `cyclotomic (l^k)` は `L` 上既約。

★★★★これが葉 4 で実際に使う形である。周囲の数体は証明の中で
`L ⊔ ℚ(ζ)` として作る。 -/
theorem cyclotomic_irreducible_of_not_dvd_discr' (B : Type*) [Field B] [NumberField B]
    (l k : ℕ) [Fact l.Prime] [NeZero (l ^ k)]
    (hL : ¬ (l : ℤ) ∣ NumberField.discr B) :
    Irreducible (Polynomial.cyclotomic (l ^ k) B) := by
  classical
  set Ω := AlgebraicClosure B with hΩ
  haveI : NeZero ((l ^ k : ℕ) : Ω) := ⟨by
    simpa using (Nat.cast_ne_zero (R := Ω)).2 (NeZero.ne (l ^ k))⟩
  -- ★原始 `l^k` 乗根を取る
  obtain ⟨ζ, hζ0⟩ : ∃ ζ : Ω, (Polynomial.cyclotomic (l ^ k) Ω).IsRoot ζ := by
    refine IsAlgClosed.exists_root _ ?_
    rw [Polynomial.degree_cyclotomic]
    exact_mod_cast (Nat.totient_pos.2 (NeZero.pos _)).ne'
  have hζ : IsPrimitiveRoot ζ (l ^ k) :=
    (Polynomial.isRoot_cyclotomic_iff (n := l ^ k) (R := Ω)).1 hζ0
  -- ★`B` の像と `ℚ(ζ)` を `Ω` の中で取り、その合成体を周囲の数体にする
  set Bim : IntermediateField ℚ Ω := (IsScalarTower.toAlgHom ℚ B Ω).fieldRange with hBim
  set Z : IntermediateField ℚ Ω := IntermediateField.adjoin ℚ ({ζ} : Set Ω) with hZ
  set M : IntermediateField ℚ Ω := Bim ⊔ Z with hMdef
  have eB : B ≃ₐ[ℚ] ↥Bim := AlgEquiv.ofInjectiveField (IsScalarTower.toAlgHom ℚ B Ω)
  haveI : FiniteDimensional ℚ ↥Bim := eB.toLinearEquiv.finiteDimensional
  haveI : Algebra.IsIntegral ℚ Ω := Algebra.IsIntegral.trans B
  haveI : FiniteDimensional ℚ ↥Z :=
    IntermediateField.finiteDimensional_adjoin
      (fun x _ => Algebra.IsIntegral.isIntegral (R := ℚ) x)
  haveI : FiniteDimensional ℚ ↥M := IntermediateField.finiteDimensional_sup Bim Z
  haveI : NumberField ↥M := ⟨⟩
  -- ★`M` の中で見た `B` と `ζ`
  set B'' : IntermediateField ℚ ↥M := Bim.restrict le_sup_left with hB''
  have eB'' : ↥Bim ≃ₐ[ℚ] ↥B'' := IntermediateField.restrict_algEquiv _
  have hζZ : ζ ∈ Z := IntermediateField.subset_adjoin ℚ _ rfl
  have hζM : ζ ∈ M := (le_sup_right : Z ≤ M) hζZ
  have hζ' : IsPrimitiveRoot (⟨ζ, hζM⟩ : ↥M) (l ^ k) := by
    refine IsPrimitiveRoot.of_map_of_injective (f := M.val.toMonoidHom) ?_ M.val.injective
    exact hζ
  have hdiscr : NumberField.discr ↥B'' = NumberField.discr B :=
    (NumberField.discr_eq_discr_of_algEquiv B (eB.trans eB'')).symm
  have hirr : Irreducible (Polynomial.cyclotomic (l ^ k) ↥B'') :=
    cyclotomic_irreducible_of_not_dvd_discr B'' l k hζ' (by rw [hdiscr]; exact hL)
  -- ★`B ≃ₐ[ℚ] B''` で既約性を移す
  have e : B ≃+* ↥B'' := (eB.trans eB'').toRingEquiv
  have hmap : Polynomial.map (e : B →+* ↥B'') (Polynomial.cyclotomic (l ^ k) B)
      = Polynomial.cyclotomic (l ^ k) ↥B'' := Polynomial.map_cyclotomic _ _
  have hiff := MulEquiv.irreducible_iff (M := Polynomial B) (N := Polynomial ↥B'')
    (Polynomial.mapEquiv e) (x := Polynomial.cyclotomic (l ^ k) B)
  rw [show (Polynomial.mapEquiv e) (Polynomial.cyclotomic (l ^ k) B)
      = Polynomial.cyclotomic (l ^ k) ↥B'' from hmap] at hiff
  exact hiff.1 hirr

end BaseFree

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

def finrank_adjoin_zeta.src : Source :=
  { paper := "GenEll", pdfPage := 20,
    item := "Theorem 3.8([L(ζ) : L] = φ(l^k)——原始根の形)",
    sectionId := "genell-thm-3-8" }

def cyclotomic_irreducible_of_not_dvd_discr.src : Source :=
  { paper := "GenEll", pdfPage := 20,
    item := "Theorem 3.8(cyclotomic (l^k) は L 上既約)",
    sectionId := "genell-thm-3-8" }

def exists_algEquiv_pow.src : Source :=
  { paper := "GenEll", pdfPage := 20,
    item := "Theorem 3.8(ζ ↦ ζ^c を実現する K-自己同型)",
    sectionId := "genell-thm-3-8" }

def exists_algEquiv_pow_adjoin.src : Source :=
  { paper := "GenEll", pdfPage := 20,
    item := "Theorem 3.8(l が L で不分岐なら ζ ↦ ζ^c が L 上実現できる)",
    sectionId := "genell-thm-3-8" }

def cyclotomic_irreducible_of_not_dvd_discr'.src : Source :=
  { paper := "GenEll", pdfPage := 20,
    item := "Theorem 3.8(cyclotomic (l^k) は L 上既約——周囲体を要らない形)",
    sectionId := "genell-thm-3-8" }

def cyclo_linearDisjoint.src : Source :=
  { paper := "GenEll", pdfPage := 20,
    item := "Theorem 3.8(ℚ(ζ_{l^k}) と L は ℚ 上線型無関連)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GenEll
