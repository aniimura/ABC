/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NumberField.DecompGroup
import Mathlib.FieldTheory.Galois.Basic

/-!
# 分解群の制限(鎖 `cheb` の `cheb-spl-det` の第 2 段)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.116。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

## ★★★筋 —— 共通の `Ω` を固定して部分群で話す

合成体 `L ⊔ M` へ制限を降ろす形より、**共通の Galois 拡大 `Ω` を固定して
`G = Gal(Ω/ℚ)` の部分群で話す形**のほうが短い。

`𝔔` を `𝓞 Ω` の `p` の上の素、`D = stabilizer G 𝔔` とすると

    p ∈ Spl(K)  ⟺  D の `K` への制限が自明  ⟺  D ⊆ K.fixingSubgroup

であり、Galois 対応 `fixingSubgroup (L ⊔ M) = fixingSubgroup L ⊓ fixingSubgroup M`
(mathlib の `IntermediateField.fixingSubgroup_sup`)から**ただちに**

    Spl(L ⊔ M) = Spl(L) ∩ Spl(M)

が出る。★合成体ごとの素点の突き合わせが要らない。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `algebraMap_restrictNormal` | 制限した自己同型は `𝓞 K → 𝓞 Ω` と両立 |
| `under_smul` | ★★**素点の制限は Galois 作用と両立**(`𝔔 ↦ 𝔔 ∩ 𝓞K`) |
| `smul_restrictScalars_ideal` | `Gal(Ω/K)` と `Gal(Ω/ℚ)` の作用は一致 |
| `restrictNormalHom_restrictScalars_eq_one` | `Gal(Ω/K)` の元は `K` へ制限すると `1` |
| `exists_stabilizer_restrict` | ★★★★**分解群の制限は全射**(`D_𝔔 ↠ D_{𝔔∩K}`) |
| `mem_SplQ_iff_restrictNormalHom` | ★★★★★**`p ∈ Spl(K)` ⟺ 分解群の `K` への制限が自明** |

★★全射性の証明は 3 行 —— 持ち上げ(`AlgEquiv.restrictNormalHom_surjective`)、
`σ₀(𝔔)` も `𝔔∩K` の上にある、推移性(`Ideal.exists_smul_eq_of_isGaloisGroup`)で戻す。
-/

namespace ABC3.Found.NF

open _root_.NumberField IntermediateField Ideal
open scoped _root_.NumberField Pointwise

variable {Ω : Type*} [Field Ω] [NumberField Ω]

/-! ## ★1. 制限と `𝓞` の両立 -/

/-- ★制限した自己同型は `𝓞 K → 𝓞 Ω` と両立する。 -/
theorem algebraMap_restrictNormal (K : IntermediateField ℚ Ω) [Normal ℚ ↥K]
    (σ : Ω ≃ₐ[ℚ] Ω) (y : 𝓞 ↥K) :
    algebraMap (𝓞 ↥K) (𝓞 Ω) ((AlgEquiv.restrictNormalHom ↥K σ) • y)
      = σ • algebraMap (𝓞 ↥K) (𝓞 Ω) y := by
  refine FaithfulSMul.algebraMap_injective (𝓞 Ω) Ω ?_
  have h1 : ∀ z : 𝓞 ↥K, algebraMap (𝓞 Ω) Ω (algebraMap (𝓞 ↥K) (𝓞 Ω) z)
      = algebraMap ↥K Ω (algebraMap (𝓞 ↥K) ↥K z) := fun z => rfl
  have h2 : ∀ w : 𝓞 Ω, algebraMap (𝓞 Ω) Ω (σ • w) = σ (algebraMap (𝓞 Ω) Ω w) :=
    fun w => rfl
  rw [h1, h2, h1]
  have h3 : algebraMap (𝓞 ↥K) ↥K ((AlgEquiv.restrictNormalHom ↥K σ) • y)
      = (AlgEquiv.restrictNormalHom ↥K σ) (algebraMap (𝓞 ↥K) ↥K y) := rfl
  rw [h3]
  exact AlgEquiv.restrictNormal_commutes σ ↥K (algebraMap (𝓞 ↥K) ↥K y)

/-- ★★**素点の制限は Galois 作用と両立する**。 -/
theorem under_smul (K : IntermediateField ℚ Ω) [Normal ℚ ↥K] (σ : Ω ≃ₐ[ℚ] Ω)
    (Q : Ideal (𝓞 Ω)) :
    Ideal.under (𝓞 ↥K) (σ • Q)
      = (AlgEquiv.restrictNormalHom ↥K σ) • Ideal.under (𝓞 ↥K) Q := by
  ext y
  rw [Ideal.mem_pointwise_smul_iff_inv_smul_mem]
  show algebraMap (𝓞 ↥K) (𝓞 Ω) y ∈ σ • Q ↔ algebraMap (𝓞 ↥K) (𝓞 Ω) _ ∈ Q
  rw [Ideal.mem_pointwise_smul_iff_inv_smul_mem]
  have hinv : ((AlgEquiv.restrictNormalHom ↥K σ)⁻¹ • y)
      = (AlgEquiv.restrictNormalHom ↥K σ⁻¹) • y := by
    rw [map_inv]
  rw [hinv, algebraMap_restrictNormal K σ⁻¹ y]

/-! ## ★2. `Gal(Ω/K)` の元 -/

theorem restrictNormalHom_restrictScalars_eq_one (K : IntermediateField ℚ Ω) [Normal ℚ ↥K]
    (ρ : Ω ≃ₐ[↥K] Ω) : AlgEquiv.restrictNormalHom ↥K (ρ.restrictScalars ℚ) = 1 := by
  refine AlgEquiv.ext fun x => ?_
  have h := AlgEquiv.restrictNormal_commutes (ρ.restrictScalars ℚ) ↥K x
  have h2 : (ρ.restrictScalars ℚ) (algebraMap ↥K Ω x) = algebraMap ↥K Ω x := ρ.commutes x
  rw [h2] at h
  have h3 : ((AlgEquiv.restrictScalars ℚ ρ).restrictNormal ↥K) x = x :=
    (FaithfulSMul.algebraMap_injective ↥K Ω) h
  exact h3

theorem smul_restrictScalars_ideal (K : IntermediateField ℚ Ω) (ρ : Ω ≃ₐ[↥K] Ω)
    (I : Ideal (𝓞 Ω)) : (ρ.restrictScalars ℚ) • I = ρ • I := by
  ext x
  rw [Ideal.mem_pointwise_smul_iff_inv_smul_mem, Ideal.mem_pointwise_smul_iff_inv_smul_mem]
  rfl

/-! ## ★3. 分解群の制限は全射 -/

variable [IsGalois ℚ Ω]

/-- ★★★★**分解群の制限は全射** —— `D_𝔔 ↠ D_{𝔔∩K}`。

★★中身は 3 行:
1. `τ` を `σ₀ ∈ Gal(Ω/ℚ)` へ持ち上げる(`AlgEquiv.restrictNormalHom_surjective`)
2. `σ₀(𝔔)` も `𝔔∩K` の上にある(`under_smul`)
3. 推移性(`Ideal.exists_smul_eq_of_isGaloisGroup`)で `ρ ∈ Gal(Ω/K)` が戻す -/
theorem exists_stabilizer_restrict (K : IntermediateField ℚ Ω) [Normal ℚ ↥K]
    (Q : Ideal (𝓞 Ω)) (hQ : Q.IsPrime)
    (τ : ↥K ≃ₐ[ℚ] ↥K) (hτ : τ • (Ideal.under (𝓞 ↥K) Q) = Ideal.under (𝓞 ↥K) Q) :
    ∃ σ : Ω ≃ₐ[ℚ] Ω, σ • Q = Q ∧ AlgEquiv.restrictNormalHom ↥K σ = τ := by
  haveI := hQ
  obtain ⟨σ₀, hσ₀⟩ :=
    AlgEquiv.restrictNormalHom_surjective (F := ℚ) (K₁ := ↥K) Ω τ
  set pK : Ideal (𝓞 ↥K) := Ideal.under (𝓞 ↥K) Q with hpK
  haveI : Q.LiesOver pK := ⟨rfl⟩
  haveI : (σ₀ • Q).IsPrime := (Ideal.IsPrime.smul_iff σ₀).mpr hQ
  haveI : (σ₀ • Q).LiesOver pK := by
    refine ⟨?_⟩
    rw [hpK, under_smul K σ₀ Q, hσ₀, hτ]
  obtain ⟨ρ, hρ⟩ := Ideal.exists_smul_eq_of_isGaloisGroup pK (σ₀ • Q) Q (Ω ≃ₐ[↥K] Ω)
  refine ⟨(ρ.restrictScalars ℚ) * σ₀, ?_, ?_⟩
  · rw [mul_smul, smul_restrictScalars_ideal, hρ]
  · rw [map_mul, restrictNormalHom_restrictScalars_eq_one, one_mul, hσ₀]

/-! ## ★4. `p ∈ Spl(K)` の分解群による特徴づけ -/

theorem isPrime_under_int (K : IntermediateField ℚ Ω) (Q : Ideal (𝓞 Ω)) (hQ : Q.IsPrime) :
    (Ideal.under (𝓞 ↥K) Q).IsPrime := by
  haveI := hQ
  exact Ideal.IsPrime.under (𝓞 ↥K) Q

theorem liesOver_under_int (K : IntermediateField ℚ Ω) (Q : Ideal (𝓞 Ω)) {p : ℕ}
    (hQlo : Q.LiesOver (Ideal.span {(p : ℤ)})) :
    (Ideal.under (𝓞 ↥K) Q).LiesOver (Ideal.span {(p : ℤ)}) := by
  haveI := hQlo
  refine ⟨?_⟩
  have h : Ideal.under ℤ (Ideal.under (𝓞 ↥K) Q) = Ideal.under ℤ Q := Ideal.under_under Q
  rw [h]
  exact hQlo.over

/-- ★★★★★**`p ∈ Spl(K)` ⟺ 分解群の `K` への制限が自明**。

★→ は `under_smul`、← は `exists_stabilizer_restrict`(全射性)。 -/
theorem mem_SplQ_iff_restrictNormalHom (K : IntermediateField ℚ Ω) [IsGalois ℚ ↥K]
    {p : Nat.Primes} (Q : Ideal (𝓞 Ω)) (hQ : Q.IsPrime)
    (hQlo : Q.LiesOver (Ideal.span {((p : ℕ) : ℤ)})) :
    p ∈ SplQ ↥K ↔ ∀ σ : Ω ≃ₐ[ℚ] Ω, σ • Q = Q → AlgEquiv.restrictNormalHom ↥K σ = 1 := by
  have hPu := isPrime_under_int K Q hQ
  have hLOu := liesOver_under_int K Q hQlo
  rw [mem_SplQ_iff_stabilizer_trivial (Ideal.under (𝓞 ↥K) Q) hPu hLOu]
  constructor
  · intro h σ hσ
    refine h _ ?_
    rw [← under_smul K σ Q, hσ]
  · intro h τ hτ
    obtain ⟨σ, hσQ, hσres⟩ := exists_stabilizer_restrict K Q hQ τ hτ
    rw [← hσres]
    exact h σ hσQ

/-! ## ★5. `Spl(L ⊔ M) = Spl(L) ∩ Spl(M)` -/

/-- ★★`K` への制限が自明 ⟺ `K` を各点固定する。 -/
theorem restrictNormalHom_eq_one_iff (K : IntermediateField ℚ Ω) [Normal ℚ ↥K]
    (σ : Ω ≃ₐ[ℚ] Ω) :
    AlgEquiv.restrictNormalHom ↥K σ = 1 ↔ σ ∈ K.fixingSubgroup := by
  rw [IntermediateField.mem_fixingSubgroup_iff]
  constructor
  · intro h x hx
    have hh := AlgEquiv.restrictNormal_commutes σ ↥K ⟨x, hx⟩
    have h2 : (σ.restrictNormal ↥K) ⟨x, hx⟩ = ⟨x, hx⟩ := DFunLike.congr_fun h ⟨x, hx⟩
    rw [h2] at hh
    exact hh.symm
  · intro h
    refine AlgEquiv.ext fun y => ?_
    refine FaithfulSMul.algebraMap_injective ↥K Ω ?_
    have hh := AlgEquiv.restrictNormal_commutes σ ↥K y
    show algebraMap ↥K Ω ((σ.restrictNormal ↥K) y) = algebraMap ↥K Ω y
    rw [hh]
    exact h (algebraMap ↥K Ω y) y.2

theorem exists_prime_over {p : ℕ} (hp : p.Prime) :
    ∃ Q : Ideal (𝓞 Ω), Q.IsPrime ∧ Q.LiesOver (Ideal.span {(p : ℤ)}) := by
  haveI : (Ideal.span {(p : ℤ)}).IsMaximal := span_intCast_isMaximal hp
  obtain ⟨Q, hQ1, hQ2⟩ := Ideal.exists_isPrime_liesOver_of_faithfullyFlat
    (A := ℤ) (B := 𝓞 Ω) (Ideal.span {(p : ℤ)})
  exact ⟨Q, hQ1, hQ2⟩

/-- ★逸脱の記録: `IntermediateField.normal_sup` に `[Normal ℚ ↥L]` を
インスタンス探索で渡すと落ちる(2026-08-25 実測)。明示に渡す。 -/
theorem isGalois_sup (L M : IntermediateField ℚ Ω) (hL : Normal ℚ ↥L) (hM : Normal ℚ ↥M) :
    IsGalois ℚ ↥(L ⊔ M) := by
  haveI : Normal ℚ ↥(L ⊔ M) :=
    @IntermediateField.normal_sup ℚ Ω _ _ _ L M hL hM
  exact IsGalois.mk

/-- ★★★★★★**`Spl(L ⊔ M) = Spl(L) ∩ Spl(M)`** —— Galois 対応
(`IntermediateField.fixingSubgroup_sup`)1 行。

★★これが原文の「[again by Tchebotarev's density theorem] `L₁ ⊆ L₂`」の要である。 -/
theorem SplQ_sup (L M : IntermediateField ℚ Ω) (hL : IsGalois ℚ ↥L) (hM : IsGalois ℚ ↥M) :
    SplQ ↥(L ⊔ M) = SplQ ↥L ∩ SplQ ↥M := by
  haveI := hL
  haveI := hM
  haveI hLn : Normal ℚ ↥L := hL.to_normal
  haveI hMn : Normal ℚ ↥M := hM.to_normal
  haveI hsup : IsGalois ℚ ↥(L ⊔ M) := isGalois_sup L M hLn hMn
  haveI hsupn : Normal ℚ ↥(L ⊔ M) := hsup.to_normal
  ext p
  obtain ⟨Q, hQ, hQlo⟩ := exists_prime_over (Ω := Ω) p.2
  rw [mem_SplQ_iff_restrictNormalHom (L ⊔ M) Q hQ hQlo, Set.mem_inter_iff,
    mem_SplQ_iff_restrictNormalHom L Q hQ hQlo, mem_SplQ_iff_restrictNormalHom M Q hQ hQlo]
  constructor
  · intro h
    constructor
    · intro σ hσ
      rw [restrictNormalHom_eq_one_iff]
      have hmem := (restrictNormalHom_eq_one_iff (L ⊔ M) σ).mp (h σ hσ)
      rw [IntermediateField.fixingSubgroup_sup] at hmem
      exact hmem.1
    · intro σ hσ
      rw [restrictNormalHom_eq_one_iff]
      have hmem := (restrictNormalHom_eq_one_iff (L ⊔ M) σ).mp (h σ hσ)
      rw [IntermediateField.fixingSubgroup_sup] at hmem
      exact hmem.2
  · rintro ⟨h1, h2⟩ σ hσ
    rw [restrictNormalHom_eq_one_iff, IntermediateField.fixingSubgroup_sup]
    exact ⟨(restrictNormalHom_eq_one_iff L σ).mp (h1 σ hσ),
      (restrictNormalHom_eq_one_iff M σ).mp (h2 σ hσ)⟩

/-! ## ★6. 完全分解する素数の集合が体を決める -/

open Filter Topology in
/-- ★★★★★★★**[cheb-spl-det] 完全分解する素数の集合が体を決める**
([MilneCFT] Chapter V, Theorem 3.25 の底が `ℚ` の場合)。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

★★中身は 3 行 ——
`Spl(L ⊔ M) = Spl(L) ∩ Spl(M) = Spl(M)`、
密度(`tendsto_splQ_div_log`)を両方に当てて**極限の一意性**で `[L⊔M:ℚ] = [M:ℚ]`、
`IntermediateField.eq_of_le_of_finrank_eq` で `M = L ⊔ M`。 -/
theorem le_of_SplQ_subset (L M : IntermediateField ℚ Ω)
    (hL : IsGalois ℚ ↥L) (hM : IsGalois ℚ ↥M) (h : SplQ ↥M ⊆ SplQ ↥L) : L ≤ M := by
  haveI := hL
  haveI := hM
  haveI hsup : IsGalois ℚ ↥(L ⊔ M) := isGalois_sup L M hL.to_normal hM.to_normal
  have hset : SplQ ↥(L ⊔ M) = SplQ ↥M := by
    rw [SplQ_sup L M hL hM]
    exact Set.inter_eq_right.mpr h
  have h1 := tendsto_splQ_div_log (L := ↥(L ⊔ M))
  have h2 := tendsto_splQ_div_log (L := ↥M)
  rw [hset] at h1
  have heq : (1 : ℝ) / (Module.finrank ℚ ↥(L ⊔ M) : ℝ) = 1 / (Module.finrank ℚ ↥M : ℝ) :=
    tendsto_nhds_unique h1 h2
  have hfr : Module.finrank ℚ ↥M = Module.finrank ℚ ↥(L ⊔ M) := by
    have hpos1 : (0 : ℝ) < (Module.finrank ℚ ↥(L ⊔ M) : ℝ) := by
      have hh : 0 < Module.finrank ℚ ↥(L ⊔ M) := Module.finrank_pos
      positivity
    have hpos2 : (0 : ℝ) < (Module.finrank ℚ ↥M : ℝ) := by
      have hh : 0 < Module.finrank ℚ ↥M := Module.finrank_pos
      positivity
    field_simp at heq
    exact_mod_cast heq
  have hMeq : M = L ⊔ M := IntermediateField.eq_of_le_of_finrank_eq le_sup_right hfr
  rw [hMeq]
  exact le_sup_left

/-- ★★★★★★★**完全分解する素数の集合が体を決める**(等号版)——
`Theorem 6.4, (iv)` の「`L₁ ≅ L₂`」の中身。 -/
theorem eq_of_SplQ_eq (L M : IntermediateField ℚ Ω)
    (hL : IsGalois ℚ ↥L) (hM : IsGalois ℚ ↥M) (h : SplQ ↥L = SplQ ↥M) : L = M :=
  le_antisymm (le_of_SplQ_subset L M hL hM (h ▸ subset_rfl))
    (le_of_SplQ_subset M L hM hL (h ▸ subset_rfl))

/-! ## ★7. `Spl` は体の同型で不変 -/

section Congr

/-- ★環同型が誘導するイデアルの全単射。 -/
def idealEquivOfRingEquiv {R S : Type*} [CommRing R] [CommRing S] (e : R ≃+* S) :
    Ideal R ≃ Ideal S where
  toFun I := Ideal.map (e : R →+* S) I
  invFun J := Ideal.map (e.symm : S →+* R) J
  left_inv I := by
    show Ideal.map (e.symm : S →+* R) (Ideal.map (e : R →+* S) I) = I
    rw [Ideal.map_map]
    have h : (e.symm : S →+* R).comp (e : R →+* S) = RingHom.id R := by ext x; simp
    rw [h, Ideal.map_id]
  right_inv J := by
    show Ideal.map (e : R →+* S) (Ideal.map (e.symm : S →+* R) J) = J
    rw [Ideal.map_map]
    have h : (e : R →+* S).comp (e.symm : S →+* R) = RingHom.id S := by ext x; simp
    rw [h, Ideal.map_id]

/-- ★ノルムは環同型で不変(`absNorm I = #(S/I)` だから)。 -/
theorem absNorm_map_ringEquiv {L M : Type*} [Field L] [Field M] [NumberField L] [NumberField M]
    (e : L ≃+* M) (I : Ideal (𝓞 L)) :
    Ideal.absNorm (Ideal.map (RingOfIntegers.mapRingEquiv e : 𝓞 L →+* 𝓞 M) I)
      = Ideal.absNorm I := by
  rw [Ideal.absNorm_apply, Ideal.absNorm_apply, Submodule.cardQuot_apply,
    Submodule.cardQuot_apply]
  refine Nat.card_congr ?_
  exact (Ideal.quotientEquiv I (Ideal.map (RingOfIntegers.mapRingEquiv e : 𝓞 L →+* 𝓞 M) I)
    (RingOfIntegers.mapRingEquiv e) rfl).toEquiv.symm

/-- ★★イデアル計数は体の同型で不変。 -/
theorem idealCount_congr {L M : Type*} [Field L] [Field M] [NumberField L] [NumberField M]
    (e : L ≃+* M) (n : ℕ) : idealCount L n = idealCount M n := by
  refine Nat.card_congr ?_
  refine Equiv.subtypeEquiv (idealEquivOfRingEquiv (RingOfIntegers.mapRingEquiv e)) ?_
  intro I
  show Ideal.absNorm I = n ↔ Ideal.absNorm (Ideal.map
    (RingOfIntegers.mapRingEquiv e : 𝓞 L →+* 𝓞 M) I) = n
  rw [absNorm_map_ringEquiv e I]

/-- ★★`p ∈ Spl(L)` ⟺ `a_L(p) = [L:ℚ]`。 -/
theorem mem_SplQ_iff_idealCount {L : Type*} [Field L] [NumberField L] [IsGalois ℚ L]
    {p : Nat.Primes} : p ∈ SplQ L ↔ idealCount L (p : ℕ) = Module.finrank ℚ L := by
  constructor
  · exact fun h => idealCount_eq_finrank_of_splitsCompletely p.2 h
  · intro h
    have h1 : idealCount L (p : ℕ)
        ≤ (primesOver (Ideal.span {((p : ℕ) : ℤ)}) (𝓞 L)).ncard := by
      rw [idealCount_eq_ncard p.2]
      exact Set.ncard_le_ncard (fun P hP => hP.1) (finite_primesOver_int p.2)
    have h2 := ncard_primesOver_le_finrank_int (L := L) p.2
    show (primesOver (Ideal.span {((p : ℕ) : ℤ)}) (𝓞 L)).ncard = Module.finrank ℚ L
    omega

/-- ★★★**`Spl` は体の同型で不変**。 -/
theorem SplQ_congr {L M : Type*} [Field L] [Field M] [NumberField L] [NumberField M]
    [IsGalois ℚ L] [IsGalois ℚ M] (e : L ≃ₐ[ℚ] M) : SplQ L = SplQ M := by
  have hfr : Module.finrank ℚ L = Module.finrank ℚ M :=
    LinearEquiv.finrank_eq (e.toLinearEquiv)
  ext p
  rw [mem_SplQ_iff_idealCount, mem_SplQ_iff_idealCount,
    idealCount_congr (e.toRingEquiv) (p : ℕ), hfr]

/-! ## ★8. 抽象な 2 つの数体の版 -/

theorem isGalois_sup_ac (L M : IntermediateField ℚ (AlgebraicClosure ℚ))
    (hL : Normal ℚ ↥L) (hM : Normal ℚ ↥M) : IsGalois ℚ ↥(L ⊔ M) := by
  haveI : Normal ℚ ↥(L ⊔ M) :=
    @IntermediateField.normal_sup ℚ (AlgebraicClosure ℚ) _ _ _ L M hL hM
  exact IsGalois.mk

/-- ★★★★★★★**完全分解する素数の集合が体を決める**(抽象な数体の版)——
`Theorem 6.4, (iv)` の結論「`L₁ ≅ L₂`」の体論の側。

★★★中身は「両方を `ℚ̄` へ埋め込んで合成体の中で `eq_of_SplQ_eq` を当てる」だけ。
`Spl` が同型で不変(`SplQ_congr`)なので埋め込みの取り方によらない。 -/
theorem nonempty_algEquiv_of_SplQ_eq (L₁ L₂ : Type) [Field L₁] [Field L₂]
    [NumberField L₁] [NumberField L₂] [IsGalois ℚ L₁] [IsGalois ℚ L₂]
    (h : SplQ L₁ = SplQ L₂) : Nonempty (L₁ ≃ₐ[ℚ] L₂) := by
  classical
  let i₁ : L₁ →ₐ[ℚ] AlgebraicClosure ℚ := IsAlgClosed.lift
  let i₂ : L₂ →ₐ[ℚ] AlgebraicClosure ℚ := IsAlgClosed.lift
  let F₁ : IntermediateField ℚ (AlgebraicClosure ℚ) := i₁.fieldRange
  let F₂ : IntermediateField ℚ (AlgebraicClosure ℚ) := i₂.fieldRange
  let e₁ : L₁ ≃ₐ[ℚ] ↥F₁ := AlgHom.equivFieldRange
  let e₂ : L₂ ≃ₐ[ℚ] ↥F₂ := AlgHom.equivFieldRange
  haveI hf1 : FiniteDimensional ℚ ↥F₁ := LinearEquiv.finiteDimensional e₁.toLinearEquiv
  haveI hf2 : FiniteDimensional ℚ ↥F₂ := LinearEquiv.finiteDimensional e₂.toLinearEquiv
  haveI : NumberField ↥F₁ := ⟨⟩
  haveI : NumberField ↥F₂ := ⟨⟩
  haveI hg1 : IsGalois ℚ ↥F₁ := (AlgEquiv.transfer_galois e₁).mp ‹_›
  haveI hg2 : IsGalois ℚ ↥F₂ := (AlgEquiv.transfer_galois e₂).mp ‹_›
  haveI hfsup : FiniteDimensional ℚ ↥(F₁ ⊔ F₂) :=
    IntermediateField.finiteDimensional_sup F₁ F₂
  haveI : NumberField ↥(F₁ ⊔ F₂) := ⟨⟩
  haveI : IsGalois ℚ ↥(F₁ ⊔ F₂) := isGalois_sup_ac F₁ F₂ hg1.to_normal hg2.to_normal
  let G₁ : IntermediateField ℚ ↥(F₁ ⊔ F₂) :=
    (IntermediateField.inclusion (le_sup_left : F₁ ≤ F₁ ⊔ F₂)).fieldRange
  let G₂ : IntermediateField ℚ ↥(F₁ ⊔ F₂) :=
    (IntermediateField.inclusion (le_sup_right : F₂ ≤ F₁ ⊔ F₂)).fieldRange
  let f₁ : ↥F₁ ≃ₐ[ℚ] ↥G₁ := AlgHom.equivFieldRange
  let f₂ : ↥F₂ ≃ₐ[ℚ] ↥G₂ := AlgHom.equivFieldRange
  haveI : FiniteDimensional ℚ ↥G₁ := LinearEquiv.finiteDimensional f₁.toLinearEquiv
  haveI : FiniteDimensional ℚ ↥G₂ := LinearEquiv.finiteDimensional f₂.toLinearEquiv
  haveI : NumberField ↥G₁ := ⟨⟩
  haveI : NumberField ↥G₂ := ⟨⟩
  haveI hgg1 : IsGalois ℚ ↥G₁ := (AlgEquiv.transfer_galois f₁).mp hg1
  haveI hgg2 : IsGalois ℚ ↥G₂ := (AlgEquiv.transfer_galois f₂).mp hg2
  have hsplG : SplQ ↥G₁ = SplQ ↥G₂ := by
    rw [← SplQ_congr f₁, ← SplQ_congr e₁, h, SplQ_congr e₂, SplQ_congr f₂]
  have hGeq : G₁ = G₂ := eq_of_SplQ_eq G₁ G₂ hgg1 hgg2 hsplG
  exact ⟨e₁.trans (f₁.trans ((IntermediateField.equivOfEq hGeq).trans (f₂.symm.trans e₂.symm)))⟩

/-! ## ★9. `deg(Ψ^rlf) = 1` の算術の核 -/

/-- ★★★**次数 1 の素点が両側にあれば `deg = 1`**。

原文 (FrdI p.116):
> §4, Theorem 10], it follows that [Li : Q] is equal to the maximum of the deg(Li, vi),

★`d₂(f p₀) = c·d₁(p₀) = c` は正整数、`1 = d₂(q₀) = c·d₁(f⁻¹ q₀)` から
`c = 1/d₁(f⁻¹ q₀)`。積が `1` なのでどちらも `1`。
★★次数 1 の素点は「完全分解する素数」であり、無限個ある(`cheb-spl-infinite`)。 -/
theorem deg_eq_one_of_degOne_both (c : ℝ) (hc : 0 < c)
    (d₁ d₂ : Nat.Primes → ℕ) (f : Nat.Primes ≃ Nat.Primes)
    (hd : ∀ p, (d₂ (f p) : ℝ) = c * (d₁ p : ℝ))
    (h1 : ∃ p, d₁ p = 1) (h2 : ∃ q, d₂ q = 1)
    (hpos1 : ∀ p, 0 < d₁ p) : c = 1 := by
  obtain ⟨p₀, hp₀⟩ := h1
  obtain ⟨q₀, hq₀⟩ := h2
  have hA : (d₂ (f p₀) : ℝ) = c := by rw [hd p₀, hp₀]; norm_num
  have hB : c * (d₁ (f.symm q₀) : ℝ) = 1 := by
    have hh := hd (f.symm q₀)
    rw [Equiv.apply_symm_apply, hq₀] at hh
    simpa using hh.symm
  set n := d₂ (f p₀) with hn
  set m := d₁ (f.symm q₀) with hm
  have hmpos : 0 < m := hpos1 _
  have hnm : (n : ℝ) * (m : ℝ) = 1 := by rw [hA]; linarith
  have hnm' : n * m = 1 := by exact_mod_cast hnm
  have hn1 : n = 1 := Nat.eq_one_of_mul_eq_one_right (by rwa [mul_comm] at hnm')
  rw [← hA, hn1]
  norm_num

end Congr

/-! ### ★出典の紐付け -/

/-- ★★★★★locator —— `cheb-spl-det` の第 2 段。 -/
def mem_SplQ_iff_restrictNormalHom.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — p ∈ Spl(K) ⟺ 分解群の K への制限が自明",
    sectionId := "frdi-thm-6-4" }

/-- ★★★★★★★locator —— `cheb-spl-det` そのもの。 -/
def le_of_SplQ_subset.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — 完全分解する素数の集合が体を決める",
    sectionId := "frdi-thm-6-4" }

def le_of_SplQ_subset.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "tendsto_splQ_div_log(密度は 1/[L:ℚ])"
      (.inProject "ABC3" "ABC3.Found.NF.tendsto_splQ_div_log") 116,
    .citation "[ABC3]" "SplQ_sup(Spl(L ⊔ M) = Spl(L) ∩ Spl(M))"
      (.inProject "ABC3" "ABC3.Found.NF.SplQ_sup") 116,
    .citation "[mathlib]" "IntermediateField.fixingSubgroup_sup(Galois 対応)"
      (.inMathlib "IntermediateField.fixingSubgroup_sup") 116 ]

end ABC3.Found.NF
