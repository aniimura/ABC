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

/-! ### ★出典の紐付け -/

/-- ★★★★★locator —— `cheb-spl-det` の第 2 段。 -/
def mem_SplQ_iff_restrictNormalHom.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — p ∈ Spl(K) ⟺ 分解群の K への制限が自明",
    sectionId := "frdi-thm-6-4" }

end ABC3.Found.NF
