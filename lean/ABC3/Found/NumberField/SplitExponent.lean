/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NumberField.SplitSeparable
import Mathlib.RingTheory.Discriminant

/-!
# 完全分解する素数の存在は**無条件**になる（`exponent θ ≠ 0` の回収）

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.116。

## ★★何が残っていたか

`SplitSeparable.lean` の `infinite_splitsCompletely'` は
「`ℚ` 上 Galois な数体で完全分解する素数は無限個」を **Chebotarev を使わずに**与えるが、

    hexp : RingOfIntegers.exponent θ ≠ 0

を仮定として残していた。`exponent θ = absNorm (under ℤ (conductor ℤ θ))` なので、
これは「**`ℤ[θ]` が `𝓞 N` の中で有限指数**」＝「`θ` が原始元」に相当する。

## ★★★回収の道筋 —— 判別式

| 段 | 主張 | 根拠 |
|---|---|---|
| 1 | 原始元は**整**に取れる | `exists_integral_multiple`(分母を払う) |
| 2 | 整な原始元は `PowerBasis ℚ N` を与える | `adjoin.powerBasis` ＋ `topEquiv` |
| 3 | `disc(θ) · 𝓞 N ⊆ ℤ[θ]` | ★mathlib の `discr_mul_isIntegral_mem_adjoin` |
| 4 | `disc(θ)` は `0` でない**整数** | `discr_isIntegral` ＋ `discr_isUnit_of_basis` |
| 5 | よって `exponent θ ≠ 0` | `absNorm_eq_zero_iff` |

★★段 3 が要点で、これは古典的な「`disc(θ)` は導手を割る」である。
mathlib に `Algebra.discr_mul_isIntegral_mem_adjoin` として入っていた。

## ★★★★これで無条件になるもの

* `exists_splitsCompletely_prime` —— 完全分解する有理素数の存在（`Theorem 6.4, (i)` の数論）
* `finrank_isGreatest_deg'` —— `[L:ℚ]` は `deg(L,v)` の最大値（`Theorem 6.4, (iv)` の (a)）

★どちらも **Chebotarev の密度定理を 1 度も使わない**。
-/

namespace ABC3.Found.NF

open _root_.NumberField Polynomial Ideal Algebra
open scoped _root_.NumberField

/-! ## ★段 1・2. 整な原始元と `PowerBasis` -/

open IntermediateField in
/-- ★★★**数体には整な生成元をもつ `PowerBasis ℚ N` がある**。

★原始元 `α`（`Field.exists_primitive_element`）は整とは限らないが、
`exists_integral_multiple` で分母を払った `β = m • α` は整で、
`m` が `ℚ` で可逆なので `ℚ⟮β⟯ = ℚ⟮α⟯ = ⊤` である。 -/
theorem exists_integral_powerBasis (N : Type*) [Field N] [NumberField N] :
    ∃ B : PowerBasis ℚ N, IsIntegral ℤ B.gen := by
  classical
  obtain ⟨α, hα⟩ := Field.exists_primitive_element ℚ N
  have halg : IsAlgebraic ℤ α :=
    (IsFractionRing.isAlgebraic_iff ℤ ℚ N).mpr (Algebra.IsAlgebraic.isAlgebraic α)
  obtain ⟨m, hm0, hmi⟩ := halg.exists_integral_multiple
  set β : N := m • α with hβ
  have hmQ : ((m : ℚ)) ≠ 0 := Int.cast_ne_zero.mpr hm0
  have hβα : β = ((m : ℚ)) • α := by
    rw [hβ, ← Int.cast_smul_eq_zsmul ℚ m α]
  have hge : ℚ⟮α⟯ ≤ ℚ⟮β⟯ := by
    rw [adjoin_simple_le_iff]
    have hαe : α = ((m : ℚ))⁻¹ • β := by
      rw [hβα, smul_smul, inv_mul_cancel₀ hmQ, one_smul]
    rw [hαe, Algebra.smul_def]
    exact mul_mem (IntermediateField.algebraMap_mem _ _) (mem_adjoin_simple_self ℚ β)
  have htop : ℚ⟮β⟯ = ⊤ := le_antisymm le_top (hα ▸ hge)
  have hβint : IsIntegral ℚ β := Algebra.IsIntegral.isIntegral β
  refine ⟨(adjoin.powerBasis (K := ℚ) (L := N) hβint).map
    ((IntermediateField.equivOfEq htop).trans IntermediateField.topEquiv), ?_⟩
  show IsIntegral ℤ (((IntermediateField.equivOfEq htop).trans IntermediateField.topEquiv)
    (adjoin.powerBasis (K := ℚ) (L := N) hβint).gen)
  have hg : ((IntermediateField.equivOfEq htop).trans IntermediateField.topEquiv)
      (adjoin.powerBasis (K := ℚ) (L := N) hβint).gen = β := rfl
  rw [hg, hβ]
  exact hmi

/-! ## ★段 3〜5. 判別式が導手に入る整数を与える -/

/-- ★**`𝓞 N` の中の `ℤ[θ]` は `N` の中の `ℤ[θ]` から読み取れる**（包含が単射なので）。 -/
theorem mem_adjoin_of_coe_mem (N : Type*) [Field N] [NumberField N] (θ x : 𝓞 N)
    (h : (x : N) ∈ Algebra.adjoin ℤ ({(θ : N)} : Set N)) :
    x ∈ Algebra.adjoin ℤ ({θ} : Set (𝓞 N)) := by
  have hinj : Function.Injective (algebraMap (𝓞 N) N) := RingOfIntegers.coe_injective
  have hmap := AlgHom.map_adjoin_singleton (IsScalarTower.toAlgHom ℤ (𝓞 N) N) θ
  rw [show ((IsScalarTower.toAlgHom ℤ (𝓞 N) N) θ) = ((θ : N)) from rfl] at hmap
  rw [← hmap] at h
  obtain ⟨y, hy, hyx⟩ := h
  exact hinj hyx ▸ hy

/-- ★★★★★**`PowerBasis` の生成元は `exponent ≠ 0`**。

★中身は `disc(θ) · 𝓞 N ⊆ ℤ[θ]`（mathlib の `discr_mul_isIntegral_mem_adjoin`）で、
`disc(θ)` は `0` でない整数（`discr_isIntegral` ＋ `discr_isUnit_of_basis`）だから
導手 `conductor ℤ θ` が `ℤ` と自明でなく交わる。 -/
theorem exponent_ne_zero_of_powerBasis (N : Type*) [Field N] [NumberField N]
    (B : PowerBasis ℚ N) (hgen : IsIntegral ℤ B.gen) :
    _root_.RingOfIntegers.exponent (⟨B.gen, hgen⟩ : 𝓞 N) ≠ 0 := by
  classical
  set θ : 𝓞 N := ⟨B.gen, hgen⟩ with hθ
  have hθc : ((θ : N)) = B.gen := rfl
  have hbi : ∀ i, IsIntegral ℤ (B.basis i) := by
    intro i; rw [PowerBasis.coe_basis]; exact hgen.pow _
  have hdint : IsIntegral ℤ (discr ℚ ⇑B.basis) := discr_isIntegral ℚ hbi
  obtain ⟨n, hn⟩ := IsIntegrallyClosed.isIntegral_iff.mp hdint
  have hu : IsUnit (discr ℚ ⇑B.basis) := discr_isUnit_of_basis ℚ B.basis
  have hn0 : n ≠ 0 := by
    intro h0
    apply hu.ne_zero
    rw [← hn, h0]
    simp
  have hmem : (n : 𝓞 N) ∈ conductor ℤ θ := by
    rw [mem_conductor_iff]
    intro b
    refine mem_adjoin_of_coe_mem N θ _ ?_
    rw [hθc]
    have hz : IsIntegral ℤ ((b : N)) := b.isIntegral_coe
    have hd := discr_mul_isIntegral_mem_adjoin (K := ℚ) (R := ℤ) (B := B) hgen hz
    refine Eq.mpr (congrArg (fun t : N => t ∈ Algebra.adjoin ℤ ({B.gen} : Set N)) ?_) hd
    rw [Algebra.smul_def, ← hn]
    simp
  rw [_root_.RingOfIntegers.exponent, Ne, Ideal.absNorm_eq_zero_iff]
  intro hbot
  have hmem2 : n ∈ Ideal.under ℤ (conductor ℤ θ) := by
    show (algebraMap ℤ (𝓞 N)) n ∈ conductor ℤ θ
    simpa using hmem
  rw [hbot] at hmem2
  exact hn0 (Ideal.mem_bot.mp hmem2)

/-- ★★★**どの数体にも `exponent ≠ 0` の代数的整数がある**。 -/
theorem exists_exponent_ne_zero (N : Type*) [Field N] [NumberField N] :
    ∃ θ : 𝓞 N, _root_.RingOfIntegers.exponent θ ≠ 0 := by
  obtain ⟨B, hgen⟩ := exists_integral_powerBasis N
  exact ⟨⟨B.gen, hgen⟩, exponent_ne_zero_of_powerBasis N B hgen⟩

/-! ## ★★★★★無条件になった主張 -/

/-- ★★★★★★**`ℚ` 上 Galois な数体で完全分解する素数は無限個**（★仮定なし）。

★★**Chebotarev の密度定理は 1 度も使わない**。 -/
theorem infinite_splitsCompletely_of_isGalois (N : Type*) [Field N] [NumberField N]
    [IsGalois ℚ N] :
    {p : ℕ | p.Prime ∧
      (primesOver (span {(p : ℤ)}) (𝓞 N)).ncard = Module.finrank ℚ N}.Infinite := by
  obtain ⟨θ, hexp⟩ := exists_exponent_ne_zero N
  exact infinite_splitsCompletely' θ hexp

/-- ★★★★★★**完全分解する有理素数が存在する**（★仮定なし）。

★これが `Theorem 6.4, (i)` の「数体の自己同型で全素点を固定するものは恒等」の
唯一残っていた入力である。 -/
theorem exists_splitsCompletely_prime (N : Type*) [Field N] [NumberField N] [IsGalois ℚ N] :
    ∃ p : ℕ, p.Prime ∧
      (primesOver (span {(p : ℤ)}) (𝓞 N)).ncard = Module.finrank ℚ N :=
  (infinite_splitsCompletely_of_isGalois N).nonempty

/-! ### ★出典の紐付け -/

def exists_integral_powerBasis.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4 — 数体には整な生成元をもつ PowerBasis がある",
    sectionId := "frdi-thm-6-4" }

def exponent_ne_zero_of_powerBasis.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4 — 判別式が導手に入る整数を与える(exponent ≠ 0)",
    sectionId := "frdi-thm-6-4" }

def exponent_ne_zero_of_powerBasis.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "discr_mul_isIntegral_mem_adjoin(disc(θ)·𝓞 N ⊆ ℤ[θ])"
      (.inMathlib "Algebra.discr_mul_isIntegral_mem_adjoin") 116,
    .citation "[mathlib]" "discr_isUnit_of_basis(判別式は 0 でない)"
      (.inMathlib "Algebra.discr_isUnit_of_basis") 116 ]

def exists_splitsCompletely_prime.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4 — 完全分解する有理素数の存在(仮定なし・Chebotarev 不使用)",
    sectionId := "frdi-thm-6-4" }

def exists_splitsCompletely_prime.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "infinite_splitsCompletely'(Chebotarev を使わない迂回路)"
      (.inProject "ABC3" "ABC3.Found.NF.infinite_splitsCompletely'") 116,
    .citation "[ABC3]" "exists_exponent_ne_zero(残っていた仮定 hexp の回収)"
      (.inProject "ABC3" "ABC3.Found.NF.exists_exponent_ne_zero") 116 ]

end ABC3.Found.NF
