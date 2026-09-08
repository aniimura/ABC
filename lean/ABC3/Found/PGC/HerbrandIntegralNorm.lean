import ABC3.Found.PGC.IntegerMiscInstances
import ABC3.Found.PGC.HasseArfStrongInduction

/-!
# [pGC] ★★★★★★Hasse–Arf の整数性を**ノルムの言葉だけ**で言う

`HasseArfStrongInduction.lean:447-459`
`exists_natCast_herbrandPhiGroup_of_lowerRamificationGroup_one_eq_top` は
`A`(= `𝒪_K`) `B`(= `𝒪_M`) `G` の 3 つの型と **13 個の仮説**を要求していた。
直前の波までに 9 項目が個別に埋まったので、本ファイルは

★★**それらを実際に組み立て、`K`・`M`・ノルムだけの仮説から結論を出す**。

## 何が出たか

`exists_herbrandPhiGroup_natCast_of_norm` ——
`φ_G(m) ∈ ℕ`(Hasse–Arf の整数性)が、
★**環・イデアル・分岐群の語を 1 つも仮説に持たない形**で出る。
仮説はすべて `K`・`M`・`π`・`g` とノルムの不等式である。

## ★★測定（本ファイルの主目的）

★★**組み立ては通った。** すなわち `[Fintype G]` `[IsDomain A]` `[IsDomain B]`
`[SMulCommClass G A B]` `[FaithfulSMul G B]` `[CharP (ResidueField B) p]`
`[IsNoetherian A B]` `[IsDiscreteValuationRing A]` `[IsDiscreteValuationRing B]`
`hπ'` `hresA` `hadj` `hAinj` `habel` `h1` の**全部**が
`letI` / `have` で供給でき、`exact` が通った(実測 8.3 秒)。

★前波までの表で「未測定」「高い」としていた項目はすべて閉じた。
★**`hne`(`G_m ≠ G_{m+1}`)だけが仮説として残る**が、これは原典でも仮説である
(跳びが起きる `m` を取る、という意味であり、空虚回避の条件)。

## ★★★残り(正確に、`file:line` つき)

`harith` (3) `p^{m+1} ∣ u(m+1) − u m` を閉じるのに要るのは、
`HasseArfCongruence.lean` の抽象核

    dvd_sub_of_phi_intCast (hint : ∀ m, φ m = (j m : ℝ))
      (hrec : ∀ m, φ (m+1) = φ m + (u(m+1) − u m)/p^{m+1}) : p^{m+1} ∣ u(m+1) − u m

の 2 入力である。本ファイルは **`hint` を供給する**(`φ_G(m) ∈ ℕ`)。
残るのは `hrec`(`φ` の漸化式)で、証拠は
`HerbrandComposition.lean:455 herbrandPhiGroup_natCast`
(`φ_G(n) = (1/|G|) Σ_{i=1}^{n} |G_i|`)と
`RamificationSubgroupCard.lean` の `card_eq_pow_of_mem_iff`(`|G_i| = p^{k+1−m}`)。

## 逸脱の記録

1. `G := M ≃ₐ[K] M` に固定している(`RamNormBridge.mem_ramification_iff` と
   `AlgEquiv.fintype` がこの形だから)。
2. `hiso` を全 `σ` について仮説で受ける。木の `norm_algEquiv_eq` が供給する形。
3. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**(import のみ)。
-/

namespace ABC3.Found.PGC

open IsLocalRing IsDiscreteValuationRing

namespace IntegerNorm

section Assemble

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]
  [FiniteDimensional K M]

/-- ★★★★★★**Hasse–Arf の整数性、ノルムの言葉だけの版**。

`φ_G(m)` が自然数になる。仮説は `K`・`M`・`π`・`π_K`・`g` とノルムの不等式のみ。 -/
theorem exists_herbrandPhiGroup_natCast_of_norm
    {π : M} {πK : K} {n t p m : ℕ} {g : M ≃ₐ[K] M} (hn1 : 1 < n)
    (hiso : ∀ (σ : M ≃ₐ[K] M) (z : M), ‖σ • z‖ = ‖z‖)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ i : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * i))
    (hval : ∀ z : M, z ≠ 0 → ∃ i : ℤ, ‖z‖ = ‖π‖ ^ i)
    (hπmem : π ∈ integerSubring M)
    (hbr : ‖g π - π‖ = ‖π‖ ^ (t + 1)) (ht : 1 ≤ t)
    (hgen : ∀ h : M ≃ₐ[K] M, ∃ j : ℤ, h = g ^ j)
    (hK0 : 0 < ‖algebraMap K M πK‖) (hK1 : ‖algebraMap K M πK‖ < 1)
    (hvalKK : ∀ a : K, a ≠ 0 → ∃ i : ℤ, ‖algebraMap K M a‖ = ‖algebraMap K M πK‖ ^ i)
    (hKmem : πK ∈ baseIntegerSubring K M)
    (hp : p.Prime) (hpM : ‖((p : ℕ) : M)‖ < 1)
    (hne : letI := isLocalRing_integerSubring (M := M)
      letI := integerMulSemiringAction hiso
      lowerRamificationGroup ↥(integerSubring M) (M ≃ₐ[K] M) m ≠
        lowerRamificationGroup ↥(integerSubring M) (M ≃ₐ[K] M) (m + 1)) :
    letI := integerMulSemiringAction hiso
    letI := isDiscreteValuationRing_integerSubring hπ0 hπ1 hval hπmem
    ∃ j : ℕ, herbrandPhiGroup (M ≃ₐ[K] M) (⟨π, hπmem⟩ : ↥(integerSubring M)) (m : ℝ) = (j : ℝ) := by
  letI := isLocalRing_integerSubring (M := M)
  letI := integerMulSemiringAction hiso
  letI := baseAlgebra K M
  letI := isDiscreteValuationRing_baseIntegerSubring hK0 hK1 hvalKK hKmem
  letI := isDiscreteValuationRing_integerSubring hπ0 hπ1 hval hπmem
  letI := isNoetherian_baseIntegerSubring hπ0 hπ1 hn hvalK hπmem hK0 hK1 hvalKK hKmem
  letI := smulCommClass_base (K := K) hiso (fun σ a => σ.commutes a)
  letI := faithfulSMul_integer hiso hπ0 hπ1 hπmem
  letI := charP_residueField (M := M) hp hpM
  exact exists_natCast_herbrandPhiGroup_of_lowerRamificationGroup_one_eq_top
    (A := ↥(baseIntegerSubring K M)) (B := ↥(integerSubring M)) (G := M ≃ₐ[K] M)
    p hp (irreducible_pi hπ0 hπ1 hval hπmem)
    (fun b => exists_sub_mem_maximalIdeal hn1 hπ0 hπ1 hn hvalK hval hπmem b)
    (adjoin_pi_eq_top hπ0 hπ1 hn hvalK hπmem)
    (injective_baseRingHom (K := K) (M := M))
    (mul_comm_of_forall_zpow hgen)
    (lowerRamificationGroup_one_eq_top hiso hπ0 hπ1 hn hvalK hval hπmem hbr ht hgen)
    hne

end Assemble

/-! ## `.src` と 公理 -/

def exists_herbrandPhiGroup_natCast_of_norm.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Theorem 6.11", sectionId := "thm-6-11" }

#print axioms exists_herbrandPhiGroup_natCast_of_norm

end IntegerNorm

end ABC3.Found.PGC
