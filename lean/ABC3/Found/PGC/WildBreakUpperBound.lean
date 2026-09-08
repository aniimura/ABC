import ABC3.Found.PGC.WildBreakLowerBound
import ABC3.Found.PGC.RamificationJumpBound

/-!
# 骨組み(作業中) —— 跳びの上界 `(p−1)t ≤ p·e`
-/

namespace ABC3.Found.PGC

namespace WildBreakUpper

open Finset WildBreak

/-! ## §1 抽象核 —— 共役はすべて等距離 -/

section Equidistant

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-- ★**抽象核** —— `g` が等長なら `‖g^j π − π‖ ≤ ‖g π − π‖`(超距離＋群の性質だけ)。 -/
theorem norm_pow_apply_sub_le (g : M ≃ₐ[K] M) (hiso : ∀ z : M, ‖g z‖ = ‖z‖) (π : M) :
    ∀ j : ℕ, ‖(g ^ j) π - π‖ ≤ ‖g π - π‖ := by
  intro j
  induction j with
  | zero => simp
  | succ j ih =>
      have hrw : (g ^ (j + 1)) π - π = g ((g ^ j) π - π) + (g π - π) := by
        rw [map_sub, ← AlgEquiv.mul_apply, ← pow_succ']
        ring
      rw [hrw]
      refine (IsUltrametricDist.norm_add_le_max _ _).trans (max_le ?_ le_rfl)
      rw [hiso]
      exact ih

/-- ★★**抽象核** —— `orderOf g = p`(素数)で `p ∤ j` なら `‖g^j π − π‖ = ‖g π − π‖`。

★`g^j` も生成元なので上の不等式を**両向き**に使える。★分岐・付値の語彙は 1 語も出ない。 -/
theorem norm_pow_apply_sub_eq {p : ℕ} (hp : p.Prime) (g : M ≃ₐ[K] M) (hg : orderOf g = p)
    (hiso : ∀ z : M, ‖g z‖ = ‖z‖) (π : M) {j : ℕ} (hj : ¬ (p ∣ j)) :
    ‖(g ^ j) π - π‖ = ‖g π - π‖ := by
  haveI : Fact p.Prime := ⟨hp⟩
  refine le_antisymm (norm_pow_apply_sub_le g hiso π j) ?_
  -- `j` の法 `p` 逆元 `k` を取る
  have hjne : ((j : ZMod p)) ≠ 0 := by
    rw [Ne, ZMod.natCast_eq_zero_iff]
    exact hj
  set k : ℕ := ((j : ZMod p)⁻¹).val with hkdef
  have hmod : j * k ≡ 1 [MOD p] := by
    rw [← ZMod.natCast_eq_natCast_iff]
    push_cast
    rw [hkdef, ZMod.natCast_val, ZMod.cast_id]
    field_simp
  have hgjk : (g ^ j) ^ k = g := by
    rw [← pow_mul, ← pow_one g]
    exact (pow_eq_pow_iff_modEq).mpr (by simpa [hg] using hmod)
  have hisoj : ∀ z : M, ‖(g ^ j) z‖ = ‖z‖ := norm_pow_apply g hiso j
  have := norm_pow_apply_sub_le (g ^ j) hisoj π k
  rwa [hgjk] at this

end Equidistant

/-! ## §2 Galois の段 —— ★**ノルムが 1 語も出ない**(#323 の薬: `[Field M]` だけで書く) -/

section GaloisPart

/-- `M = K(π)` なら `deg minpoly = [M:K]`。★`NormedField` を**使わない**形で書く(#323)。 -/
theorem natDegree_minpoly_of_adjoin_eq_top {K M : Type*} [Field K] [Field M] [Algebra K M]
    [FiniteDimensional K M] {π : M} (htop : Algebra.adjoin K ({π} : Set M) = ⊤) :
    (minpoly K π).natDegree = Module.finrank K M := by
  have hint : IsIntegral K π := IsIntegral.of_finite K π
  have htop' : IntermediateField.adjoin K ({π} : Set M) = ⊤ :=
    IntermediateField.adjoin_eq_top_of_algebra K _ htop
  have h1 : Module.finrank K (IntermediateField.adjoin K ({π} : Set M))
      = (minpoly K π).natDegree := IntermediateField.adjoin.finrank hint
  have h2 : Module.finrank K (IntermediateField.adjoin K ({π} : Set M))
      = Module.finrank K M := by
    rw [htop']
    exact LinearEquiv.finrank_eq (IntermediateField.topEquiv (F := K) (E := M)).toLinearEquiv
  rw [← h1, h2]

/-- ★★**`minpoly K π` の根はすべて `g` の軌道の中**。

`IsConjRoot.exists_algEquiv`(`Normal` が要る)で共役を与える `σ` を取り、
`Nat.card ⟨g⟩ = orderOf g = p = [M:K] = Nat.card Gal(M/K)` から `⟨g⟩ = ⊤` を出す。
★ノルムが 1 語も出ない(#323 の薬)。 -/
theorem exists_pow_of_isRoot {K M : Type*} [Field K] [Field M] [Algebra K M]
    [FiniteDimensional K M] [IsGalois K M] {p : ℕ} {π : M} (g : M ≃ₐ[K] M)
    (hg : orderOf g = p) (hnK : Module.finrank K M = p)
    {a : M} (haroot : Polynomial.aeval a (minpoly K π) = 0) :
    ∃ j : ℕ, (g ^ j) π = a := by
  have hint : IsIntegral K π := IsIntegral.of_finite K π
  have hconj : IsConjRoot K π a := (isConjRoot_iff_aeval_eq_zero hint).mpr haroot
  obtain ⟨σ, hσ⟩ := IsConjRoot.exists_algEquiv hconj
  have hcard : Nat.card (Subgroup.zpowers g) = Nat.card (M ≃ₐ[K] M) := by
    rw [Nat.card_zpowers, hg, IsGalois.card_aut_eq_finrank, hnK]
  have htop : Subgroup.zpowers g = ⊤ := Subgroup.eq_top_of_card_eq _ hcard
  have hmem : σ⁻¹ ∈ Subgroup.zpowers g := by rw [htop]; trivial
  obtain ⟨j, hj⟩ := mem_powers_iff_mem_zpowers.mpr hmem
  have hj' : g ^ j = σ⁻¹ := hj
  refine ⟨j, ?_⟩
  rw [hj', ← hσ]
  simp

end GaloisPart
/-! ## §3 ★★★★上界 `(p−1)·i ≤ p·e` の供給 -/

section Upper

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-- ★★★★★**跳びの上界** —— 次数 `p` の全分岐巡回拡大で `(p−1)·i ≤ p·e`。

`RamificationJumpBound.norm_natCast_le_pow_of_splits`(monic だけで足りる。★Eisenstein 性も
係数の整性も要らない)に、§1(共役は等距離)と §2(根は軌道の中)を差し込む。

★★**前波(＝自分)の記述の訂正**: 第 1105 の報告は
「上界には★**最小多項式の係数が整**であることが要る」と書いたが、★**これは誤りである。**
`norm_natCast_le_pow_of_splits` が要求するのは `Monic` / `Separable` / `Splits` /
`natDegree = n` だけで、係数の整性(Eisenstein 性)は**一度も使わない**
(`RamificationJumpBound.lean` 冒頭 docstring の「monic だけで足りた」が正しい)。 -/
theorem sub_one_mul_le_of_totallyRamified [FiniteDimensional K M] [IsGalois K M]
    {p e i : ℕ} (hp : p.Prime) {π : M}
    (g : M ≃ₐ[K] M) (hg : orderOf g = p) (hiso : ∀ z : M, ‖g z‖ = ‖z‖)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hnK : Module.finrank K M = p)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((p : ℤ) * m))
    (heM : ‖(p : M)‖ = ‖π‖ ^ (p * e))
    (hbr : ‖g π - π‖ = ‖π‖ ^ (i + 1)) :
    (p - 1) * i ≤ p * e := by
  have hint : IsIntegral K π := IsIntegral.of_finite K π
  have htop : Algebra.adjoin K ({π} : Set M) = ⊤ :=
    TotallyRamifiedLayer.adjoin_eq_top_of_valK hπ0 (ne_of_lt hπ1) hnK hvalK
  have hdeg : (minpoly K π).natDegree = p := by
    rw [natDegree_minpoly_of_adjoin_eq_top htop, hnK]
  have hgp : g ^ p = 1 := by rw [← hg]; exact pow_orderOf_eq_one g
  have hbreak : ∀ a : M, ((minpoly K π).map (algebraMap K M)).IsRoot a → a ≠ π →
      ‖π - a‖ = ‖π‖ ^ (i + 1) := by
    intro a haroot hane
    have haeval : Polynomial.aeval a (minpoly K π) = 0 := by
      rw [Polynomial.IsRoot, Polynomial.eval_map, ← Polynomial.aeval_def] at haroot
      exact haroot
    obtain ⟨j, hj⟩ := exists_pow_of_isRoot (p := p) g hg hnK haeval
    have hnd : ¬ (p ∣ j) := by
      intro ⟨c, hc⟩
      apply hane
      rw [← hj, hc, pow_mul, hgp, one_pow]
      rfl
    rw [norm_sub_rev, ← hj, norm_pow_apply_sub_eq hp g hg hiso π hnd, hbr]
  have hle : ‖(p : M)‖ ≤ ‖π‖ ^ ((p - 1) * i) :=
    norm_natCast_le_pow_of_splits hp.pos hπ0 hπ1 hvalK (minpoly.monic hint) hdeg
      (Algebra.IsSeparable.isSeparable K π) (minpoly.aeval K π)
      (Normal.splits (IsGalois.to_normal) π) hbreak
  exact sub_one_mul_le_of_norm_natCast_eq_pow hπ0 hπ1 heM hle

end Upper

end WildBreakUpper

end ABC3.Found.PGC
