import ABC3.Found.PGC.WildBreakPowDegree

/-!
# [pGC] ★★★★★★★★`harith` が組み上がった —— 出口から `harith` と `htop` が消えた

`JumpFromValueGroup.lean:265` の
`exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_jumpFree` は
`harith`（4 条件）を仮説で受けていた。本ファイルはそれを**内部で供給**する。

## ★★持ち場の 3 つの問いに対する実測

1. **`hsg`（`s j` と `g^{p^j}` の同一視）は本当に 1 行か** → ★**1 行**。
   `rw [← hsg j π]` だけ（本ファイル `hbr`）。
2. **`hiso` を全 `σ` で受ける形が `PureStepSetup.norm_algEquiv_eq` から出るか** → ★**出ない**。
   `PureStepSetup.lean:283` は
   `[NontriviallyNormedField k] [IsUltrametricDist k] [CompleteSpace k] [NormedAlgebra k M]
    [Algebra.IsAlgebraic k M]` を要求する。
   ★**素の `[Field K]`（ノルムなし）の設定では供給できない**。
   ⇒ `hiso : ∀ z, ‖g z‖ = ‖z‖`（**生成元だけ**）を仮説で受ける。逸脱に記録。
   ★全 `σ` 版は `htop`（`Gal = ⟨g⟩`）と `norm_zpow_apply` で**内部で作れる**ので、
   仮説は生成元 1 つで済む。
3. **(2)(3) が `∀ m < k` に限られている件** → ★(2) は**ちょうど合う**
   （`jump_lt_succ` は `m` と `m+1` の跡びを使うので `m < k` で `m+1 ≤ k`）。
   ★★(3) は**合わなかった**: `HasseArfCongruenceNorm.dvd_sub_jump_of_norm` は
   `hsmono : ∀ s t, s < t → u s < u t` を**全域**で要求するが、
   `harith` からは `k` 以下しか出ない。
   ⇒ ★**跡びの列を `k` の先へ延ばして回避した**:

       uN j := if j ≤ k then (u j).toNat else (u k).toNat + (j − k)

   これは全域で狭義単調で、`j ≤ k` では `u j` に一致するので
   `hjump`（`s ≤ k` のみ）も結論（`m+1 ≤ k`）もそのまま通る。

## ★内部で作れたもの（仮説で受けなかったもの）

| 必要だったもの | 供給元 |
|---|---|
| `hπ0 : 0 < ‖π‖` | `TotallyRamified.norm_pos_of_normp`（`TotallyRamifiedValueGroup.lean:391`） |
| `hπ1 : ‖π‖ < 1` | 本ファイル `norm_lt_one_of_normp` |
| `hval`（`M` の値群） | `TotallyRamified.exists_zpow_norm`（同 `:163`） |
| `htop`（`Gal = ⟨g⟩`） | 本ファイル `forall_mem_zpowers_of_orderOf_eq_finrank` |
| ★`πK`（`K` の素元） | 本ファイル `exists_base_uniformizer` |
| `htop`（`adjoin K {π} = ⊤`） | `TotallyRamifiedLayer.adjoin_eq_top_of_valK` |

★★**`πK` が一番危うかった**: `hvalK` は「`K` の値群が `‖π‖^{qℤ}` に**入る**」しか
言わないので、生成元を**実現する元**があるかは別問題である。
`TotallyRamified.exists_valSub_gen`（同 `:252`）の同値式

    (∃ d : K, d ≠ 0 ∧ ‖algebraMap K M d‖ = ‖π‖^m) ↔ (c : ℤ) ∣ m

を `m := c` で使うと★**元そのものが取れる**。
非自明性（`c ≠ 0`）は `‖(p:M)‖ = ‖π‖^{q·e}` と `q·e ≠ 0` から。

## 出したもの

| 宣言 | 内容 |
|---|---|
| `norm_lt_one_of_normp` | `‖π‖ < 1` |
| `norm_zpow_apply` | 等長性は整数冪に伝わる |
| `forall_mem_zpowers_of_orderOf_eq_finrank` | `Gal(M/K) = ⟨g⟩` |
| ★★`exists_base_uniformizer` | ★**`K` の素元の存在** |
| ★★★★★★★`harith_of_norm` | **`harith` の 4 条件をノルムだけから** |
| ★★★★★★★★`exists_norm_sub_algebraMap_le_prod_axDecay_of_norm` | **`harith` と `htop` が消えた出口** |

## 逸脱の記録

1. ★★`hiso : ∀ z, ‖g z‖ = ‖z‖` を**仮説で受ける**（上の問い 2）。
   具体層（`ℚ_p` を下に敷いた場合）では `PureStepSetup.lean:283 norm_algEquiv_eq` が供給する。
2. 跡びの列を `k` の先へ延ばした（上の問い 3）。結論は `m+1 ≤ k` でのみ使うので影響なし。
3. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
-/

namespace ABC3.Found.PGC

namespace HarithAssembly

/-! ## §1 補助 —— `harith` の文脈から取れるもの -/

section Helpers

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

omit [IsUltrametricDist M] in
/-- `‖p‖ = p⁻¹` と `‖p‖ = ‖π‖^{q·e}`（`q·e ≥ 1`）から `‖π‖ < 1`。 -/
theorem norm_lt_one_of_normp {p N : ℕ} (hp1 : 1 < p) {π : M} (_hN : 0 < N)
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹) (heM : ‖(p : M)‖ = ‖π‖ ^ N) : ‖π‖ < 1 := by
  by_contra hcon
  rw [not_lt] at hcon
  have h1 : (1 : ℝ) ≤ ‖π‖ ^ N := one_le_pow₀ hcon
  rw [← heM, hnormp] at h1
  have hp0 : (0 : ℝ) < (p : ℝ) := by
    have : (0 : ℕ) < p := by omega
    exact_mod_cast this
  have hp1' : (1 : ℝ) < (p : ℝ) := by exact_mod_cast hp1
  rw [le_inv_comm₀ one_pos hp0] at h1
  linarith

omit [IsUltrametricDist M] in
/-- 等長性は整数冪に伝わる。 -/
theorem norm_zpow_apply (g : M ≃ₐ[K] M) (hiso : ∀ z : M, ‖g z‖ = ‖z‖) :
    ∀ (j : ℤ) (z : M), ‖(g ^ j) z‖ = ‖z‖ := by
  have hnat : ∀ (n : ℕ) (z : M), ‖(g ^ n) z‖ = ‖z‖ := WildBreak.norm_pow_apply g hiso
  intro j z
  rcases j with n | n
  · rw [Int.ofNat_eq_natCast, zpow_natCast]; exact hnat n z
  · rw [zpow_negSucc]
    have h := hnat (n + 1) (((g ^ (n + 1))⁻¹ : M ≃ₐ[K] M) z)
    rw [show (g ^ (n + 1)) (((g ^ (n + 1))⁻¹ : M ≃ₐ[K] M) z) = z from
      AlgEquiv.apply_symm_apply _ z] at h
    exact h.symm

omit [IsUltrametricDist M] in
/-- ★`g` が位数 `= [M:K]` なら `Gal(M/K) = ⟨g⟩`。

`WildBreakUpper.isGalois_of_orderOf_eq_finrank` ＋ `IsGalois.card_aut_eq_finrank` ＋
`Subgroup.eq_top_of_card_eq`。 -/
theorem forall_mem_zpowers_of_orderOf_eq_finrank [FiniteDimensional K M] {n : ℕ}
    (g : M ≃ₐ[K] M) (hg : orderOf g = n) (hnK : Module.finrank K M = n) :
    ∀ σ : M ≃ₐ[K] M, σ ∈ Subgroup.zpowers g := by
  haveI := WildBreakUpper.isGalois_of_orderOf_eq_finrank g hg hnK
  have hcard : Nat.card ↥(Subgroup.zpowers g) = Nat.card (M ≃ₐ[K] M) := by
    rw [Nat.card_zpowers, hg, IsGalois.card_aut_eq_finrank, hnK]
  have htop : Subgroup.zpowers g = ⊤ := Subgroup.eq_top_of_card_eq _ hcard
  intro σ
  rw [htop]
  trivial

omit [IsUltrametricDist M] in
/-- ★★**底 `K` の素元** —— `K` の値群の生成元を実現する元が取れる。

`TotallyRamified.exists_valSub_gen`（`TotallyRamifiedValueGroup.lean:252`）で
値群の生成元 `c` を取り、`c ∣ c` の向きで**元そのもの**を取り出す。 -/
theorem exists_base_uniformizer {π : M} {n N p : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hN : 0 < N)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (hpK : (p : K) ≠ 0) (heM : ‖(p : M)‖ = ‖π‖ ^ N) :
    ∃ πK : K, 0 < ‖algebraMap K M πK‖ ∧ ‖algebraMap K M πK‖ < 1 ∧
      (∀ a : K, a ≠ 0 → ∃ m : ℤ,
        ‖algebraMap K M a‖ = ‖algebraMap K M πK‖ ^ m) := by
  obtain ⟨c, hc⟩ := TotallyRamified.exists_valSub_gen (E := K) (M := M) hπ0
  have hpmem : ∃ d : K, d ≠ 0 ∧ ‖algebraMap K M d‖ = ‖π‖ ^ ((N : ℕ) : ℤ) := by
    refine ⟨(p : K), hpK, ?_⟩
    rw [map_natCast, heM, zpow_natCast]
  have hcdvd : (c : ℤ) ∣ ((N : ℕ) : ℤ) := (hc _).mp hpmem
  have hc0 : c ≠ 0 := by
    intro h
    rw [h] at hcdvd
    simp only [Nat.cast_zero, zero_dvd_iff, Nat.cast_eq_zero] at hcdvd
    omega
  obtain ⟨πK, hπK0, hπKn⟩ := (hc ((c : ℕ) : ℤ)).mpr dvd_rfl
  have hc1 : (1 : ℤ) ≤ (c : ℤ) := by
    have : 1 ≤ c := Nat.one_le_iff_ne_zero.mpr hc0
    exact_mod_cast this
  refine ⟨πK, ?_, ?_, ?_⟩
  · rw [hπKn]
    exact zpow_pos hπ0 _
  · rw [hπKn]
    calc ‖π‖ ^ ((c : ℕ) : ℤ) < ‖π‖ ^ (0 : ℤ) :=
          zpow_lt_zpow_right_of_lt_one₀ hπ0 hπ1 (by omega)
      _ = 1 := zpow_zero _
  · intro a ha
    obtain ⟨m, hm⟩ := hvalK a ha
    have hmem : ∃ d : K, d ≠ 0 ∧ ‖algebraMap K M d‖ = ‖π‖ ^ ((n : ℤ) * m) := ⟨a, ha, hm⟩
    obtain ⟨m', hm'⟩ := (hc _).mp hmem
    refine ⟨m', ?_⟩
    rw [hm, hπKn, ← zpow_mul, hm']

end Helpers

/-! ## §2 ★★★★★★★`harith` の組み立て -/

section Assemble

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-- ★★★★★★★**`harith` をノルムの仮説だけから供給する**。

`JumpFromValueGroup.lean:272-276` の `harith` の 4 条件そのものを出す。

★★**逸脱 1 本だけ**: `hiso : ∀ z, ‖g z‖ = ‖z‖`を**仮説で受ける**。
`PureStepSetup.lean:283 norm_algEquiv_eq` は
`[NontriviallyNormedField k] [CompleteSpace k] [NormedAlgebra k M]` を要求するので、
★**素の `[Field K]` の設定では出ない**（具体層で `ℚ_p` を下に敷いてから供給する）。 -/
theorem harith_of_norm [FiniteDimensional K M] {p k e : ℕ} (hp : p.Prime) {π : M}
    (g : M ≃ₐ[K] M) (hg : orderOf g = p ^ (k + 1))
    (hiso : ∀ z : M, ‖g z‖ = ‖z‖) (he : 0 < e)
    (hnK : Module.finrank K M = p ^ (k + 1))
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ,
      ‖algebraMap K M a‖ = ‖π‖ ^ (((p ^ (k + 1) : ℕ) : ℤ) * m))
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹)
    (heM : ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e))
    (s : ℕ → (M →+* M)) (hsg : ∀ j, ∀ z : M, s j z = (g ^ p ^ j) z) :
    ∀ u : ℕ → ℤ, (∀ j, j ≤ k → ‖s j π - π‖ = ‖π‖ ^ (u j + 1)) →
      1 ≤ u 0 ∧ (∀ m, m < k → u m < u (m + 1)) ∧
        (∀ m, m < k → (p : ℤ) ^ (m + 1) ∣ u (m + 1) - u m) ∧
        ((p : ℤ) - 1) * u k ≤ (p : ℤ) ^ (k + 1) * (e : ℤ) := by
  classical
  haveI : Fact p.Prime := ⟨hp⟩
  intro u hu
  have hq0 : 0 < p ^ (k + 1) := pow_pos hp.pos _
  have hqe : 0 < p ^ (k + 1) * e := Nat.mul_pos hq0 he
  have hq1 : 1 < p ^ (k + 1) := Nat.one_lt_pow (by omega) hp.one_lt
  have hπ0 : 0 < ‖π‖ := TotallyRamified.norm_pos_of_normp (p := p) (k := k) he hnormp heM
  have hπ1 : ‖π‖ < 1 := norm_lt_one_of_normp hp.one_lt hqe hnormp heM
  have hpM : ‖(p : M)‖ < 1 := by
    rw [hnormp]
    have h1 : (1 : ℝ) < (p : ℝ) := by exact_mod_cast hp.one_lt
    rw [inv_lt_one_iff₀]
    right; exact h1
  have hbr : ∀ j, j ≤ k → ‖(g ^ p ^ j) π - π‖ = ‖π‖ ^ (u j + 1) := by
    intro j hj
    rw [← hsg j π]
    exact hu j hj
  -- (1)
  have hbr0 : ‖g π - π‖ = ‖π‖ ^ (u 0 + 1) := by
    have h := hbr 0 (Nat.zero_le k)
    simpa using h
  have h1 : 1 ≤ u 0 :=
    WildBreakPow.one_le_jump_zero_pow hp g hg hiso hπ0 hπ1 hnK hvalK he heM hbr0
  -- (2)
  have hstep : ∀ m, m < k → 1 ≤ u m → u m < u (m + 1) := fun m hm hum =>
    JumpMono.jump_lt_succ g hiso hπ0 hπ1 hnK hvalK hpM hum (hbr m (le_of_lt hm)) (hbr (m + 1) hm)
  have hone : ∀ m, m ≤ k → 1 ≤ u m := by
    intro m
    induction m with
    | zero => intro _; exact h1
    | succ r ih =>
        intro hr
        have hrk : r < k := by omega
        have := hstep r hrk (ih (by omega))
        omega
  have h2 : ∀ m, m < k → u m < u (m + 1) := fun m hm => hstep m hm (hone m (le_of_lt hm))
  have hmono_lt : ∀ b : ℕ, b ≤ k → ∀ a : ℕ, a < b → u a < u b := by
    intro b
    induction b with
    | zero => intro _ a ha; omega
    | succ r ih =>
        intro hr a ha
        rcases Nat.lt_or_ge a r with har | har
        · exact lt_trans (ih (by omega) a har) (h2 r (by omega))
        · have hae : a = r := by omega
          subst hae
          exact h2 a (by omega)
  have hmono_le : ∀ a b : ℕ, a ≤ b → b ≤ k → u a ≤ u b := by
    intro a b hab hbk
    rcases Nat.eq_or_lt_of_le hab with hh | hh
    · rw [hh]
    · exact le_of_lt (hmono_lt b hbk a hh)
  -- 共通の道具
  have htop : ∀ σ : M ≃ₐ[K] M, σ ∈ Subgroup.zpowers g :=
    forall_mem_zpowers_of_orderOf_eq_finrank g hg hnK
  have hisoAll : ∀ (σ : M ≃ₐ[K] M) (z : M), ‖σ z‖ = ‖z‖ := by
    intro σ z
    obtain ⟨j, hj⟩ := htop σ
    rw [← hj]
    exact norm_zpow_apply g hiso j z
  have hval : ∀ y : M, y ≠ 0 → ∃ m : ℤ, ‖y‖ = ‖π‖ ^ m :=
    fun y hy => TotallyRamified.exists_zpow_norm hπ0 (ne_of_lt hπ1) hnK hvalK hy
  -- (4)
  have h4 : ((p : ℤ) - 1) * u k ≤ (p : ℤ) ^ (k + 1) * (e : ℤ) :=
    WildBreakPow.hupper_of_totallyRamified_pow hp g hg hisoAll hπ0 hπ1 hnK hvalK hval heM
      (hbr k le_rfl)
  -- (3)
  have h3 : ∀ m, m < k → (p : ℤ) ^ (m + 1) ∣ u (m + 1) - u m := by
    intro m hm
    set uN : ℕ → ℕ := fun j => if j ≤ k then (u j).toNat else (u k).toNat + (j - k) with huN
    have huNle : ∀ j, j ≤ k → uN j = (u j).toNat := by intro j hj; simp [huN, hj]
    have huNcast : ∀ j, j ≤ k → ((uN j : ℕ) : ℤ) = u j := by
      intro j hj
      rw [huNle j hj]
      exact Int.toNat_of_nonneg (by have := hone j hj; omega)
    have hsmono : ∀ a b : ℕ, a < b → uN a < uN b := by
      intro a b hab
      by_cases hbk : b ≤ k
      · have hak : a ≤ k := by omega
        rw [huNle a hak, huNle b hbk]
        have hlt : u a < u b := hmono_lt b hbk a hab
        have hua : 0 ≤ u a := by have := hone a hak; omega
        omega
      · have hbig : uN b = (u k).toNat + (b - k) := by simp [huN, hbk]
        by_cases hak : a ≤ k
        · rw [huNle a hak, hbig]
          have hak' : u a ≤ u k := hmono_le a k hak le_rfl
          have hua : 0 ≤ u a := by have := hone a hak; omega
          omega
        · have hbiga : uN a = (u k).toNat + (a - k) := by simp [huN, hak]
          rw [hbiga, hbig]
          omega
    have hjump : ∀ j : ℕ, j ≤ k → ‖(g ^ p ^ j) π - π‖ = ‖π‖ ^ (uN j + 1) := by
      intro j hj
      rw [hbr j hj, ← zpow_natCast ‖π‖ (uN j + 1)]
      congr 1
      have := huNcast j hj
      push_cast
      omega
    have hu0N : 1 ≤ uN 0 := by
      have := huNcast 0 (Nat.zero_le k)
      omega
    have hπmem : π ∈ IntegerNorm.integerSubring M := le_of_lt hπ1
    have hpKne : (p : K) ≠ 0 := by
      intro hc
      have : ((p : ℕ) : M) = 0 := by
        have := congrArg (algebraMap K M) hc
        rwa [map_natCast, map_zero] at this
      rw [this, norm_zero] at hpM
      rw [this, norm_zero] at hnormp
      have hp0 : (0 : ℝ) < (p : ℝ) := by exact_mod_cast hp.pos
      have : (0 : ℝ) < ((p : ℝ))⁻¹ := inv_pos.mpr hp0
      linarith [hnormp]
    obtain ⟨πK, hK0, hK1, hvalKK⟩ :=
      exists_base_uniformizer (n := p ^ (k + 1)) (N := p ^ (k + 1) * e) hπ0 hπ1 hqe hvalK
        hpKne heM
    have hKmem : πK ∈ IntegerNorm.baseIntegerSubring K M := le_of_lt hK1
    have hdvd := IntegerNorm.dvd_sub_jump_of_norm (n := p ^ (k + 1)) (m := m)
      hq1 hp (fun σ z => hisoAll σ z) hπ0 hπ1 hnK hvalK hval hπmem hjump hsmono hu0N
      hg htop hK0 hK1 hvalKK hKmem hpM (by omega)
    rw [huNcast (m + 1) (by omega), huNcast m (by omega)] at hdvd
    exact hdvd
  exact ⟨h1, h2, h3, h4⟩

end Assemble

/-! ## §3 ★★★★★★★★出口から `harith` と `htop` が消えた -/

section Exit

open Finset GainedTowerModel GainedTowerModel.GaloisTower TotallyRamified

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

/-- ★★★★★★★★**`harith` が消えた出口**。

`JumpFromValueGroup.exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_jumpFree`
（`JumpFromValueGroup.lean:265`）の `harith` と `htop` を両方とも内部で供給する。

* `harith` ← §2 `harith_of_norm`
* `htop`  ← `TotallyRamifiedLayer.adjoin_eq_top_of_valK`

★★**残ったノルムの仮説は `hiso` / `hnormp` / `heM` / `hvalK` / `hnK` の 5 本**。
うち `hiso` だけが本波で**新たに必要になった**もので、
具体層（`ℚ_p` を下に敷いた場合）では
`PureStepSetup.lean:283 norm_algEquiv_eq` が供給する。 -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay_of_norm
    {K : Type*} [Field K] [Algebra K M] [FiniteDimensional K M]
    {p : ℕ} [Fact p.Prime] {e k : ℕ} {π : M}
    (g : M ≃ₐ[K] M) (hg : orderOf g = p ^ (k + 1))
    (hiso : ∀ z : M, ‖g z‖ = ‖z‖)
    (τ : M →+ M) (hτ : ∀ z : M, τ z = g z)
    (s : ℕ → (M →+* M)) (hsg : ∀ j, ∀ z : M, s j z = (g ^ p ^ j) z)
    (he : 0 < e)
    (hnK : Module.finrank K M = p ^ (k + 1))
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ,
      ‖algebraMap K M a‖ = ‖π‖ ^ (((p ^ (k + 1) : ℕ) : ℤ) * m))
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹)
    (heM : ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e)) (x : M) :
    ∃ y : twr g p k, ‖x - algebraMap (twr g p k) M y‖
      ≤ (∏ j ∈ Finset.Icc 1 (k + 1), axDecay p j) * ‖τ x - x‖ := by
  have hp : p.Prime := Fact.out
  have hπ0 : 0 < ‖π‖ := TotallyRamified.norm_pos_of_normp (p := p) (k := k) he hnormp heM
  have hqe : 0 < p ^ (k + 1) * e := Nat.mul_pos (pow_pos hp.pos _) he
  have hπ1 : ‖π‖ < 1 := norm_lt_one_of_normp hp.one_lt hqe hnormp heM
  have htop : Algebra.adjoin K ({π} : Set M) = ⊤ :=
    TotallyRamifiedLayer.adjoin_eq_top_of_valK hπ0 (ne_of_lt hπ1) hnK hvalK
  exact JumpFromValueGroup.exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_jumpFree
    g hg τ hτ s hsg he
    (harith_of_norm hp g hg hiso he hnK hvalK hnormp heM s hsg)
    hnK hvalK hnormp heM htop x

end Exit



/-! ## `.src` と 公理 -/

def harith_of_norm.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_norm_sub_algebraMap_le_prod_axDecay_of_norm.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

#print axioms norm_lt_one_of_normp
#print axioms norm_zpow_apply
#print axioms forall_mem_zpowers_of_orderOf_eq_finrank
#print axioms exists_base_uniformizer
#print axioms harith_of_norm
#print axioms exists_norm_sub_algebraMap_le_prod_axDecay_of_norm

end HarithAssembly

end ABC3.Found.PGC
