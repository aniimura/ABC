import ABC3.Found.PGC.WildBreakUpperBound

/-!
# [pGC] `harith` の **`hult`(狭義単調 `u m < u (m+1)`)は定理である**

`GainedTowerModel.exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic` /
`JumpFromValueGroup.…_of_cyclic_jumpFree` が要求する `harith` は 4 条件

    (1) 1 ≤ u 0                         (2) ∀ m < k, u m < u (m+1)
    (3) ∀ m < k, p^{m+1} ∣ u(m+1) − u m  (4) (p−1)·u k ≤ p^{k+1}·e

からなる。★**本ファイルは (2) を証明して仮説から落とす。**

## ★★測ったこと(先に結論)

* (2) は★**剰余体も different も Hasse–Arf も使わずに**、超距離と基底展開だけで出る。
  中身は Serre の `v(σ z − z) ≥ v(z) + i_σ`(§2)と、
  `h^pπ/π = ∏_{j<p}(1 + h^j b)` の 1 次の項の評価(§1・§3)である。
* ★`1 ≤ u m` が要る(`t ≥ 1` を使って `‖π‖^t < 1` にする)。
  ★したがって (2) は **(1) から従う**(`one_le_jump_of_zero`)。
  ⇒ ★★`harith` の 4 条件は実質 **(1) (3) (4) の 3 本**に減った。
* ★残る本丸は **(3) Hasse–Arf の合同**である。本ファイルでは**触っていない**。

## 証明の筋(★圧縮で消えても復元できるようここに書く)

`h := g^{p^m}`、`b := hπ/π − 1`(`‖b‖ = ‖π‖^t`, `t = u m ≥ 1`)と置く。telescoping で

    h^p π / π = ∏_{j<p} h^j(hπ/π) = ∏_{j<p} (1 + h^j b).

§1 `norm_prod_one_add_sub_one_sub_sum_le` で 1 次の項を取り出すと残りは `‖b‖²`。
1 次の項は `Σ_{j<p} h^j b = p·b + Σ_{j<p}(h^j b − b)` と分け、
§2 `norm_pow_sub_apply_le_mul`(`‖h^j z − z‖ ≤ ‖z‖·‖π‖^t`)で
`‖Σ_j h^j b‖ ≤ ‖b‖·max(‖p‖, ‖π‖^t)`。`‖b‖² = ‖b‖·‖π‖^t` なので合わせて

    ‖h^pπ/π − 1‖ ≤ ‖b‖·max(‖p‖, ‖π‖^t) < ‖b‖      (‖p‖ < 1 かつ t ≥ 1)

両辺 `‖π‖` 倍して `‖h^pπ − π‖ < ‖hπ − π‖`。★これが `u m < u(m+1)`。

§2 の中身は「`z = Σ_{l<n} c_l π^l`(素元冪が基底)で各項のノルムが `‖z‖` 以下」
(`exists_coeff_norm_le`、`WildBreak.norm_le_norm_sum_of_pairwise_ne'` を使う)と、
`‖(hπ)^l − π^l‖ ≤ ‖π‖^{l−1}·‖hπ − π‖` だけである。

## 節の構成

* §1 抽象核(超距離) —— `‖∏(1+η) − 1 − Σ η‖ ≤ r²`
* §2 `exists_coeff_norm_le` / ★`norm_sub_apply_le_mul` / `norm_pow_sub_apply_le_mul`
* §3 ★★`norm_pow_prime_sub_lt`(`‖h^pπ − π‖ < ‖hπ − π‖`)/ ★★★`jump_lt_succ`(＝ `hult`)/
  `one_le_jump_of_zero`((1) から (2) が全段で従う)

## ★在庫の測定(コマンドを残す)

```
grep -rn "theorem norm_pow_sub_pow_le" lean/ABC3/Found/PGC/*.lean
  → ★3 本ある。`GainedTowerStep.lean:159` の
    `norm_pow_sub_pow_le (hπ0 : 0 < ‖π‖) (hw : ‖w‖ ≤ ‖π‖) (m) : ‖w^m − π^m‖ ≤ ‖π‖^{m−1}·‖w − π‖`
    が本ファイルの要。★自作しかけて `has already been declared` で気づいた(#158 の手)。
  → `norm_sum_le_of_forall_le` も同ファイルに在る。★私が `WildBreakLowerBound` §1 で
    自作したものは**重複**であった(引数の並びが違うだけ)。★「無い」と判定する前に
    `grep -rn "theorem <名前>" lean/ABC3/Found/PGC/*.lean` を打つこと。
#check @Finset.sum_sub_distrib / @Commute.mul_geom_sum₂ / @pow_le_pow_of_le_one → 在る
```

## 逸脱の記録

1. `hpM : ‖(p : M)‖ < 1` を仮説に置いた(`hnormp` から出るが、本ファイルは正規化を要らない形にした)。
2. `hiso`(等長)は仮説。`ConcreteDegPFree.hiso_g2` の道(スペクトルノルムなら `rfl` 1 行)で供給される。
3. `GainedTowerStep` / `WildBreak*` は**読むだけ**で 1 行も書き換えていない。
-/

namespace ABC3.Found.PGC

namespace JumpMono

open Finset WildBreak

/-! ## §1 抽象核(超距離) —— 積の 1 次の項を取り出す -/

section Ultra

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

/-- ★★**抽象核** —— `‖η i‖ ≤ r ≤ 1` なら
`‖∏(1 + η i) − 1 − Σ η i‖ ≤ r²`(★1 次の項を取り出した残りは 2 次)。 -/
theorem norm_prod_one_add_sub_one_sub_sum_le {ι : Type*} (s : Finset ι) (η : ι → M) {r : ℝ}
    (hr0 : 0 ≤ r) (hr1 : r ≤ 1) (h : ∀ i ∈ s, ‖η i‖ ≤ r) :
    ‖(∏ i ∈ s, (1 + η i)) - 1 - ∑ i ∈ s, η i‖ ≤ r * r := by
  classical
  induction s using Finset.induction_on with
  | empty => simpa using mul_nonneg hr0 hr0
  | insert a s ha ih =>
      have hha : ‖η a‖ ≤ r := h a (Finset.mem_insert_self a s)
      have hsub : ∀ i ∈ s, ‖η i‖ ≤ r := fun i hi => h i (Finset.mem_insert_of_mem hi)
      have hP : ‖(∏ i ∈ s, (1 + η i)) - 1‖ ≤ r :=
        norm_prod_one_add_sub_one_le s η hr0 hr1 hsub
      have hIH : ‖(∏ i ∈ s, (1 + η i)) - 1 - ∑ i ∈ s, η i‖ ≤ r * r := ih hsub
      rw [Finset.prod_insert ha, Finset.sum_insert ha,
        show (1 + η a) * (∏ i ∈ s, (1 + η i)) - 1 - (η a + ∑ i ∈ s, η i)
          = ((∏ i ∈ s, (1 + η i)) - 1 - ∑ i ∈ s, η i)
            + η a * ((∏ i ∈ s, (1 + η i)) - 1) by ring]
      refine (IsUltrametricDist.norm_add_le_max _ _).trans (max_le hIH ?_)
      rw [norm_mul]
      exact mul_le_mul hha hP (norm_nonneg _) hr0

end Ultra

/-! ## §2 「`h` は元を `‖π‖^t` だけしか動かさない」 -/

section Expansion

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-- 素元冪の基底展開。★各項のノルムは全体のノルム以下(§1 の核)。

★`WildBreak.exists_sub_algebraMap_norm_le` の前半と同じ構成だが、
そちらは `c 0` だけを取り出す形なので、ここでは**係数の族ごと**返す形に切り直した。 -/
theorem exists_coeff_norm_le [FiniteDimensional K M] {π : M} {n : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (z : M) :
    ∃ c : Fin n → K, (∑ l, c l • π ^ (l : ℕ)) = z ∧ ∀ l, ‖c l • π ^ (l : ℕ)‖ ≤ ‖z‖ := by
  classical
  have hπne : ‖π‖ ≠ 1 := ne_of_lt hπ1
  have hli : LinearIndependent K (fun l : Fin n => π ^ (l : ℕ)) := by
    refine TotallyRamified.linearIndependent_of_ne_mod hπ0 hπne hvalK _
      (fun l => ((l : ℕ) : ℤ)) ?_ ?_
    · intro l; simp [zpow_natCast]
    · intro i j hij
      exact TotallyRamifiedLayer.not_dvd_sub_of_lt i.isLt j.isLt (fun h => hij (Fin.ext h))
  have hzmem : z ∈ Submodule.span K (Set.range (fun l : Fin n => π ^ (l : ℕ))) := by
    rw [hli.span_eq_top_of_card_eq_finrank' (by simp [hn])]; exact Submodule.mem_top
  obtain ⟨c, hc⟩ := (Submodule.mem_span_range_iff_exists_fun K).mp hzmem
  refine ⟨c, hc, ?_⟩
  set x : Fin n → M := fun l => c l • π ^ (l : ℕ) with hxdef
  have hxval : ∀ l : Fin n, x l = c l • π ^ (l : ℕ) := fun _ => rfl
  have hsum : ∑ l, x l = z := hc
  have hexp : ∀ l : Fin n, c l ≠ 0 →
      ∃ m : ℤ, ‖x l‖ = ‖π‖ ^ ((n : ℤ) * m + ((l : ℕ) : ℤ)) := by
    intro l hl
    obtain ⟨m, hm⟩ := hvalK (c l) hl
    exact ⟨m, by rw [hxval, Algebra.smul_def, norm_mul, hm, norm_pow,
      ← zpow_natCast ‖π‖ (l : ℕ), ← zpow_add₀ (ne_of_gt hπ0)]⟩
  have hczero : ∀ l : Fin n, x l ≠ 0 → c l ≠ 0 := by
    intro l hx hcl
    exact hx (by rw [hxval, hcl, zero_smul])
  have hpair : ∀ i j : Fin n, i ≠ j → x i ≠ 0 → x j ≠ 0 → ‖x i‖ ≠ ‖x j‖ := by
    intro i j hij hi hj
    obtain ⟨mi, hmi⟩ := hexp i (hczero i hi)
    obtain ⟨mj, hmj⟩ := hexp j (hczero j hj)
    rw [hmi, hmj]
    intro hcon
    have h2 : (n : ℤ) * mi + ((i : ℕ) : ℤ) = (n : ℤ) * mj + ((j : ℕ) : ℤ) :=
      (zpow_right_inj₀ hπ0 hπne).mp hcon
    refine TotallyRamifiedLayer.not_dvd_sub_of_lt i.isLt j.isLt
      (fun h => hij (Fin.ext h)) ⟨mj - mi, ?_⟩
    linarith [mul_sub (n : ℤ) mj mi]
  intro l
  have := norm_le_norm_sum_of_pairwise_ne' x hpair l
  rwa [hsum] at this

/-- ★★★★**`h` は元を `‖π‖^t` だけしか動かさない** —— `‖h z − z‖ ≤ ‖z‖·‖π‖^t`。

★これが Serre の `v(σz − z) ≥ v(z) + i_σ` のノルム版であり、★狭義単調(§3)の要である。
★証明は基底展開(`exists_coeff_norm_le`)と超距離だけ。

★在庫の測定: `‖w^m − π^m‖ ≤ ‖π‖^{m−1}·‖w − π‖` は
★**`GainedTowerStep.norm_pow_sub_pow_le`(159 行)に在った**(自作しかけて衝突で気づいた)。 -/
theorem norm_sub_apply_le_mul [FiniteDimensional K M] {π : M} {n t : ℕ}
    (h : M ≃ₐ[K] M) (hiso : ∀ w : M, ‖h w‖ = ‖w‖)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (hbr : ‖h π - π‖ = ‖π‖ ^ (t + 1)) (z : M) :
    ‖h z - z‖ ≤ ‖z‖ * ‖π‖ ^ t := by
  classical
  obtain ⟨c, hc, hnorm⟩ := exists_coeff_norm_le hπ0 hπ1 hn hvalK z
  have hzz : h z - z = ∑ l : Fin n, (c l • ((h π) ^ (l : ℕ) - π ^ (l : ℕ))) := by
    rw [← hc, map_sum, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl (fun l _ => ?_)
    rw [Algebra.smul_def, Algebra.smul_def, map_mul, AlgEquiv.commutes,
      map_pow, mul_sub]
  rw [hzz]
  refine norm_sum_le_of_forall_le _ _ (by positivity) (fun l _ => ?_)
  rcases Nat.eq_zero_or_pos (l : ℕ) with hl0 | hl1
  · rw [hl0]
    simp [mul_nonneg (norm_nonneg z) (by positivity : (0:ℝ) ≤ ‖π‖ ^ t)]
  · have hstep : ‖(h π) ^ (l : ℕ) - π ^ (l : ℕ)‖ ≤ ‖π‖ ^ ((l : ℕ) - 1) * ‖π‖ ^ (t + 1) := by
      rw [← hbr]
      exact ABC3.Found.PGC.GainedTowerStep.norm_pow_sub_pow_le hπ0 (le_of_eq (hiso π)) _
    have hcomb : ‖π‖ ^ ((l : ℕ) - 1) * ‖π‖ ^ (t + 1) = ‖π‖ ^ (l : ℕ) * ‖π‖ ^ t := by
      rw [← pow_add, ← pow_add]
      congr 1
      omega
    rw [Algebra.smul_def, norm_mul]
    calc ‖algebraMap K M (c l)‖ * ‖(h π) ^ (l : ℕ) - π ^ (l : ℕ)‖
        ≤ ‖algebraMap K M (c l)‖ * (‖π‖ ^ (l : ℕ) * ‖π‖ ^ t) := by
          rw [← hcomb]; exact mul_le_mul_of_nonneg_left hstep (norm_nonneg _)
      _ = (‖algebraMap K M (c l)‖ * ‖π‖ ^ (l : ℕ)) * ‖π‖ ^ t := by ring
      _ ≤ ‖z‖ * ‖π‖ ^ t := by
          refine mul_le_mul_of_nonneg_right ?_ (by positivity)
          have := hnorm l
          rwa [Algebra.smul_def, norm_mul, norm_pow] at this


/-- 冪でも同じ: `‖h^j z − z‖ ≤ ‖z‖·‖π‖^t`。 -/
theorem norm_pow_sub_apply_le_mul [FiniteDimensional K M] {π : M} {n t : ℕ}
    (h : M ≃ₐ[K] M) (hiso : ∀ w : M, ‖h w‖ = ‖w‖)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (hbr : ‖h π - π‖ = ‖π‖ ^ (t + 1)) (j : ℕ) (z : M) :
    ‖(h ^ j) z - z‖ ≤ ‖z‖ * ‖π‖ ^ t := by
  induction j with
  | zero => simpa using mul_nonneg (norm_nonneg z) (by positivity : (0:ℝ) ≤ ‖π‖ ^ t)
  | succ j ih =>
      have hrw : (h ^ (j + 1)) z - z = h ((h ^ j) z - z) + (h z - z) := by
        rw [map_sub, ← AlgEquiv.mul_apply, ← pow_succ']
        ring
      rw [hrw]
      refine (IsUltrametricDist.norm_add_le_max _ _).trans (max_le ?_ ?_)
      · rw [hiso]; exact ih
      · exact norm_sub_apply_le_mul h hiso hπ0 hπ1 hn hvalK hbr z

end Expansion

/-! ## §3 ★★★★★★狭義単調 `u m < u (m+1)` -/

section Mono

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-- ★★★★★★**跳びは `p` 乗で真に伸びる** —— `‖h^p π − π‖ < ‖h π − π‖`。

★これが `harith` の `hult`(狭義単調 `u m < u (m+1)`)の中身である。

筋: `u := hπ/π = 1 + b`(`‖b‖ = ‖π‖^t`)と telescoping で
`h^pπ/π = ∏_{j<p}(1 + h^j b)`。§1 で 1 次の項を取り出すと残りは `‖b‖²`、
1 次の項 `Σ_j h^j b = p·b + Σ_j (h^j b − b)` は §2 で `‖b‖·max(‖p‖, ‖π‖^t)` 以下。
`t ≥ 1` と `‖p‖ < 1` から両方 `< ‖b‖`。★剰余体も different も出てこない。 -/
theorem norm_pow_prime_sub_lt [FiniteDimensional K M] {p n t : ℕ} {π : M}
    (h : M ≃ₐ[K] M) (hiso : ∀ w : M, ‖h w‖ = ‖w‖)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (hpM : ‖(p : M)‖ < 1) (ht : 1 ≤ t)
    (hbr : ‖h π - π‖ = ‖π‖ ^ (t + 1)) :
    ‖(h ^ p) π - π‖ < ‖h π - π‖ := by
  classical
  have hπne : π ≠ 0 := norm_pos_iff.mp hπ0
  set b : M := h π / π - 1 with hbdef
  have hb : π * b = h π - π := by rw [hbdef]; field_simp
  have hbnorm : ‖b‖ = ‖π‖ ^ t := by
    have h1 : ‖π‖ * ‖b‖ = ‖π‖ ^ t * ‖π‖ := by
      rw [← norm_mul, hb, hbr, pow_succ]
    have h2 : ‖π‖ * ‖b‖ = ‖π‖ * ‖π‖ ^ t := by rw [h1]; ring
    exact mul_left_cancel₀ (ne_of_gt hπ0) h2
  have hbt : ‖π‖ ^ t ≤ ‖π‖ := by
    calc ‖π‖ ^ t ≤ ‖π‖ ^ 1 := pow_le_pow_of_le_one (le_of_lt hπ0) (le_of_lt hπ1) ht
      _ = ‖π‖ := pow_one _
  have hb1 : ‖b‖ ≤ 1 := by rw [hbnorm]; linarith
  have hgj : ∀ (j : ℕ) (w : M), ‖(h ^ j) w‖ = ‖w‖ := WildBreak.norm_pow_apply h hiso
  have hfne : ∀ j : ℕ, (h ^ j) π ≠ 0 := fun j => norm_pos_iff.mp (by rw [hgj]; exact hπ0)
  have hgu : ∀ j : ℕ, (h ^ j) (h π / π) = 1 + (h ^ j) b := by
    intro j
    rw [hbdef, map_sub, map_one]
    ring
  have hprod : ∏ j ∈ Finset.range p, (1 + (h ^ j) b) = (h ^ p) π / π := by
    rw [← WildBreak.prod_telescope h hπne hfne p]
    refine Finset.prod_congr rfl (fun j _ => ?_)
    rw [← hgu j, map_div₀, ← AlgEquiv.mul_apply, ← pow_succ]
  set c : ℝ := max ‖(p : M)‖ (‖π‖ ^ t) with hcdef
  have hc0 : 0 ≤ c := le_trans (norm_nonneg _) (le_max_left _ _)
  have hclt : c < 1 := max_lt hpM (lt_of_le_of_lt hbt hπ1)
  have hηnorm : ∀ j ∈ Finset.range p, ‖(h ^ j) b‖ ≤ ‖b‖ := fun j _ => le_of_eq (hgj j b)
  have hlin : ‖(∏ j ∈ Finset.range p, (1 + (h ^ j) b)) - 1
      - ∑ j ∈ Finset.range p, (h ^ j) b‖ ≤ ‖b‖ * ‖b‖ :=
    norm_prod_one_add_sub_one_sub_sum_le _ _ (norm_nonneg b) hb1 hηnorm
  have hsumeq : ∑ j ∈ Finset.range p, (h ^ j) b
      = (p : M) * b + ∑ j ∈ Finset.range p, ((h ^ j) b - b) := by
    rw [Finset.sum_sub_distrib, Finset.sum_const, Finset.card_range, nsmul_eq_mul]
    ring
  have hsumbd : ‖∑ j ∈ Finset.range p, (h ^ j) b‖ ≤ ‖b‖ * c := by
    rw [hsumeq]
    refine (IsUltrametricDist.norm_add_le_max _ _).trans (max_le ?_ ?_)
    · rw [norm_mul, mul_comm]
      exact mul_le_mul_of_nonneg_left (le_max_left _ _) (norm_nonneg b)
    · refine le_trans (norm_sum_le_of_forall_le _ _ (by positivity) (fun j _ =>
        norm_pow_sub_apply_le_mul h hiso hπ0 hπ1 hn hvalK hbr j b)) ?_
      exact mul_le_mul_of_nonneg_left (le_max_right _ _) (norm_nonneg b)
  have hPbd : ‖(h ^ p) π / π - 1‖ ≤ ‖b‖ * c := by
    rw [← hprod,
      show (∏ j ∈ Finset.range p, (1 + (h ^ j) b)) - 1
        = ((∏ j ∈ Finset.range p, (1 + (h ^ j) b)) - 1 - ∑ j ∈ Finset.range p, (h ^ j) b)
          + ∑ j ∈ Finset.range p, (h ^ j) b by ring]
    refine (IsUltrametricDist.norm_add_le_max _ _).trans (max_le ?_ hsumbd)
    refine le_trans hlin ?_
    refine mul_le_mul_of_nonneg_left ?_ (norm_nonneg b)
    rw [hbnorm]
    exact le_max_right _ _
  have hbpos : 0 < ‖b‖ := by rw [hbnorm]; exact pow_pos hπ0 t
  have hfin : (h ^ p) π - π = π * ((h ^ p) π / π - 1) := by field_simp
  rw [hfin, norm_mul, ← hb, norm_mul]
  refine mul_lt_mul_of_pos_left ?_ hπ0
  calc ‖(h ^ p) π / π - 1‖ ≤ ‖b‖ * c := hPbd
    _ < ‖b‖ * 1 := by exact mul_lt_mul_of_pos_left hclt hbpos
    _ = ‖b‖ := mul_one _


/-- ★★★★★★★**`harith` の `hult`(狭義単調)は定理である**。

`u m` を `‖(g^{p^m}) π − π‖ = ‖π‖^{u m + 1}` で定めるとき、`1 ≤ u m` なら `u m < u (m+1)`。 -/
theorem jump_lt_succ [FiniteDimensional K M] {p n : ℕ} {π : M} {u : ℕ → ℤ} {m : ℕ}
    (g : M ≃ₐ[K] M) (hiso : ∀ w : M, ‖g w‖ = ‖w‖)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (hpM : ‖(p : M)‖ < 1) (hum : 1 ≤ u m)
    (hbm : ‖(g ^ p ^ m) π - π‖ = ‖π‖ ^ (u m + 1))
    (hbm1 : ‖(g ^ p ^ (m + 1)) π - π‖ = ‖π‖ ^ (u (m + 1) + 1)) :
    u m < u (m + 1) := by
  have hisoh : ∀ w : M, ‖(g ^ p ^ m) w‖ = ‖w‖ := WildBreak.norm_pow_apply g hiso (p ^ m)
  have hti : (((u m).toNat : ℕ) : ℤ) = u m := Int.toNat_of_nonneg (by omega)
  have hbr : ‖(g ^ p ^ m) π - π‖ = ‖π‖ ^ ((u m).toNat + 1) := by
    rw [hbm, ← zpow_natCast ‖π‖ ((u m).toNat + 1)]
    congr 1
    push_cast [hti]
    ring
  have hpow : (g ^ p ^ m) ^ p = g ^ p ^ (m + 1) := by
    rw [← pow_mul, ← pow_succ]
  have hlt := norm_pow_prime_sub_lt (p := p) (g ^ p ^ m) hisoh hπ0 hπ1 hn hvalK hpM
    (by omega : 1 ≤ (u m).toNat) hbr
  rw [hpow, hbm1, hbm] at hlt
  by_contra hcon
  rw [not_lt] at hcon
  exact absurd hlt (not_lt.mpr (zpow_le_zpow_right_of_le_one₀ hπ0 (le_of_lt hπ1) (by omega)))

/-- ★★系: `1 ≤ u 0` から**すべての段で** `1 ≤ u m` かつ狭義単調。 -/
theorem one_le_jump_of_zero [FiniteDimensional K M] {p n : ℕ} {π : M} {u : ℕ → ℤ}
    (g : M ≃ₐ[K] M) (hiso : ∀ w : M, ‖g w‖ = ‖w‖)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (hpM : ‖(p : M)‖ < 1) (hu0 : 1 ≤ u 0)
    (hb : ∀ j : ℕ, ‖(g ^ p ^ j) π - π‖ = ‖π‖ ^ (u j + 1)) :
    ∀ m : ℕ, 1 ≤ u m ∧ u m < u (m + 1) := by
  intro m
  induction m with
  | zero =>
      exact ⟨hu0, jump_lt_succ g hiso hπ0 hπ1 hn hvalK hpM hu0 (hb 0) (hb 1)⟩
  | succ m ih =>
      have h1 : 1 ≤ u (m + 1) := by have ha := ih.1; have hb2 := ih.2; omega
      exact ⟨h1, jump_lt_succ g hiso hπ0 hπ1 hn hvalK hpM h1 (hb (m + 1)) (hb (m + 2))⟩
end Mono

/-! ## §4 `.src` と 公理 -/

def norm_sub_apply_le_mul.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_pow_prime_sub_lt.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def jump_lt_succ.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

#print axioms norm_prod_one_add_sub_one_sub_sum_le
#print axioms exists_coeff_norm_le
#print axioms norm_sub_apply_le_mul
#print axioms norm_pow_sub_apply_le_mul
#print axioms norm_pow_prime_sub_lt
#print axioms jump_lt_succ
#print axioms one_le_jump_of_zero

end JumpMono

end ABC3.Found.PGC
