/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.SixExp.LatticeHouse

/-!
# 数え上げ —— 外挿の `hgap` はいつ満たされるか

★`ResearchPaper/frdi-decomposition.json` の `sixexp` チェーン。
最終目標は [FrdI] `Lemma 6.5, (ii)`(six exponentials theorem)。

## ★超越性証明が閉じる理由

`Extrapolation.lean` の `extrapolation_induction` が要求する `hgap` は

  `M · (1/2)^{|T k|} < (H k)^{-(d-1)}`

である。ここで `|T k|` は**それまでに得た零点の個数**、`H k` は**格子点での値の
house の上界**(`LatticeHouse.lean` で `|s|·A^{L·N}·C` と評価済み)。

★同値な形は `M · H^e < 2^{|T|}`(`gap_iff`)。つまり

* **左辺**は house の増大 —— `A^{L·N}` すなわち `N` について**指数の底が固定で 1 乗**
* **右辺**は零点の個数 —— 格子の箱を `[0,N)³` に取れば `2^{N³}`

★★**3 乗は 1 乗に勝つ**(`exists_gap_index`)。これが「外挿がいつか閉じる」ことの
理由であり、six exponentials theorem で `y` が **3 個**必要な理由でもある
(`y` が 2 個だと `2^{N²}` と `A^{L·N}` の競争になり、定数の取り方に依存してしまう)。

★`exists_gap_of_growth` が `hgap` をそのままの形で供給する。
-/

namespace ABC3.Found.SixExp

/-- ★★外挿の `hgap` の**同値な書き換え** —— 「解析側 < Liouville 側」は
`M · H^e < 2^{|T|}` と同じ。★零点の個数 `|T|` が指数で効き、house は `H^e` で効く。 -/
theorem gap_iff {M H : ℝ} (hH : 0 < H) (e T : ℕ) :
    M * (1/2 : ℝ) ^ T < (H ^ e)⁻¹ ↔ M * H ^ e < 2 ^ T := by
  have hHe : (0:ℝ) < H ^ e := by positivity
  have h2 : (0:ℝ) < 2 ^ T := by positivity
  rw [one_div, inv_pow, ← div_eq_mul_inv, ← one_div (H ^ e), div_lt_div_iff₀ h2 hHe, one_mul]

/-- ★★★**数え上げの核** —— 零点の個数が `N³`、house が `A^{bN}` で増えるとき、
`N` が大きければ **`2^{N³}` が `K·A^{bN}` を追い越す**。

★これが「外挿がいつか閉じる」ことの理由である —— **3 乗は 1 乗に勝つ**。

★証明は 3 段:
`2^N > A^b`(`N` 大)から `A^{bN} < (2^N)^N = 2^{N²}`、
`K < 2^{N²}`、そして `N ≥ 2` なら `2·N² ≤ N³` なので `2^{N²}·2^{N²} ≤ 2^{N³}`。 -/
theorem exists_gap_index (K A : ℝ) (hK : 0 < K) (hA1 : 1 ≤ A) (b : ℕ) :
    ∃ N₀ : ℕ, ∀ N ≥ N₀, K * A ^ (b * N) < 2 ^ (N ^ 3) := by
  obtain ⟨n, hn⟩ := pow_unbounded_of_one_lt (A ^ b) (by norm_num : (1:ℝ) < 2)
  obtain ⟨mK, hmK⟩ := pow_unbounded_of_one_lt K (by norm_num : (1:ℝ) < 2)
  refine ⟨max 2 (max n mK), fun N hN => ?_⟩
  have hN2 : 2 ≤ N := le_trans (le_max_left _ _) hN
  have hNn : n ≤ N := le_trans (le_trans (le_max_left _ _) (le_max_right 2 _)) hN
  have hNm : mK ≤ N := le_trans (le_trans (le_max_right _ _) (le_max_right 2 _)) hN
  have hA0 : (0:ℝ) ≤ A := le_trans zero_le_one hA1
  have h1 : (A:ℝ) ^ b < 2 ^ N := lt_of_lt_of_le hn (pow_le_pow_right₀ (by norm_num) hNn)
  have h2 : (A:ℝ) ^ (b * N) < 2 ^ (N * N) := by
    have hlt : ((A:ℝ) ^ b) ^ N < ((2:ℝ) ^ N) ^ N :=
      pow_lt_pow_left₀ h1 (by positivity) (by omega)
    calc (A:ℝ) ^ (b * N) = ((A:ℝ) ^ b) ^ N := by rw [← pow_mul]
      _ < ((2:ℝ) ^ N) ^ N := hlt
      _ = 2 ^ (N * N) := by rw [← pow_mul]
  have h3 : K < 2 ^ (N * N) := by
    refine lt_of_lt_of_le hmK (pow_le_pow_right₀ (by norm_num) ?_)
    calc mK ≤ N := hNm
      _ ≤ N * N := Nat.le_mul_of_pos_left N (by omega)
  have h4 : (2:ℝ) ^ (N * N) * 2 ^ (N * N) ≤ 2 ^ (N ^ 3) := by
    rw [← pow_add]
    refine pow_le_pow_right₀ (by norm_num) ?_
    have hrw : N * N + N * N = 2 * (N * N) := by ring
    rw [hrw]
    calc 2 * (N * N) ≤ N * (N * N) := Nat.mul_le_mul_right _ hN2
      _ = N ^ 3 := by ring
  calc K * A ^ (b * N) < 2 ^ (N * N) * 2 ^ (N * N) :=
        mul_lt_mul'' h3 h2 (le_of_lt hK) (by positivity)
    _ ≤ 2 ^ (N ^ 3) := h4

/-- ★★★★**外挿の `hgap` は `N` を大きく取れば必ず満たされる**。

`|T k| = N³`(格子の箱)、`H k = |s|·A^{L·N}·C`(house の上界)とすると、
`hgap` は `M·H^d < 2^{N³}` と同値(`gap_iff`)で、
★**3 乗は 1 乗に勝つ**から `N` を大きくすれば成り立つ(`exists_gap_index`)。 -/
theorem exists_gap_of_growth {M A C cs : ℝ} (hM : 0 < M) (hA1 : 1 ≤ A) (hC : 0 < C)
    (hcs : 0 < cs) (d L : ℕ) :
    ∃ N₀ : ℕ, ∀ N ≥ N₀,
      M * (1/2 : ℝ) ^ (N ^ 3) < ((cs * (A ^ (L * N) * C)) ^ d)⁻¹ := by
  obtain ⟨N₀, hN₀⟩ := exists_gap_index (M * (cs ^ d * C ^ d)) A (by positivity) hA1 (L * d)
  refine ⟨N₀, fun N hN => ?_⟩
  have hH : 0 < cs * (A ^ (L * N) * C) := by positivity
  rw [gap_iff hH]
  have hexp : (cs * (A ^ (L * N) * C)) ^ d = (cs ^ d * C ^ d) * A ^ ((L * d) * N) := by
    have h1 : (L * d) * N = L * N * d := by ring
    rw [h1, mul_pow, mul_pow, ← pow_mul]
    ring
  rw [hexp, ← mul_assoc]
  exact hN₀ N hN

end ABC3.Found.SixExp
