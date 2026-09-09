import ABC3.Found.PGC.LayerDegreeIsP

/-!
# [pGC] ★★`p ∣ j` の桁で等号が破れる量を測った —— 損失因子 `max ‖(j:L)‖ (‖u−π‖/‖π‖)`

## 持ち場（前波で「要るのは」とした点）

前波の私の言葉（逐語）:

> ⇒ ★★**要るのは「`p` 進単数でない桁（`p ∣ j`）を持つ展開でも使えるノルムの下界」**です。

## ★本波が足したもの

木が持っているのは 2 つで、どちらも `j` に依らない:

* `CyclicJumpNorm.lean:189 norm_pow_sub_pow_eq` —— ★**等号**、ただし `‖(k:L)‖ = 1` が要る
* `CyclicJumpNorm.lean:162 norm_pow_succ_sub_pow_succ_le` —— `‖u^{i+1} − π^{i+1}‖ ≤ ‖u−π‖·‖π‖^i`
  （★`j` に依らない粗い上界。`GainedTowerStep.lean:159 norm_pow_sub_pow_le` も同型）

★**その中間が無い**（`grep -rn "norm_pow_sub_pow" lean/ABC3/Found/PGC/*.lean` で確認、
`max (‖((j` を含む行は **0 件**）。§1 がそれである:

  `‖u^{j+2} − π^{j+2}‖ ≤ ‖u − π‖ · max (‖(j+2 : L)‖·‖π‖^{j+1}) (‖u−π‖·‖π‖^{j})`

★`‖(j+2:L)‖ = 1` なら右辺は `‖u−π‖·‖π‖^{j+1}` に潰れ、木の粗い上界と一致する。
★★`p ∣ j+2` なら `‖(j+2:L)‖ ≤ 1/p` なので ★**真に小さくなる**（§2）。

## ★★測って分かったこと —— 損失因子は `max ‖(n:L)‖ (‖u−π‖/‖π‖)`

`n := j+2` として、桁 `n` の寄与は木の粗い上界 `‖u−π‖‖π‖^{n−1}` の

  `max (‖(n:L)‖) (‖u−π‖/‖π‖)`

倍以下である。★2 つの項の意味:

* `‖(n:L)‖` —— ★**形式微分 `n π^{n−1}` の係数**（`p ∣ n` で `≤ 1/p`）
* `‖u−π‖/‖π‖` —— ★**2 次以上の項**（`< 1` は `‖u−π‖ < ‖π‖` から）

⇒ ★★`‖(n:L)‖ < 1` なら**必ず**等号が破れる（§2 `norm_pow_sub_pow_lt_of_natCast_lt_one`）。
★これが「`hchar` が `n > p` で偽」（前波 `LayerDegreeIsP`）の**量的な姿**である。

## ★これが `j₀` の議論とどう繋がるか（★本波では繋いでいない）

第 1120 波台の `MaxMinIndex` / `SlotResidue` / `MainPartCoeffs` は
「桁のノルムが相異なるので最大を与える `j₀` が**一意**」を使っていた。
★`p ∤ j₀` なら §1 の右辺が `‖u−π‖‖π‖^{j₀−1}` に潰れ、`j₀` の桁が主役のまま残る。
★★しかし `p ∣ j₀` の場合（例: `x = π^p`）は §2 により主役の寄与が真に小さくなり、
★**下界が出ない**。⇒ ★残っているのは「`p ∣ j₀` の場合の扱い」ちょうど 1 つ。

★★**本ファイルはその扱いを与えていない。** ★与えられるとも書いていない。
★**破れる量を測っただけ**である。

## 逸脱の記録

- §1・§2 は超距離ノルム体だけで、★分岐・付値・Galois・桁展開の語彙が 1 語も出ない。
  ★`p` すら出てこない（`‖(n:L)‖ < 1` としか言わない）。
- §1 の証明は `geom_sum₂_mul`（mathlib）＋ 木の `norm_pow_succ_sub_pow_succ_le` ＋
  超距離の和の評価だけ。★新しい道具は使っていない。
- ★指数は `j + 2`（`n ≥ 2`）で書いた。★`ℕ` の切り捨て引き算（`n−1`, `n−2`）を
  避けるためで、`n ≤ 1` では主張が退化する（`n = 1` は木の等号がそのまま使える）。
-/

namespace ABC3.Found.PGC

namespace DigitLossFactor

/-! ## §1 ★抽象核 —— 桁 `n` の寄与の `n` 依存の上界 -/

section Kernel

variable {L : Type*} [NormedField L] [IsUltrametricDist L]

/-- ★★★**桁 `n = j+2` の寄与の、`n` に依存した上界**（木には無い形）。

  `‖u^{j+2} − π^{j+2}‖ ≤ ‖u − π‖ · max (‖(j+2 : L)‖·‖π‖^{j+1}) (‖u−π‖·‖π‖^{j})`

第 1 項は**形式微分** `n π^{n−1}`、第 2 項は **2 次以上**の寄与である。

★`‖(j+2:L)‖ = 1` なら第 1 項が勝ち、木の
`CyclicJumpNorm.norm_pow_succ_sub_pow_succ_le` と一致する。
★★`p ∣ j+2` なら第 1 項が `1/p` 倍以下になり、★**真に小さくなる**（§2）。 -/
theorem norm_pow_sub_pow_le_max {u π : L} (h : ‖u - π‖ < ‖π‖) (j : ℕ) :
    ‖u ^ (j + 2) - π ^ (j + 2)‖
      ≤ ‖u - π‖ * max (‖((j + 2 : ℕ) : L)‖ * ‖π‖ ^ (j + 1)) (‖u - π‖ * ‖π‖ ^ j) := by
  classical
  set S : L := ∑ m ∈ Finset.range (j + 2), u ^ m * π ^ (j + 2 - 1 - m) with hS
  have hgeom : S * (u - π) = u ^ (j + 2) - π ^ (j + 2) := geom_sum₂_mul u π (j + 2)
  -- `S` を「形式微分の項」と「残り」に分ける
  have hpi : ∀ m ∈ Finset.range (j + 2),
      π ^ (j + 2 - 1 - m) * (u ^ m - π ^ m) = u ^ m * π ^ (j + 2 - 1 - m) - π ^ (j + 1) := by
    intro m hm
    have hmlt : m < j + 2 := Finset.mem_range.mp hm
    have hpow : π ^ (j + 2 - 1 - m) * π ^ m = π ^ (j + 1) := by
      rw [← pow_add]
      congr 1
      omega
    rw [mul_sub, hpow]
    ring
  have hsplit : S = ((j + 2 : ℕ) : L) * π ^ (j + 1)
      + ∑ m ∈ Finset.range (j + 2), π ^ (j + 2 - 1 - m) * (u ^ m - π ^ m) := by
    rw [Finset.sum_congr rfl hpi, Finset.sum_sub_distrib, Finset.sum_const,
      Finset.card_range, hS]
    simp [nsmul_eq_mul]
  -- 残りの評価
  have hrest : ‖∑ m ∈ Finset.range (j + 2), π ^ (j + 2 - 1 - m) * (u ^ m - π ^ m)‖
      ≤ ‖u - π‖ * ‖π‖ ^ j := by
    refine IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg (by positivity) ?_
    intro m hm
    have hmlt : m < j + 2 := Finset.mem_range.mp hm
    rcases Nat.eq_zero_or_pos m with hm0 | hm1
    · subst hm0
      simp
      positivity
    · obtain ⟨i, rfl⟩ : ∃ i, m = i + 1 := ⟨m - 1, by omega⟩
      have hstep := norm_pow_succ_sub_pow_succ_le h i
      rw [norm_mul, norm_pow]
      have hmul : ‖π‖ ^ (j + 2 - 1 - (i + 1)) * (‖u - π‖ * ‖π‖ ^ i)
          = ‖u - π‖ * ‖π‖ ^ j := by
        rw [show ‖π‖ ^ (j + 2 - 1 - (i + 1)) * (‖u - π‖ * ‖π‖ ^ i)
            = ‖u - π‖ * (‖π‖ ^ (j + 2 - 1 - (i + 1)) * ‖π‖ ^ i) by ring, ← pow_add]
        congr 2
        omega
      calc ‖π‖ ^ (j + 2 - 1 - (i + 1)) * ‖u ^ (i + 1) - π ^ (i + 1)‖
          ≤ ‖π‖ ^ (j + 2 - 1 - (i + 1)) * (‖u - π‖ * ‖π‖ ^ i) :=
            mul_le_mul_of_nonneg_left hstep (by positivity)
        _ = ‖u - π‖ * ‖π‖ ^ j := hmul
  -- 合成
  have hSle : ‖S‖ ≤ max (‖((j + 2 : ℕ) : L)‖ * ‖π‖ ^ (j + 1)) (‖u - π‖ * ‖π‖ ^ j) := by
    rw [hsplit]
    refine (IsUltrametricDist.norm_add_le_max _ _).trans (max_le_max ?_ hrest)
    rw [norm_mul, norm_pow]
  calc ‖u ^ (j + 2) - π ^ (j + 2)‖ = ‖S‖ * ‖u - π‖ := by rw [← hgeom, norm_mul]
    _ ≤ max (‖((j + 2 : ℕ) : L)‖ * ‖π‖ ^ (j + 1)) (‖u - π‖ * ‖π‖ ^ j) * ‖u - π‖ :=
        mul_le_mul_of_nonneg_right hSle (norm_nonneg _)
    _ = ‖u - π‖ * max (‖((j + 2 : ℕ) : L)‖ * ‖π‖ ^ (j + 1)) (‖u - π‖ * ‖π‖ ^ j) := by ring

end Kernel

/-! ## §2 ★★`‖(n:L)‖ < 1` なら等号は**必ず**破れる -/

section Strict

variable {L : Type*} [NormedField L] [IsUltrametricDist L]

/-- ★★★**`p ∣ n` の桁では木の上界が真に達成されない**。

`‖(j+2 : L)‖ < 1`（具体層では `p ∣ j+2`）なら

  `‖u^{j+2} − π^{j+2}‖ < ‖u − π‖ · ‖π‖^{j+1}`

★これが `LayerDegreeIsP`（前波）の「`hchar` は `n > p` で偽」の**量的な姿**である。
★★`x` の主役の桁 `j₀` が `p ∣ j₀` だと、この分だけ `‖x^u − x‖` の下界が落ちる。 -/
theorem norm_pow_sub_pow_lt_of_natCast_lt_one {u π : L} (hne : u ≠ π)
    (h : ‖u - π‖ < ‖π‖) (j : ℕ) (hn : ‖((j + 2 : ℕ) : L)‖ < 1) :
    ‖u ^ (j + 2) - π ^ (j + 2)‖ < ‖u - π‖ * ‖π‖ ^ (j + 1) := by
  have hd0 : 0 < ‖u - π‖ := by
    rw [norm_pos_iff]
    exact sub_ne_zero_of_ne hne
  have hπ0 : 0 < ‖π‖ := lt_of_le_of_lt (norm_nonneg _) h
  have hmax : max (‖((j + 2 : ℕ) : L)‖ * ‖π‖ ^ (j + 1)) (‖u - π‖ * ‖π‖ ^ j)
      < ‖π‖ ^ (j + 1) := by
    refine max_lt ?_ ?_
    · have : (0 : ℝ) < ‖π‖ ^ (j + 1) := by positivity
      nlinarith
    · have hstep : ‖u - π‖ * ‖π‖ ^ j < ‖π‖ * ‖π‖ ^ j :=
        mul_lt_mul_of_pos_right h (by positivity)
      calc ‖u - π‖ * ‖π‖ ^ j < ‖π‖ * ‖π‖ ^ j := hstep
        _ = ‖π‖ ^ (j + 1) := by ring
  calc ‖u ^ (j + 2) - π ^ (j + 2)‖
      ≤ ‖u - π‖ * max (‖((j + 2 : ℕ) : L)‖ * ‖π‖ ^ (j + 1)) (‖u - π‖ * ‖π‖ ^ j) :=
        norm_pow_sub_pow_le_max h j
    _ < ‖u - π‖ * ‖π‖ ^ (j + 1) := mul_lt_mul_of_pos_left hmax hd0

end Strict

/-! ## §3 使っている公理の一覧 -/

#print axioms norm_pow_sub_pow_le_max
#print axioms norm_pow_sub_pow_lt_of_natCast_lt_one

end DigitLossFactor

end ABC3.Found.PGC
