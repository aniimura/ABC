import ABC3.Found.PGC.HasseArfCongruenceNorm

/-!
# [pGC] ★★★★★★★`harith` (1) を `p^{k+1}` 次に上げた —— `1 ≤ u 0`

## ★★まず測ったこと（持ち場の問い：(1) と (4) のどちらが安いか）

★**(1) が安い**。理由は本ファイルの§1 にある通り、
`k = 0` の証明（`WildBreakLowerBound.lean:367`）で
★**素数性を使っているのは最後の 1 行だけ**だったからである。

★★持ち場の指摘は当たっていた: `JumpStrictMono.lean:333 one_le_jump_of_zero` は
`hu0 : 1 ≤ u 0` を ★**仮説のまま**受けており、(1) を証明していない。
同ファイルの docstring の断定「(1) からすべての段で `1 ≤ u m`」は真だが、
★**その入口である (1) 自体は `k = 0`（次数 `p`）しか無かった**。

## どこが変わったか（★**1 行だけ**）

`WildBreak.norm_sub_le_sq_of_totallyRamified`（`WildBreakLowerBound.lean:367`）の段取りは

1. `w := gπ/π`、`‖w‖ = 1`。telescoping で `∏_{j<q} g^j w = 1`。
2. `w` を `K` の元 `A` で `‖w − A‖ ≤ ‖π‖` まで近似（`exists_sub_algebraMap_norm_le`）。
3. `A^q·∏(1+η_j) = 1` から `‖A^q − 1‖ ≤ ‖π‖`、値群で `≤ ‖π‖^q`。
4. ★**`‖A − 1‖^q ≤ ‖π‖^q`**。
5. `‖w − 1‖ ≤ ‖π‖` ⇒ `‖gπ − π‖ ≤ ‖π‖²`。

★素数性を使っているのは **4 だけ**（`norm_sub_one_pow_le`、二項係数 `p ∣ C(p,l)`）。
1・2・3・5 は `p` を `q = p^{k+1}` に置き換えるだけで通る。

★★**4 を `p^{k+1}` に上げるのに二項係数は要らない**。
`Nat.Prime.dvd_choose_pow`（`Data/Nat/Multiplicity.lean:260`、実在する）を使う道もあるが、
★素数の場合を **`r` 回反復する方が安い**：

    ‖A−1‖^{p^{r+1}} = (‖A−1‖^p)^{p^r} ≤ max(‖A^p−1‖, ‖p‖)^{p^r} ≤ c

（`‖A^p−1‖^{p^r} ≤ c` は帰納法、`‖p‖^{p^r} ≤ ‖p‖ ≤ c` は `‖p‖ ≤ 1`）。

## 出したもの

| 宣言 | 内容 |
|---|---|
| ★★★`norm_sub_one_pow_pow_le` | ★抽象核（純ノルム）`‖A^{p^r}−1‖ ≤ c ⇒ ‖A−1‖^{p^r} ≤ c` |
| ★★`norm_sub_le_sq_of_totallyRamified_pow` | `‖gπ − π‖ ≤ ‖π‖²`（次数 `p^{k+1}`） |
| ★★★`one_le_jump_zero_pow` | **`harith` (1)** `1 ≤ u 0` |

## ★★残りは (4) だけ。どこで止まるかを測った（`file:line`）

(4) `(p−1)·u k ≤ p^{k+1}·e` の `k = 0` 版は
`WildBreakUpperBound.lean:228 sub_one_mul_le_of_totallyRamified`。
その中身は `RamificationJumpBound.norm_natCast_le_pow_of_splits` であり、
★**決定的な入力は `hbreak`**（同 `:242-256`）:

    hbreak : ∀ a, (minpoly K π).map(algebraMap).IsRoot a → a ≠ π → ‖π − a‖ = ‖π‖^{i+1}

すなわち★**共役がすべて同じ距離**であることを使っている。
これは `|G| = p`（素数位数）だから成り立つのであって、
★★**`q = p^{k+1}` では偽である**——共役 `g^j π` の距離は `‖π‖^{u_{v_p(j)}+1}` で
`j` ごとに違う。実際 `ℤ₃(ζ₈₁)` では `u = (2, 8, 26)` で 3 通りの距離が出る。

★検算した代替経路と、それが足りない理由:
`h := g^{p^k}`（位数 `p`）に対して `∏_{j<p}(1 + h^j b) = 1`（`b := hπ/π − 1`）を
展開すると `‖p‖ ≤ ‖π‖^{u k}`、すなわち `u k ≤ p^{k+1}·e` は出るが、
★**係数 `p−1` が出ない**。`p ≥ 3` では不十分である。
係数 `p−1` は「差式（全共役の積）」からしか出ない。

⇒ ★★★**(4) の一般化に残るのは 1 ノードだけ**:
★**中間体 `L := M^{⟨g^{p^k}⟩}` を立て、`M/L` に `k = 0` 版を適用すること**。
そのとき `e' = v_L(p) = p^k·e` なので
`(p−1)·u k ≤ p·e' = p^{k+1}·e` となりちょうど目標に合う（★手で検算済み）。
要るのは `finrank L M = p`・`IsGalois L M`・`L` の値群（`‖algebraMap L M a‖ = ‖π‖^{p·m}`）であり、
既知の危険は `lean-idioms.md #59`（中間体の 2 層をまたぐ `rfl`）と `#69`、`#296`。
★**本波はそこに降りていない。**

★数値の検算（本波で手で）: `ℚ₃(ζ₈₁)/ℚ₃(ζ₃)` は `p=3, k=2, q=27, e=2, u=(2,8,26)`。
(4) は `2·26 = 52 ≤ 27·2 = 54` ★真。`k=1` の例は `1·8 = 8 ≤ 4·2 = 8` ★等号。
★上で述べた弱い経路が与える `u k ≤ q·e` は `26 ≤ 54`、`8 ≤ 8` で、
★**`k=1` の例では偶然一致するが `p=3` の例では弱い**。

## 逸脱の記録

1. §2 は `WildBreakLowerBound.lean:367` の段取りを `p → q` で写している。
   ★既存ファイルは読み取り専用なので一般化版を新しく書いた（元は消していない）。
   ★差分は `norm_sub_one_pow_le` → `norm_sub_one_pow_pow_le` の 1 行である。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
-/

namespace ABC3.Found.PGC

namespace WildBreakPow

open Finset WildBreak

/-! ## §1 抽象核（純ノルム）—— `x^{p^r} = 1` のノルム版 -/

section Core

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

/-- ★★★**抽象核** —— `‖A‖ = 1`・`‖A^{p^r} − 1‖ ≤ c`・`‖(p:M)‖ ≤ c ≤ 1` なら
`‖A − 1‖^{p^r} ≤ c`。

★`WildBreak.norm_sub_one_pow_le`（素数 `p` の場合、`WildBreakLowerBound.lean:319`）の
**`r` 段の反復**。★二項係数を `p^r` に持ち上げる必要はない
（`Nat.Prime.dvd_choose_pow` は使わない）—— 素数の場合を `r` 回使えばよい。

★これが「標数 `p` で `x^{p^r} = 1 ⇒ x = 1`」のノルム版である。
★体の分岐も剰余体も出てこない。 -/
theorem norm_sub_one_pow_pow_le {p : ℕ} (hp : p.Prime) {c : ℝ} (hc1 : c ≤ 1)
    (hc : ‖(p : M)‖ ≤ c) :
    ∀ (r : ℕ) (A : M), ‖A‖ = 1 → ‖A ^ p ^ r - 1‖ ≤ c → ‖A - 1‖ ^ p ^ r ≤ c := by
  have hpM1 : ‖(p : M)‖ ≤ 1 := le_trans hc hc1
  intro r
  induction r with
  | zero => intro A _ h; simpa using h
  | succ r ih =>
      intro A hA h
      have hpow : ∀ z : M, (z ^ p) ^ p ^ r = z ^ p ^ (r + 1) := by
        intro z; rw [← pow_mul, ← pow_succ']
      have hApn : ‖A ^ p‖ = 1 := by rw [norm_pow, hA, one_pow]
      have hh : ‖(A ^ p) ^ p ^ r - 1‖ ≤ c := by rw [hpow]; exact h
      have hIH : ‖A ^ p - 1‖ ^ p ^ r ≤ c := ih (A ^ p) hApn hh
      have hstep : ‖A - 1‖ ^ p ≤ max ‖A ^ p - 1‖ ‖(p : M)‖ :=
        norm_sub_one_pow_le hp hA (le_max_right _ _) (le_max_left _ _)
      have hd0 : (0 : ℝ) ≤ max ‖A ^ p - 1‖ ‖(p : M)‖ := le_trans (norm_nonneg _) (le_max_left _ _)
      have hr1 : 1 ≤ p ^ r := Nat.one_le_pow _ _ hp.pos
      have hdr : (max ‖A ^ p - 1‖ ‖(p : M)‖) ^ p ^ r ≤ c := by
        rcases max_cases ‖A ^ p - 1‖ ‖(p : M)‖ with ⟨he, _⟩ | ⟨he, _⟩
        · rw [he]; exact hIH
        · rw [he]
          calc ‖(p : M)‖ ^ p ^ r ≤ ‖(p : M)‖ ^ 1 :=
                pow_le_pow_of_le_one (norm_nonneg _) hpM1 hr1
            _ = ‖(p : M)‖ := pow_one _
            _ ≤ c := hc
      calc ‖A - 1‖ ^ p ^ (r + 1) = (‖A - 1‖ ^ p) ^ p ^ r := by rw [← pow_mul, ← pow_succ']
        _ ≤ (max ‖A ^ p - 1‖ ‖(p : M)‖) ^ p ^ r :=
            pow_le_pow_left₀ (by positivity) hstep _
        _ ≤ c := hdr

end Core

/-! ## §2 ★★★`harith` (1) の `p^{k+1}` 次版 -/

section Assembly

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-- ★★★★★★**次数 `p^{k+1}` の全分岐巡回拡大は暴分岐である** —— `harith` (1)。

`‖g π − π‖ ≤ ‖π‖²`、すなわち `‖g π − π‖ = ‖π‖^{u₀+1}` と書いたときの `1 ≤ u₀`。

★★`WildBreak.norm_sub_le_sq_of_totallyRamified`（`WildBreakLowerBound.lean:367`）は
`g^p = 1` と `[M:K] = p` で **`k = 0` に固定**されていた。
本定理は `q := p^{k+1}` に上げる。

★**何が変わったかは 1 行だけ**: 最後の段で
`WildBreak.norm_sub_one_pow_le`（素数 `p`）の代わりに
`norm_sub_one_pow_pow_le`（§1、`p^r`）を使う。
★他の段はすべて `p` を `q` に置き換えるだけで通る（素数性を使っていない）。
★二項係数を `p^{k+1}` に持ち上げる必要はない。 -/
theorem norm_sub_le_sq_of_totallyRamified_pow [FiniteDimensional K M] {p k e : ℕ}
    (hp : p.Prime) {π : M} (g : M ≃ₐ[K] M) (hgq : g ^ p ^ (k + 1) = 1)
    (hiso : ∀ z : M, ‖g z‖ = ‖z‖)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = p ^ (k + 1))
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ,
      ‖algebraMap K M a‖ = ‖π‖ ^ (((p ^ (k + 1) : ℕ) : ℤ) * m))
    (he : 0 < e) (heM : ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e)) :
    ‖g π - π‖ ≤ ‖π‖ ^ 2 := by
  classical
  set q : ℕ := p ^ (k + 1) with hq
  have hq1 : 1 < q := by
    rw [hq]
    exact Nat.one_lt_pow (by omega) hp.one_lt
  have hq0 : 0 < q := by omega
  have hπne : π ≠ 0 := norm_pos_iff.mp hπ0
  have hgj : ∀ (j : ℕ) (z : M), ‖(g ^ j) z‖ = ‖z‖ := norm_pow_apply g hiso
  have hfne : ∀ j : ℕ, (g ^ j) π ≠ 0 := fun j => norm_pos_iff.mp (by rw [hgj]; exact hπ0)
  set w : M := g π / π with hwdef
  have hw1 : ‖w‖ = 1 := by rw [hwdef, norm_div, hiso]; exact div_self (ne_of_gt hπ0)
  have hgw : ∀ j : ℕ, (g ^ j) w = (g ^ (j + 1)) π / (g ^ j) π := by
    intro j
    rw [hwdef, map_div₀, ← AlgEquiv.mul_apply, ← pow_succ]
  have hprod : ∏ j ∈ Finset.range q, (g ^ j) w = 1 := by
    simp_rw [hgw]
    rw [prod_telescope g hπne hfne q, hq, hgq]
    simpa using div_self hπne
  obtain ⟨a, ha⟩ := exists_sub_algebraMap_norm_le (n := q) hq1 hπ0 hπ1 (by rw [hn, hq])
    (by rw [hq]; exact hvalK) (le_of_eq hw1)
  set A : M := algebraMap K M a with hAdef
  have hAnorm : ‖A‖ = 1 := by
    have hlt : ‖A - w‖ < ‖w‖ := by
      rw [norm_sub_rev, hw1]
      exact lt_of_le_of_lt ha hπ1
    have hmax := IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm (x := w) (y := A - w)
      (ne_of_gt hlt)
    rw [show w + (A - w) = A by ring] at hmax
    rw [hmax, max_eq_left (le_of_lt hlt), hw1]
  have hAne : A ≠ 0 := norm_pos_iff.mp (by rw [hAnorm]; norm_num)
  have hgA : ∀ j : ℕ, (g ^ j) A = A := fun j => (g ^ j).commutes a
  set η : ℕ → M := fun j => (g ^ j) w / A - 1 with hηdef
  have hfactor : ∀ j : ℕ, (g ^ j) w = A * (1 + η j) := by
    intro j
    rw [hηdef]
    field_simp
    ring
  have hηnorm : ∀ j : ℕ, ‖η j‖ ≤ ‖π‖ := by
    intro j
    have h1 : η j = ((g ^ j) w - A) / A := by rw [hηdef]; field_simp
    have h2 : (g ^ j) w - A = (g ^ j) (w - A) := by rw [map_sub, hgA j]
    rw [h1, norm_div, hAnorm, div_one, h2, hgj]
    exact ha
  have hprodA : A ^ q * ∏ j ∈ Finset.range q, (1 + η j) = 1 := by
    have hcong : ∏ j ∈ Finset.range q, (g ^ j) w = ∏ j ∈ Finset.range q, (A * (1 + η j)) :=
      Finset.prod_congr rfl (fun j _ => hfactor j)
    rw [Finset.prod_mul_distrib, Finset.prod_const, Finset.card_range] at hcong
    rw [← hcong]
    exact hprod
  have hP : ‖(∏ j ∈ Finset.range q, (1 + η j)) - 1‖ ≤ ‖π‖ :=
    norm_prod_one_add_sub_one_le _ _ (le_of_lt hπ0) (le_of_lt hπ1) (fun j _ => hηnorm j)
  have hAqm1 : ‖A ^ q - 1‖ ≤ ‖π‖ := by
    have hAq : ‖A ^ q‖ = 1 := by rw [norm_pow, hAnorm, one_pow]
    have hid : A ^ q - 1 = A ^ q * (1 - ∏ j ∈ Finset.range q, (1 + η j)) := by
      rw [mul_sub, mul_one, hprodA]
    rw [hid, norm_mul, hAq, one_mul, ← norm_neg]
    simpa using hP
  have hAqK : A ^ q - 1 = algebraMap K M (a ^ q - 1) := by
    rw [map_sub, map_pow, map_one, hAdef]
  have hAqm1' : ‖A ^ q - 1‖ ≤ ‖π‖ ^ q := by
    rcases eq_or_ne (a ^ q - 1) 0 with h0 | h0
    · rw [hAqK, h0, map_zero, norm_zero]
      exact pow_nonneg (norm_nonneg π) q
    · obtain ⟨m, hm⟩ := hvalK (a ^ q - 1) h0
      rw [hAqK, hm]
      have hlt1 : ‖π‖ ^ ((q : ℤ) * m) < 1 := by
        rw [← hm, ← hAqK]
        exact lt_of_le_of_lt hAqm1 hπ1
      have hpos : 0 < (q : ℤ) * m := (zpow_lt_one_iff_right_of_lt_one₀ hπ0 hπ1).mp hlt1
      have hqZ : (0 : ℤ) < (q : ℤ) := by exact_mod_cast hq0
      have hm1 : 1 ≤ m := by nlinarith
      have hqle : (q : ℤ) ≤ (q : ℤ) * m := by nlinarith
      calc ‖π‖ ^ ((q : ℤ) * m) ≤ ‖π‖ ^ ((q : ℕ) : ℤ) :=
            zpow_le_zpow_right_of_le_one₀ hπ0 (le_of_lt hπ1) hqle
        _ = ‖π‖ ^ q := zpow_natCast _ _
  have hcp : ‖(p : M)‖ ≤ ‖π‖ ^ q := by
    rw [heM]
    exact pow_le_pow_of_le_one (norm_nonneg π) (le_of_lt hπ1) (Nat.le_mul_of_pos_right q he)
  have hc1 : ‖π‖ ^ q ≤ 1 := pow_le_one₀ (norm_nonneg π) (le_of_lt hπ1)
  have hAB : ‖A - 1‖ ^ q ≤ ‖π‖ ^ q :=
    norm_sub_one_pow_pow_le hp hc1 hcp (k + 1) A hAnorm hAqm1'
  have hA1 : ‖A - 1‖ ≤ ‖π‖ := le_of_pow_le_pow_left₀ (by omega) (norm_nonneg π) hAB
  have hw1' : ‖w - 1‖ ≤ ‖π‖ := by
    have hsplit : w - 1 = (w - A) + (A - 1) := by ring
    rw [hsplit]
    exact (IsUltrametricDist.norm_add_le_max _ _).trans (max_le ha hA1)
  have hfin : g π - π = π * (w - 1) := by
    rw [hwdef]
    field_simp
  rw [hfin, norm_mul, sq]
  exact mul_le_mul_of_nonneg_left hw1' (norm_nonneg π)

end Assembly

/-! ## §3 `harith` (1) の字面そのもの -/

section Harith

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-- ★★★★★★★**`harith` (1)** —— `1 ≤ u 0`。

`JumpFromValueGroup.lean:273` の第 1 項そのものの形である
（`u : ℕ → ℤ`、`‖s 0 π − π‖ = ‖π‖^{u 0 + 1}`、`s 0 = g^{p^0} = g`）。

★★**前波までの状態**: `JumpStrictMono.lean:333 one_le_jump_of_zero` は
`hu0 : 1 ≤ u 0` を ★**仮説のまま**受けていた（そこから全段の `1 ≤ u m` と狭義単調を出していた）。
`k = 0` では `WildBreakLowerBound.lean:509` が
`norm_sub_le_sq_of_totallyRamified` で埋めていたが、
★**`p^{k+1}` 次版は無かった**。本定理がそれである。 -/
theorem one_le_jump_zero_pow [FiniteDimensional K M] {p k e : ℕ} (hp : p.Prime)
    {π : M} {u0 : ℤ} (g : M ≃ₐ[K] M) (hg : orderOf g = p ^ (k + 1))
    (hiso : ∀ z : M, ‖g z‖ = ‖z‖)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = p ^ (k + 1))
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ,
      ‖algebraMap K M a‖ = ‖π‖ ^ (((p ^ (k + 1) : ℕ) : ℤ) * m))
    (he : 0 < e) (heM : ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e))
    (hbr : ‖g π - π‖ = ‖π‖ ^ (u0 + 1)) : 1 ≤ u0 := by
  have hgq : g ^ p ^ (k + 1) = 1 := by rw [← hg]; exact pow_orderOf_eq_one g
  have hle := norm_sub_le_sq_of_totallyRamified_pow hp g hgq hiso hπ0 hπ1 hn hvalK he heM
  by_contra hcon
  rw [not_le] at hcon
  have hz2 : ‖π‖ ^ (2 : ℤ) = ‖π‖ ^ (2 : ℕ) := by
    rw [show (2 : ℤ) = ((2 : ℕ) : ℤ) by norm_num, zpow_natCast]
  have hlt2 : ‖π‖ ^ (2 : ℤ) < ‖π‖ ^ (u0 + 1) :=
    zpow_lt_zpow_right_of_lt_one₀ hπ0 hπ1 (by omega)
  rw [hz2, ← hbr] at hlt2
  linarith

end Harith



/-! ## `.src` と 公理 -/

def norm_sub_le_sq_of_totallyRamified_pow.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def one_le_jump_zero_pow.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

#print axioms norm_sub_one_pow_pow_le
#print axioms norm_sub_le_sq_of_totallyRamified_pow
#print axioms one_le_jump_zero_pow

end WildBreakPow

end ABC3.Found.PGC
