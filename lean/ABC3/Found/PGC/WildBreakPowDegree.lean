import ABC3.Found.PGC.HasseArfCongruenceNorm

/-!
# [pGC] ★★★★★★★`harith` の (1) と (4) を `p^{k+1}` 次に上げた

★これで `JumpFromValueGroup.lean:272-276` の `harith` の 4 条件はすべて
`p^{k+1}` 次で揃った（(2) は `JumpStrictMono.lean:308`、
(3) は `HasseArfCongruenceNorm.lean dvd_sub_jump_of_norm`）。

## ★★持ち場の問い：(1) と (4) のどちらが安いか

★**(1) が安かった**。`k = 0` の証明（`WildBreakLowerBound.lean:367`）で
★**素数性を使っているのは 5 段のうち 4 段目の 1 行だけ**だったからである。

★★持ち場の指摘は当たっていた: `JumpStrictMono.lean:333 one_le_jump_of_zero` は
`hu0 : 1 ≤ u 0` を ★**仮説のまま**受けており、(1) を証明していない。
同ファイルの docstring の断定「(1) からすべての段で `1 ≤ u m`」は真だが、
★**その入口である (1) 自体は `k = 0`（次数 `p`）しか無かった**。

## (1) で変わったのは★**1 行だけ**

`norm_sub_le_sq_of_totallyRamified` の段取り:
1. `w := gπ/π`、telescoping で `∏_{j<q} g^j w = 1`。
2. `w` を `K` の元 `A` で `‖w − A‖ ≤ ‖π‖` まで近似。
3. `A^q·∏(1+η_j) = 1` から `‖A^q − 1‖ ≤ ‖π‖`、値群で `≤ ‖π‖^q`。
4. ★**`‖A − 1‖^q ≤ ‖π‖^q`**　← ここだけが素数性を使う。
5. `‖w − 1‖ ≤ ‖π‖` ⇒ `‖gπ − π‖ ≤ ‖π‖²`。

★★**4 を `p^{k+1}` に上げるのに二項係数は要らない**。
`Nat.Prime.dvd_choose_pow`（`Data/Nat/Multiplicity.lean:260`、実在する）を使う道もあるが、
★素数の場合を **`r` 回反復する方が安い**：

    ‖A−1‖^{p^{r+1}} = (‖A−1‖^p)^{p^r} ≤ max(‖A^p−1‖, ‖p‖)^{p^r} ≤ c.

## ★★★(4) も閉じた —— 前波の自分の見積もりを覚す

★★**偽だった見積もり**: 本波の途中まで私は
「(4) には中間体が要り、`lean-idioms.md` #59/#69 の危険区間に降りることになる」
と書いていた。★**これは偽である。**
中間体を `IntermediateField` ではなく
★**`FixedPoints.subfield ↥(Subgroup.zpowers h) M`（単なる `Subfield`）**で建てれば、
危険区間には入らない（#332「重いと記録された罠は層に固有」のもう 1 例）。

★★在庫の測定（すべて `inferInstance` で通った。これが決め手）:

```lean
example (h : M ≃ₐ[K] M) : MulSemiringAction ↥(Subgroup.zpowers h) M := inferInstance  -- ★在る
example (h : M ≃ₐ[K] M) : FaithfulSMul     ↥(Subgroup.zpowers h) M := inferInstance  -- ★在る
example (h) : Algebra          ↥(FixedPoints.subfield ↥(Subgroup.zpowers h) M) M := inferInstance
example (h) : FiniteDimensional ↥(FixedPoints.subfield ↥(Subgroup.zpowers h) M) M := inferInstance
example (h) : IsGalois          ↥(FixedPoints.subfield ↥(Subgroup.zpowers h) M) M := inferInstance
example (h) (a) : algebraMap ↥(FixedPoints.subfield …) M a = (a : M) := rfl              -- ★rfl
example (h) (z) : (FixedPoints.toAlgAutMulEquiv … ⟨h, …⟩) z = h z := rfl                -- ★rfl
```

★最初 `MulSemiringAction.compHom` で `letI` して `FaithfulSMul` の合成に失敗したが
（`failed to synthesize FaithfulSMul (↥(Subgroup.zpowers h)) M`）、
★**そもそも両方とも instance として存在していた**ので `letI` 自体が不要だった。

## ★★(4) の鍵：中間体の値群をノルムだけで出す

`valuation_of_fixed`（§4）—— `h z = z` なら `‖z‖ = ‖π‖^{p·m}`。
★本波で新しく見つけた道: `p ∤ m` とし、`N := ∏_{j<p} h^jπ`（`h` 不変、`‖N‖ = ‖π‖^p`）と
Bézout `a·m + b·p = 1` から `w := z^a N^b` を作ると `h w = w` かつ `‖w‖ = ‖π‖`。
`v := w/π` は `‖v‖ = 1` で `‖h v − v‖ = ‖π‖^t` だが、
`RamNormBridge.norm_sub_apply_le_of_norm_le_one` は `≤ ‖π‖^{t+1}` を要求する。矛盾。

## ★代替経路と、それが足りない理由（測定済み、他の波のために残す）

* `WildBreakUpperBound.lean:242-256` の `hbreak`「共役がすべて同じ距離」は
  ★**`q = p^{k+1}` では偽**（`ℚ₃(ζ₈₁)` では `u = (2,8,26)` で 3 通りの距離）。
  だから `K` の上で直接 `norm_natCast_le_pow_of_splits` を使う道は閉じている。
* `h := g^{p^k}` の telescoping だけでは `‖p‖ ≤ ‖π‖^{u k}`、すなわち `u k ≤ p^{k+1}e` しか出ず、
  ★**係数 `p−1` が出ない**（`p ≥ 3` で不十分）。係数 `p−1` は差式からしか出ない。

★数値の検算: `ℚ₃(ζ₈₁)/ℚ₃(ζ₃)` は `p=3, k=2, q=27, e=2, u=(2,8,26)`。
(4) は `2·26 = 52 ≤ 27·2 = 54` ★真。`k=1` の例は `1·8 = 8 ≤ 4·2 = 8` ★等号。
弱い経路が与える `u k ≤ q·e` は `26 ≤ 54`、`8 ≤ 8` ——
★`k=1` の例では偶然一致するが `p=3` の例では弱い。

## 出したもの

| 宣言 | 内容 |
|---|---|
| ★★★`norm_sub_one_pow_pow_le` | ★抽象核（純ノルム）`‖A^{p^r}−1‖ ≤ c ⇒ ‖A−1‖^{p^r} ≤ c` |
| ★★`norm_sub_le_sq_of_totallyRamified_pow` | `‖gπ − π‖ ≤ ‖π‖²`（次数 `p^{k+1}`） |
| ★★★`one_le_jump_zero_pow` | **`harith` (1)** `1 ≤ u 0` |
| `apply_prod_orbit` | `∏_{j<p} h^jπ` は `h` 不変 |
| ★★★★★`valuation_of_fixed` | ★**中間体の値群は `‖π‖^{pℤ}`**（体を建てずに） |
| `finrank_fixedSubfield_zpowers` | `[M : M^{⟨h⟩}] = orderOf h` |
| ★★★`sub_one_mul_le_of_totallyRamified_pow` | `(p−1)·t ≤ p^{k+1}·e`（ℕ 版） |
| ★★★★★★★`hupper_of_totallyRamified_pow` | **`harith` (4)**（ℤ 版、字面そのもの） |

## 逸脱の記録

1. §2 は `WildBreakLowerBound.lean:367` の段取りを `p → q` で写している。
   ★既存ファイルは読み取り専用なので一般化版を新しく書いた（元は消していない）。
2. §5 は `hiso` を「すべての `σ`」について受ける（`g^{p^k}` に使うため）。
   木の `PureStepSetup.norm_algEquiv_eq`（スペクトルノルム）がこの形で供給する。
3. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
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

/-! ## §4 ★★★(4) の残り 1 ノードのうち、**中間体を建てずに**取れる部分 -/

section FixedValue

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

omit [IsUltrametricDist M] in
/-- `N(π) := ∏_{j<p} h^j π` は `h` で不変（`h^p = 1`）。

★`WildBreak.prod_telescope` をそのまま使う（帰納法を書き直さない）。 -/
theorem apply_prod_orbit {p : ℕ} {π : M} {h : M ≃ₐ[K] M} (hhp : h ^ p = 1)
    (hiso : ∀ z : M, ‖h z‖ = ‖z‖) (hπ0 : 0 < ‖π‖) :
    h (∏ j ∈ Finset.range p, (h ^ j) π) = ∏ j ∈ Finset.range p, (h ^ j) π := by
  classical
  have hπne : π ≠ 0 := norm_pos_iff.mp hπ0
  have hgj : ∀ (j : ℕ) (z : M), ‖(h ^ j) z‖ = ‖z‖ := norm_pow_apply h hiso
  have hfne : ∀ j : ℕ, (h ^ j) π ≠ 0 := fun j => norm_pos_iff.mp (by rw [hgj]; exact hπ0)
  have hP0 : (∏ j ∈ Finset.range p, (h ^ j) π) ≠ 0 :=
    Finset.prod_ne_zero_iff.mpr (fun j _ => hfne j)
  have hmap : h (∏ j ∈ Finset.range p, (h ^ j) π)
      = ∏ j ∈ Finset.range p, (h ^ (j + 1)) π := by
    rw [map_prod]
    exact Finset.prod_congr rfl (fun j _ => by rw [← AlgEquiv.mul_apply, ← pow_succ'])
  have hdiv : (∏ j ∈ Finset.range p, (h ^ (j + 1)) π)
      / (∏ j ∈ Finset.range p, (h ^ j) π) = 1 := by
    rw [← Finset.prod_div_distrib, prod_telescope h hπne hfne p, hhp]
    simpa using div_self hπne
  rw [hmap, ← div_eq_one_iff_eq hP0]
  exact hdiv

/-- ★★★★★**中間体の値群は `‖π‖^{pℤ}`** —— `h` で不変な元のノルムは `‖π‖^{p·m}`。

★★**`IntermediateField` を一切建てない**（`lean-idioms.md` #59/#69/#296 の危険区間に入らない）。

証明（★本波で新しく見つけた道）: `‖z‖ = ‖π‖^m` で `p ∤ m` とする。
`N := ∏_{j<p} h^jπ` は `h` 不変で `‖N‖ = ‖π‖^p`。Bézout で `a·m + b·p = 1` を取り
`w := z^a N^b` とおくと `h w = w` かつ `‖w‖ = ‖π‖`。`v := w/π` は `‖v‖ = 1` で

    h v − v = w(π − hπ)/((hπ)π),　‖h v − v‖ = ‖π‖·‖π‖^{t+1}/‖π‖² = ‖π‖^t.

一方 `RamNormBridge.norm_sub_apply_le_of_norm_le_one` は `‖h v − v‖ ≤ ‖π‖^{t+1}` を与える。
`‖π‖ < 1` なので矛盾。★すなわち `p ∣ m`。 -/
theorem valuation_of_fixed [FiniteDimensional K M] {p n t : ℕ} (hp : p.Prime)
    {π : M} {h : M ≃ₐ[K] M} (hhp : h ^ p = 1) (hiso : ∀ z : M, ‖h z‖ = ‖z‖)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (hval : ∀ y : M, y ≠ 0 → ∃ m : ℤ, ‖y‖ = ‖π‖ ^ m)
    (hbr : ‖h π - π‖ = ‖π‖ ^ (t + 1))
    {z : M} (hz : h z = z) (hz0 : z ≠ 0) :
    ∃ m : ℤ, ‖z‖ = ‖π‖ ^ ((p : ℤ) * m) := by
  classical
  obtain ⟨m, hm⟩ := hval z hz0
  by_cases hdvd : (p : ℤ) ∣ m
  · obtain ⟨c, rfl⟩ := hdvd
    exact ⟨c, hm⟩
  exfalso
  have hπne : π ≠ 0 := norm_pos_iff.mp hπ0
  have hgj : ∀ (j : ℕ) (y : M), ‖(h ^ j) y‖ = ‖y‖ := norm_pow_apply h hiso
  have hfne : ∀ j : ℕ, (h ^ j) π ≠ 0 := fun j => norm_pos_iff.mp (by rw [hgj]; exact hπ0)
  have hpprime : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  obtain ⟨b, a, hab⟩ := (hpprime.coprime_iff_not_dvd).mpr hdvd
  set N : M := ∏ j ∈ Finset.range p, (h ^ j) π with hNdef
  have hN0 : N ≠ 0 := Finset.prod_ne_zero_iff.mpr (fun j _ => hfne j)
  have hNfix : h N = N := apply_prod_orbit hhp hiso hπ0
  have hNnorm : ‖N‖ = ‖π‖ ^ ((p : ℕ) : ℤ) := by
    rw [hNdef, norm_prod]
    rw [Finset.prod_congr rfl (fun j _ => hgj j π), Finset.prod_const, Finset.card_range,
      zpow_natCast]
  set w : M := z ^ a * N ^ b with hwdef
  have hw0 : w ≠ 0 := by
    rw [hwdef]
    exact mul_ne_zero (zpow_ne_zero _ hz0) (zpow_ne_zero _ hN0)
  have hwfix : h w = w := by
    rw [hwdef, map_mul, map_zpow₀, map_zpow₀, hz, hNfix]
  have hwnorm : ‖w‖ = ‖π‖ := by
    rw [hwdef, norm_mul, norm_zpow, norm_zpow, hm, hNnorm, ← zpow_mul, ← zpow_mul,
      ← zpow_add₀ (ne_of_gt hπ0)]
    have hexp : m * a + ((p : ℕ) : ℤ) * b = 1 := by linarith [hab]
    rw [hexp, zpow_one]
  set v : M := w / π with hvdef
  have hv1 : ‖v‖ = 1 := by rw [hvdef, norm_div, hwnorm, div_self (ne_of_gt hπ0)]
  have hhπ0 : h π ≠ 0 := by
    refine norm_pos_iff.mp ?_
    rw [hiso]; exact hπ0
  have hvsub : h v - v = w * (π - h π) / ((h π) * π) := by
    rw [hvdef, map_div₀, hwfix]
    field_simp
  have hnv : ‖h v - v‖ = ‖π‖ ^ t := by
    rw [hvsub, norm_div, norm_mul, hwnorm, norm_mul, hiso, norm_sub_rev, hbr]
    rw [pow_succ]
    field_simp
  have hle := RamNormBridge.norm_sub_apply_le_of_norm_le_one (n := n) (t := t) h hiso hπ0 hπ1
    hn hvalK hbr (z := v) (le_of_eq hv1)
  rw [hnv] at hle
  have hlt : ‖π‖ ^ (t + 1) < ‖π‖ ^ t := pow_lt_pow_right_of_lt_one₀ hπ0 hπ1 (by omega)
  linarith

end FixedValue

/-! ## §5 ★★★★★★★`harith` (4) を `p^{k+1}` 次に上げる

★中間体を `FixedPoints.subfield ↥(Subgroup.zpowers h) M` として建てる。
★★`IntermediateField` ではなく `Subfield` なので `lean-idioms.md` #59/#69 の層に入らない。 -/

section UpperPow

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]
  [FiniteDimensional K M]

omit [IsUltrametricDist M] in
/-- ★中間体の次数は `orderOf h`。 -/
theorem finrank_fixedSubfield_zpowers (h : M ≃ₐ[K] M) {p : ℕ} (hg : orderOf h = p) :
    Module.finrank ↥(FixedPoints.subfield ↥(Subgroup.zpowers h) M) M = p := by
  classical
  haveI : Fintype ↥(Subgroup.zpowers h) := Fintype.ofFinite _
  rw [FixedPoints.finrank_eq_card, ← Nat.card_eq_fintype_card, Nat.card_zpowers, hg]

/-- ★★★★★★★**`harith` (4)** —— `(p−1)·t ≤ p^{k+1}·e`、`t` は `g^{p^k}` の跡び。

★★**段取り**: `h := g^{p^k}` は位数 `p`。`L := M^{⟨h⟩}` に対し
`WildBreakUpper.sub_one_mul_le_of_totallyRamified`（`k = 0` 版）を当てる。
`[M:L] = p`、`v_L(p) = p^k·e` なので結論は `(p−1)t ≤ p·(p^k e) = p^{k+1} e`。

★鍵は `L` の値群が `‖π‖^{pℤ}` であること（§4 `valuation_of_fixed`）。 -/
theorem sub_one_mul_le_of_totallyRamified_pow {p k e t : ℕ} (hp : p.Prime) {π : M}
    (g : M ≃ₐ[K] M) (hg : orderOf g = p ^ (k + 1))
    (hiso : ∀ (σ : M ≃ₐ[K] M) (z : M), ‖σ z‖ = ‖z‖)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hnK : Module.finrank K M = p ^ (k + 1))
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ,
      ‖algebraMap K M a‖ = ‖π‖ ^ (((p ^ (k + 1) : ℕ) : ℤ) * m))
    (hval : ∀ y : M, y ≠ 0 → ∃ m : ℤ, ‖y‖ = ‖π‖ ^ m)
    (heM : ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e))
    (hbr : ‖(g ^ p ^ k) π - π‖ = ‖π‖ ^ (t + 1)) :
    (p - 1) * t ≤ p ^ (k + 1) * e := by
  classical
  set h : M ≃ₐ[K] M := g ^ p ^ k with hhdef
  have hordh : orderOf h = p := by
    have hkk : k + 1 - k = 1 := by omega
    rw [hhdef, RamCard.orderOf_pow_prime_pow hp hg (by omega), hkk, pow_one]
  have hhp : h ^ p = 1 := by rw [← hordh]; exact pow_orderOf_eq_one h
  haveI : Fintype ↥(Subgroup.zpowers h) := Fintype.ofFinite _
  set L := FixedPoints.subfield ↥(Subgroup.zpowers h) M with hLdef
  have hnL : Module.finrank ↥L M = p := finrank_fixedSubfield_zpowers h hordh
  have hvalL : ∀ a : ↥L, a ≠ 0 → ∃ m : ℤ,
      ‖algebraMap ↥L M a‖ = ‖π‖ ^ ((p : ℤ) * m) := by
    intro a ha
    have hfix : h (a : M) = (a : M) := a.2 ⟨h, Subgroup.mem_zpowers h⟩
    have ha0 : (a : M) ≠ 0 := fun hc => ha (Subtype.ext hc)
    exact valuation_of_fixed hp hhp (fun z => hiso h z) hπ0 hπ1 hnK hvalK hval hbr hfix ha0
  have heM' : ‖(p : M)‖ = ‖π‖ ^ (p * (p ^ k * e)) := by
    rw [heM]
    congr 1
    rw [← mul_assoc, ← pow_succ']
  set h' : M ≃ₐ[↥L] M :=
    FixedPoints.toAlgAutMulEquiv ↥(Subgroup.zpowers h) M ⟨h, Subgroup.mem_zpowers h⟩ with hh'def
  have hordh' : orderOf h' = p := by
    have h1 : orderOf h' = orderOf (⟨h, Subgroup.mem_zpowers h⟩ : ↥(Subgroup.zpowers h)) :=
      orderOf_injective (FixedPoints.toAlgAutMulEquiv ↥(Subgroup.zpowers h) M).toMonoidHom
        (MulEquiv.injective _) _
    rw [h1, Subgroup.orderOf_mk, hordh]
  have happ : ∀ z : M, h' z = h z := fun z => rfl
  have hiso' : ∀ z : M, ‖h' z‖ = ‖z‖ := fun z => by rw [happ]; exact hiso h z
  have hbr' : ‖h' π - π‖ = ‖π‖ ^ (t + 1) := by rw [happ]; exact hbr
  have hmain := WildBreakUpper.sub_one_mul_le_of_totallyRamified (K := ↥L) (e := p ^ k * e)
    hp h' hordh' hiso' hπ0 hπ1 hnL hvalL heM' hbr'
  calc (p - 1) * t ≤ p * (p ^ k * e) := hmain
    _ = p ^ (k + 1) * e := by rw [← mul_assoc, ← pow_succ']

end UpperPow

/-! ## §6 `harith` (4) の字面そのもの（`ℤ` 版） -/

section HarithFour

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]
  [FiniteDimensional K M]

/-- ★★★★★★★**`harith` (4)** —— `((p:ℤ)−1)·u k ≤ (p:ℤ)^{k+1}·e`。

`JumpFromValueGroup.lean:276` の第 4 項そのものの形。
`t ≤ 0` なら左辺 ≤ 0 ≤ 右辺で自明、`0 < t` なら §5 を使う。 -/
theorem hupper_of_totallyRamified_pow {p k e : ℕ} (hp : p.Prime) {π : M} {t : ℤ}
    (g : M ≃ₐ[K] M) (hg : orderOf g = p ^ (k + 1))
    (hiso : ∀ (σ : M ≃ₐ[K] M) (z : M), ‖σ z‖ = ‖z‖)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hnK : Module.finrank K M = p ^ (k + 1))
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ,
      ‖algebraMap K M a‖ = ‖π‖ ^ (((p ^ (k + 1) : ℕ) : ℤ) * m))
    (hval : ∀ y : M, y ≠ 0 → ∃ m : ℤ, ‖y‖ = ‖π‖ ^ m)
    (heM : ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e))
    (hbr : ‖(g ^ p ^ k) π - π‖ = ‖π‖ ^ (t + 1)) :
    ((p : ℤ) - 1) * t ≤ (p : ℤ) ^ (k + 1) * (e : ℤ) := by
  have hp1 : (1 : ℤ) ≤ (p : ℤ) := by exact_mod_cast hp.one_le
  have hpe : (0 : ℤ) ≤ (p : ℤ) ^ (k + 1) * (e : ℤ) := by positivity
  by_cases ht : 0 < t
  · have hi : ((t.toNat : ℕ) : ℤ) = t := Int.toNat_of_nonneg (le_of_lt ht)
    have hbr' : ‖(g ^ p ^ k) π - π‖ = ‖π‖ ^ (t.toNat + 1) := by
      have hcast1 : ((t.toNat + 1 : ℕ) : ℤ) = t + 1 := by omega
      rw [hbr, ← zpow_natCast ‖π‖ (t.toNat + 1), hcast1]
    have hmain := sub_one_mul_le_of_totallyRamified_pow hp g hg hiso hπ0 hπ1 hnK hvalK hval heM hbr'
    have hcast : ((p - 1 : ℕ) : ℤ) * ((t.toNat : ℕ) : ℤ)
        ≤ ((p ^ (k + 1) : ℕ) : ℤ) * ((e : ℕ) : ℤ) := by exact_mod_cast hmain
    rw [Nat.cast_sub hp.one_le, hi] at hcast
    push_cast at hcast
    linarith
  · rw [not_lt] at ht
    nlinarith

end HarithFour








/-! ## `.src` と 公理 -/

def norm_sub_le_sq_of_totallyRamified_pow.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def one_le_jump_zero_pow.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def sub_one_mul_le_of_totallyRamified_pow.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

#print axioms norm_sub_one_pow_pow_le
#print axioms norm_sub_le_sq_of_totallyRamified_pow
#print axioms one_le_jump_zero_pow
#print axioms apply_prod_orbit
#print axioms valuation_of_fixed
#print axioms finrank_fixedSubfield_zpowers
#print axioms sub_one_mul_le_of_totallyRamified_pow
#print axioms hupper_of_totallyRamified_pow

end WildBreakPow

end ABC3.Found.PGC
