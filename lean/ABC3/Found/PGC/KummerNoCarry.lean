import ABC3.Found.PGC.KummerReduction

/-!
# [pGC] ★★★Kummer は mathlib に**在った** —— `Nat.factorization_choose'` が繰り上がりの数

## 持ち場（前波で「残っているちょうど 1 点」とした点）

前波の私の言葉（逐語）:

> `p ∤ C(p^a·m − 1, p^a − 1)`。★mathlib に Lucas が無いので**自分で書く**か、
> `Nat.factorization_choose` 系で回るかを測る。★**どちらもまだ測っていません**。

## ★★★測った結果 —— **回る。しかも一発で。**

`Nat.factorization_choose'`（`Data/Nat/Choose/Factorization.lean:114`、逐語）:

```lean
theorem factorization_choose' {p n k b : ℕ} (hp : p.Prime) (hnb : log p (n + k) < b) :
    (choose (n + k) k).factorization p = #{i ∈ Ico 1 b | p ^ i ≤ k % p ^ i + n % p ^ i}
```

★★**これは Kummer の定理そのもの**である（右辺が「基数 `p` で `k` と `n` を足すときの
繰り上がりの回数」）。★前波の私は `grep -in "lucas"` しか叩かず `LucasLehmer` しか出なかったので
「無い」と書いた。★★**名前で引いて外した**（`Kummer` でも `Lucas` でも出ない。
出るのは `factorization_choose'` という**結論の形**の名前である）。

★★★**予防 5 本目**: ★**定理の「名前」ではなく「結論の形」で引く**。
本波は `grep -n "Nat.factorization_choose" .cache/mathlib-index.txt` で当てた。
（既存 4 本: 数は `wc -l` ／ grep は実行してから書く ／ 名前を引くときは定義も開く ／
抽象核を書く前に「抽象核」と自称する節を先に見る。）

## ★★本ファイルの結果 —— 前波の目標を**追い越した**

前波は `p ∤ C(p^a·m, p^a)` を `p ∤ m` と `p ∤ C(p^a·m−1, p^a−1)` に**分けた**が、
★`factorization_choose'` は ★**分けずに直接**出す（§2 `not_dvd_choose_prime_pow_mul`）:

  `p` 素数, `0 < m`, `p ∤ m` ⇒ `¬ p ∣ C(p^a·m, p^a)`

⇒ ★前波の `KummerReduction.not_dvd_choose_prime_pow_of_not_dvd`（分ける形）は
★**使う必要が無くなった**（誤りではない）。

### 繰り上がりが 0 である理由（§1 が型にしたもの）

`k = p^a`, `n = p^a(m−1)` として、各 `i ≥ 1` で `k % p^i + n % p^i < p^i`:

* `i ≤ a`: `p^i ∣ p^a` なので両方 `0`。和は `0 < p^i`。
* `i = a + c`（`c ≥ 1`）: `p^a % p^{a+c} = p^a`、
  `p^a(m−1) % p^{a+c} = p^a·((m−1) % p^c)`（`Nat.mul_mod_mul_left`）。
  和は `p^a·(1 + (m−1) % p^c)` で、★`p ∤ m` から `(m−1) % p^c ≠ p^c − 1`
  （さもなくば `p^c ∣ m`）なので `1 + (m−1) % p^c < p^c`。よって和 `< p^{a+c}`。

★★**`p ∤ m` が効くのは 2 番目の場合だけ**である。

## ★これで `p ∣ j₀` の桁の「主役」の候補が使えるようになった

第 1137 の `DominantTermLowerBound.norm_pow_sub_pow_eq_of_dominant` に
`k₀ = j − p^a` を入れるとき、係数 `C(j, p^a)` が単数であることが要る。
★★**それが本ファイルで埋まった**（`j = p^a·m`, `p ∤ m`）。

★ただし「その項が**主役である**」（他をすべて上回る）ことは**別の条件**であり、
第 1138 の `term_lt_of_index_lt`（`‖δ‖ < ‖π‖` なら大きい `k` ほど大きい）から
`k = j − p^a < j − 1`（`a ≥ 1`）なので ★**単調性だけでは主役にならない**。
★形式微分の項（`k = j−1`、係数 `‖(j:L)‖ = p^{−a}`）との比較が要る。
★★**本ファイルはその比較を与えていない。** ★与えられるとも書いていない。

## 逸脱の記録

- §1・§2 は `ℕ` だけ。★ノルムも体も分岐も出ない。
- ★前波の `KummerReduction` は**訂正しない**（誤りではない）。★本ファイルで
  「分けなくても直接出る」と記録した。
-/

namespace ABC3.Found.PGC

namespace KummerNoCarry

open Finset

/-! ## §1 ★繰り上がりが起きないこと（`ℕ` の剰余だけ） -/

section NoCarry

/-- ★★`p ∤ m`（`c ≥ 1`）なら `(m − 1) % p^c ≠ p^c − 1`。

★さもなくば `p^c ∣ m`、したがって `p ∣ m`。 -/
theorem mod_ne_of_not_dvd {p m c : ℕ} (hp : 1 < p) (hc : 0 < c) (hm : 0 < m)
    (hpm : ¬ p ∣ m) : (m - 1) % p ^ c ≠ p ^ c - 1 := by
  intro hcon
  have hpc : 0 < p ^ c := pow_pos (by omega) c
  have h := Nat.div_add_mod (m - 1) (p ^ c)
  rw [hcon] at h
  have hmul : p ^ c * ((m - 1) / p ^ c + 1) = p ^ c * ((m - 1) / p ^ c) + p ^ c := by ring
  have hdvd : p ^ c ∣ m := ⟨(m - 1) / p ^ c + 1, by omega⟩
  exact hpm (dvd_trans (dvd_pow_self p (by omega)) hdvd)

/-- ★★**繰り上がりは起きない** —— `k = p^a`, `n = p^a(m−1)` について
`∀ i ≥ 1, k % p^i + n % p^i < p^i`。

★これが Kummer の判定条件（`p ^ i ≤ k % p ^ i + n % p ^ i`）の否定である。 -/
theorem no_carry {p a m : ℕ} (hp : 1 < p) (hm : 0 < m) (hpm : ¬ p ∣ m) {i : ℕ} (hi : 0 < i) :
    p ^ a % p ^ i + (p ^ a * (m - 1)) % p ^ i < p ^ i := by
  have hp0 : 0 < p := by omega
  have hpi : 0 < p ^ i := pow_pos hp0 i
  -- ★`le_or_lt` は無い（#302）ので `Nat.lt_or_ge` を使う
  rcases Nat.lt_or_ge a i with hai | hia
  · -- `a < i`: `i = a + c`, `c ≥ 1`
    obtain ⟨c, hc, rfl⟩ : ∃ c, 0 < c ∧ i = a + c := ⟨i - a, by omega, by omega⟩
    have hpc : 0 < p ^ c := pow_pos hp0 c
    have hpa : 0 < p ^ a := pow_pos hp0 a
    have hsplit : p ^ (a + c) = p ^ a * p ^ c := pow_add p a c
    have h1 : p ^ a % p ^ (a + c) = p ^ a := by
      refine Nat.mod_eq_of_lt ?_
      rw [hsplit]
      have hc2 : 1 < p ^ c := Nat.one_lt_pow (by omega) hp
      nlinarith
    have h2 : (p ^ a * (m - 1)) % p ^ (a + c) = p ^ a * ((m - 1) % p ^ c) := by
      rw [hsplit]
      exact Nat.mul_mod_mul_left _ _ _
    have hne : (m - 1) % p ^ c ≠ p ^ c - 1 := mod_ne_of_not_dvd hp hc hm hpm
    have hlt : (m - 1) % p ^ c < p ^ c := Nat.mod_lt _ hpc
    have hkey : 1 + (m - 1) % p ^ c < p ^ c := by omega
    rw [h1, h2, hsplit]
    calc p ^ a + p ^ a * ((m - 1) % p ^ c) = p ^ a * (1 + (m - 1) % p ^ c) := by ring
      _ < p ^ a * p ^ c := mul_lt_mul_of_pos_left hkey hpa
  · -- `i ≤ a`: どちらも `0`
    obtain ⟨d, hd⟩ : p ^ i ∣ p ^ a := pow_dvd_pow p hia
    have h1 : p ^ a % p ^ i = 0 := by rw [hd]; exact Nat.mul_mod_right _ _
    have h2 : (p ^ a * (m - 1)) % p ^ i = 0 := by
      rw [hd, mul_assoc]
      exact Nat.mul_mod_right _ _
    omega

end NoCarry

/-! ## §2 ★★★Kummer を当てる -/

section Kummer

/-- ★★★**`p ∤ m` なら `p ∤ C(p^a·m, p^a)`**（Kummer、繰り上がり 0）。

★前波は `p ∤ m` と `p ∤ C(p^a·m−1, p^a−1)` に**分けた**が、
★`Nat.factorization_choose'` は**分けずに直接**出す。 -/
theorem not_dvd_choose_prime_pow_mul {p a m : ℕ} (hp : p.Prime) (hm : 0 < m)
    (hpm : ¬ p ∣ m) : ¬ p ∣ (p ^ a * m).choose (p ^ a) := by
  classical
  have hp1 : 1 < p := hp.one_lt
  have hpa : 0 < p ^ a := pow_pos hp.pos a
  have hsum : p ^ a * (m - 1) + p ^ a = p ^ a * m := by
    have : p ^ a * (m - 1) + p ^ a * 1 = p ^ a * (m - 1 + 1) := by ring
    rw [mul_one] at this
    rw [this]
    congr 1
    omega
  set b := Nat.log p (p ^ a * m) + 1 with hb
  have hnb : Nat.log p (p ^ a * (m - 1) + p ^ a) < b := by
    rw [hsum, hb]
    omega
  have hfac := Nat.factorization_choose' (p := p) (n := p ^ a * (m - 1)) (k := p ^ a) hp hnb
  rw [hsum] at hfac
  have hempty : ((Finset.Ico 1 b).filter
      (fun i => p ^ i ≤ p ^ a % p ^ i + (p ^ a * (m - 1)) % p ^ i)) = ∅ := by
    refine Finset.filter_eq_empty_iff.mpr ?_
    intro i hi
    have hi1 : 0 < i := (Finset.mem_Ico.mp hi).1
    exact not_le_of_gt (no_carry hp1 hm hpm hi1)
  rw [hempty, Finset.card_empty] at hfac
  have hne : (p ^ a * m).choose (p ^ a) ≠ 0 := by
    refine Nat.ne_of_gt (Nat.choose_pos ?_)
    exact Nat.le_mul_of_pos_right _ hm
  intro hcon
  have := (hp.dvd_iff_one_le_factorization hne).mp hcon
  omega

end Kummer

/-! ## §3 使っている公理の一覧 -/

#print axioms mod_ne_of_not_dvd
#print axioms no_carry
#print axioms not_dvd_choose_prime_pow_mul

end KummerNoCarry

end ABC3.Found.PGC
