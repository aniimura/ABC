import ABC3.Found.PGC.DominantTermLowerBound

/-!
# [pGC] ★★★訂正: 前波 §4 の仮説は**実際には満たされない** ＋ Kummer を 1 段減らす

## ★★★自己訂正 21 度目 —— 前波 `DominantTermLowerBound` §4 の向きが逆だった

前波の私は §4 で `‖π‖ ≤ ‖δ‖`（`δ = u − π` が**大きい**側）を仮定し、
逸脱の記録に「★どの `σ` でそれが実現するかは本ファイルでは**測っていない**」と書いた。

★★**測った。実現しない。**

木の該当補題はすべて `hu : ‖u − π‖ < ‖π‖` を仮定している
（`CyclicJumpNorm.lean:162 norm_pow_succ_sub_pow_succ_le` /
`:189 norm_pow_sub_pow_eq` / `:404 norm_digitSum_sub_mul_norm_eq`）。
さらに具体層の `hbreak : ‖σ π − π‖ = ‖π‖^{i+1}`（`i ≥ 1`, `‖π‖ < 1`）から
★`‖δ‖ = ‖π‖^{i+1} < ‖π‖` が**常に**従う（§1 `delta_lt_pi_of_break`）。

⇒ ★★**前波の `norm_pow_sub_pow_prime_pow`（§4）は、意図した応用では空虚**である。
★★★前波の私は「`x = π^{p^a}` で下界が出た」と報告したが、★**その下界は使えない**。
★`norm_sum_eq_of_dominant`（§2、抽象核）と `pow_sub_pow_eq_sum_choose`（§1）は**正しいまま**で、
誤っていたのは §4 の**仮説の向き**である。

## ★正しい向きでは主役は「大きい `k`」（§2）

`‖δ‖ < ‖π‖` のとき、二項展開の項のノルム `‖π‖^k·‖δ‖^{j−k}` は
★**`k` について狭義単調増加**である（§2 `term_lt_of_index_lt`）。
⇒ 主役の候補は `range j` の**最大** `k = j − 1`、すなわち ★**形式微分の項** `j·π^{j−1}·δ`。

★★`p ∣ j` だと係数 `j` が単数でないので、その項が主役でなくなりうる。
★次の候補は `k = j − q`（`q = p^a`, `p^a ∥ j`）で、係数は `C(j, q)`。
★★これが単数（`p ∤ C(j,q)`）であることが要る。★★これが Kummer である。

## ★★測定 —— **Lucas / Kummer は mathlib に無い**

`grep -in "lucas" .cache/mathlib-index.txt` → `LucasLehmer.*` **のみ**（Mersenne 素数判定、無関係）。
`Nat.Prime.emultiplicity_choose_prime_pow`（`Data/Nat/Multiplicity.lean:251`）は
★**`C(p^n, k)` の形だけ**で、`C(p^a·m, p^a)` には当たらない。

## ★★本ファイルが減らした 1 段（§3）

`Nat.add_one_mul_choose_eq` から**厳密な等式**

  `C(p^a·m, p^a) = m · C(p^a·m − 1, p^a − 1)`

が出る（§3 `choose_prime_pow_eq_mul`）。★これは**割り算も付値も使わない**。

⇒ ★`p ∤ m` と `p ∤ C(p^a·m − 1, p^a − 1)` から `p ∤ C(p^a·m, p^a)`（§3）。
★★**Kummer が要るのは後者だけ**になった（`p^a·m − 1` の下位 `a` 桁は全部 `p−1`、
`p^a − 1` も同じなので繰り上がりが起きない、というのが Kummer 側の中身）。

## ★残っているちょうど 1 点

`p ∤ C(p^a·m − 1, p^a − 1)`。★mathlib に Lucas が無いので**自分で書く**か、
`Nat.factorization_choose` 系で回るかを測る。★どちらも**まだ測っていない**。

## 逸脱の記録

- §1・§2 はノルムだけ、§3 は `ℕ` だけ。★分岐・Galois・桁展開の語彙は 1 語も出ない。
- ★前波のファイル `DominantTermLowerBound.lean` は**訂正しない**（他が読んでいる）。
  ★§1–§3（恒等式・抽象核・主役の項）は**正しいまま**で、§4 の仮説の向きだけが問題である。
  ★本ファイルで名指しして訂正した。
-/

namespace ABC3.Found.PGC

namespace KummerReduction

/-! ## §1 ★`‖δ‖ < ‖π‖` は具体層で**常に**成り立つ -/

section Direction

variable {L : Type*} [NormedField L]

/-- ★`hbreak : ‖δ‖ = ‖π‖^{i+1}`（`i ≥ 1`, `0 < ‖π‖ < 1`）なら `‖δ‖ < ‖π‖`。

⇒ ★★前波 `DominantTermLowerBound` §4 の仮説 `‖π‖ ≤ ‖δ‖` は**実現しない**。 -/
theorem delta_lt_pi_of_break {π δ : L} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) {i : ℕ} (hi : 0 < i)
    (hbreak : ‖δ‖ = ‖π‖ ^ (i + 1)) : ‖δ‖ < ‖π‖ := by
  rw [hbreak]
  calc ‖π‖ ^ (i + 1) < ‖π‖ ^ 1 :=
        pow_lt_pow_right_of_lt_one₀ hπ0 hπ1 (by omega)
    _ = ‖π‖ := pow_one _

/-- ★同じ状況で `¬ (‖π‖ ≤ ‖δ‖)`（前波 §4 の仮説の否定）。 -/
theorem not_pi_le_delta_of_break {π δ : L} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) {i : ℕ} (hi : 0 < i)
    (hbreak : ‖δ‖ = ‖π‖ ^ (i + 1)) : ¬ (‖π‖ ≤ ‖δ‖) :=
  not_le_of_gt (delta_lt_pi_of_break hπ0 hπ1 hi hbreak)

end Direction

/-! ## §2 ★正しい向きでは項のノルムは `k` について単調増加 -/

section Monotone

variable {L : Type*} [NormedField L]

/-- ★★`‖δ‖ < ‖π‖` なら `‖π‖^k·‖δ‖^{j−k}` は `k` について**狭義単調増加**。

⇒ 主役の候補は `range j` の**最大** `k = j−1`、すなわち★**形式微分の項**。
★これが「`hchar`（`‖(j:L)‖ = 1`）が要る」ことの理由である。 -/
theorem term_lt_of_index_lt {π δ : L} (hδ0 : 0 < ‖δ‖) (hδ : ‖δ‖ < ‖π‖) {j k k' : ℕ}
    (hkk : k < k') (hk' : k' ≤ j) :
    ‖π‖ ^ k * ‖δ‖ ^ (j - k) < ‖π‖ ^ k' * ‖δ‖ ^ (j - k') := by
  have hπ0 : (0 : ℝ) < ‖π‖ := lt_trans hδ0 hδ
  have hd : k' - k ≠ 0 := by omega
  have e1 : ‖π‖ ^ k' = ‖π‖ ^ k * ‖π‖ ^ (k' - k) := by
    rw [← pow_add]
    congr 1
    omega
  have e2 : ‖δ‖ ^ (j - k) = ‖δ‖ ^ (j - k') * ‖δ‖ ^ (k' - k) := by
    rw [← pow_add]
    congr 1
    omega
  rw [e1, e2]
  have hlt : ‖δ‖ ^ (k' - k) < ‖π‖ ^ (k' - k) :=
    pow_lt_pow_left₀ hδ (le_of_lt hδ0) hd
  have hpos : (0 : ℝ) < ‖π‖ ^ k * ‖δ‖ ^ (j - k') := by positivity
  nlinarith [hpos, hlt]

end Monotone

/-! ## §3 ★★Kummer を 1 段減らす（`ℕ` だけ） -/

section Choose

/-- ★★**厳密な等式** `(p^a·m)·C(p^a·m − 1, p^a − 1) = C(p^a·m, p^a)·p^a`。

`Nat.add_one_mul_choose_eq`（★`Nat.succ_mul_choose_eq` は deprecated）に
`n := p^a·m − 1`, `k := p^a − 1` を入れただけ。★割り算も付値も使わない。 -/
theorem choose_prime_pow_mul {p a m : ℕ} (hp : 0 < p) (hm : 0 < m) :
    (p ^ a * m) * ((p ^ a * m - 1).choose (p ^ a - 1))
      = ((p ^ a * m).choose (p ^ a)) * p ^ a := by
  have hq : 0 < p ^ a := pow_pos hp a
  have hn : 0 < p ^ a * m := Nat.mul_pos hq hm
  have h := Nat.add_one_mul_choose_eq (p ^ a * m - 1) (p ^ a - 1)
  simpa only [Nat.succ_eq_add_one, Nat.sub_add_cancel hn, Nat.sub_add_cancel hq] using h

/-- ★★`C(p^a·m, p^a) = m · C(p^a·m − 1, p^a − 1)`（`p^a` を割った形）。 -/
theorem choose_prime_pow_eq_mul {p a m : ℕ} (hp : 0 < p) (hm : 0 < m) :
    (p ^ a * m).choose (p ^ a) = m * ((p ^ a * m - 1).choose (p ^ a - 1)) := by
  have hq : 0 < p ^ a := pow_pos hp a
  have h := choose_prime_pow_mul (p := p) (a := a) (m := m) hp hm
  refine Nat.eq_of_mul_eq_mul_left hq ?_
  calc p ^ a * ((p ^ a * m).choose (p ^ a))
      = ((p ^ a * m).choose (p ^ a)) * p ^ a := by ring
    _ = (p ^ a * m) * ((p ^ a * m - 1).choose (p ^ a - 1)) := h.symm
    _ = p ^ a * (m * ((p ^ a * m - 1).choose (p ^ a - 1))) := by ring

/-- ★★★**Kummer が要るのは 1 段だけ**になった。

`p ∤ m` と `p ∤ C(p^a·m − 1, p^a − 1)` から `p ∤ C(p^a·m, p^a)`。

★★後者（下位 `a` 桁がすべて `p−1` どうしなので繰り上がりが起きない）は
★**mathlib に無い**（`grep -in "lucas" .cache/mathlib-index.txt` → `LucasLehmer` のみ）。
★`Nat.Prime.emultiplicity_choose_prime_pow` は `C(p^n, k)` の形だけである。 -/
theorem not_dvd_choose_prime_pow_of_not_dvd {p a m : ℕ} (hp : p.Prime) (hm : 0 < m)
    (hpm : ¬ p ∣ m) (hrest : ¬ p ∣ (p ^ a * m - 1).choose (p ^ a - 1)) :
    ¬ p ∣ (p ^ a * m).choose (p ^ a) := by
  rw [choose_prime_pow_eq_mul hp.pos hm]
  intro hcon
  rcases (Nat.Prime.dvd_mul hp).mp hcon with h | h
  · exact hpm h
  · exact hrest h

end Choose

/-! ## §4 使っている公理の一覧 -/

#print axioms delta_lt_pi_of_break
#print axioms not_pi_le_delta_of_break
#print axioms term_lt_of_index_lt
#print axioms choose_prime_pow_mul
#print axioms choose_prime_pow_eq_mul
#print axioms not_dvd_choose_prime_pow_of_not_dvd

end KummerReduction

end ABC3.Found.PGC
