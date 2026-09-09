import ABC3.Found.PGC.GainedBudgetDeficit

/-!
# [pGC] ★★★答: 桁展開は底まで届くが**ノルムの評価が届かない** —— 塞いでいるのは `hchar` 1 本

## 持ち場（前波で「次の 1 点」とした点）

前波の私の言葉（逐語）:

> ★**`Gained*` の中身（`exists_digitSum_of_adjoin_eq_top` による桁展開）が、
> 底の `K` まで一度に展開できるか**を読むこと。★まだ読んでいません。

## ★★★読んだ結果 —— 桁展開は**次数に依らない**（底まで届く）

`MinpolyOrbitSplit.lean:528`（名前空間 `ABC3.Found.PGC`、`grep -n '^namespace'` で確認）:

```lean
theorem exists_digitSum_of_adjoin_eq_top {p : ℕ} {π : L} (hint : IsIntegral K π)
    (hdeg : (minpoly K π).natDegree = p) (htop : Algebra.adjoin K ({π} : Set L) = ⊤)
    (x : L) : ∃ a : ℕ → K, x = digitSum p π a
```

★**`p` は「素数」である必要がない** —— 単に `minpoly K π` の次数である。
⇒ `M = K(π)` で `[M:K] = p^{k+1}` なら、`x = Σ_{j < p^{k+1}} a_j π^j`（★係数は**底の `K`**）。
★★**桁展開そのものは底まで一度に届く。**

## ★★★ではなぜ層ごと（`[M:F] = p`）なのか —— `hchar` ちょうど 1 本

値段を与えるのは `CyclicJumpNorm.lean:404 norm_digitSum_sub_mul_norm_eq`（★**等号**）:

  `‖x^u − x‖ · ‖π‖ = ‖u − π‖ · ‖x − a_0‖`

その仮説に

  `hchar : ∀ j, 0 < j → j < n → ‖(j : L)‖ = 1`

がある。理由は同ファイル `:188` の「★★★**抽象核(本ファイルの心臓)**」
`norm_pow_sub_pow_eq`（`‖(k:L)‖ = 1` ならば `‖u^k − π^k‖ = ‖u−π‖·‖π‖^{k−1}` が**等号**）で、
docstring が自分でこう言っている（逐語）:

> ★`‖(k:L)‖ = 1` は具体層では「`0 < k < p` は `p` 進単数」の言い換えである。

★★**`hchar` は `n > p` で偽である**（`j := p` を取れば `‖(p:L)‖ < 1 ≠ 1`）。§1 で型にした。

⇒ ★★★**答: 桁展開は底まで届くが、ノルムの等式が届かない。
層の次数が `p` でなければならない理由は `hchar` ちょうど 1 本である。**

## ★測定（コマンドと出力）

`grep -rhn "hchar : ∀ j, 0 < j → j < " lean/ABC3/Found/PGC/*.lean | sed 's/.*j < //' | sort | uniq -c | sort -rn`

```
     14 p → ‖(j : M)‖ = 1)
      5 p → ‖(j : M)‖ = 1 :=
      5 n → ‖(j : L)‖ = 1)
      1 p → ‖(j : L)‖ = 1)
```

⇒ ★**具体層 20 箇所すべてが `j < p`**。`j < n` の 5 箇所は `CyclicJumpNorm` の抽象核
（`n` は変数）。★`hchar` を出現ファイル数で数えると **9** ファイル。

★予防 4 本目（前波で決めたもの）を書く前に実行した:
`grep -rn "抽象核" lean/ABC3/Found/PGC/CyclicJumpNorm.lean …` →
`CyclicJumpNorm.lean:179` が「★★★抽象核(本ファイルの心臓)」だと分かり、
★**そこが `hchar` の出どころ**だと特定できた。★重複を作らずに済んだ。

## ★これで「残っている 1 点」が言い換わった

* 桁展開: 底まで届く（`exists_digitSum_of_adjoin_eq_top`、次数に依らない）
* ノルムの等式: `n ≤ p` でしか使えない（`hchar`）
* 過払い: 着地の深さ `d'` に対して因子 `∏_{Icc 1 d'}`（前波 `gained_overpay_factor`、等式）

⇒ ★★**要るのは「`p` 進単数でない桁（`p | j`）を持つ展開でも使えるノルムの下界」**である。
★`p | j` の桁では `‖u^j − π^j‖` が**真に小さくなる**ので、
★★等号は諦めて「**どの桁が主役か**」を選ぶ形（私の第 1120 波台の `j₀` の議論と同じ形）になる。
★これは本ファイルでは**やっていない**。

## 逸脱の記録

- §1・§2 は `NormedField` の `‖(j : L)‖` だけで、★分岐も Galois も桁展開も出ない。
- ★本ファイルは「`hchar` を回避する評価」を**与えていない**。★与えられるとも書いていない。
  ★塞いでいるものを 1 本に**特定した**だけである。
-/

namespace ABC3.Found.PGC

namespace LayerDegreeIsP

/-! ## §1 ★抽象核 —— `hchar` は `n > p` で偽 -/

section Kernel

variable {L : Type*} [NormedField L]

/-- ★★**`hchar` は `n > p` で成り立たない**。

`j := p` を取れば `0 < p < n` かつ `‖(p : L)‖ < 1 ≠ 1`。

⇒ `CyclicJumpNorm.norm_digitSum_sub_mul_norm_eq`（★ノルムの**等号**）は
★**桁数 `n` が `p` を超えると使えない**。 -/
theorem hchar_fails_of_prime_lt {p n : ℕ} (hp : 0 < p) (hpn : p < n)
    (hple : ‖(p : L)‖ < 1) :
    ¬ (∀ j : ℕ, 0 < j → j < n → ‖(j : L)‖ = 1) :=
  fun h => absurd (h p hp hpn) (ne_of_lt hple)

/-- ★★**`hchar` が成り立つ `n` はちょうど `n ≤ p`**。

`hlt`（`0 < j < p` は単数）と `hple`（`‖(p:L)‖ < 1`）は `p` 進局所体の言い換えである。 -/
theorem hchar_iff_le {p n : ℕ} (hp : 0 < p)
    (hlt : ∀ j : ℕ, 0 < j → j < p → ‖(j : L)‖ = 1) (hple : ‖(p : L)‖ < 1) :
    (∀ j : ℕ, 0 < j → j < n → ‖(j : L)‖ = 1) ↔ n ≤ p := by
  constructor
  · intro h
    by_contra hn
    exact hchar_fails_of_prime_lt hp (by omega) hple h
  · intro hn j hj0 hjn
    exact hlt j hj0 (by omega)

/-- ★★★**層の次数の上限はちょうど `p`**。

`{n | hchar が n で成り立つ}` の最大元が `p` である。
⇒ ★`Gained*` / `CyclicJumpNorm` が `[M:F] = p` の 1 層しか扱えない理由は**これ**であり、
★予算（`axDecay` の積）の側の問題ではない。 -/
theorem max_layer_degree {p : ℕ} (hp : 0 < p)
    (hlt : ∀ j : ℕ, 0 < j → j < p → ‖(j : L)‖ = 1) (hple : ‖(p : L)‖ < 1) :
    IsGreatest {n : ℕ | ∀ j : ℕ, 0 < j → j < n → ‖(j : L)‖ = 1} p := by
  constructor
  · exact fun j hj0 hjp => hlt j hj0 hjp
  · intro n hn
    exact (hchar_iff_le hp hlt hple).mp hn

end Kernel

/-! ## §2 ★`p` 進局所体での具体化 -/

section Concrete

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-- `‖(p : K.carrier)‖ < 1`（`PureStepSetup.norm_natCast_p_eq_inv`（`:455`）から）。 -/
theorem norm_natCast_p_lt_one_carrier (K : PAdicLocalField p) :
    ‖((p : ℕ) : K.carrier)‖ < 1 := by
  have h1 : (1 : ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  rw [PureStepSetup.norm_natCast_p_eq_inv K]
  rw [inv_lt_one_iff₀]
  exact Or.inr h1

/-- ★★**`K.carrier` では `hchar` は `n > p` で偽**。

⇒ ★桁展開を `p` 桁より深く取ると `CyclicJumpNorm` の等号が使えない。 -/
theorem hchar_fails_carrier (K : PAdicLocalField p) {n : ℕ} (hpn : p < n) :
    ¬ (∀ j : ℕ, 0 < j → j < n → ‖((j : ℕ) : K.carrier)‖ = 1) :=
  hchar_fails_of_prime_lt (Fact.out : p.Prime).pos hpn (norm_natCast_p_lt_one_carrier K)

end Concrete

/-! ## §3 使っている公理の一覧 -/

#print axioms hchar_fails_of_prime_lt
#print axioms hchar_iff_le
#print axioms max_layer_degree
#print axioms norm_natCast_p_lt_one_carrier
#print axioms hchar_fails_carrier

end LayerDegreeIsP

end ABC3.Found.PGC
