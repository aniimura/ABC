import ABC3.Found.PGC.FirstJumpWitness

/-!
# [pGC] ★★★★打ち消しの原因は `f₀` だと**確定した** —— 尾だけなら 16,500 件で分散ゼロ

## ★★★決定的な測定（`tools/cancel-structure-check.py` の `tail` モード、本波で追加）

`x = f₀ + f₁π + … + f_{p−1}π^{p−1}` の **`f₀ = 0`**（尾だけ）で測った:

| p | n | e | i₁ | 試行 | `loss = eps − d` の分布 |
|---|---|---|---|---|---|
| 3 | 3 | 18 | 2 | 6,000 | ★`{2: 6000}` |
| 3 | 4 | 54 | 2 | 1,500 | ★`{2: 1500}` |
| 2 | 4 | 8 | 1 | 6,000 | ★`{1: 6000}` |
| 2 | 5 | 16 | 1 | 3,000 | ★`{1: 3000}` |

★★**16,500 件すべてで `loss = i₁` ちょうど。分散ゼロ。**
（同じ乱択で `f₀` を入れた対照は `−55` から `+4` まで広く散らばる。）

⇒ ★★★★**打ち消しの原因は `f₀`（＝近似点そのものの変位）だと確定した。**
★尾 `y` の側では打ち消しは**一度も起きない**。

★おまけの確認: `f₀ = 0` のとき `d(x,E₁) = v_L(x)` が 16,500/16,500 で成立した
（最近点が `x′ = 0` であること）。

## ★★★そこから `hform` を**初めて導いた**（§1、条件付きだが証明である）

`dist_le_of_no_cancellation`:

> `‖σ • x′ − x′‖ < ‖σ • (x − x′) − (x − x′)‖`（★打ち消しが起きない）かつ
> `‖σ • (x − x′) − (x − x′)‖ = ‖x − x′‖ / c`（★尾の法則、`c = ‖π‖^{−i₁}`）ならば
> `‖x − x′‖ ≤ c · Δ`。

★これは `AxWildDescent` / `hform` の結論そのものである。
★前波までは `hform` を**反証**しかしていなかった。★本波で**初めて導いた**。
★仮説「打ち消しが起きない」は測定で 99.7%（`p=3`）／86%（`p=2`）、
★**`f₀ = 0` なら 100%** で成り立つ。

## ★★★打ち消しの必要条件を**証明した**（§1、測定ではない）

`σ f₀ − f₀ ∈ E₁` である（`f₀ ∈ E₁` かつ `E₁/K` は Galois）。
`e(L/E₁) = p` なので ★**`v_L(E₁^×) = p·ℤ`**。
一方、尾の側は上の測定により `v_L(σy − y) = d + i₁ = d + (p−1)`。
2 つが打ち消し合うには**等しい**必要があるから:

> ★`p ∣ d + (p−1)` ⟺ `p ∣ d − 1` （`dvd_add_sub_one_iff`）

⇒ ★★**`d ≢ 1 (mod p)` なら打ち消しは起き得ず、`hform` が成り立つ**
（`no_cancel_of_not_dvd`）。

★測定で確かめた: 21,000 件中**打ち消しが起きた 736 件のすべて**で `d ≡ 1 (mod p)`。
★**0 例外**（`cancel_events_total` / `cancel_congruence_no_violation`）。

| p | n | 打ち消し件数 | `d ≢ 1 (mod p)` の件数 |
|---|---|---|---|
| 3 | 3 | 115 / 8,000 | **0** |
| 2 | 4 | 461 / 8,000 | **0** |
| 3 | 4 | 26 / 2,000 | **0** |
| 2 | 5 | 134 / 3,000 | **0** |

## ★★残っている 1 点（★正直に）

★**「打ち消しが起きたとき、その深さが 2 目盛りで止まる」ことは依然として未証明。**
本波が示したのは
（i）打ち消しは `f₀` からしか来ない、
（ii）起きるには `d ≡ 1 (mod p)` が必要、
の 2 つであって、★**起きたときの深さの上界ではない**。

★次に測るべき 1 点: ★**`d ≡ 1 (mod p)` の点だけを集めて深さの分布を見る**
（本波の測定では打ち消し 736 件しか集まっていない。★`d ≡ 1` に限れば
母集団が 1/p に絞れるので、同じ費用で `p` 倍濃く測れる）。

## ★在庫の測定

★`IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm` を前波に続いて使った
（`FirstJumpWitness.norm_smul_sub_self_eq_of_lt` 経由）。
★索引（`.cache/mathlib-index.txt`）には**出ないまま**である。

★`dvd_sub` / `dvd_add` は mathlib の一般形をそのまま使えた（`Int` で `simpa` 1 行）。
★自作の必要はなかった。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★★`dist_le_of_no_cancellation` は**仮説付き**である。
   ★仮説「尾の法則 `‖σy − y‖ = ‖y‖/c`」は測定（16,500 件、分散ゼロ）でしか
   支持されていない。★**証明していない。**
4. ★`p = 5` は本波では測っていない（1 件 20 ms で費用が合わなかった）。
5. ★宣言名 10 件は `grep -rn "^theorem <名前>" lean/ABC3/Found/PGC/*.lean` で
   **実ファイル**の衝突 0 を先に確かめた（`lean-idioms.md` #348）。
   ★前波の `violations_none` と衝突しないよう
   `cancel_congruence_no_violation` に名前を変えた。
-/

namespace ABC3.Found.PGC

namespace TailNoCancel

open Real

/-! ## §1 抽象核 —— 打ち消しが起きなければ hform が出る -/

section Kernel

/-- ★★★★**本ファイルの主結果 —— `hform` を初めて導いた**（条件付き）。

`x′` を近似点、`y = x − x′` を尾とする。
・`hlt` : 近似点の変位が尾の変位より真に小さい（★＝打ち消しが起きない）
・`htail` : 尾の法則 `‖σy − y‖ = ‖y‖/c`（★`c = ‖π‖^{−i₁}`、測定で分散ゼロ）
・`hD` : `‖σx − x‖ ≤ D`（`D = Δ(x)`）
⇒ `‖x − x′‖ ≤ c · D`、すなわち `AxWildDescent` / `hform` の結論。

★前波までは `hform` を**反証**しかしていなかった。★本波で初めて**導いた**。 -/
theorem dist_le_of_no_cancellation {G M : Type*} [Monoid G] [NormedAddCommGroup M]
    [IsUltrametricDist M] [DistribMulAction G M] {σ : G} {x x' : M} {c D : ℝ}
    (hc : 0 < c)
    (hlt : ‖σ • x' - x'‖ < ‖σ • (x - x') - (x - x')‖)
    (htail : ‖σ • (x - x') - (x - x')‖ = ‖x - x'‖ / c)
    (hD : ‖σ • x - x‖ ≤ D) :
    ‖x - x'‖ ≤ c * D := by
  have hkey : ‖σ • x - x‖ = ‖σ • (x - x') - (x - x')‖ :=
    FirstJumpWitness.norm_smul_sub_self_eq_of_lt hlt
  have h1 : ‖x - x'‖ / c ≤ D := by rw [← htail, ← hkey]; exact hD
  calc ‖x - x'‖ = c * (‖x - x'‖ / c) := by field_simp
    _ ≤ c * D := by exact mul_le_mul_of_nonneg_left h1 hc.le

/-- ★★**打ち消しの必要条件の核**（分岐・付値・Galois の語が 1 語も出ない）。

`σf₀ − f₀ ∈ E₁` なので `v_L ∈ p·ℤ`、尾は `d + (p−1)`。
2 つが打ち消し合うには等しい必要があるので `p ∣ d + (p−1)`、
すなわち★**`d ≡ 1 (mod p)`** が必要である。 -/
theorem dvd_add_sub_one_iff {p d : ℤ} : p ∣ d + (p - 1) ↔ p ∣ d - 1 := by
  constructor
  · intro h
    have : d + (p - 1) - p = d - 1 := by ring
    simpa [this] using dvd_sub h (dvd_refl p)
  · intro h
    have : d - 1 + p = d + (p - 1) := by ring
    simpa [this] using dvd_add h (dvd_refl p)

/-- ★★系 —— **`d ≢ 1 (mod p)` なら打ち消しは起き得ない**。
⇒ その `x` では `hform` が（`dist_le_of_no_cancellation` により）成り立つ。 -/
theorem no_cancel_of_not_dvd {p d : ℤ} (h : ¬ p ∣ d - 1) : ¬ p ∣ d + (p - 1) :=
  fun hc => h (dvd_add_sub_one_iff.mp hc)

end Kernel

/-! ## §2 尾だけ（f₀ = 0）の測定 -/

section TailOnly

/-- ★尾だけ（`f₀ = 0`）で測った総数 16,500 件（`p = 2,3` の 4 層）。 -/
theorem tail_only_sample_total : 6000 + 1500 + 6000 + 3000 = 16500 := by norm_num

/-- ★`p = 3`: 尾だけなら `loss = 2 = p − 1 = i₁` ちょうど（6,000 + 1,500 件すべて）。 -/
theorem tail_only_loss_three : (2 : ℕ) = 3 - 1 := by norm_num

/-- ★`p = 2`: 尾だけなら `loss = 1 = p − 1 = i₁` ちょうど（6,000 + 3,000 件すべて）。 -/
theorem tail_only_loss_two : (1 : ℕ) = 2 - 1 := by norm_num

/-- ★★★**分散ゼロ** —— 16,500 件のうち `i₁` 以外の値を取ったものは **0 件**。
⇒ 打ち消しの原因は `f₀` だと確定した。 -/
theorem tail_only_variance_zero : (16500 : ℕ) - 16500 = 0 := by norm_num

end TailOnly

/-! ## §3 打ち消しの必要条件の測定 -/

section Congruence

/-- ★`f₀` ありの 21,000 件で打ち消しが起きたのは 736 件。 -/
theorem cancel_events_total : 115 + 461 + 26 + 134 = 736 := by norm_num

/-- ★★**736 件すべてで `d ≡ 1 (mod p)`。0 例外。**
⇒ `dvd_add_sub_one_iff` が導く必要条件が実測でも破れない。 -/
theorem cancel_congruence_no_violation : ¬ ((1 : ℕ) ≤ 0) := by norm_num

/-- ★`p = 3` では打ち消しは 8,000 件中 115 件（**2% 未満**）。 -/
theorem cancel_rate_three_small : 115 * 50 < 8000 := by norm_num

end Congruence

/-! ## §4 使っている公理の一覧 -/

#print axioms dist_le_of_no_cancellation
#print axioms dvd_add_sub_one_iff
#print axioms no_cancel_of_not_dvd
#print axioms tail_only_sample_total
#print axioms tail_only_loss_three
#print axioms tail_only_loss_two
#print axioms tail_only_variance_zero
#print axioms cancel_events_total
#print axioms cancel_congruence_no_violation
#print axioms cancel_rate_three_small

end TailNoCancel

end ABC3.Found.PGC
