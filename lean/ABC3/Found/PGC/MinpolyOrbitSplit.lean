import ABC3.Found.PGC.RamificationJumpBound
import Mathlib.RingTheory.Adjoin.Polynomial.Basic

/-!
# [pGC] `minpoly K π` の根 = `σ`-軌道、および `AxWildDescent K (axDecay p)` の `k = 1` の供給

`Found/PGC/RamificationJumpBound.lean` は sharp な跳びの上界
`‖(p:L)‖ ≤ ‖π‖^{(p−1)i}` を閉じ、その出口

```
norm_sub_digit_zero_le_rpow_mul :  ‖x − a₀‖ ≤ p^{1/(p−1)} · ‖σ x − x‖
```

の右辺の定数を `AxEpsilonDecay.axDecay p 1` に一致させた。ただしその定理は
4 つの入力を仮定として受け取ったままだった:

1. `hsplit`: `minpoly K π` が `L` 上で `∏_k (X − σ^k π)` と分解すること、
2. `hchar`: `0 < j < p` で `‖(j:L)‖ = 1`、
3. `hval`: `K^×` のノルムが `‖π‖` の `p` 乗部分群に入ること、
4. 桁展開 `x = Σ_{j<p} a_j π^j` が実際に取れること。

本ファイルは 1・2・4 を証明し、3 を「`K` が離散付値で一意化元 `ϖ` を持ち
`‖ϖ‖ = ‖π‖^p`」という 2 行の形に落とす。

## ★★★まず: 配られた字面「`minpoly` の根 = `σ`-軌道」は 偽であった

★配られた持ち場は「`minpoly K π` の根が `L` 上でちょうど `σ`-軌道になる」と書いていた。
★★この字面は `L/K` が Galois でないと偽である。★反例(手計算):

`p` を奇素数、`K = ℚ_p`、`π = p^{1/p}`、`L = K(π)`。`minpoly K π = X^p − p`
(Eisenstein なので既約、次数 `p`)で `L/K` は全分岐 `p` 次。しかし

* `X^p − p` の根は `ζ_p^j π` (`j < p`)、
* `ζ_p ∈ L` なら `p − 1 ∣ [L:K] = p` だが `gcd(p−1, p) = 1` なので `p − 1 = 1`、
  すなわち `p = 2`。★奇素数では `ζ_p ∉ L`。
* よって `L` の中に在る `X^p − p` の根は `π` ただ 1 つで、`Aut(L/K) = 1`。

★`σ` は恒等写像しか取れず、軌道は `{π}`(長さ 1)、根は `p` 個。
★★全分岐 `p` 次でも Galois とは限らない、というのがここで効いている
(`L/K` の Galois 閉包は `K(π, ζ_p)` で `p(p−1)` 次)。

★mathlib 側の裏づけ: `Normal.minpoly_eq_iff_mem_orbit`
(`FieldTheory/Normal/Basic.lean:238`)は
`minpoly F x = minpoly F y ↔ x ∈ MulAction.orbit Gal(E/F) y` を
★`[Normal F E]` の下でのみ主張している。★「根 = 軌道」は正規性の言い換えである。

### ★正しい形 —— Galois より弱い 3 つの仮定で足りる

本ファイルが証明したのは次である(`map_minpoly_eq_prod_iterate`):

```
p 素数、σ : L →ₐ[K] L、σ^[p] π = π、σ π ≠ π、(minpoly K π).natDegree = p
  ⟹  (minpoly K π).map (algebraMap K L) = ∏_{k<p} (X − σ^[k] π)
```

★`Normal` も `IsGalois` も `Separable` も `Finite` も要らない。要るのは
「`σ` が `π` の上で位数ちょうど `p`」と「`[K(π):K] = p`」だけである。
★具体層ではこの 3 つが「`L/K` は `σ` が生成する巡回 `p` 次拡大で `L = K(π)`」から出る。

### ★検算 3 手の結果

1. どの仮定が要るか: 上のとおり。★Galois そのものは要らないが、
   「`σ` の `π` 上の位数が `p`」は要る(上の反例はここが 1 になる)。
2. 重根が無いこと(`Separable`): ★仮定にも結論にも要らなかった。
   `p` 個の相異なる根を持つ次数 `p` の monic 多項式は自動的に分離的であり、
   本ファイルの道(`eq_prod_X_sub_C_of_nodup_of_card`)は
   `Polynomial.splits_iff_card_roots` を通るので分離性を結論として得る。
   ★`RamificationJumpBound.exists_multiset_of_splits` は `Separable` を要求するが、
   ★本ファイルはそれを使わない(仮定が 1 つ減る)。
3. 軌道の長さが `p`: `injOn_iterate_of_prime` が保証する。
   ★退化(`σ π = π`)は仮定で排除し、★組み立ての段では `hjump` と `0 < ‖π‖` から
   導出しているので、呼ぶ側が余分に仮定を書く必要は無い。

## ★★点 2(桁展開の基底の食い違い)を測った —— 食い違いは回避できる

`CyclicJumpNorm.lean` は「`UniformizerExpansion.exists_digits` の基底は `∏_{i<n} σ^i π`
であって `π^j` ではない。貼り合わせるときはここが食い違う」と警告していた。
★測った結果、★★食い違いは重くない。`exists_digits` を使わなければよいからである。

| | `exists_digits`(既存) | 本ファイル `exists_digitSum_of_mem_adjoin` |
|---|---|---|
| 基底 | `uniformizerProd σ π n = ∏_{i<n} σ^i π` | `π^j` |
| 係数 | `𝒪[K]`(整数環) | `K`(体) |
| 項数 | `Finset.Ico N M`(無限級数の切り詰め) | `Finset.range p`(有限、ちょうど `p` 項) |
| 主張の型 | `x − Σ ∈ 𝔪^M`(イデアル) | `x = Σ`(★等式) |

★★つまり `exists_digits` は「任意精度の近似」であって、`CyclicJumpNorm` が要求する
「`K`-ベクトル空間 `L = ⊕_{j<p} K·π^j` の座標」ではない。★後者は
★分岐も付値も使わない: `x = aeval π g` を `minpoly K π` で割った余りに置き換え
(`Polynomial.modByMonic_add_div`)、`Polynomial.aeval_eq_sum_range'` で
`Σ_{j<p} a_j π^j` に直すだけである(★20 行)。
★`uniformizerProd σ π n` と `π^n` は `B` の単数倍しか違わないが、
その単数は `K` に入らないので 2 つの展開は本当に別物である。
★★しかし要るのは後者だけなので、食い違いは障害にならない。

## ★★継ぎ目 —— `WildDepthFieldDescent` との関係(2026-09-08 に地形が変わった)

同じ日に別の波が `WildDepthFieldDescent.axWildDescent_prime` で
★`AxWildDescent K (fun _ => (p:ℝ))` を無条件に証明した。その道は
★`p`-Sylow `P ≤ Q`、`[Q:P] = p` を取って `y := (1/p)·Σ_{c ∈ Q/P} c•x` と平均する形で、
損失は `‖(p:E)‖⁻¹ = p` ちょうど(`WildDepthDescent.exists_natDegree_minpoly_descent_div`)。

★★本ファイルはその `p` を `k = 1` の段で `p^{1/(p−1)} = axDecay p 1` に絞る側である。
★平均の道と本ファイルの道は数学が別物なので、そのまま差し替えることはできない。
★★継ぎ目に何が要るかを正確に書いておく:

| `exists_natDegree_minpoly_descent_div`(平均) | 本ファイル(桁展開) |
|---|---|
| 仮定は `E/F` 有限次 Galois だけ | `L = K(π)`・`[L:K] = p`・`σ` が `π` 上で位数 `p` |
| 分岐を一切見ない | 全分岐(`hval`)と跳び `i ≥ 1` を見る |
| 損失 `‖(p:E)‖⁻¹ = p` | 損失 `p^{1/(p−1)}`(★`p ≥ 3` で真に良い) |
| `x` はどこにあってもよい | ★`x ∈ K(π)` が要る |

★★つまり差し替えの障害は定数ではなく **`x` が全分岐巡回 `p` 次の層に入ること**である。
★これが「点 1(塔の分解)」であり、本ファイルは**そこには手を付けていない**。
★★平均の道が「中間体で 1 段下がる」枠組みを回避したのと同じ理由で、
点 1 も「`x` を含む全分岐巡回 `p` 次の層を取る」という形では一般には取れない
(順分岐の部分と非可換な部分が残る)。★次の波はここを測ること。

## ★★抽象核(分岐・付値・Galois・p 進の語彙が 1 語も出ない)

| 宣言 | 内容 | 語彙 |
|---|---|---|
| ★`eq_of_iterate_eq_of_coprime` | `f^[n]a = a`, `gcd(k,n)=1`, `f^[k]a = a` ⟹ `f a = a` | 反復と `Nat` だけ |
| `iterate_ne_of_lt_of_prime` | `n` 素数・`f a ≠ a` なら `f^[j]a ≠ f^[k]a` (`j<k<n`) | 同上 |
| ★`injOn_iterate_of_prime` | 軌道 `{f^[k]a : k < n}` はちょうど `n` 元 | 同上 |
| ★`eq_prod_X_sub_C_of_nodup_of_card` | monic で「相異なる根が次数だけある」⟹ `∏(X−a)` | 多項式だけ |
| ★`norm_natCast_eq_one_of_coprime` | `‖(p:L)‖<1`・`gcd(j,p)=1` ⟹ `‖(j:L)‖=1` | 超距離ノルム体 + Bézout |
| `norm_iterate_sub_self_eq_of_coprime_of_fix` | 跳びが生成元の取り方に依らない(★`f^[p]π = π` だけで) | 超距離 + 反復 |

★1 本目が心臓である。「位数 `p` の巡回群には非自明な部分群が無い」を
★群も部分群も使わず、`Function.iterate` と `Nat.ModEq.pow_totient` だけで書いた。

★★`RamificationJumpBound` §7 の `norm_iterate_sub_self_eq_of_coprime` は
`hid : f^[p] = id`(`L` 全体で恒等)を要求していたが、本ファイルは
★`f^[p] π = π`(`π` 1 点でだけ)に弱めた。★`σ` が `L` 全体で位数 `p` である必要は無い。

## 在庫の測定(自分で測った。★MCP は 1 度も呼んでいない。コマンドを残す)

```
grep -nE "splits_iff_card_roots|card_roots'|eq_of_le_of_card_le|Multiset\.le_iff_subset" \
  .cache/mathlib-index.txt
  → ★Polynomial.splits_iff_card_roots (Splits.lean:376) —— これが要だった。
     「根の個数 = 次数」から Splits が出るので、分離性を仮定せずに済む。
     Multiset.le_iff_subset / eq_of_le_of_card_le / Polynomial.card_roots' も在る。
grep -nE "prodXSubSMul|FixedPoints\.minpoly|Normal\.splits|minpoly.*orbit" .cache/mathlib-index.txt
  → ★Normal.minpoly_eq_iff_mem_orbit (Normal/Basic.lean:238) が「根 = 軌道」そのもの。
     ★ただし [Normal F E] を要求する。★これが「配られた字面が偽」の裏づけになった。
  → prodXSubSMul (Algebra/Polynomial/GroupRingAction.lean:82) も在るが、
     ★こちらは G が R に作用する有限群であることを要求し、
     結論も「係数が固定される」までで minpoly との一致は別に要る。
     ★本ファイルの道(根の個数だけ数える)の方が仮定が少ないので使わなかった。
grep -nE "norm_natCast_le_one|norm_intCast_le_one" .cache/mathlib-index.txt
  → ★IsUltrametricDist.norm_natCast_le_one / norm_intCast_le_one (Ring/Ultra.lean)。
     ★★索引の行は `lemma norm_natCast_le_one (n : ℕ) : ‖(n : R)‖ ≤ 1` と表示するが、
     ★実際には section variable `(R)` が明示引数として先頭に付く。
     索引は宣言の字面しか持たないので section variable が見えない(下の逸脱 6 を見よ)。
grep -nE "modByMonic_lt|modByMonic_add_div|adjoin_singleton_eq_range_aeval" .cache/mathlib-index.txt
  → ★Polynomial.natDegree_modByMonic_lt / modByMonic_add_div /
     Algebra.adjoin_singleton_eq_range_aeval。★点 2 はこの 3 本で 20 行で閉じた。
grep -cE "\.<名前>(\.|\t)" .cache/decl-index.txt (本ファイルの全宣言名について)
  → すべて 0(衝突なし)。
```

★使わなかった在庫とその理由:

* `RamificationJumpBound.exists_multiset_of_splits` —— ★`Separable` と `Splits` を
  要求する。本ファイルは根の個数から `Splits` を結論するので、両方要らない。
  ★在庫を使わない方が仮定が減る例である。
* `Polynomial.Splits.eq_prod_roots_of_monic` —— ★こちらは使った
  (`eq_prod_X_sub_C_of_nodup_of_card` の最後の 1 行)。配られた受け口の名前は当たっていた。
* `UniformizerExpansion.exists_digits` —— 上の表のとおり別物なので使わなかった。

## ★逸脱の記録(CLAUDE.md「逸脱」)

1. ★原典(Serre / Ax)は Galois 理論の言葉で「根 = 共役 = 軌道」と書く。
   本ファイルは Galois を仮定しない形(`σ^[p]π = π` と `[K(π):K] = p`)に分解した。
   ★★配られた字面(「全分岐 `p` 次なら根 = `σ`-軌道」)は偽なので、
   これは単なる一般化ではなく訂正である。理由は本ファイル冒頭に書いた。
2. ★`Separable` を仮定しない(標数 0 なので自動だが、そもそも要らなかった)。
   ★結果として本ファイルは標数の仮定を 1 つも置いていない。
3. ★ノルム言語で書いた(`CyclicJumpNorm` / `RamificationJumpBound` と同じ理由)。
4. ★`hval`(値群が `p` 乗部分群)は仮定のままである。本ファイルはそれを
   「`K` の離散性 + `‖ϖ‖ = ‖π‖^p`」に落としただけで、
   ★`K` が離散付値体であること自体は具体層(`PAdicLocalField`)の仕事である。
5. ★`‖(p:L)‖ = (p:ℝ)⁻¹` は仮定のままである。具体層では
   `Padic.norm_p`(idiom #291)と「ノルムが `ℚ_p` のそれを延長する」から出る。
6. ★`.src` は `RamificationJumpBound.lean` と同じ項目(pGC 物理 p.6 Corollary 3.1)を指す。
   本ファイルの中身は `AxSenTate` の入力であって、原典が独立に立てた項目ではない。

## ★退化の自己検査

* `p = 2`, `K = ℚ₂`, `L = ℚ₂(√2)`, `π = √2`, `σπ = −√2`: `σ^[2]π = π`、`σπ ≠ π`、
  `minpoly = X²−2` で次数 2 = `p`。★仮定はすべて満たされ、結論は
  `X²−2 = (X−√2)(X+√2)`。★非空虚である。
* `σ π ≠ π` を落とすと偽: 軌道が縮んで根の個数が足りない(上の `ℚ_p(p^{1/p})` の反例)。
* `(minpoly K π).natDegree = p` を落とすと偽: `π ∈ K` なら次数 1 で軌道も 1 点。
* `p` が素数でないと偽: `p = 4` で `σ^[2] π = π` の場合、軌道は 2 点しかないのに
  次数は 4 なので根が足りない。★`injOn_iterate_of_prime` の
  「`0 < k < p` なら `gcd(k,p) = 1`」がここで効いている。
-/

namespace ABC3.Found.PGC

open Polynomial

/-! ## §1 抽象核(反復と `Nat.Coprime` だけ) —— 「位数 `p` の巡回群に非自明な部分群は無い」 -/

section IterateCore

variable {α : Type*}

/-- ★★抽象核: `f^[n] a = a` かつ `gcd(k, n) = 1` かつ `f^[k] a = a` ならば `f a = a`。

★「`a` の固定部分群が `ℤ/n` の部分群で、`k` を含むなら全体」の群を使わない版。
★`k^{φ(n)} ≡ 1 (mod n)`(`Nat.ModEq.pow_totient`)を使って
`f^[k^{φ(n)}] a = f^[1 + n·q] a = f (f^[n·q] a) = f a` と `f^[k^{φ(n)}] a = a` を突き合わせる。
★群論・体論・分岐・付値の語彙が 1 語も出てこない。 -/
theorem eq_of_iterate_eq_of_coprime {f : α → α} {a : α} {n k : ℕ}
    (hn : 1 < n) (hfix : f^[n] a = a) (hk : Nat.Coprime k n) (hka : f^[k] a = a) :
    f a = a := by
  have ht : 0 < Nat.totient n := Nat.totient_pos.mpr (by omega)
  have hkt : k ^ Nat.totient n ≡ 1 [MOD n] := Nat.ModEq.pow_totient hk
  have h1 : k ^ Nat.totient n % n = 1 := by
    have h := hkt
    rwa [Nat.ModEq, Nat.mod_eq_of_lt hn] at h
  have hq : k ^ Nat.totient n = 1 + n * (k ^ Nat.totient n / n) := by
    conv_lhs => rw [← Nat.div_add_mod (k ^ Nat.totient n) n]
    rw [h1]; ring
  have hsplit : k ^ Nat.totient n = k * k ^ (Nat.totient n - 1) := by
    rw [← pow_succ']
    congr 1
    omega
  have hA : f^[k ^ Nat.totient n] a = a := by
    rw [hsplit, Function.iterate_mul]
    exact Function.iterate_fixed hka _
  rw [hq, Function.iterate_add_apply, Function.iterate_mul,
    Function.iterate_fixed hfix] at hA
  simpa using hA

/-- 抽象核: `n` が素数で `f^[n] a = a`・`f a ≠ a` なら、`j < k < n` について
`f^[j] a ≠ f^[k] a`。★軌道は `n` 個の相異なる点からなる。 -/
theorem iterate_ne_of_lt_of_prime {f : α → α} {a : α} {n : ℕ} (hn : n.Prime)
    (hfix : f^[n] a = a) (hne : f a ≠ a) {j k : ℕ} (hjk : j < k) (hkn : k < n) :
    f^[j] a ≠ f^[k] a := by
  intro h
  have h1 : f^[n - k + j] a = a := by
    have h2 := congrArg (fun y => f^[n - k] y) h
    simp only [← Function.iterate_add_apply] at h2
    rwa [show n - k + k = n by omega, hfix] at h2
  have hnd : ¬ n ∣ (n - k + j) :=
    fun hd => absurd (Nat.le_of_dvd (by omega) hd) (by omega)
  have hcop : Nat.Coprime (n - k + j) n :=
    ((Nat.Prime.coprime_iff_not_dvd hn).mpr hnd).symm
  exact hne (eq_of_iterate_eq_of_coprime hn.one_lt hfix hcop h1)

/-- ★抽象核: 軌道 `k ↦ f^[k] a`(`k < n`)は単射。★「軌道の長さがちょうど `n`」。 -/
theorem injOn_iterate_of_prime {f : α → α} {a : α} {n : ℕ} (hn : n.Prime)
    (hfix : f^[n] a = a) (hne : f a ≠ a) :
    Set.InjOn (fun k => f^[k] a) (Set.Iio n) := by
  intro j hj k hk h
  simp only [Set.mem_Iio] at hj hk
  rcases lt_trichotomy j k with hlt | heq | hgt
  · exact absurd h (iterate_ne_of_lt_of_prime hn hfix hne hlt hk)
  · exact heq
  · exact absurd h.symm (iterate_ne_of_lt_of_prime hn hfix hne hgt hj)

end IterateCore

/-! ## §2 抽象核(多項式だけ / 超距離ノルム体だけ) -/

/-- ★★抽象核: monic `f` が「相異なる根を次数だけ」持つなら `f = ∏(X − a)`。

★`Polynomial.splits_iff_card_roots`(根の個数 = 次数 ⟺ `Splits`)を通るので、
★分離性を仮定せずに結論として得る。`RamificationJumpBound.exists_multiset_of_splits`
より仮定が 1 つ少ない。 -/
theorem eq_prod_X_sub_C_of_nodup_of_card {L : Type*} [Field L] {f : L[X]} (hf : f.Monic)
    {S : Multiset L} (hnd : S.Nodup) (hroot : ∀ a ∈ S, f.IsRoot a)
    (hcard : Multiset.card S = f.natDegree) :
    f = (S.map fun a => X - C a).prod := by
  have hsub : S ⊆ f.roots := fun a ha => Polynomial.mem_roots'.mpr ⟨hf.ne_zero, hroot a ha⟩
  have hle : S ≤ f.roots := (Multiset.le_iff_subset hnd).mpr hsub
  have hcard2 : Multiset.card f.roots ≤ Multiset.card S := by
    rw [hcard]; exact Polynomial.card_roots' f
  have heq : S = f.roots := Multiset.eq_of_le_of_card_le hle hcard2
  have hsp : f.Splits := Polynomial.splits_iff_card_roots.mpr (by rw [← heq, hcard])
  rw [heq]
  exact hsp.eq_prod_roots_of_monic hf

/-- ★抽象核: 超距離ノルム体で `‖(p:L)‖ < 1` かつ `gcd(j,p) = 1` なら `‖(j:L)‖ = 1`。

★Bézout `1 = j·u + p·v`(`Nat.gcd_eq_gcd_ab`)を超距離で潰すだけ。
★p 進も剰余体も出てこない。 -/
theorem norm_natCast_eq_one_of_coprime {L : Type*} [NormedField L] [IsUltrametricDist L]
    {p j : ℕ} (hp : ‖(p : L)‖ < 1) (h : Nat.Coprime j p) : ‖(j : L)‖ = 1 := by
  have hbez : ((Nat.gcd j p : ℕ) : ℤ) = (j : ℤ) * Nat.gcdA j p + (p : ℤ) * Nat.gcdB j p :=
    Nat.gcd_eq_gcd_ab j p
  rw [h] at hbez
  have hL : (1 : L)
      = (j : L) * ((Nat.gcdA j p : ℤ) : L) + (p : L) * ((Nat.gcdB j p : ℤ) : L) := by
    have h2 := congrArg (fun z : ℤ => (z : L)) hbez
    push_cast at h2
    simpa using h2
  rcases eq_or_lt_of_le (IsUltrametricDist.norm_natCast_le_one L j) with heq | hlt
  · exact heq
  · exfalso
    have h1 : ‖(1 : L)‖ < 1 := by
      rw [hL]
      refine lt_of_le_of_lt (IsUltrametricDist.norm_add_le_max _ _) (max_lt ?_ ?_)
      · rw [norm_mul]
        calc ‖(j : L)‖ * ‖((Nat.gcdA j p : ℤ) : L)‖ ≤ ‖(j : L)‖ * 1 :=
              mul_le_mul_of_nonneg_left (IsUltrametricDist.norm_intCast_le_one L _)
                (norm_nonneg _)
          _ < 1 := by rw [mul_one]; exact hlt
      · rw [norm_mul]
        calc ‖(p : L)‖ * ‖((Nat.gcdB j p : ℤ) : L)‖ ≤ ‖(p : L)‖ * 1 :=
              mul_le_mul_of_nonneg_left (IsUltrametricDist.norm_intCast_le_one L _)
                (norm_nonneg _)
          _ < 1 := by rw [mul_one]; exact hp
    simp at h1

/-- ★`RamificationJumpBound.norm_sub_digit_zero_le_rpow_mul` の `hchar` そのもの。 -/
theorem norm_natCast_eq_one_of_lt_prime {L : Type*} [NormedField L] [IsUltrametricDist L]
    {p : ℕ} (hp : p.Prime) (hnorm : ‖(p : L)‖ < 1) {j : ℕ} (hj0 : 0 < j) (hjp : j < p) :
    ‖(j : L)‖ = 1 :=
  norm_natCast_eq_one_of_coprime hnorm
    (((Nat.Prime.coprime_iff_not_dvd hp).mpr
      (fun hd => absurd (Nat.le_of_dvd hj0 hd) (by omega))).symm)

/-! ## §3 具体層 —— `minpoly K π` の根がちょうど `σ`-軌道 -/

section Orbit

variable {K L : Type*} [Field K] [Field L] [Algebra K L]

/-- `K`-代数準同型の反復は `minpoly K π` の根を根に写す。★`σ` に全単射性も等長性も要らない。 -/
theorem aeval_iterate_minpoly_eq_zero (σ : L →ₐ[K] L) (π : L) (k : ℕ) :
    Polynomial.aeval ((σ : L → L)^[k] π) (minpoly K π) = 0 := by
  induction k with
  | zero => rw [Function.iterate_zero_apply]; exact minpoly.aeval K π
  | succ k ih =>
    rw [Function.iterate_succ_apply', Polynomial.aeval_algHom_apply σ, ih, map_zero]

/-- ★★★本ファイルの主結果 1 —— `minpoly K π` は `L` 上で `σ`-軌道に分解する。

```
p 素数、σ : L →ₐ[K] L、σ^[p] π = π、σ π ≠ π、(minpoly K π).natDegree = p
  ⟹  (minpoly K π).map (algebraMap K L) = ∏_{k<p} (X − σ^[k] π)
```

★★配られた字面は「全分岐 `p` 次なら根 = `σ`-軌道」だったが、
★それは `L/K` が Galois でないと偽である(ファイル冒頭の反例)。
★正しい形は上のとおりで、`Normal` も `IsGalois` も `Separable` も `Finite` も使わない。

★`σ` が `π` の上で位数ちょうど `p` であること(`σ^[p]π = π` かつ `σπ ≠ π`)から
軌道が `p` 個の相異なる点になり(`injOn_iterate_of_prime`)、
それが次数 `p` の monic 多項式の根をすべて尽くす(`eq_prod_X_sub_C_of_nodup_of_card`)。 -/
theorem map_minpoly_eq_prod_iterate {p : ℕ} (hp : p.Prime) (σ : L →ₐ[K] L) {π : L}
    (hfix : (σ : L → L)^[p] π = π) (hne : σ π ≠ π)
    (hdeg : (minpoly K π).natDegree = p) :
    (minpoly K π).map (algebraMap K L)
      = (((Multiset.range p).map fun k => (σ : L → L)^[k] π).map fun a => X - C a).prod := by
  have hint : IsIntegral K π := by
    by_contra hc
    rw [minpoly.eq_zero hc, Polynomial.natDegree_zero] at hdeg
    exact hp.ne_zero hdeg.symm
  have hmonic : (minpoly K π).Monic := minpoly.monic hint
  have hM : ((minpoly K π).map (algebraMap K L)).Monic := hmonic.map _
  have hinj : ∀ x ∈ Multiset.range p, ∀ y ∈ Multiset.range p,
      (σ : L → L)^[x] π = (σ : L → L)^[y] π → x = y := by
    intro x hx y hy h
    exact injOn_iterate_of_prime hp hfix hne (Multiset.mem_range.mp hx)
      (Multiset.mem_range.mp hy) h
  refine eq_prod_X_sub_C_of_nodup_of_card hM
    (Multiset.Nodup.map_on hinj (Multiset.nodup_range p)) ?_ ?_
  · intro a ha
    obtain ⟨k, _, rfl⟩ := Multiset.mem_map.mp ha
    show ((minpoly K π).map (algebraMap K L)).IsRoot _
    rw [Polynomial.IsRoot, Polynomial.eval_map, ← Polynomial.aeval_def]
    exact aeval_iterate_minpoly_eq_zero σ π k
  · rw [Multiset.card_map, Multiset.card_range, hmonic.natDegree_map, hdeg]

/-- ★`RamificationJumpBound.norm_natCast_le_pow_of_prod_X_sub_C` が要求する
`π ::ₘ T` の形にほどいたもの。`T` の元は `0 < k < p` の `σ^[k] π` である。 -/
theorem exists_multiset_iterate_of_minpoly {p : ℕ} (hp : p.Prime) (σ : L →ₐ[K] L) {π : L}
    (hfix : (σ : L → L)^[p] π = π) (hne : σ π ≠ π)
    (hdeg : (minpoly K π).natDegree = p) :
    ∃ T : Multiset L,
      (minpoly K π).map (algebraMap K L) = ((π ::ₘ T).map fun b => X - C b).prod ∧
      ∀ b ∈ T, ∃ k : ℕ, 0 < k ∧ k < p ∧ b = (σ : L → L)^[k] π := by
  classical
  set S : Multiset L := (Multiset.range p).map fun k => (σ : L → L)^[k] π with hS
  have hinj : ∀ x ∈ Multiset.range p, ∀ y ∈ Multiset.range p,
      (σ : L → L)^[x] π = (σ : L → L)^[y] π → x = y := by
    intro x hx y hy h
    exact injOn_iterate_of_prime hp hfix hne (Multiset.mem_range.mp hx)
      (Multiset.mem_range.mp hy) h
  have hnd : S.Nodup := Multiset.Nodup.map_on hinj (Multiset.nodup_range p)
  have hmemπ : π ∈ S := by
    rw [hS]
    exact Multiset.mem_map.mpr ⟨0, Multiset.mem_range.mpr hp.pos, by simp⟩
  have hcount : Multiset.count π (S.erase π) = 0 := by
    rw [Multiset.count_erase_self]
    have := Multiset.nodup_iff_count_le_one.mp hnd π
    omega
  have hnot : π ∉ S.erase π := Multiset.count_eq_zero.mp hcount
  refine ⟨S.erase π, ?_, ?_⟩
  · rw [Multiset.cons_erase hmemπ, hS]
    exact map_minpoly_eq_prod_iterate hp σ hfix hne hdeg
  · intro b hb
    have hbS : b ∈ S := Multiset.mem_of_mem_erase hb
    rw [hS] at hbS
    obtain ⟨k, hk, rfl⟩ := Multiset.mem_map.mp hbS
    refine ⟨k, ?_, Multiset.mem_range.mp hk, rfl⟩
    rcases Nat.eq_zero_or_pos k with rfl | hpos
    · rw [Function.iterate_zero_apply] at hb
      exact absurd hb hnot
    · exact hpos

end Orbit

/-! ## §4 抽象核 —— 跳びの一様性を「`f^[p] π = π`(1 点でだけ)」から出す

★`RamificationJumpBound` §7 は `hid : f^[p] = id`(`L` 全体で恒等)を要求していた。
★ここではそれを `π` 1 点での固定に弱める。証明の骨は同じ望遠鏡である。 -/

section IterateNorm

variable {L : Type*} [NormedField L] [IsUltrametricDist L]

/-- ★抽象核: 等長で差を保つ `f` が `f^[p] π = π` を満たし `gcd(k,p) = 1` なら
`‖f^[k]π − π‖ = ‖fπ − π‖`。★`f^[p] = id` は要らない。 -/
theorem norm_iterate_sub_self_eq_of_coprime_of_fix {f : L → L}
    (hsub : ∀ x y, f (x - y) = f x - f y) (hiso : ∀ y, ‖f y‖ = ‖y‖)
    {p : ℕ} (hp : 1 < p) {π : L} (hfix : f^[p] π = π) {k : ℕ} (hk : Nat.Coprime k p) :
    ‖f^[k] π - π‖ = ‖f π - π‖ := by
  have ht : 0 < Nat.totient p := Nat.totient_pos.mpr (by omega)
  have hkt : k ^ Nat.totient p ≡ 1 [MOD p] := Nat.ModEq.pow_totient hk
  have h1 : k ^ Nat.totient p % p = 1 := by
    have h := hkt
    rwa [Nat.ModEq, Nat.mod_eq_of_lt hp] at h
  have hq : k ^ Nat.totient p = 1 + p * (k ^ Nat.totient p / p) := by
    conv_lhs => rw [← Nat.div_add_mod (k ^ Nat.totient p) p]
    rw [h1]; ring
  have hiterπ : f^[k ^ Nat.totient p] π = f π := by
    rw [hq, Function.iterate_add_apply, Function.iterate_mul, Function.iterate_fixed hfix]
    simp
  have hsplit : k ^ Nat.totient p = k * k ^ (Nat.totient p - 1) := by
    rw [← pow_succ']
    congr 1
    omega
  refine le_antisymm (norm_iterate_sub_self_le_of_isometry hsub hiso π k) ?_
  have h2 : ‖(f^[k])^[k ^ (Nat.totient p - 1)] π - π‖ ≤ ‖f^[k] π - π‖ :=
    norm_iterate_sub_self_le_of_isometry (iterate_map_sub_of_map_sub hsub k)
      (norm_iterate_of_norm_eq hiso k) π _
  rwa [← Function.iterate_mul, ← hsplit, hiterπ] at h2

/-- ★`p` 素数・`0 < k < p` 版。`norm_natCast_le_pow_of_prod_X_sub_C` の `hbreak` を作る。 -/
theorem norm_iterate_sub_self_eq_of_lt_of_fix {f : L → L}
    (hsub : ∀ x y, f (x - y) = f x - f y) (hiso : ∀ y, ‖f y‖ = ‖y‖)
    {p : ℕ} (hp : p.Prime) {π : L} (hfix : f^[p] π = π) {i k : ℕ}
    (hk0 : 0 < k) (hkp : k < p) (hjump : ‖f π - π‖ = ‖π‖ ^ (i + 1)) :
    ‖f^[k] π - π‖ = ‖π‖ ^ (i + 1) := by
  have hnd : ¬ p ∣ k := fun h => absurd (Nat.le_of_dvd hk0 h) (by omega)
  have hcop : Nat.Coprime k p := ((Nat.Prime.coprime_iff_not_dvd hp).mpr hnd).symm
  rw [norm_iterate_sub_self_eq_of_coprime_of_fix hsub hiso hp.one_lt hfix hcop, hjump]

end IterateNorm

/-! ## §5 供給 —— `hval` を「`K` の一意化元 `ϖ` と `‖ϖ‖ = ‖π‖^n`」に落とす -/

section Supply

variable {K L : Type*} [Field K] [NormedField L] [Algebra K L]

/-- ★`CyclicJumpNorm` / `RamificationJumpBound` の `hval` の供給。

`K` の値群が `‖ϖ‖` で生成され(`hK`)、`‖ϖ‖ = ‖π‖^n`(= 分岐指数が `n`)なら、
`K^×` のノルムはすべて `‖π‖` の `n` 乗部分群に入る。
★これが「`L/K` は全分岐で `e = n`」のノルム言語での中身の全部である。 -/
theorem norm_algebraMap_mem_zpow_of_uniformizer {n : ℕ} {π : L} {ϖ : K}
    (hϖ : ‖algebraMap K L ϖ‖ = ‖π‖ ^ (n : ℤ))
    (hK : ∀ c : K, c ≠ 0 → ∃ m : ℤ, ‖algebraMap K L c‖ = ‖algebraMap K L ϖ‖ ^ m) :
    ∀ c : K, c ≠ 0 → ∃ m : ℤ, ‖algebraMap K L c‖ = ‖π‖ ^ ((n : ℤ) * m) := by
  intro c hc
  obtain ⟨m, hm⟩ := hK c hc
  exact ⟨m, by rw [hm, hϖ, ← zpow_mul]⟩

end Supply

/-! ## §6 桁展開の供給(★点 2) —— `x ∈ K(π)` なら `x = Σ_{j<p} a_j π^j`

★`UniformizerExpansion.exists_digits`(基底 `∏σ^iπ`、係数 `𝒪[K]`、イデアルでの近似)
とは別物である。要るのは `K`-ベクトル空間としての座標の方で、
★こちらは分岐も付値も使わずに `modByMonic` だけで出る。 -/

section Digits

variable {K L : Type*} [Field K] [Field L] [Algebra K L]

/-- ★★点 2 の解決: `x ∈ K[π]` なら `x = digitSum p π a`(`a : ℕ → K`)。

★`x = aeval π g` を `minpoly K π` で割った余り `g %ₘ minpoly K π`(次数 `< p`)に
置き換え、`Polynomial.aeval_eq_sum_range'` で `Σ_{j<p} a_j π^j` に直すだけ。
★★分岐・付値・p 進の語彙が 1 語も出てこない。 -/
theorem exists_digitSum_of_mem_adjoin {p : ℕ} {π : L} (hint : IsIntegral K π)
    (hdeg : (minpoly K π).natDegree = p) {x : L} (hx : x ∈ Algebra.adjoin K ({π} : Set L)) :
    ∃ a : ℕ → K, x = digitSum p π a := by
  rw [Algebra.adjoin_singleton_eq_range_aeval] at hx
  obtain ⟨g, rfl⟩ := hx
  have hm : (minpoly K π).Monic := minpoly.monic hint
  have hne1 : minpoly K π ≠ 1 := by
    intro h
    have hpos : 0 < (minpoly K π).natDegree := minpoly.natDegree_pos hint
    rw [h, Polynomial.natDegree_one] at hpos
    omega
  have hr : Polynomial.aeval π (g %ₘ minpoly K π) = Polynomial.aeval π g := by
    conv_rhs => rw [← Polynomial.modByMonic_add_div g (minpoly K π)]
    simp [minpoly.aeval]
  have hdlt : (g %ₘ minpoly K π).natDegree < p := by
    rw [← hdeg]; exact Polynomial.natDegree_modByMonic_lt g hm hne1
  refine ⟨fun j => (g %ₘ minpoly K π).coeff j, ?_⟩
  show (Polynomial.aeval π) g = _
  rw [digitSum, ← hr, Polynomial.aeval_eq_sum_range' hdlt]
  exact Finset.sum_congr rfl fun j _ => by rw [Algebra.smul_def]

/-- ★`L = K(π)`(`Algebra.adjoin K {π} = ⊤`)なら `L` のすべての元が桁展開を持つ。 -/
theorem exists_digitSum_of_adjoin_eq_top {p : ℕ} {π : L} (hint : IsIntegral K π)
    (hdeg : (minpoly K π).natDegree = p) (htop : Algebra.adjoin K ({π} : Set L) = ⊤)
    (x : L) : ∃ a : ℕ → K, x = digitSum p π a :=
  exists_digitSum_of_mem_adjoin hint hdeg (by rw [htop]; trivial)

end Digits

/-! ## §7 組み立て —— `AxWildDescent K (axDecay p)` の `k = 1` の 1 段 -/

section Assemble

variable {K L : Type*} [Field K] [NormedField L] [IsUltrametricDist L] [Algebra K L]

/-- ★★`RamificationJumpBound.norm_natCast_le_pow_of_prod_X_sub_C` に
`hsplit` と `hbreak` を供給した形。

`(p−1)·i ≤ e_L`(ノルム言語で `‖(p:L)‖ ≤ ‖π‖^{(p−1)i}`)が
★「`σ` が `π` 上で位数 `p`」「`σ` が等長」「`[K(π):K] = p`」「`hval`」だけから出る。
★`hne : σπ ≠ π` は `hjump` と `0 < ‖π‖` から導出しているので呼ぶ側は書かなくてよい。 -/
theorem norm_natCast_le_pow_of_orbit {p i : ℕ} (hp : p.Prime) (σ : L →ₐ[K] L) {π : L}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hiso : ∀ y : L, ‖σ y‖ = ‖y‖)
    (hfix : (σ : L → L)^[p] π = π) (hjump : ‖σ π - π‖ = ‖π‖ ^ (i + 1))
    (hdeg : (minpoly K π).natDegree = p)
    (hval : ∀ c : K, c ≠ 0 → ∃ m : ℤ, ‖algebraMap K L c‖ = ‖π‖ ^ ((p : ℤ) * m)) :
    ‖(p : L)‖ ≤ ‖π‖ ^ ((p - 1) * i) := by
  have hne : σ π ≠ π := by
    intro h
    rw [h, sub_self, norm_zero] at hjump
    exact absurd hjump.symm (ne_of_gt (pow_pos hπ0 _))
  have hint : IsIntegral K π := by
    by_contra hc
    rw [minpoly.eq_zero hc, Polynomial.natDegree_zero] at hdeg
    exact hp.ne_zero hdeg.symm
  obtain ⟨T, hT, hTmem⟩ := exists_multiset_iterate_of_minpoly hp σ hfix hne hdeg
  refine norm_natCast_le_pow_of_prod_X_sub_C hp.pos hπ0 hπ1 hval (minpoly.monic hint) hdeg
    hT ?_
  intro b hb
  obtain ⟨k, hk0, hkp, rfl⟩ := hTmem b hb
  rw [norm_sub_rev]
  exact norm_iterate_sub_self_eq_of_lt_of_fix (f := (σ : L → L))
    (fun x y => map_sub σ x y) hiso hp hfix hk0 hkp hjump

/-- ★★★本ファイルの主結果 2 —— 桁展開の元について 1 段の損失が `axDecay p 1` 以下。

```
‖x − a₀‖ ≤ p^{1/(p−1)} · ‖σ x − x‖        (x = Σ_{j<p} a_j π^j)
```

★右辺の定数 `p^{1/(p−1)}` は `AxEpsilonDecay.axDecay p 1` そのもの
(`AxEpsilonDecay.axDecay_one`)。
★`RamificationJumpBound.norm_sub_digit_zero_le_rpow_mul` の 4 つの仮定のうち
`hsplit`(点 4)と `hchar`(点 3)を本ファイルが埋めたので、
★残る仮定は `hval` と `‖(p:L)‖ = (p:ℝ)⁻¹` の 2 つだけになった。 -/
theorem norm_sub_digit_zero_le_axDecay_of_orbit {p i : ℕ} (hp : p.Prime) (σ : L →ₐ[K] L)
    {π : L} {a : ℕ → K} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hi : 0 < i)
    (hiso : ∀ y : L, ‖σ y‖ = ‖y‖)
    (hfix : (σ : L → L)^[p] π = π) (hjump : ‖σ π - π‖ = ‖π‖ ^ (i + 1))
    (hdeg : (minpoly K π).natDegree = p)
    (hval : ∀ c : K, c ≠ 0 → ∃ m : ℤ, ‖algebraMap K L c‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hnormp : ‖(p : L)‖ = ((p : ℝ))⁻¹) :
    ‖digitSum p π a - algebraMap K L (a 0)‖
      ≤ (p : ℝ) ^ (1 / ((p : ℝ) - 1)) * ‖σ (digitSum p π a) - digitSum p π a‖ := by
  have hp1 : 1 < p := hp.one_lt
  have hplt : ‖(p : L)‖ < 1 := by
    rw [hnormp]
    have h2 : (1 : ℝ) < (p : ℝ) := by exact_mod_cast hp1
    rw [inv_lt_one_iff₀]
    exact Or.inr h2
  have hchar : ∀ j : ℕ, 0 < j → j < p → ‖(j : L)‖ = 1 :=
    fun j hj0 hjp => norm_natCast_eq_one_of_lt_prime hp hplt hj0 hjp
  have hne : σ π ≠ π := by
    intro h
    rw [h, sub_self, norm_zero] at hjump
    exact absurd hjump.symm (ne_of_gt (pow_pos hπ0 _))
  have hint : IsIntegral K π := by
    by_contra hc
    rw [minpoly.eq_zero hc, Polynomial.natDegree_zero] at hdeg
    exact hp.ne_zero hdeg.symm
  obtain ⟨T, hT, hTmem⟩ := exists_multiset_iterate_of_minpoly hp σ hfix hne hdeg
  refine norm_sub_digit_zero_le_rpow_mul hp1 σ hπ0 hπ1 hi hjump hchar hval hnormp
    (minpoly.monic hint) hdeg hT ?_
  intro b hb
  obtain ⟨k, hk0, hkp, rfl⟩ := hTmem b hb
  rw [norm_sub_rev]
  exact norm_iterate_sub_self_eq_of_lt_of_fix (f := (σ : L → L))
    (fun x y => map_sub σ x y) hiso hp hfix hk0 hkp hjump

/-- ★★★`k = 1` の 1 段、`x` を仮定しない形。

`L = K(π)` が「`σ` が `π` 上で位数 `p`・等長・`[K(π):K] = p`」を満たすとき、
★`L` の任意の元 `x` について `K` の元 `c` が取れて

```
‖x − c‖ ≤ p^{1/(p−1)} · ‖σ x − x‖ = axDecay p 1 · ‖σ x − x‖.
```

★★これが `AxWildDescent K (axDecay p)` の `k = 1` の段が要求する評価そのものである
(残るのは `K.closure` / `K.absGal` 側の配管:
「wild 深さ 1 の `x` から巡回 `p` 次の 1 段を切り出す塔の分解」)。 -/
theorem exists_norm_sub_algebraMap_le_axDecay_of_orbit {p i : ℕ} (hp : p.Prime)
    (σ : L →ₐ[K] L) {π : L} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hi : 0 < i)
    (hiso : ∀ y : L, ‖σ y‖ = ‖y‖)
    (hfix : (σ : L → L)^[p] π = π) (hjump : ‖σ π - π‖ = ‖π‖ ^ (i + 1))
    (hdeg : (minpoly K π).natDegree = p)
    (htop : Algebra.adjoin K ({π} : Set L) = ⊤)
    (hval : ∀ c : K, c ≠ 0 → ∃ m : ℤ, ‖algebraMap K L c‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hnormp : ‖(p : L)‖ = ((p : ℝ))⁻¹) (x : L) :
    ∃ c : K, ‖x - algebraMap K L c‖
      ≤ (p : ℝ) ^ (1 / ((p : ℝ) - 1)) * ‖σ x - x‖ := by
  have hint : IsIntegral K π := by
    by_contra hc
    rw [minpoly.eq_zero hc, Polynomial.natDegree_zero] at hdeg
    exact hp.ne_zero hdeg.symm
  obtain ⟨a, rfl⟩ := exists_digitSum_of_adjoin_eq_top hint hdeg htop x
  exact ⟨a 0, norm_sub_digit_zero_le_axDecay_of_orbit hp σ hπ0 hπ1 hi hiso hfix hjump
    hdeg hval hnormp⟩

end Assemble

/-! ## §8 `.src` -/

def map_minpoly_eq_prod_iterate.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_multiset_iterate_of_minpoly.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_natCast_eq_one_of_lt_prime.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_algebraMap_mem_zpow_of_uniformizer.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_digitSum_of_mem_adjoin.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_natCast_le_pow_of_orbit.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_sub_digit_zero_le_axDecay_of_orbit.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_norm_sub_algebraMap_le_axDecay_of_orbit.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §9 使っている公理の一覧 -/

#print axioms eq_of_iterate_eq_of_coprime
#print axioms iterate_ne_of_lt_of_prime
#print axioms injOn_iterate_of_prime
#print axioms eq_prod_X_sub_C_of_nodup_of_card
#print axioms norm_natCast_eq_one_of_coprime
#print axioms norm_natCast_eq_one_of_lt_prime
#print axioms aeval_iterate_minpoly_eq_zero
#print axioms map_minpoly_eq_prod_iterate
#print axioms exists_multiset_iterate_of_minpoly
#print axioms norm_iterate_sub_self_eq_of_coprime_of_fix
#print axioms norm_iterate_sub_self_eq_of_lt_of_fix
#print axioms norm_algebraMap_mem_zpow_of_uniformizer
#print axioms exists_digitSum_of_mem_adjoin
#print axioms exists_digitSum_of_adjoin_eq_top
#print axioms norm_natCast_le_pow_of_orbit
#print axioms norm_sub_digit_zero_le_axDecay_of_orbit
#print axioms exists_norm_sub_algebraMap_le_axDecay_of_orbit

end ABC3.Found.PGC
