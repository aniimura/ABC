import ABC3.Found.PGC.CyclicJumpNorm
import Mathlib.Algebra.Polynomial.Derivative
import Mathlib.Algebra.Polynomial.AlgebraMap
import Mathlib.Algebra.Polynomial.BigOperators
import Mathlib.Algebra.Polynomial.Splits
import Mathlib.FieldTheory.Separable
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Data.Nat.Totient
import Mathlib.FieldTheory.Finite.Basic

/-!
# [pGC] sharp な跳びの上界 `(p−1)·i ≤ e_L`

`L/K` を p 進局所体の全分岐な巡回 `p` 次拡大、`i` をその(下付き番号付けの)跳び
(`G_i = G`, `G_{i+1} = 1`、同値に `v_L(σπ_L − π_L) = i+1`)とするとき

```
(p − 1) · i  ≤  e_L        (e_L = v_L(p) = p·e_K、絶対分岐指数)
```

これは Serre, Corps Locaux III §6 Prop. 13(different の上界 `d ≤ e − 1 + v_L(e)`)を
`e = p` の場合に当てたものである。`Found/PGC/CyclicJumpNorm.lean` の
`norm_sub_digit_zero_eq_zpow_mul` が出す 1 段の最良定数 `‖π_L‖^{−i}` を
`axDecay p 1 = p^{1/(p−1)}` で押さえるのがこの不等式の役目である。

## ★★まず: 主張は本当か(自分で検算した。★配られた形はそのまま真であった)

| 例 | `p` | `e_K` | `e_L` | 跳び `i` | `(p−1)i` | 判定 |
|---|---|---|---|---|---|---|
| `ℚ₂(√2)/ℚ₂` | 2 | 1 | 2 | 2 | 2 | ★等号 |
| `ℚ₂(√−1)/ℚ₂` | 2 | 1 | 2 | 1 | 1 | 真(狭義) |
| `ℚ₂(√3)/ℚ₂` | 2 | 1 | 2 | 1 | 1 | 真(狭義) |
| `ℚ_p(ζ_{p²})/ℚ_p(ζ_p)` | p | p−1 | p(p−1) | p−1 | (p−1)² | 真(狭義、差 p−1) |

`ℚ₂(√2)`: `π = √2`、`σπ = −√2`、`v_L(σπ−π) = v_L(2√2) = 3 = i+1` なので `i = 2`、
`(p−1)i = 2 = e_L`。★等号が実際に起きるので、この不等式はこれ以上強くできない。

`ℚ_p(ζ_{p²})/ℚ_p(ζ_p)`: `Gal(ℚ_p(ζ_{p^n})/ℚ_p)` の下付き分岐群は
`G_u = Gal(·/ℚ_p(ζ_p))` for `1 ≤ u ≤ p−1` なので、部分群への制限(`H_u = H ∩ G_u`)から
跳びは `i = p−1`。`(p−1)² < p(p−1) = e_L`。

★等号条件: `(p−1)i = e_L = p·e_K` なら `(p−1) ∣ p·e_K` と `gcd(p−1,p)=1` から
`(p−1) ∣ e_K`、よって `i = p·(e_K/(p−1))` は `p` の倍数。★対偶を取ると
`p ∤ i` ならば狭義の不等式(`sub_one_mul_lt_of_not_dvd`、本ファイルで証明した)。

## ★通った道 —— (A) different 経由。ただし `differentIdeal` を一切通らない

配られた候補は (A) different / (B) 単数の norm / (C) 直接展開の 3 つだった。
★(B) は測って外れたので記録しておく:

`u = σπ/π = 1 + c`(`v_L(c) = i`)は `N_{L/K}(u) = σ^p(π)/π = 1` を満たすので
`Σ_{m=1}^{p} e_m(c, σc, …, σ^{p−1}c) = 0`。ここで各 `e_m ∈ K` だから
`‖e_m‖ ∈ ‖π‖^{pℤ}` で、`‖e_p‖ = ‖N(c)‖ = ‖π‖^{pi}`。
★しかし `e_m ≈ C(p,m)c^m` の誤差項が `‖σc − c‖ ≤ ‖π‖^{2i}` としか押さえられず、
`m = 1` の項が `‖π‖^{2i}` まで大きくなりうる。`2i < pi`(`p ≥ 3`)なので
`e_p` が最小になることが言えず、望遠鏡和では最小項の分離が起きない。
★同じ理由で「`σ^p(π) = π` を足し合わせる」素朴な (C) も
`Tr(σπ−π) = 0` が恒等的に真になるだけで `i ≤ e_L` しか出ない
(それは `CyclicJumpNorm.lean` が既に「既存と同じ定数」と記録している弱い形)。

★実際に効いたのは (A) だが、原典 (Serre III §6 Prop 13) の道とは違う。
原典は different イデアル `𝔡 = 𝔡_{L/K}` を trace 双対で定義し、
`𝔡 = (f'(π))`(単項生成)と `d ≤ e−1+v_L(e)` を別々に証明する。
★本ファイルは different を一度も定義せず、`f'(π)` を 2 通りに測るだけで閉じる:

1. 根の側: `f(X) = ∏_{k<p}(X − σ^k π)` から `f'(π) = ∏_{k≠0}(π − σ^kπ)`、
   よって `‖f'(π)‖ = (‖π‖^{i+1})^{p−1}`。
2. 係数の側: `f` は `K` 係数の monic なので `f'(π) = Σ_{l<p} f'_l π^l` は
   `K` 係数の桁展開であり、`CyclicJumpNorm.nnnorm_sum_digit_eq_sup`(桁の分離)から
   `‖f'(π)‖ = max_l ‖f'_l π^l‖ ≥ ‖f'_{p−1} π^{p−1}‖ = ‖(p:L)‖·‖π‖^{p−1}`
   (`f'_{p−1} = p·f_p = p`、monic だから)。

`(‖π‖^{i+1})^{p−1} ≥ ‖(p:L)‖·‖π‖^{p−1}` の両辺を `‖π‖^{p−1}` で割ると
`‖(p:L)‖ ≤ ‖π‖^{(p−1)i}` ——これが主結果である。

★★原典より短くなった点: 原典の `d ≤ e−1+v_L(e)` は「相異なる剰余類 mod e の項は
相殺しない」という補題を Eisenstein 多項式に当てて示すが、★その補題は
この木では `nnnorm_sum_digit_eq_sup` として既に在庫にあった(前波 CyclicJumpNorm)。
つまり「Eisenstein であること」すら要らず、★monic であることだけで足りる
(`f` の下位係数が `K` に在れば桁は分離する)。★実際に本ファイルは
`Eisenstein`・`differentIdeal`・`Ideal.ramificationIdx`・付値環のどれも import していない。

## ★★抽象核(分岐・付値・Galois・p 進の語彙が 1 語も出ない)

| 宣言 | 内容 | 語彙 |
|---|---|---|
| `eval_derivative_of_eq_X_sub_C_mul` | `F = (X−a)·g ⟹ F'(a) = g(a)` | 可換環だけ |
| `norm_multiset_prod_map_eq_pow` | ノルムが一定な族の積のノルム = `c^{card}` | ノルム体だけ |
| ★`norm_natCast_mul_pow_le_norm_aeval_derivative` | `‖(n:L)‖·‖π‖^{n−1} ≤ ‖f'(π)‖` | 超距離 + 桁の分離 |
| `rpow_inv_natCast_le_of_le_pow` | `c ≤ a^m ⟹ c^{1/m} ≤ a` | 実数だけ |
| `sub_one_mul_lt_of_not_dvd` | `(p−1)i ≤ pe` かつ `p ∤ i` ⟹ 狭義 | ℕ だけ |
| ★`norm_iterate_sub_self_eq_of_coprime` | `f^[p] = id`, `gcd(k,p)=1` ⟹ `‖f^[k]π−π‖ = ‖fπ−π‖` | 超距離 + 反復だけ |

★3 本目が心臓である。仮定は `CyclicJumpNorm` と同じ 1 つ
(`hval`: `K^×` のノルムが `‖π‖` の `n` 乗部分群に入る = 「全分岐で `e = n`」)だけで、
★`f` が Eisenstein である必要も、`n` が素数である必要も無い。

## 在庫の測定(自分で測った。★MCP は 0 回。コマンドを残す)

```
grep -n "derivative_mul\b\|derivative_prod\b\|derivative_X_sub_C\b" .cache/mathlib-index.txt
  → Polynomial.derivative_mul / derivative_X_sub_C / derivative_prod(Multiset 版)
grep -n "derivative.*prod_X_sub_C\|prod_X_sub_C.*derivative" .cache/mathlib-index.txt
  → ★Polynomial.eval_multiset_prod_X_sub_C_derivative(Derivative.lean:680)
     `eval r (derivative (∏(X−a))) = ∏_{a ∈ S.erase r}(r−a)` ——★在る。
     ただし `Multiset.erase` は `DecidableEq L` を要求するので、
     本ファイルは使わず `F = (X−π)·g` の 3 行(`eval_derivative_of_eq_X_sub_C_mul`)に
     置き換えた。★その方が仮定が減る(可換環で成り立ち、`DecidableEq` が要らない)。
grep -n "aeval_eq_sum_range|coeff_derivative\b|natDegree_derivative_lt" .cache/mathlib-index.txt
  → ★Polynomial.aeval_eq_sum_range'(natDegree < n なら aeval = Σ_{i<n} coeff i • x^i)
     ——★これが「桁展開に直す」作業を丸ごと消した(再添字が要らない)。
     Polynomial.coeff_derivative / natDegree_derivative_le
grep -n "natDegree_multiset_prod_X_sub_C_eq_card" .cache/mathlib-index.txt
  → ★在る(BigOperators.lean:328)。根の個数 = 次数がこれ 1 行で出る。
grep -n "eq_prod_roots" .cache/mathlib-index.txt
  → ★Polynomial.Splits.eq_prod_roots_of_monic(Splits.lean:310)。
     ★索引に `eq_prod_roots_of_monic_of_splits_id` は無い(名前が変わっている)。
grep -n "normHom|Multiset.prod_hom" .cache/mathlib-index.txt
  → ★normHom : α →*₀ ℝ(Normed/Ring/Basic.lean:720)+ Multiset.prod_hom'。
     「積のノルム = ノルムの積」を Multiset で使うにはこの 2 つ。
grep -n "differentIdeal|Ideal.different" .cache/mathlib-index.txt
  → ★differentIdeal は在る(DedekindDomain/Different.lean:476)。
     `aeval_derivative_mem_differentIdeal`・`pow_sub_one_dvd_differentIdeal` も在る。
     ★しかし `d ≤ e − 1 + v_L(e)`(Serre III §6 Prop 13 の上界)は
     索引に見当たらない(`pow_sub_one_dvd_differentIdeal` は下からの評価)。
     ★★本ファイルはこの不在を回避した——different を通らない道で閉じたので、
     ★「上界が無い」ことは障害にならない。
grep -n "RamificationJumpBound|norm_natCast_le_pow" .cache/decl-index.txt
  → 0 件(名前の衝突なし)。
```

## ★逸脱の記録(CLAUDE.md「逸脱」)

1. ★原典(Serre)は付値言語(`v_L(𝔡) = d`)で書く。本ファイルはノルム言語で書いた。
   `CyclicJumpNorm.lean` と同じ理由(`AxDescentStep` がノルム言語)。
   跳び `i` は `‖π − σ^kπ‖ = ‖π‖^{i+1}` として現れる。
2. ★「`L/K` は全分岐な巡回 `p` 次拡大」という分岐理論の仮定を一切置いていない。
   代わりに `hval`(値群が `n` 乗部分群)と「`f` の根がすべて `π` から距離 `‖π‖^{i+1}`」
   というノルムと代数だけの仮定に分解した。★その結果:
   * `n` が素数である必要が無い、
   * `f` が Eisenstein である必要が無い(monic で十分)、
   * `σ` も `Gal` も現れない(根の多重集合 `T` だけ)。
   ★原典より弱い仮定である。
3. ★主結果は「すべての非自明な `σ^k` が同じ跳びを持つ」(= `G_{i+1} = 1` が部分群で
   あること)を仮定として受け取る(`hbreak`)。★ただし §7 でその仮定を
   別に証明してある(`norm_iterate_sub_self_eq_of_lt`)ので、宣言としては分けたが
   穴は開いていない。★§7 は「`f` が等長で差を保ち `f^[p] = id`」だけを使い、
   分岐・付値・Galois・多項式の語彙が 1 語も出ない。
   ★残る接続点は「`minpoly` の根の多重集合が `{σ^k π}` である」という Galois の側だけ。
4. `.src` は `CyclicJumpNorm.lean` と同じ項目(pGC 物理 p.6 Corollary 3.1)を指す。
   本ファイルの中身は `AxSenTate` の入力であって、原典が独立に立てた項目ではない
   (数学的な典拠は Serre, Corps Locaux III §6 Prop. 13 / IV §1 Prop. 4)。

## ★退化の自己検査

* `n = 1` のとき: `T = 0`、結論は `‖(1:L)‖ ≤ ‖π‖^0`、すなわち `1 ≤ 1`。★空虚に真。
* `i = 0`(順分岐)のとき: 結論は `‖(n:L)‖ ≤ 1`。★真だが情報が無い。
  暴分岐(`i ≥ 1`)でこそ意味がある。
* `hbreak` を落とすと偽: 根が `π` から遠いと `‖f'(π)‖` が大きくなり、
  `‖(n:L)‖` の上界が得られない。実際 `f = X^n − a`(`a ∈ K` で `‖a‖` が大きい)を
  考えれば分かる。
* `hval` を落とすと偽: 桁の分離が壊れ、`f'(π)` の中で相殺が起きうる。
  ★これが「全分岐」を使っている唯一の場所である。
-/

namespace ABC3.Found.PGC

open Polynomial

/-! ## §1 抽象核(可換環) —— `F = (X−a)·g` なら `F'(a) = g(a)` -/

/-- 抽象核: `F = (X − a)·g` ならば `F'(a) = g(a)`。

★分岐も付値もノルムも出てこない。証明は `derivative_mul` 1 回。
★mathlib の `Polynomial.eval_multiset_prod_X_sub_C_derivative` は同じ内容を
根の多重集合の言葉で持っているが、`Multiset.erase` のために `DecidableEq` を要求する。
こちらは要求しない。 -/
theorem eval_derivative_of_eq_X_sub_C_mul {R : Type*} [CommRing R] {F g : R[X]} {a : R}
    (h : F = (X - C a) * g) : eval a (derivative F) = eval a g := by
  subst h
  rw [derivative_mul, derivative_X_sub_C]
  simp

/-- 抽象核: ノルムが一定 `c` の族の積のノルムは `c ^ (要素数)`。 -/
theorem norm_multiset_prod_map_eq_pow {L : Type*} [NormedField L] {S : Multiset L}
    {g : L → L} {c : ℝ} (h : ∀ a ∈ S, ‖g a‖ = c) :
    ‖(S.map g).prod‖ = c ^ (Multiset.card S) := by
  have h1 : (S.map fun a => ‖g a‖).prod = ‖(S.map g).prod‖ :=
    Multiset.prod_hom' (f := normHom (α := L)) S g
  rw [← h1]
  have h2 : (S.map fun a => ‖g a‖) = Multiset.replicate (Multiset.card S) c :=
    Multiset.eq_replicate.mpr ⟨by simp, fun b hb => by
      obtain ⟨a, ha, rfl⟩ := Multiset.mem_map.mp hb; exact h a ha⟩
  rw [h2, Multiset.prod_replicate]

/-! ## §2 抽象核(超距離 + 桁の分離) —— `‖(n:L)‖·‖π‖^{n−1} ≤ ‖f'(π)‖` -/

variable {K L : Type*} [Field K] [NormedField L] [IsUltrametricDist L] [Algebra K L]

/-- ★★心臓。`f ∈ K[X]` が monic で `deg f = n` なら

`‖(n : L)‖ · ‖π‖^{n−1} ≤ ‖f'(π)‖`。

`f'(π) = Σ_{l<n} f'_l π^l` は `K` 係数の桁展開で、`hval` により桁のノルムは
相異なる(`nnnorm_sum_digit_eq_sup`)。したがって和のノルムは最大値であり、
特に `l = n−1` の桁 `f'_{n−1} = n·f_n = n` の項以上である。

★仮定は `hval`(`K^×` のノルムが `‖π‖` の `n` 乗部分群に入る)だけ。
★Eisenstein である必要も `n` が素数である必要も無い。 -/
theorem norm_natCast_mul_pow_le_norm_aeval_derivative {n : ℕ} (hn : 0 < n) {π : L}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hval : ∀ c : K, c ≠ 0 → ∃ m : ℤ, ‖algebraMap K L c‖ = ‖π‖ ^ ((n : ℤ) * m))
    {f : K[X]} (hf : f.Monic) (hdeg : f.natDegree = n) :
    ‖(n : L)‖ * ‖π‖ ^ (n - 1) ≤ ‖aeval π (derivative f)‖ := by
  classical
  have hd : (derivative f).natDegree < n := by
    have := Polynomial.natDegree_derivative_le f
    omega
  have hsum : aeval π (derivative f)
      = ∑ j ∈ Finset.range n, algebraMap K L ((derivative f).coeff j) * π ^ j := by
    rw [Polynomial.aeval_eq_sum_range' hd]
    exact Finset.sum_congr rfl fun j _ => by rw [Algebra.smul_def]
  have hsup := nnnorm_sum_digit_eq_sup (K := K) (n := n) hπ0 hπ1 hval
    (fun j => (derivative f).coeff j) (s := Finset.range n)
    (fun j hj => Finset.mem_range.mp hj)
  have h1 : n - 1 + 1 = n := by omega
  have hc1 : f.coeff n = 1 := by rw [← hdeg]; exact hf.coeff_natDegree
  have h2 : ((n - 1 : ℕ) : K) + 1 = (n : K) := by
    have hx : ((n - 1 + 1 : ℕ) : K) = (n : K) := by rw [h1]
    simpa using hx
  have hcoeff : (derivative f).coeff (n - 1) = (n : K) := by
    rw [Polynomial.coeff_derivative, h1, hc1, one_mul, h2]
  have hmem : n - 1 ∈ Finset.range n := Finset.mem_range.mpr (by omega)
  have hle : ‖algebraMap K L ((derivative f).coeff (n - 1)) * π ^ (n - 1)‖₊
      ≤ ‖∑ j ∈ Finset.range n, algebraMap K L ((derivative f).coeff j) * π ^ j‖₊ := by
    rw [hsup]
    exact Finset.le_sup (f := fun j => ‖algebraMap K L ((derivative f).coeff j) * π ^ j‖₊) hmem
  rw [hcoeff] at hle
  have h3 := (NNReal.coe_le_coe).mpr hle
  simpa [hsum, norm_mul, norm_pow] using h3

/-! ## §3 根の側 —— `f'(π) = ∏_{a ∈ T}(π − a)` -/

omit [IsUltrametricDist L] in
/-- `f` が `L` 上で `∏_{a ∈ π ::ₘ T}(X − a)` と分解するなら `f'(π) = ∏_{a ∈ T}(π − a)`。 -/
theorem aeval_derivative_eq_prod_of_split {π : L} {f : K[X]} {T : Multiset L}
    (hsplit : f.map (algebraMap K L) = ((π ::ₘ T).map fun a => X - C a).prod) :
    aeval π (derivative f) = (T.map fun a => π - a).prod := by
  have h1 : aeval π (derivative f) = eval π (derivative (f.map (algebraMap K L))) := by
    rw [derivative_map, eval_map, aeval_def]
  rw [h1, hsplit, Multiset.map_cons, Multiset.prod_cons]
  rw [eval_derivative_of_eq_X_sub_C_mul (a := π) rfl, eval_multiset_prod, Multiset.map_map]
  simp

/-! ## §4 主結果 —— `‖(n:L)‖ ≤ ‖π‖^{(n−1)i}` -/

/-- ★★★主結果(ノルム言語)。

`f ∈ K[X]` は monic・`deg f = n`、`L` 上で `∏_{a ∈ π ::ₘ T}(X − a)` と分解し、
`π` 以外の根 `a ∈ T` がすべて `‖π − a‖ = ‖π‖^{i+1}`(= 跳びが `i`)を満たすとき

`‖(n : L)‖ ≤ ‖π‖^{(n−1)·i}`。

★`n = p`、`‖(p:L)‖ = ‖π‖^{e_L}` を入れると `(p−1)·i ≤ e_L`
(`sub_one_mul_le_of_norm_natCast_eq_pow`)。 -/
theorem norm_natCast_le_pow_of_prod_X_sub_C {n i : ℕ} (hn : 0 < n) {π : L}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hval : ∀ c : K, c ≠ 0 → ∃ m : ℤ, ‖algebraMap K L c‖ = ‖π‖ ^ ((n : ℤ) * m))
    {f : K[X]} (hf : f.Monic) (hdeg : f.natDegree = n) {T : Multiset L}
    (hsplit : f.map (algebraMap K L) = ((π ::ₘ T).map fun a => X - C a).prod)
    (hbreak : ∀ a ∈ T, ‖π - a‖ = ‖π‖ ^ (i + 1)) :
    ‖(n : L)‖ ≤ ‖π‖ ^ ((n - 1) * i) := by
  have hcard : Multiset.card T = n - 1 := by
    have h := congrArg Polynomial.natDegree hsplit
    rw [hf.natDegree_map, hdeg, Polynomial.natDegree_multiset_prod_X_sub_C_eq_card] at h
    simp at h
    omega
  have heval := aeval_derivative_eq_prod_of_split hsplit
  have hnorm : ‖aeval π (derivative f)‖ = (‖π‖ ^ (i + 1)) ^ (n - 1) := by
    rw [heval, norm_multiset_prod_map_eq_pow hbreak, hcard]
  have hmain := norm_natCast_mul_pow_le_norm_aeval_derivative (K := K) hn hπ0 hπ1 hval hf hdeg
  rw [hnorm] at hmain
  have hrw : (‖π‖ ^ (i + 1)) ^ (n - 1) = ‖π‖ ^ ((n - 1) * i) * ‖π‖ ^ (n - 1) := by
    rw [← pow_mul, ← pow_add]
    ring_nf
  rw [hrw] at hmain
  exact le_of_mul_le_mul_right hmain (pow_pos hπ0 _)

omit [IsUltrametricDist L] in
/-- monic・separable・`L` 上で分解する `f` の根の多重集合を `π ::ₘ T` の形にほどく。
`T` の元は `f` の根であって `π` とは異なる。 -/
theorem exists_multiset_of_splits {f : K[X]} (hf : f.Monic) (hsep : f.Separable) {π : L}
    (hroot : aeval π f = 0) (hsp : (f.map (algebraMap K L)).Splits) :
    ∃ T : Multiset L, f.map (algebraMap K L) = ((π ::ₘ T).map fun a => X - C a).prod ∧
      ∀ a ∈ T, (f.map (algebraMap K L)).IsRoot a ∧ a ≠ π := by
  classical
  have hFm : (f.map (algebraMap K L)).Monic := hf.map _
  have hFroot : (f.map (algebraMap K L)).IsRoot π := by
    rw [IsRoot, eval_map, ← aeval_def, hroot]
  have hmem : π ∈ (f.map (algebraMap K L)).roots :=
    Polynomial.mem_roots'.mpr ⟨hFm.ne_zero, hFroot⟩
  have hnd : (f.map (algebraMap K L)).roots.Nodup := Polynomial.nodup_roots hsep.map
  have hcount : Multiset.count π ((f.map (algebraMap K L)).roots.erase π) = 0 := by
    rw [Multiset.count_erase_self]
    have := Multiset.nodup_iff_count_le_one.mp hnd π
    omega
  have hnot : π ∉ (f.map (algebraMap K L)).roots.erase π := Multiset.count_eq_zero.mp hcount
  refine ⟨(f.map (algebraMap K L)).roots.erase π, ?_, ?_⟩
  · rw [Multiset.cons_erase hmem]
    exact hsp.eq_prod_roots_of_monic hFm
  · intro a ha
    exact ⟨(Polynomial.mem_roots'.mp (Multiset.mem_of_mem_erase ha)).2,
      fun h => hnot (h ▸ ha)⟩

/-- ★主結果の `Splits` 版。根の多重集合を自分で作らなくてよい。 -/
theorem norm_natCast_le_pow_of_splits {n i : ℕ} (hn : 0 < n) {π : L}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hval : ∀ c : K, c ≠ 0 → ∃ m : ℤ, ‖algebraMap K L c‖ = ‖π‖ ^ ((n : ℤ) * m))
    {f : K[X]} (hf : f.Monic) (hdeg : f.natDegree = n) (hsep : f.Separable)
    (hroot : aeval π f = 0) (hsp : (f.map (algebraMap K L)).Splits)
    (hbreak : ∀ a : L, (f.map (algebraMap K L)).IsRoot a → a ≠ π →
      ‖π - a‖ = ‖π‖ ^ (i + 1)) :
    ‖(n : L)‖ ≤ ‖π‖ ^ ((n - 1) * i) := by
  obtain ⟨T, hT, hTmem⟩ := exists_multiset_of_splits hf hsep hroot hsp
  exact norm_natCast_le_pow_of_prod_X_sub_C hn hπ0 hπ1 hval hf hdeg hT
    (fun a ha => hbreak a (hTmem a ha).1 (hTmem a ha).2)

/-! ## §5 付値言語への翻訳と、`AxWildDescent` が使う実数の形 -/

omit [IsUltrametricDist L] in
/-- ★★配られた形そのもの: `‖(n:L)‖ = ‖π‖^e`(= `e = v_L(n)`)なら `(n−1)·i ≤ e`。

`n = p` で `e = e_L = v_L(p)` を入れると `(p−1)·i ≤ e_L`。 -/
theorem sub_one_mul_le_of_norm_natCast_eq_pow {n i e : ℕ} {π : L}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (he : ‖(n : L)‖ = ‖π‖ ^ e)
    (h : ‖(n : L)‖ ≤ ‖π‖ ^ ((n - 1) * i)) : (n - 1) * i ≤ e := by
  rw [he] at h
  exact (pow_le_pow_iff_right_of_lt_one₀ hπ0 hπ1).mp h

/-- 抽象核(ℕ だけ): `(p−1)·i ≤ p·e` かつ `p ∤ i` なら狭義の不等式。

★等号 `(p−1)i = p·e` は `p ∣ (p−1)i` を強いるが `gcd(p, p−1) = 1` なので `p ∣ i`。
★つまり等号が起きるのは `p ∣ i` のときだけである
(`ℚ₂(√2)/ℚ₂` は `p = 2`, `i = 2` で確かに `p ∣ i`)。 -/
theorem sub_one_mul_lt_of_not_dvd {p i e : ℕ} (hp : 0 < p) (hnd : ¬ p ∣ i)
    (h : (p - 1) * i ≤ p * e) : (p - 1) * i < p * e := by
  rcases lt_or_eq_of_le h with h' | h'
  · exact h'
  · exfalso
    have hdvd : p ∣ (p - 1) * i := ⟨e, h'⟩
    have hcop : Nat.Coprime p (p - 1) := by
      have h1 : p - (p - 1) = 1 := by omega
      rw [← Nat.coprime_sub_self_left (Nat.sub_le p 1), h1]
      exact Nat.coprime_one_left _
    exact hnd (hcop.dvd_of_dvd_mul_left hdvd)

/-- 抽象核(実数だけ): `0 ≤ c ≤ a^m`(`m ≥ 1`)なら `c^{1/m} ≤ a`。 -/
theorem rpow_inv_natCast_le_of_le_pow {a c : ℝ} (ha : 0 ≤ a) (hc : 0 ≤ c) {m : ℕ}
    (hm : 0 < m) (h : c ≤ a ^ m) : c ^ ((m : ℝ))⁻¹ ≤ a := by
  have hm0 : (0 : ℝ) < (m : ℝ) := by exact_mod_cast hm
  have h1 : c ^ ((m : ℝ))⁻¹ ≤ (a ^ m) ^ ((m : ℝ))⁻¹ :=
    Real.rpow_le_rpow hc h (by positivity)
  rwa [← Real.rpow_natCast a m, ← Real.rpow_mul ha, mul_inv_cancel₀ (ne_of_gt hm0),
    Real.rpow_one] at h1

omit [IsUltrametricDist L] in
/-- ★★★`AxWildDescent` の `k = 1` が使う形。

`‖(p:L)‖ = (p:ℝ)⁻¹`(p 進の正規化)のもとで、主結果 `‖(p:L)‖ ≤ ‖π‖^{(p−1)i}` は

`‖π‖^{−i} ≤ p^{1/(p−1)} = axDecay p 1`

を与える。★`CyclicJumpNorm.norm_sub_digit_zero_eq_zpow_mul` の右辺の係数が
ちょうど `‖π‖^{−i}` なので、この 1 本で 1 段の損失が `axDecay p 1` で押さえられる。 -/
theorem zpow_neg_jump_le_rpow {p i : ℕ} (hp : 1 < p) {π : L} (hπ0 : 0 < ‖π‖)
    (hnorm : ‖(p : L)‖ = ((p : ℝ))⁻¹)
    (h : ‖(p : L)‖ ≤ ‖π‖ ^ ((p - 1) * i)) :
    ‖π‖ ^ (-(i : ℤ)) ≤ (p : ℝ) ^ (1 / ((p : ℝ) - 1)) := by
  have hp0 : (0 : ℝ) < (p : ℝ) := by positivity
  have hcast : ((p - 1 : ℕ) : ℝ) = (p : ℝ) - 1 := by
    have h1 : (1 : ℕ) ≤ p := le_of_lt hp
    push_cast [Nat.cast_sub h1]
    ring
  have hai : (0 : ℝ) < ‖π‖ ^ i := pow_pos hπ0 i
  have hpow : ‖π‖ ^ ((p - 1) * i) = (‖π‖ ^ i) ^ (p - 1) := by
    rw [← pow_mul, mul_comm]
  rw [hnorm, hpow] at h
  have hstep := rpow_inv_natCast_le_of_le_pow (le_of_lt hai) (by positivity)
    (m := p - 1) (by omega) h
  rw [Real.inv_rpow (le_of_lt hp0), hcast] at hstep
  have hrpos : (0 : ℝ) < (p : ℝ) ^ ((p : ℝ) - 1)⁻¹ := Real.rpow_pos_of_pos hp0 _
  have hinv : (‖π‖ ^ i)⁻¹ ≤ (p : ℝ) ^ ((p : ℝ) - 1)⁻¹ := by
    have h2 := inv_anti₀ (inv_pos.mpr hrpos) hstep
    rwa [inv_inv] at h2
  rw [zpow_neg, zpow_natCast, one_div]
  exact hinv

/-! ## §6 `CyclicJumpNorm` との接続 —— 1 段の損失が `axDecay p 1` で押さえられる -/

/-- ★★★`CyclicJumpNorm.norm_sub_digit_zero_eq_zpow_mul`(1 段の最良定数 `‖π‖^{−i}`)と
本ファイルの `‖(p:L)‖ ≤ ‖π‖^{(p−1)i}` を掛け合わせた形。

`d(x, K) = ‖π‖^{−i}·‖σx − x‖ ≤ p^{1/(p−1)}·‖σx − x‖`

★右辺の定数 `p^{1/(p−1)}` は `AxEpsilonDecay.axDecay p 1` そのものである
(`AxEpsilonDecay.axDecay_one`)。★これが `AxWildDescent K (axDecay p)` の
`k = 1` の段で必要な評価の全部である(残りは具体層の配管:
`x = Σ_{j<p} a_j π^j` の桁展開と、`f = minpoly` の分解を作ること)。 -/
theorem norm_sub_digit_zero_le_rpow_mul {p i : ℕ} (hp : 1 < p) {π : L} {a : ℕ → K}
    (σ : L →ₐ[K] L) (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hi : 0 < i)
    (hjump : ‖σ π - π‖ = ‖π‖ ^ (i + 1))
    (hchar : ∀ j, 0 < j → j < p → ‖(j : L)‖ = 1)
    (hval : ∀ c : K, c ≠ 0 → ∃ m : ℤ, ‖algebraMap K L c‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hnormp : ‖(p : L)‖ = ((p : ℝ))⁻¹)
    {f : K[X]} (hf : f.Monic) (hdeg : f.natDegree = p) {T : Multiset L}
    (hsplit : f.map (algebraMap K L) = ((π ::ₘ T).map fun b => X - C b).prod)
    (hroots : ∀ b ∈ T, ‖π - b‖ = ‖π‖ ^ (i + 1)) :
    ‖digitSum p π a - algebraMap K L (a 0)‖
      ≤ (p : ℝ) ^ (1 / ((p : ℝ) - 1)) * ‖σ (digitSum p π a) - digitSum p π a‖ := by
  have hbound : ‖(p : L)‖ ≤ ‖π‖ ^ ((p - 1) * i) :=
    norm_natCast_le_pow_of_prod_X_sub_C (by omega) hπ0 hπ1 hval hf hdeg hsplit hroots
  have hconst : ‖π‖ ^ (-(i : ℤ)) ≤ (p : ℝ) ^ (1 / ((p : ℝ) - 1)) :=
    zpow_neg_jump_le_rpow hp hπ0 hnormp hbound
  rw [norm_sub_digit_zero_eq_zpow_mul σ (by omega) hπ0 hπ1 hi hjump hchar hval]
  exact mul_le_mul_of_nonneg_right hconst (norm_nonneg _)

/-! ## §7 接続点のための抽象核 —— 「跳びは `σ` の取り方に依らない」

★逸脱の記録 3 が言っている `hbreak`(すべての非自明な `σ^k` が同じ跳びを持つ)を、
★分岐・付値・Galois・多項式のどれも使わずに証明する。
必要なのは「`σ` が等長で差を保ち、`p` 回で恒等になる」ことだけである。 -/

section Iterate

variable {L : Type*} [NormedField L] [IsUltrametricDist L]

omit [IsUltrametricDist L] in
/-- 差を保つ写像の反復も差を保つ。 -/
theorem iterate_map_sub_of_map_sub {f : L → L} (hsub : ∀ x y, f (x - y) = f x - f y) :
    ∀ (k : ℕ) (x y : L), f^[k] (x - y) = f^[k] x - f^[k] y := by
  intro k
  induction k with
  | zero => simp
  | succ k ih =>
    intro x y
    rw [Function.iterate_succ_apply, Function.iterate_succ_apply, Function.iterate_succ_apply,
      hsub, ih]

omit [IsUltrametricDist L] in
/-- 等長写像の反復も等長。 -/
theorem norm_iterate_of_norm_eq {f : L → L} (hiso : ∀ y, ‖f y‖ = ‖y‖) :
    ∀ (k : ℕ) (y : L), ‖f^[k] y‖ = ‖y‖ := by
  intro k
  induction k with
  | zero => simp
  | succ k ih => intro y; rw [Function.iterate_succ_apply, ih (f y), hiso]

/-- 抽象核: 等長で差を保つ `f` について `‖f^[k]π − π‖ ≤ ‖fπ − π‖`。

★望遠鏡 `f^[k+1]π − π = f(f^[k]π − π) + (fπ − π)` を超距離で `max` に潰すだけ。
★これは「`G_j` が合成で閉じている(部分群である)」ことのノルム版である。 -/
theorem norm_iterate_sub_self_le_of_isometry {f : L → L}
    (hsub : ∀ x y, f (x - y) = f x - f y) (hiso : ∀ y, ‖f y‖ = ‖y‖) (π : L) :
    ∀ k : ℕ, ‖f^[k] π - π‖ ≤ ‖f π - π‖ := by
  intro k
  induction k with
  | zero => simp
  | succ k ih =>
    have h1 : f^[k + 1] π - π = f (f^[k] π - π) + (f π - π) := by
      rw [hsub, Function.iterate_succ_apply']
      ring
    rw [h1]
    refine (IsUltrametricDist.norm_add_le_max _ _).trans ?_
    rw [hiso]
    exact max_le ih le_rfl

/-- ★★抽象核: `f^[p] = id` かつ `gcd(k,p) = 1` なら `‖f^[k]π − π‖ = ‖fπ − π‖`。

★片側は上の望遠鏡、もう片側は `f` 自身が `f^[k]` の反復であること
(`k·k^{φ(p)−1} ≡ 1 mod p`)から出る。
★★これが「位数 `p` の巡回群では跳びが生成元の取り方に依らない」ことの全部で、
★分岐・付値・Galois・多項式の語彙が 1 語も出てこない。 -/
theorem norm_iterate_sub_self_eq_of_coprime {f : L → L}
    (hsub : ∀ x y, f (x - y) = f x - f y) (hiso : ∀ y, ‖f y‖ = ‖y‖)
    {p : ℕ} (hp : 1 < p) (hid : f^[p] = id) (π : L) {k : ℕ} (hk : Nat.Coprime k p) :
    ‖f^[k] π - π‖ = ‖f π - π‖ := by
  have ht : 0 < Nat.totient p := Nat.totient_pos.mpr (by omega)
  have hkt : k ^ Nat.totient p ≡ 1 [MOD p] := Nat.ModEq.pow_totient hk
  have h1 : k ^ Nat.totient p % p = 1 := by
    have h := hkt
    rwa [Nat.ModEq, Nat.mod_eq_of_lt hp] at h
  have hq : k ^ Nat.totient p = p * (k ^ Nat.totient p / p) + 1 := by
    conv_lhs => rw [← Nat.div_add_mod (k ^ Nat.totient p) p]
    rw [h1]
  have hiter : f^[k ^ Nat.totient p] = f := by
    rw [hq, Function.iterate_add, Function.iterate_mul, hid, Function.iterate_id]
    simp
  have hsplit : k ^ Nat.totient p = k * k ^ (Nat.totient p - 1) := by
    rw [← pow_succ']
    congr 1
    omega
  refine le_antisymm (norm_iterate_sub_self_le_of_isometry hsub hiso π k) ?_
  have h2 : ‖(f^[k])^[k ^ (Nat.totient p - 1)] π - π‖ ≤ ‖f^[k] π - π‖ :=
    norm_iterate_sub_self_le_of_isometry (iterate_map_sub_of_map_sub hsub k)
      (norm_iterate_of_norm_eq hiso k) π _
  rwa [← Function.iterate_mul, ← hsplit, hiter] at h2

/-- ★`p` が素数のとき、`0 < k < p` なら `f^[k]` の跳びは `f` の跳びと等しい。

★`norm_natCast_le_pow_of_prod_X_sub_C` の `hbreak`(根 `T` 上の一様な距離)を
`T = {f^[k] π ∣ 0 < k < p}` の形で供給するときに使う。 -/
theorem norm_iterate_sub_self_eq_of_lt {f : L → L}
    (hsub : ∀ x y, f (x - y) = f x - f y) (hiso : ∀ y, ‖f y‖ = ‖y‖)
    {p : ℕ} (hp : p.Prime) (hid : f^[p] = id) (π : L) {i k : ℕ}
    (hk0 : 0 < k) (hkp : k < p) (hjump : ‖f π - π‖ = ‖π‖ ^ (i + 1)) :
    ‖f^[k] π - π‖ = ‖π‖ ^ (i + 1) := by
  have hnd : ¬ p ∣ k := fun h => absurd (Nat.le_of_dvd hk0 h) (by omega)
  have hcop : Nat.Coprime k p := (Nat.Prime.coprime_iff_not_dvd hp).mpr hnd |>.symm
  rw [norm_iterate_sub_self_eq_of_coprime hsub hiso hp.one_lt hid π hcop, hjump]

end Iterate

/-! ## §8 `.src` -/

def norm_natCast_mul_pow_le_norm_aeval_derivative.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_natCast_le_pow_of_prod_X_sub_C.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_natCast_le_pow_of_splits.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def sub_one_mul_le_of_norm_natCast_eq_pow.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def zpow_neg_jump_le_rpow.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_sub_digit_zero_le_rpow_mul.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §9 使っている公理の一覧 -/

#print axioms eval_derivative_of_eq_X_sub_C_mul
#print axioms norm_multiset_prod_map_eq_pow
#print axioms norm_natCast_mul_pow_le_norm_aeval_derivative
#print axioms aeval_derivative_eq_prod_of_split
#print axioms norm_natCast_le_pow_of_prod_X_sub_C
#print axioms exists_multiset_of_splits
#print axioms norm_natCast_le_pow_of_splits
#print axioms sub_one_mul_le_of_norm_natCast_eq_pow
#print axioms sub_one_mul_lt_of_not_dvd
#print axioms rpow_inv_natCast_le_of_le_pow
#print axioms zpow_neg_jump_le_rpow
#print axioms norm_sub_digit_zero_le_rpow_mul
#print axioms norm_iterate_sub_self_le_of_isometry
#print axioms norm_iterate_sub_self_eq_of_coprime
#print axioms norm_iterate_sub_self_eq_of_lt

end ABC3.Found.PGC
