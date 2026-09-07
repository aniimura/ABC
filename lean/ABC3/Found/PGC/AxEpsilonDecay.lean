import ABC3.Found.PGC.AxTowerDecay

/-!
# [pGC] 「深い段では `ε` が減る」—— `AxWildDescent` に残るただ 1 点の詰め

`Found/PGC/AxTowerDecay.lean` は Ax–Sen–Tate を
**`AxWildDescent K c`(wild 深さ 1 段の降下)ただ 1 点**に落とし、
その減衰の形として

  `hdecay : ∀ k, c k ≤ p ^ ((1/(p−1)) · (1/p)^k)`   (★以下「旧仮説」)

を仮定した `axLemma_of_wildDescent_geometric` を用意した。
本ファイルは (a) **旧仮説が `k = 1` で古典的な最良定数を下回る**ことを測り、
(b) **正しい指数**(`(1/p)^{k−1}`)で同じ結論 `AxLemma K (axConstant p)` を出し、
(c) **`ε` が減る仕組みそのもの**を、分岐・付値・Galois・`p` 進の語彙を
1 語も使わない抽象核として切り出す。

## 在庫調査(自分で測った。★MCP は 0 回。★コマンドを残す)

```
grep -n "sum_range_sub" .cache/mathlib-index.txt
  → 1 件しか出ない(`Real.sum_range_sub_log_div_le`)。★**本命が出ない**
  ★`#check @Finset.sum_range_sub` を leanfile.mjs に投げると**在る**:
     `∑ i ∈ range n, (f (i+1) - f i) = f n - f 0`
grep -n "norm_sum_le_of_forall_le" .cache/mathlib-index.txt        → ★0 件
  ★しかし `#check @IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg` は**在る**:
     `0 ≤ C → (∀ i ∈ s, ‖f i‖ ≤ C) → ‖∑ i ∈ s, f i‖ ≤ C`
  ★★これで「索引に無い ⇒ 不在」の反例が **6 例目**である。
  ☆★しかも `Found/PGC/AxLemma.lean` の docstring は
     「超距離での Multiset 和の一様上界は mathlib に無い(AxLemma.lean が自前で持つ)」
     と書いている ——★**Finset 版は在る。**(Multiset 版は確かに無い。)
grep -n "sum_Ico_eq_sum_range" .cache/mathlib-index.txt            → ★0 件
  ★`#check @Finset.sum_Ico_eq_sum_range` は**在る**(`∑_{Ico m n} f = ∑_{range (n−m)} f (m+·)`)
grep -n "Ico_succ_right\|Icc.*=.*Ico" .cache/mathlib-index.txt
  → `Finset.Ico_succ_right_eq_Icc` は在るが `Order.succ` 版。ℕ では `ext; omega` の方が速い
grep -n "norm_p\b" .cache/mathlib-index.txt
  → ★Padic.norm_p (NumberTheory/Padics/PadicNumbers.lean:852) `‖(p : ℚ_p)‖ = (p:ℝ)⁻¹`
  ★`padicNormE.norm_p` は**無い**(`Invalid field \`norm_p\`: ... `AbsoluteValue.norm_p``)
grep -c "\b<名前>\b" .cache/decl-index.txt (本ファイルの全宣言名について) → すべて 0(衝突なし)
```

## ★★★測って分かったこと 1 —— 旧仮説は `k = 1` で**古典的な最良定数を下回る**

☆★★**持ち場の指示は「`c k ≤ p^{(1/(p−1))(1/p)^k}` を示せ」だった。★これは示せない。**

`AxWildDescent K c` は wild 深さ `1` の `x` について
「深さ `0` の `x′` が `‖x − x′‖ ≤ c 1 · ε` かつ `∀σ ‖σx′ − x′‖ ≤ c 1 · ε` で取れる」
を含む。深さ `0` は tame なので木の
`exists_norm_sub_algebraMap_le_of_natDegree_tame` が `‖x′ − a‖ ≤ c 1 · ε` を返し、
超距離で `d(x, K) ≤ c 1 · ε` になる。★これを定理にしたのが本ファイルの
`axWildDescent_depth_one` である。

ところが **`d(x,K) = p^{1/(p−1)} · ε` を満たす深さ 1 の `x` が実在する**(手計算):

  `K = ℚ_p(ζ_p)`(`e_K = p−1`)、`x = p^{1/p}`、`L = K(x)`。
  `T^p − p` は Eisenstein なので `[L:K] = p`、`L/K` は巡回(`ζ_p ∈ K`)、
  `σx = ζ_p x`。`e_L = v_L(p) = p(p−1)`、`v_L(x) = p−1`、`v_L(ζ_p − 1) = p` だから
  `v_L(σx − x) = 2p − 1`。`v_L(K^×) = pℤ ∌ p−1` なので `d(x,K) = |π_L|^{p−1}`。
  ⇒ `d(x,K)/Δ(x) = |π_L|^{−p} = p^{p/(p(p−1))} = p^{1/(p−1)}`。
  また `wildDepth K x = v_p([K(x):K]) = v_p(p) = 1`。

⇒ ★★**`c 1 ≥ p^{1/(p−1)}` が強制される。**
一方、旧仮説は `c 1 ≤ p^{(1/(p−1))·(1/p)} = p^{1/(p(p−1))} < p^{1/(p−1)}` を要求する。
⇒ ★★**`AxWildDescent K c ∧ 旧仮説` は `K = ℚ_p(ζ_p)` について空虚である。**
★`axLemma_of_wildDescent_geometric` は**真だが空虚**な含意だった。
★この段落は「次の波が旧仮説を証明しようとするのを止めるため」に書いてある。
★(★手計算の部分は形式化していない。形式化したのは `axWildDescent_depth_one`
  ——「深さ 1 の反例があれば `c 1` の下界が出る」という**橋**の方である。)

### ★指数がずれている、という診断

旧仮説の指数の総和は `Σ_{k≥1}(1/(p−1))p^{−k} = 1/(p−1)²` であって、
結論に書いてある `axConstant p = p^{p/(p−1)²}` **より真に小さい**。
★「仮説の総和が結論より小さい」こと自体が**仮説が強すぎる**という合図である。
正しい指数は **`(1/p)^{k−1}`** で、そのとき

  `Σ_{k≥1} (1/(p−1))·p^{−(k−1)} = (1/(p−1))·(p/(p−1)) = p/(p−1)²`

と **`axConstant` にちょうど一致する**(`axExponent_eq`)。
★本ファイルの `axDecay` / `axLemma_of_axDecay` がこの正しい形である。
★`axDecay p 1 = p^{1/(p−1)}` は上の古典的な最良定数と**ぴったり一致する**
(`axDecay_one`)。★★**指数が 1 つずれていた、というだけの話だが、
ずれていた側では仮説が満たせない。**

★旧仮説 ⇒ 新仮説(`rpow_geometric_le_axDecay`)なので、
★**本ファイルの `axLemma_of_axDecay` は `axLemma_of_wildDescent_geometric` を含む**
(機械が検査した意味で真に弱い仮説である)。

## ★★★測って分かったこと 2 —— 「`ε` が減る」の中身は**純粋に群論と超距離**である

`σ` を 1 つ固定し、`y := σx − x` と置く。`σ^n x − x` は**望遠鏡和**で

  `σ^n x − x = Σ_{j<n} σ^j y = n·y + Σ_{j<n}(σ^j y − y)`

と書ける。超距離なので

  `‖σ^n x − x‖ ≤ max( ‖n·y‖ , max_{j<n} ‖σ^j y − y‖ )`。

★★**第 1 項は `‖(n : K)‖·‖y‖`** —— `n = p` なら **`p^{−1}‖y‖`**(激減する)。
★★**第 2 項は `σ` が `y` をどれだけ動かすか** —— 深い分岐群の元なら小さい。
⇒ ★★**`σ → σ^p` と上がるたびに変位が縮む。**★これが「深い段では `ε` が減る」の中身で、
★**分岐・付値・Galois・`p` 進の語彙が 1 語も要らない**(`iterate_sub_self_eq_nsmul_add` /
`norm_iterate_sub_self_le`)。★具体層は `smulClosureAddHom` を代入するだけである
(`norm_smul_pow_sub_le` / `norm_smul_pow_prime_sub_le`)。

★★**持ち場の指示が挙げた `σ^p x − x = Tr_{M/M′}(σx − x)` は、
`σ` の位数が `p` のときは両辺 `0` になって使えない**(`σ^p = id`)。
★実際に効くのは上の**望遠鏡和と `‖(p:K)‖ = 1/p`** の方である。
★跡写像も `Tr(𝒪_M) ⊆ 𝔭^c` も**要らなかった**。

## ★★測って分かったこと 3 —— 「`ε` が増えない」なら定数は 1 段ぶんで済む

`AxTowerDecay.exists_mem_of_descent_budget` の予算関数を `F ≡ max 1 C` に取ると、
**`g ≤ 1`(`ε` が増えない)なら総損失は 1 段の損失 `C` そのもの**である
(`exists_mem_of_descent_uniform` / `axLemma_of_wildDescentDecay`)。
★★**警告**: これは `C = p^{1/(p−1)}` で `AxLemma K (p^{1/(p−1)})` を与え、
**Ax の定数 `p^{p/(p−1)²}` より真に良い**(`1/(p−1) < p/(p−1)²`)。
⇒ ★**`θ ≤ 1` は一般には成り立たないと見るべきである。**★掘る前に測ること。
★それでもこの補題を置くのは、**2 つのつまみ(1 段の損失 `C` と `ε` の伸び `θ`)の
どちらを詰めれば閉じるかを機械が言える形にする**ためである。

## 何が言えたか

### §1 抽象核(★分岐・付値・Galois・`p` 進の語彙が **1 語も出てこない**)

| 宣言 | 内容 |
|---|---|
| `iterate_addMonoidHom_sub` | `f^[j](a−b) = f^[j]a − f^[j]b` |
| `iterate_sub_self_eq_sum` | ★**望遠鏡和** `f^[n]x − x = Σ_{j<n} f^[j](fx−x)` |
| `iterate_sub_self_eq_nsmul_add` | `f^[n]x − x = n•(fx−x) + Σ_{j<n}(f^[j](fx−x) − (fx−x))` |
| ★★`norm_iterate_sub_self_le` | ★**`ε` が減る** `‖f^[n]x − x‖ ≤ max a b·‖fx − x‖` |
| ★`norm_iterate_pow_sub_self_le` | ★**塔に沿った幾何減衰** `‖f^[n^k]x − x‖ ≤ θ^k‖fx − x‖` |
| `finsetSum_geometric_shift_le` | `Σ_{k∈[1,n]} A r^{k−1} ≤ A/(1−r)`(★ずらした版) |
| `prod_le_rpow_of_geometric_shift` | `c k ≤ D^{A r^{k−1}}` ⇒ `∏_{[1,n]} c k ≤ D^{A/(1−r)}` |
| ★`exists_mem_of_descent_uniform` | `c ≡ C`・`g ≡ θ ≤ 1` ⇒ 総損失は `C` |

### §2 具体層

| 宣言 | 内容 |
|---|---|
| `smulClosureAddHom` / `iterate_smulClosureAddHom` | `Γ_K` の作用を加法準同型として反復 |
| `norm_nsmul_closure` | `‖n•y‖ = ‖(n:K)‖·‖y‖` |
| ★★`norm_smul_pow_sub_le` | `‖σ^n x − x‖ ≤ max ‖(n:K)‖ b · ‖σx − x‖` |
| ★★`norm_smul_pow_prime_sub_le` | `n = p` の場合。**第 1 のつまみが `1/p`** |
| `axLemma_of_wildDescent_Icc` | 積の条件を `Finset.Icc 1 n` だけに弱めた版 |
| ★`axDecay` / `one_le_axDecay` / `axDecay_one` | 正しい減衰列 `p^{(1/(p−1))(1/p)^{k−1}}` |
| ★`rpow_geometric_le_axDecay` | ★**旧仮説 ⇒ 新仮説**(機械が検査した包含) |
| ★★`axLemma_of_axDecay` | ★**`AxWildDescent K axDecay` ⇒ `AxLemma K (axConstant p)`** |
| ★★`axSenTate_of_axDecay` | ★そこから `AxSenTate K` |
| ★★`axWildDescent_depth_one` | ★**測定器**: 深さ 1 で `d(x,K) ≤ c 1 · ε` |
| `AxWildDescentDecay` / `axLemma_of_wildDescentDecay` | もう 1 つのつまみ(`ε` が増えない版) |

## ★残った穴(正直に書く)

★★**`AxLemma K C` の項は作れていない。**
残っているのは `AxWildDescent K axDecay`、すなわち
**「wild 深さ `k` の 1 段の損失が `p^{(1/(p−1))p^{1−k}}` で済む」**である。
### ★★★2026-09-08 追記 —— 上の見立ては**2 点とも古い**

☆★**(1) `(p−1)i ≤ e_L` は木に入った。**
`Found/PGC/RamificationJumpBound.lean`(561 行、`sorry` 0)が閉じた。
★`differentIdeal` を通らず、`f'(π)` を「根の側」と「係数の側」の 2 通りに測るだけで出る。
★等号が `ℚ₂(√2)/ℚ₂` で実現するので**これ以上強い形は無い**。
出口は `norm_sub_digit_zero_le_rpow_mul : ‖x − a₀‖ ≤ p^(1/(p−1)) · ‖σ x − x‖`。

☆★★**(2) 「`k ≥ 2` は中間体 2 層の壁」は誤りだった。**
★そもそも「深さ `k` の中間体で 1 段下がる」という**枠組み自体が偽**である
(`WildDepthDescent.lean::not_forall_exists_relIndex_padicValNat_eq`、
反例 `G = A₄`・`H` = 1 点固定群・`p = 2` を形式化した)。
★★正しい形は `P` = p-Sylow、`P ≤ Q`・`[Q:P] = p` を取り
`y := (1/p)Σ_{c∈Q/P} c•x` と**平均する**ことで、`x'` が `K(x)` に入る必要はない。
★これで `WildDepthFieldDescent.lean::axWildDescent_prime` が
**`AxWildDescent K (fun _ => (p:ℝ))` を無条件に**与えた(`k ≥ 1` 全部、`k = 1` と分ける必要も無い)。
★配管は `K(x)` を `IntermediateField` として作らず `MulAction.stabilizer` の指数で次数を測ることで
**#59 に一度も触らずに**抜けた(`lean-idioms.md` #296)。

★★**したがって残るのは定数だけ**である: `axWildDescent_prime` の `c k = p`(積 `p^n`、非有界)を
`axDecay p k = p^{(1/(p−1))p^{1−k}}`(積 `p^{p/(p−1)²}`、有界)に絞る 1 点。
★`WildDepthDescent.lean` §6 の `exists_natDegree_minpoly_descent_div` を差し替える形になっている。

☆★また、本ファイルが「望遠鏡和が `ε` 減衰の中身」と診断した点も**半分外れていた** ——
1 段の評価には効くが、深さの降下で実際に効いたのは**剰余類 `Q/P` 上の平均**である。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. ★原典 (Ax 1970) は定数の勘定を地の文で畳む。本ファイルは
   `AxTowerDecay` の予算関数の枠組みに乗せて**指数を明示**した。原典に無い定式化である。
2. ★`AxTowerDecay.axLemma_of_wildDescent_geometric` の仮説を**修正**した
   (`(1/p)^k → (1/p)^{k−1}`)。★既存の宣言は**触っていない**(消すと他が壊れる)。
   本ファイルは `rpow_geometric_le_axDecay` で「旧 ⇒ 新」を証明し、包含関係を明示した。
3. §1 の抽象核は `M` に超距離な半ノルム加法群しか要求せず、`f` はただの加法準同型でよい
   (★**全単射性も等長性も要らない**)。原典より弱い設定である。
4. `.src` は `AxLemma.lean` / `SenLemma.lean` / `AxTowerDecay.lean` と**同じ項目**
   (pGC 物理 p.6 Corollary 3.1)を指す。原典が独立に立てた項目ではない。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued

/-! ## §1 抽象核

★以下の 8 本には**分岐・付値・Galois・`p` 進の語彙が 1 つも出てこない**。 -/

section EpsCore

variable {M : Type*} [SeminormedAddCommGroup M]

/-- 加法準同型の反復は差を保つ。 -/
theorem iterate_addMonoidHom_sub (f : M →+ M) (j : ℕ) (a b : M) :
    f^[j] (a - b) = f^[j] a - f^[j] b := by
  induction j generalizing a b with
  | zero => simp
  | succ k ih => simp [Function.iterate_succ_apply, ih, map_sub]

def iterate_addMonoidHom_sub.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★**望遠鏡和** —— `f^[n] x − x = Σ_{j<n} f^[j] (f x − x)`。

★「`σ` の `n` 乗による変位は、1 乗による変位の `σ`-軌道の和」。
★群も可逆性も要らない(加法準同型 1 本)。 -/
theorem iterate_sub_self_eq_sum (f : M →+ M) (x : M) (n : ℕ) :
    f^[n] x - x = ∑ j ∈ Finset.range n, f^[j] (f x - x) := by
  have h : ∀ j : ℕ, f^[j] (f x - x) = f^[j + 1] x - f^[j] x := by
    intro j
    rw [iterate_addMonoidHom_sub, Function.iterate_succ_apply]
  simp only [h]
  rw [Finset.sum_range_sub (fun j => f^[j] x) n]
  simp

def iterate_sub_self_eq_sum.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★**望遠鏡和の「跡」型の書き換え** ——
`f^[n] x − x = n•(f x − x) + Σ_{j<n}(f^[j](f x − x) − (f x − x))`。

★右辺第 1 項が「`n` 倍」、第 2 項が「`σ` が `y = f x − x` をどれだけ動かすか」。
★★この 2 項に分けるのが「深い段では `ε` が減る」の全部である。 -/
theorem iterate_sub_self_eq_nsmul_add (f : M →+ M) (x : M) (n : ℕ) :
    f^[n] x - x =
      n • (f x - x) + ∑ j ∈ Finset.range n, (f^[j] (f x - x) - (f x - x)) := by
  have h : ∑ j ∈ Finset.range n, (f^[j] (f x - x) - (f x - x))
      = (∑ j ∈ Finset.range n, f^[j] (f x - x)) - n • (f x - x) := by
    rw [Finset.sum_sub_distrib, Finset.sum_const, Finset.card_range]
  rw [h, iterate_sub_self_eq_sum]
  abel

def iterate_sub_self_eq_nsmul_add.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

variable [IsUltrametricDist M]

/-- ★★★**「`ε` が減る」の抽象核**(本ファイルの主結果 1)。

`y := f x − x` について

* `‖n • y‖ ≤ a‖y‖`(★`M` が体で `n` が可逆なら `a = ‖(n : M)‖`。`n = p` なら `1/p`)、
* `∀ j < n, ‖f^[j] y − y‖ ≤ b‖y‖`(★`f` が `y` をどれだけ動かすか)

ならば **`‖f^[n] x − x‖ ≤ max a b · ‖f x − x‖`**。

★★`max a b < 1` なら **`f` を `f^[n]` に取り替えると変位が真に縮む**。
★超距離性が効くのは 1 箇所(和が `max` に潰れる)だけである。
★★分岐も付値も Galois も出てこない。 -/
theorem norm_iterate_sub_self_le (f : M →+ M) (x : M) (n : ℕ) {a b : ℝ}
    (ha : ‖n • (f x - x)‖ ≤ a * ‖f x - x‖)
    (hb : ∀ j < n, ‖f^[j] (f x - x) - (f x - x)‖ ≤ b * ‖f x - x‖) :
    ‖f^[n] x - x‖ ≤ max a b * ‖f x - x‖ := by
  set y := f x - x with hy
  have hy0 : (0:ℝ) ≤ ‖y‖ := norm_nonneg _
  have hay : 0 ≤ a * ‖y‖ := le_trans (norm_nonneg _) ha
  have hmax : a * ‖y‖ ≤ max a b * ‖y‖ :=
    mul_le_mul_of_nonneg_right (le_max_left a b) hy0
  have hC : 0 ≤ max a b * ‖y‖ := le_trans hay hmax
  have hsum : ‖∑ j ∈ Finset.range n, (f^[j] y - y)‖ ≤ max a b * ‖y‖ := by
    refine IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg hC ?_
    intro j hj
    exact (hb j (Finset.mem_range.mp hj)).trans
      (mul_le_mul_of_nonneg_right (le_max_right a b) hy0)
  rw [iterate_sub_self_eq_nsmul_add]
  exact (IsUltrametricDist.norm_add_le_max _ _).trans (max_le (ha.trans hmax) hsum)

def norm_iterate_sub_self_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

omit [IsUltrametricDist M] in
/-- ★**塔に沿った幾何減衰** —— 1 段で `θ` 倍に縮むなら `k` 段で `θ^k` 倍。

★`norm_iterate_sub_self_le` を `f^[n^m]` に当てたものが仮説 `H` である。
★★これが「`Σ` が収束する」ことの中身であり、`AxTowerDecay.prod_le_rpow_of_geometric`
と対になる(あちらは損失の積、こちらは `ε` の積)。 -/
theorem norm_iterate_pow_sub_self_le (f : M →+ M) (n : ℕ) {θ : ℝ} (hθ : 0 ≤ θ) (x : M)
    (H : ∀ (m : ℕ) (z : M), ‖f^[n ^ (m + 1)] z - z‖ ≤ θ * ‖f^[n ^ m] z - z‖) (k : ℕ) :
    ‖f^[n ^ k] x - x‖ ≤ θ ^ k * ‖f x - x‖ := by
  induction k with
  | zero => simp
  | succ m ih =>
      refine (H m x).trans ?_
      calc θ * ‖f^[n ^ m] x - x‖ ≤ θ * (θ ^ m * ‖f x - x‖) :=
            mul_le_mul_of_nonneg_left ih hθ
        _ = θ ^ (m + 1) * ‖f x - x‖ := by ring

def norm_iterate_pow_sub_self_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end EpsCore

/-! ### §1.2 ずらした幾何級数(添字が `1` から始まる版)

★`AxTowerDecay.finsetSum_geometric_le` は `Σ_{k∈s} A r^k` を扱う。
降下の添字は `Finset.Icc 1 n` を走るので、**指数を 1 つずらした版**が要る。
★このずれが「測って分かったこと 1」の正体である。 -/

/-- ★`Σ_{k ∈ [1,n]} A r^{k−1} ≤ A/(1−r)`。★`k = 1` の項が `A` そのものになる。 -/
theorem finsetSum_geometric_shift_le {A r : ℝ} (hA : 0 ≤ A) (hr0 : 0 ≤ r) (hr1 : r < 1)
    (n : ℕ) : ∑ k ∈ Finset.Icc 1 n, A * r ^ (k - 1) ≤ A / (1 - r) := by
  have h : Finset.Icc 1 n = Finset.Ico 1 (n + 1) := by
    ext m; simp only [Finset.mem_Icc, Finset.mem_Ico]; omega
  rw [h, Finset.sum_Ico_eq_sum_range]
  simp only [Nat.add_sub_cancel_left, Nat.add_sub_cancel]
  exact finsetSum_geometric_le hA hr0 hr1 (Finset.range n)

def finsetSum_geometric_shift_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★`c k ≤ D^{A r^{k−1}}` なら `∏_{k∈[1,n]} c k ≤ D^{A/(1−r)}`。 -/
theorem prod_le_rpow_of_geometric_shift {D A r : ℝ} (hD : 1 ≤ D) (hA : 0 ≤ A)
    (hr0 : 0 ≤ r) (hr1 : r < 1) {c : ℕ → ℝ} (hc0 : ∀ k, 0 ≤ c k)
    (hc : ∀ k, c k ≤ D ^ (A * r ^ (k - 1))) (n : ℕ) :
    ∏ k ∈ Finset.Icc 1 n, c k ≤ D ^ (A / (1 - r)) := by
  have hD0 : (0:ℝ) < D := lt_of_lt_of_le zero_lt_one hD
  calc ∏ k ∈ Finset.Icc 1 n, c k ≤ ∏ k ∈ Finset.Icc 1 n, D ^ (A * r ^ (k - 1)) :=
        Finset.prod_le_prod (fun i _ => hc0 i) (fun i _ => hc i)
    _ = D ^ (∑ k ∈ Finset.Icc 1 n, A * r ^ (k - 1)) := (Real.rpow_sum_of_pos hD0 _ _).symm
    _ ≤ D ^ (A / (1 - r)) :=
        Real.rpow_le_rpow_of_exponent_le hD (finsetSum_geometric_shift_le hA hr0 hr1 n)

def prod_le_rpow_of_geometric_shift.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★**もう 1 つのつまみ** —— 1 段の損失が一様に `C` で、しかも **`ε` が増えない**
(`g ≡ θ ≤ 1`)なら、総損失は **1 段ぶんの `C`** で済む。

★`AxTowerDecay.exists_mem_of_descent_budget` の予算関数を `F ≡ C` に取っただけ。
★★`AxTowerDecay.prod_unbounded_of_one_lt`(`θ > 1` だと積が発散する)と対になる。 -/
theorem exists_mem_of_descent_uniform {M : Type*} [SeminormedAddCommGroup M]
    [IsUltrametricDist M] {S : Set M} (deg : M → ℕ) (P : ℝ → M → Prop) {C θ : ℝ}
    (hC : 1 ≤ C) (hθ0 : 0 ≤ θ) (hθ1 : θ ≤ 1)
    (hbase : ∀ (ε : ℝ) (x : M), 0 ≤ ε → P ε x → deg x = 0 → ∃ y ∈ S, ‖x - y‖ ≤ ε)
    (hstep : ∀ (ε : ℝ) (x : M), 0 ≤ ε → P ε x → deg x ≠ 0 →
      ∃ x', deg x' < deg x ∧ ‖x - x'‖ ≤ C * ε ∧ P (θ * ε) x') :
    ∀ (x : M) (ε : ℝ), 0 ≤ ε → P ε x → ∃ y ∈ S, ‖x - y‖ ≤ C * ε := by
  intro x ε hε hP
  exact exists_mem_of_descent_budget deg P (fun _ => C) (fun _ => θ) (fun _ => C)
    (fun _ => hC) (fun _ => hθ0) (fun _ _ => le_refl C)
    (fun _ _ _ => by nlinarith) hbase hstep x ε hε hP

def exists_mem_of_descent_uniform.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §2 具体層 -/

variable {p : ℕ} [Fact p.Prime]

/-- `Γ_K` の `K̄` への作用を加法準同型として見たもの(`K̄` 版)。

★`AxLemma.smulAddHom` は `ℂ_K` 版である。本ファイルは `K̄` の上で使う。 -/
noncomputable def smulClosureAddHom (K : PAdicLocalField p) (σ : K.absGal) :
    K.closure →+ K.closure :=
  DistribSMul.toAddMonoidHom K.closure σ

def smulClosureAddHom.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

@[simp] theorem smulClosureAddHom_apply (K : PAdicLocalField p) (σ : K.absGal)
    (x : K.closure) : smulClosureAddHom K σ x = σ • x := rfl

/-- 反復は冪と一致する。 -/
theorem iterate_smulClosureAddHom (K : PAdicLocalField p) (σ : K.absGal) (n : ℕ)
    (x : K.closure) : (smulClosureAddHom K σ)^[n] x = (σ ^ n) • x := by
  induction n generalizing x with
  | zero => simp
  | succ m ih =>
      rw [Function.iterate_succ_apply, ih, smulClosureAddHom_apply, ← mul_smul, ← pow_succ]

def iterate_smulClosureAddHom.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- `‖n • y‖ = ‖(n : K)‖ · ‖y‖`(`K̄` は乗法的ノルムを持つ体)。 -/
theorem norm_nsmul_closure (K : PAdicLocalField p) (n : ℕ) (y : K.closure) :
    ‖n • y‖ = ‖((n : ℕ) : K.carrier)‖ * ‖y‖ := by
  rw [nsmul_eq_mul, norm_mul, norm_natCast_closure]

def norm_nsmul_closure.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**「深い段では `ε` が減る」の具体層**(本ファイルの主結果 2)。

`σ ∈ Γ_K` と `x ∈ K̄` について、`y := σx − x` が `σ` の冪でどれだけ動くかを
`b` で押さえると

  `‖σ^n x − x‖ ≤ max ‖(n : K)‖ b · ‖σ x − x‖`。

★★**第 1 のつまみ `‖(n:K)‖` は無条件に効く** —— `p ∣ n` なら `1` より真に小さい。
★第 2 のつまみ `b` は「`σ` が下付き分岐群のどこに居るか」で決まる。
★★分岐理論はここ(=`b` を与えるところ)にしか要らない。 -/
theorem norm_smul_pow_sub_le (K : PAdicLocalField p) (σ : K.absGal) (x : K.closure) (n : ℕ)
    {b : ℝ} (hb : ∀ j < n, ‖(σ ^ j) • (σ • x - x) - (σ • x - x)‖ ≤ b * ‖σ • x - x‖) :
    ‖(σ ^ n) • x - x‖ ≤ max ‖((n : ℕ) : K.carrier)‖ b * ‖σ • x - x‖ := by
  have key := norm_iterate_sub_self_le (smulClosureAddHom K σ) x n
    (a := ‖((n : ℕ) : K.carrier)‖) (b := b)
    (le_of_eq (norm_nsmul_closure K n _))
    (by simpa [iterate_smulClosureAddHom] using hb)
  simpa [iterate_smulClosureAddHom] using key

def norm_smul_pow_sub_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★`n = p` の場合。★**第 1 のつまみはちょうど `1/p`** である。

  `‖σ^p x − x‖ ≤ max (1/p) b · ‖σ x − x‖`。

★★`b < 1` なら **`σ` を `σ^p` に取り替えるだけで変位が真に縮む**。
★これが「wild に深い段では `ε` が減る」の 1 行の姿である。 -/
theorem norm_smul_pow_prime_sub_le (K : PAdicLocalField p) (σ : K.absGal) (x : K.closure)
    {b : ℝ} (hb : ∀ j < p, ‖(σ ^ j) • (σ • x - x) - (σ • x - x)‖ ≤ b * ‖σ • x - x‖) :
    ‖(σ ^ p) • x - x‖ ≤ max ((p : ℝ))⁻¹ b * ‖σ • x - x‖ := by
  have h := norm_smul_pow_sub_le K σ x p hb
  rwa [norm_natCast_carrier, Padic.norm_p] at h

def norm_smul_pow_prime_sub_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ### §2.1 正しい指数への配線 -/

/-- ★`AxTowerDecay.axLemma_of_wildDescent` の、積の条件を **`Finset.Icc 1 n` だけ**に
弱めた版。★ずらした減衰では `c 0` に条件が付かないので、この形が要る。 -/
theorem axLemma_of_wildDescent_Icc (K : PAdicLocalField p) {c : ℕ → ℝ} {C : ℝ}
    (hc : ∀ k, 1 ≤ c k) (hC : ∀ n : ℕ, ∏ k ∈ Finset.Icc 1 n, c k ≤ C)
    (h : AxWildDescent K c) : AxLemma K C := by
  classical
  have key : ∀ (z : K.closure) (e : ℝ), 0 ≤ e → (∀ σ : K.absGal, ‖σ • z - z‖ ≤ e) →
      ∃ y ∈ Set.range (algebraMap K.carrier K.closure),
        ‖z - y‖ ≤ (∏ k ∈ Finset.Icc 1 (wildDepth K z), c k) * e := by
    refine exists_mem_of_descent_prod (wildDepth K)
      (fun e z => ∀ σ : K.absGal, ‖σ • z - z‖ ≤ e) c hc ?_ ?_
    · intro e z he hz hz0
      obtain ⟨y, hy⟩ := exists_norm_sub_algebraMap_le_of_natDegree_tame K he
        ((wildDepth_eq_zero_iff K z).mp hz0) hz
      exact ⟨algebraMap K.carrier K.closure y, ⟨y, rfl⟩, hy⟩
    · intro e z he hz hz0
      obtain ⟨z', hlt, hd, hinv⟩ := h e he z hz
        (by
          by_contra hdvd
          exact hz0 ((wildDepth_eq_zero_iff K z).mpr hdvd))
      exact ⟨z', hlt, hd, hinv⟩
  intro ε hε x hx
  obtain ⟨y, ⟨a, ha⟩, hy⟩ := key x ε hε hx
  refine ⟨a, ?_⟩
  rw [ha]
  exact hy.trans (mul_le_mul_of_nonneg_right (hC _) hε)

def axLemma_of_wildDescent_Icc.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★**正しい減衰列** `axDecay p k = p^{(1/(p−1))·(1/p)^{k−1}}`。

★`k = 1` で `p^{1/(p−1)}`(=巡回 `p` 次 1 段の**古典的な最良定数**)、
★指数の総和は `Σ_{k≥1}(1/(p−1))p^{1−k} = p/(p−1)²` で `axConstant` に**ちょうど**一致する。
★`AxTowerDecay.axLemma_of_wildDescent_geometric` の `(1/p)^k` は
**指数が 1 つずれており**、`k = 1` で最良定数を下回るので満たせない
(モジュール docstring「測って分かったこと 1」)。 -/
noncomputable def axDecay (p : ℕ) (k : ℕ) : ℝ :=
  (p : ℝ) ^ ((1 / ((p : ℝ) - 1)) * (1 / (p : ℝ)) ^ (k - 1))

def axDecay.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- `1 ≤ axDecay p k`。 -/
theorem one_le_axDecay (p : ℕ) [Fact p.Prime] (k : ℕ) : 1 ≤ axDecay p k := by
  have h1 : (1:ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  refine Real.one_le_rpow (le_of_lt h1) ?_
  have : (0:ℝ) < 1 / ((p:ℝ) - 1) := by
    apply div_pos one_pos; linarith
  positivity

def one_le_axDecay.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★`axDecay p 1 = p^{1/(p−1)}` —— **巡回 `p` 次 1 段の古典的な最良定数**。

★`CyclicJumpNorm.norm_sub_digit_zero_eq_zpow_mul` の `‖π_L‖^{−i}` と
`(p−1)i ≤ e_L` から出るはずの値であり、
★`K = ℚ_p(ζ_p)`・`x = p^{1/p}` で**等号が実現する**(モジュール docstring)。 -/
theorem axDecay_one (p : ℕ) [Fact p.Prime] :
    axDecay p 1 = (p : ℝ) ^ (1 / ((p : ℝ) - 1)) := by
  simp [axDecay]

def axDecay_one.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★**旧仮説 ⇒ 新仮説**(機械が検査した包含)。

`AxTowerDecay.axLemma_of_wildDescent_geometric` の `hdecay` は
本ファイルの `axDecay` による条件を**含意する**ので、
★`axLemma_of_axDecay` は旧定理を**真に一般化している**。 -/
theorem rpow_geometric_le_axDecay (p : ℕ) [Fact p.Prime] (k : ℕ) :
    (p : ℝ) ^ ((1 / ((p : ℝ) - 1)) * (1 / (p : ℝ)) ^ k) ≤ axDecay p k := by
  have h1 : (1:ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  have hp0 : (0:ℝ) < (p : ℝ) := by linarith
  have hr0 : (0:ℝ) ≤ 1 / (p:ℝ) := by positivity
  have hr1 : 1 / (p:ℝ) ≤ 1 := by rw [div_le_one hp0]; linarith
  have hA : (0:ℝ) ≤ 1 / ((p:ℝ) - 1) := by
    apply le_of_lt; apply div_pos one_pos; linarith
  refine Real.rpow_le_rpow_of_exponent_le (le_of_lt h1) ?_
  refine mul_le_mul_of_nonneg_left ?_ hA
  exact pow_le_pow_of_le_one hr0 hr1 (Nat.sub_le k 1)

def rpow_geometric_le_axDecay.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**正しい減衰から Ax の定数が出る**(本ファイルの主結果 3)。

wild 深さ `k` の 1 段の損失が `axDecay p k = p^{(1/(p−1))p^{1−k}}` で押さえられるなら
**`AxLemma K (axConstant p)`**、すなわち `d(x,K) ≤ |p|^{−p/(p−1)²}·ε`。

★★`AxTowerDecay.axLemma_of_wildDescent_geometric` との違いは**指数のずれ 1 つ**だけだが、
★あちらの仮説は `k = 1` で古典的な最良定数を下回るため**満たせない**。
★`rpow_geometric_le_axDecay` により本定理はあちらを含む。 -/
theorem axLemma_of_axDecay (K : PAdicLocalField p) {c : ℕ → ℝ}
    (hc : ∀ k, 1 ≤ c k) (hdecay : ∀ k, c k ≤ axDecay p k)
    (h : AxWildDescent K c) : AxLemma K (axConstant p) := by
  have h1 : (1:ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  have hp0 : (0:ℝ) < (p : ℝ) := by linarith
  refine axLemma_of_wildDescent_Icc K hc (fun n => ?_) h
  have hbound := prod_le_rpow_of_geometric_shift (D := (p : ℝ)) (A := 1 / ((p : ℝ) - 1))
    (r := 1 / (p : ℝ)) (le_of_lt h1) (by positivity) (by positivity)
    (by rw [div_lt_one hp0]; linarith)
    (fun k => le_trans zero_le_one (hc k)) hdecay n
  rwa [axExponent_eq h1] at hbound

def axLemma_of_axDecay.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**そこから Ax–Sen–Tate**。

★これで pGC §3 の全部が **`AxWildDescent K (axDecay p)` ただ 1 点**に落ちた。 -/
theorem axSenTate_of_axDecay (K : PAdicLocalField p) {c : ℕ → ℝ}
    (hc : ∀ k, 1 ≤ c k) (hdecay : ∀ k, c k ≤ axDecay p k)
    (h : AxWildDescent K c) : AxSenTate K :=
  axSenTate_of_axLemma (le_trans zero_le_one (one_le_axConstant p))
    (axLemma_of_axDecay K hc hdecay h)

def axSenTate_of_axDecay.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ### §2.2 測定器 —— `c 1` の下界 -/

/-- ★★**測定器**(本ファイルの主結果 4)。

`AxWildDescent K c` は、**wild 深さ 1** の `x` について
**`d(x, K) ≤ c 1 · ε` を無条件に含意する**(降下 1 回 + tame な基底段)。

★★したがって「`d(x,K) = D·ε` を満たす深さ 1 の `x`」が 1 つでもあれば `D ≤ c 1`。
★`K = ℚ_p(ζ_p)`・`x = p^{1/p}` が `D = p^{1/(p−1)}` を実現する(モジュール docstring)。
⇒ ★★`AxTowerDecay` の旧仮説 `c 1 ≤ p^{1/(p(p−1))}` は**満たせない**。
★この 1 本が「旧仮説は空虚」という診断の**形式化された部分**である。 -/
theorem axWildDescent_depth_one (K : PAdicLocalField p) {c : ℕ → ℝ}
    (h : AxWildDescent K c) (hc1 : 0 ≤ c 1) {ε : ℝ} (hε : 0 ≤ ε) {x : K.closure}
    (hx : ∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) (hdepth : wildDepth K x = 1) :
    ∃ a : K.carrier, ‖x - algebraMap K.carrier K.closure a‖ ≤ c 1 * ε := by
  have hdvd : p ∣ (minpoly K.carrier x).natDegree := by
    by_contra hd
    have : wildDepth K x = 0 := (wildDepth_eq_zero_iff K x).mpr hd
    omega
  obtain ⟨x', hlt, hd, hinv⟩ := h ε hε x hx hdvd
  rw [hdepth] at hlt hd hinv
  have h0 : wildDepth K x' = 0 := by omega
  have hCe : 0 ≤ c 1 * ε := mul_nonneg hc1 hε
  obtain ⟨a, ha⟩ := exists_norm_sub_algebraMap_le_of_natDegree_tame K hCe
    ((wildDepth_eq_zero_iff K x').mp h0) hinv
  refine ⟨a, ?_⟩
  have hmax := IsUltrametricDist.norm_add_le_max (x - x') (x' - algebraMap K.carrier K.closure a)
  have heq : x - x' + (x' - algebraMap K.carrier K.closure a)
      = x - algebraMap K.carrier K.closure a := by abel
  rw [heq] at hmax
  exact hmax.trans (max_le hd ha)

def axWildDescent_depth_one.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ### §2.3 もう 1 つのつまみ —— 「`ε` が増えない」形の降下 -/

/-- ★**`ε` の伸びを明示した wild 降下**。

`AxTowerDecay.AxWildDescent K c` は「距離の損失」と「`ε` の伸び」を**同じ `c`** で
押さえていた。本定義は**分離**して、損失を一様な `C`、伸びを `θ` とする。

★★**警告(掘る前に読むこと)**: `θ ≤ 1` なら `AxLemma K C` が出る
(`axLemma_of_wildDescentDecay`)。★`C = p^{1/(p−1)}` を入れると
**Ax の定数 `p^{p/(p−1)²}` より真に良い**ので、
★`θ ≤ 1` は一般には成り立たないと見るべきである。★**非空虚性は測っていない。** -/
def AxWildDescentDecay (K : PAdicLocalField p) (C θ : ℝ) : Prop :=
  ∀ ε : ℝ, 0 ≤ ε → ∀ x : K.closure, (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) →
    p ∣ (minpoly K.carrier x).natDegree →
      ∃ x' : K.closure, wildDepth K x' < wildDepth K x ∧
        ‖x - x'‖ ≤ C * ε ∧
        ∀ σ : K.absGal, ‖σ • x' - x'‖ ≤ θ * ε

def AxWildDescentDecay.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★**`ε` が増えないなら定数は 1 段ぶんで済む**。

`AxWildDescentDecay K C θ` と `θ ≤ 1` から **`AxLemma K C`**。
★段数によらない(積が立たない)のが要点で、`exists_mem_of_descent_uniform` そのもの。 -/
theorem axLemma_of_wildDescentDecay (K : PAdicLocalField p) {C θ : ℝ}
    (hC : 1 ≤ C) (hθ0 : 0 ≤ θ) (hθ1 : θ ≤ 1) (h : AxWildDescentDecay K C θ) :
    AxLemma K C := by
  classical
  have key : ∀ (z : K.closure) (e : ℝ), 0 ≤ e → (∀ σ : K.absGal, ‖σ • z - z‖ ≤ e) →
      ∃ y ∈ Set.range (algebraMap K.carrier K.closure), ‖z - y‖ ≤ C * e := by
    refine exists_mem_of_descent_uniform (wildDepth K)
      (fun e z => ∀ σ : K.absGal, ‖σ • z - z‖ ≤ e) hC hθ0 hθ1 ?_ ?_
    · intro e z he hz hz0
      obtain ⟨y, hy⟩ := exists_norm_sub_algebraMap_le_of_natDegree_tame K he
        ((wildDepth_eq_zero_iff K z).mp hz0) hz
      exact ⟨algebraMap K.carrier K.closure y, ⟨y, rfl⟩, hy⟩
    · intro e z he hz hz0
      obtain ⟨z', hlt, hd, hinv⟩ := h e he z hz
        (by
          by_contra hdvd
          exact hz0 ((wildDepth_eq_zero_iff K z).mpr hdvd))
      exact ⟨z', hlt, hd, hinv⟩
  intro ε hε x hx
  obtain ⟨y, ⟨a, ha⟩, hy⟩ := key x ε hε hx
  exact ⟨a, by rw [ha]; exact hy⟩

def axLemma_of_wildDescentDecay.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★そこから `AxSenTate K`。 -/
theorem axSenTate_of_wildDescentDecay (K : PAdicLocalField p) {C θ : ℝ}
    (hC : 1 ≤ C) (hθ0 : 0 ≤ θ) (hθ1 : θ ≤ 1) (h : AxWildDescentDecay K C θ) :
    AxSenTate K :=
  axSenTate_of_axLemma (le_trans zero_le_one hC) (axLemma_of_wildDescentDecay K hC hθ0 hθ1 h)

def axSenTate_of_wildDescentDecay.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §3 使っている公理の一覧 -/

#print axioms iterate_addMonoidHom_sub
#print axioms iterate_sub_self_eq_sum
#print axioms iterate_sub_self_eq_nsmul_add
#print axioms norm_iterate_sub_self_le
#print axioms norm_iterate_pow_sub_self_le
#print axioms finsetSum_geometric_shift_le
#print axioms prod_le_rpow_of_geometric_shift
#print axioms exists_mem_of_descent_uniform
#print axioms smulClosureAddHom
#print axioms smulClosureAddHom_apply
#print axioms iterate_smulClosureAddHom
#print axioms norm_nsmul_closure
#print axioms norm_smul_pow_sub_le
#print axioms norm_smul_pow_prime_sub_le
#print axioms axLemma_of_wildDescent_Icc
#print axioms axDecay
#print axioms one_le_axDecay
#print axioms axDecay_one
#print axioms rpow_geometric_le_axDecay
#print axioms axLemma_of_axDecay
#print axioms axSenTate_of_axDecay
#print axioms axWildDescent_depth_one
#print axioms AxWildDescentDecay
#print axioms axLemma_of_wildDescentDecay
#print axioms axSenTate_of_wildDescentDecay

end ABC3.Found.PGC
