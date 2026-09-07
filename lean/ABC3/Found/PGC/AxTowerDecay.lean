import ABC3.Found.PGC.SenLemma
import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Algebra.Order.BigOperators.GroupWithZero.Finset
import Mathlib.Algebra.BigOperators.Intervals

/-!
# [pGC] 塔に沿った降下の「予算」—— Ax の定数 `p/(p−1)²` はどこから来るか

`Found/PGC/SenLemma.lean` は Ax–Sen–Tate を **`AxDescentStep K C`(1 段の降下、
しかも「もとの `ε` を保つ」)ただ 1 点**に落とし、
`Found/PGC/CyclicJumpNorm.lean` は **1 段の最良定数が `‖π_L‖^{−i}`(等式)**であることを出した。
本ファイルは残っていた**「段を重ねたときの勘定」**を詰める。

## 在庫調査(自分で測った。★MCP は 0 回。コマンドを残す)

```
grep -n "geom_sum_lt\|geom_sum_eq\|hasSum_geometric_of_lt_one" .cache/mathlib-index.txt
  → ★hasSum_geometric_of_lt_one (Analysis/SpecificLimits/Basic.lean:316)
  → geom_sum_lt (Algebra/Order/Field/GeomSum.lean:44) —— ただし `CanonicallyOrderedAdd` 要求で ℝ に効かない
grep -n "le_hasSum\|sum_le_tsum\b\|HasSum.sum_le" .cache/mathlib-index.txt
  → 5 件出るが `sum_le_hasSum` は**その中に無い**(出るのは `Lp_add_le_hasSum` /
    `exists_le_hasSum_of_le` / `ENNReal.sum_le_tsum` など、別物ばかり)
  ★★しかし `#check @sum_le_hasSum` を `leanfile.mjs` に投げると **在る**
     (`(s : Finset ι) → (∀ i ∉ s, 0 ≤ f i) → HasSum f a L → ∑ i ∈ s, f i ≤ a`)。
  ★「索引に無い ⇒ 不在」は言えない、の **5 例目**である(直前の波の `to_additive` に続く)。
grep -n "Real.rpow_sum_of_pos\|rpow_le_rpow_of_exponent_le" .cache/mathlib-index.txt
  → ★Real.rpow_sum_of_pos (Pow/Real.lean:242) `a^(Σ f) = Π a^f` ——本ファイルの要
  → ★Real.rpow_le_rpow_of_exponent_le (Pow/Real.lean:613)
grep -n "prod_le_prod_of_subset" .cache/mathlib-index.txt
  → ★Finset.prod_le_prod_of_subset_of_one_le (GroupWithZero/Finset.lean:70)
     ★`'` 付き(Group 版)は `MulLeftMono N` を要求するので **ℝ には効かない**(乗法が単調でない)。
     同様に `Finset.one_le_prod''` も駄目で、`Finset.one_le_prod`(GroupWithZero 版)を使う。
grep -n "Padic.valuation_natCast\|Padic.norm_eq_zpow_neg_valuation" .cache/mathlib-index.txt
  → ★両方在る。`‖(n : ℚ_p)‖ = p^{−v_p(n)}` はこの 2 本で出る。
grep -n "pow_unbounded_of_one_lt" .cache/mathlib-index.txt
  → ★在る(Algebra/Order/Archimedean/Basic.lean:158)
grep -c "\b<名前>\b" .cache/decl-index.txt (本ファイルの全宣言名について)  → すべて 0(衝突なし)
```

## ★★★測って分かったこと 1 —— 「`Σ_j i_j/e_j` が収束する」は **偽**である

☆★★**持ち場の指示は「`Σ_j i_j/e_j` が収束し、その和が `p/(p−1)²` で押さえられること」
を埋めよ、だった。★これは一般の塔については成り立たない。**
反例は**円分塔**である(手計算。本ファイルでは形式化していない。`p` は奇素数)。

`F_n = ℚ_p(μ_{p^n})` とおく。`e_n := v_{F_n}(p) = p^{n−1}(p−1)`、
`v_p(𝔡_{F_n/ℚ_p}) = n − 1/(p−1)` だから `v_p(𝔡_{F_n/F_{n−1}}) = 1`、すなわち
`v_{F_n}(𝔡_{F_n/F_{n−1}}) = e_n`。次数 `p` の全分岐拡大では `d = (p−1)(i+1)` なので

  `i_n = p^{n−1} − 1`,  `i_n/e_n = (p^{n−1}−1)/(p^{n−1}(p−1)) → 1/(p−1) > 0`。

⇒ ★★`Σ_{n≥2} i_n/e_n ≈ (N−1)/(p−1) → ∞`。**発散する。**
★この「各項が正の定数以上なら部分和は非有界」という骨は
`sum_unbounded_of_pos_le` として**形式化した**(残りは上の手計算だけである)。

★★**しかもこれは「sharp な `(p−1)i ≤ e_L` を証明すれば直る」類の話ではない。**
上の計算は `(p−1)i_n = (p−1)(p^{n−1}−1) < e_n` で **sharp な評価をほぼ等号で満たしており**、
1 段あたりの損失は `|p|^{−1/(p−1)}` が**最良**である(`CyclicJumpNorm.lean` が
等式で出した `‖π_L‖^{−i}` がそれ)。
☆★直前の波は「価値があるのは sharp な `(p−1)i ≤ e_L` で、それには different の評価が要る。
未着手」と書き残したが、★**その未着手項目を埋めても Ax の定数は出ない。**
★次の波がそこに費やすのを止めるために書いておく。

## ★★★測って分かったこと 2 —— 「一様な 1 段定数」は原理的に足りない(形式化した)

`prod_unbounded_of_one_lt` を証明した:

  `1 < D` かつ `∀k, D ≤ c k` ⇒ `∀ C, ∃ n, C < ∏_{k ∈ Icc 1 n} c k`。

★これは直前の波が地の文で書いた「超距離は**和**を `max` に潰すが**積**は潰さない」を
**定理にしたもの**である。★段ごとに `ε` を `C` 倍(`C > 1`)する勘定では、
段数が伸びると必ず破綻する。

## ★★★では Ax の定数はどこから来るのか(本ファイルの答え)

★★**`ε` が段ごとに減る**からである。1 段の降下は 2 つの量を持つ:

* `c d` —— **距離の損失**(`‖x − x′‖ ≤ c d · ε`)、
* `g d` —— **`ε` の伸び**(`∀σ ‖σx′ − x′‖ ≤ g d · ε`)。

`SenLemma.AxDescentStep` は `g ≡ 1`(`ε` を保つ)を要求して定数を `C` に固定した。
★しかし wild な段では `g d ≥ c d > 1` になりうる。逆に、**深い段では `g d < 1`**
になりうる —— 激しく分岐した拡大では**跡が小さい**ので
`σ^p x − x = Tr(σx − x)` が `σx − x` より真に小さくなるからである。
★★本ファイルの抽象核 `exists_mem_of_descent_budget` は
**`c` と `g` を分離**し、両者を束ねる「予算関数」`F` で総損失を押さえる:

  `hFc : d ≠ 0 → c d ≤ F d`,  `hFg : d′ < d → g d · F d′ ≤ F d`  ⇒  `d(x,S) ≤ F(deg x)·ε`。

★`g ≡ c` を入れると `F d = ∏_{k ≤ d} c k`(積の勘定)に退化する
(`exists_mem_of_descent_prod`)。★`g d < 1` ならば `F` は有界に取れる。
★★**Ax の定数 `p/(p−1)²` はこの `F` の値である**:

  `p/(p−1)² = Σ_{k≥0} (1/(p−1))·p^{−k}`  (`axExponent_eq` / `tsum_axExponent`)

なので、`c k ≤ |p|^{−(1/(p−1))p^{−k}}` という**幾何級数的な減衰**があれば
`F ≡ ‖p‖^{−p/(p−1)²} = axConstant p` で足りる(`prod_le_rpow_of_geometric`)。

## 何が言えたか

### §1 抽象核(★分岐・付値・Galois・`p` 進の語彙が **1 語も出てこない**)

| 宣言 | 内容 |
|---|---|
| `finsetSum_geometric_le` | 幾何級数の**任意の**有限部分和は `A/(1−r)` 以下 |
| `axExponent_eq` | `(1/(P−1))/(1−1/P) = P/(P−1)²`(Ax の指数の閉じた形) |
| `tsum_axExponent` | `Σ'_k (1/(P−1))(1/P)^k = P/(P−1)²` |
| ★`prod_le_rpow_of_geometric` | `c k ≤ D^{A r^k}` ⇒ 任意の有限積が `D^{A/(1−r)}` 以下 |
| ★`prod_unbounded_of_one_lt` | ★**一様定数 `> 1` では積が非有界**(否定的結果) |
| ★`sum_unbounded_of_pos_le` | ★**各項が `δ > 0` 以上なら部分和は非有界**(同、加法版) |
| ★★`exists_mem_of_descent_budget` | ★**予算つき降下**(`c` と `g` を分離した降下の主定理) |
| `exists_mem_of_descent_prod` | その系(`g ≡ c`、`F = ∏`) |
| `exists_mem_of_descent_of_ultrametric` | その系(`c ≡ g ≡ 1`。`SenLemma.exists_mem_of_descent` の再導出) |

★`exists_mem_of_descent_budget` は `‖·‖` と `ℕ` しか使わない。強帰納法 1 本、13 行。

### §2 具体層 —— wild 深さによる降下

* `wildDepth K x := v_p([K(x):K])` —— ★**tame ⟺ `wildDepth = 0`**。
* `AxWildDescent K c` —— 「`p ∣ [K(x):K]` なら **wild 深さが真に下がる** `x′` が
  `‖x − x′‖ ≤ c(深さ)·ε` かつ `∀σ ‖σx′ − x′‖ ≤ c(深さ)·ε` で取れる」。
* ★★`axLemma_of_wildDescent` —— **`AxWildDescent K c` + 「積が `C` 以下」⇒ `AxLemma K C`**。
  ★★**tame な基底段は仮説ではなく、木の `exists_norm_sub_algebraMap_le_of_natDegree_tame`
  で埋めてある**(`AxDescentStep` と違い、仮説は wild な段だけ)。
* ★★`axLemma_of_wildDescent_geometric` —— 減衰 `c k ≤ p^{(1/(p−1))(1/p)^k}` から
  **`AxLemma K (axConstant p)`**、`axConstant p = p^{p/(p−1)²} = |p|^{−p/(p−1)²}`
  ——★**原典 (Ax 1970) の定数そのもの**。
* ★★`axSenTate_of_wildDescent_geometric` —— そこから `AxSenTate K`。
* ★`axWildDescent_pow` —— ★★**非空虚性(無条件に証明した)**:
  `AxWildDescent K (fun k => p^k)` は**真**である(重心を 1 回取るだけ)。
  ★もちろん `∏_{k≤n} p^k = p^{n(n+1)/2}` は非有界なので `AxLemma` は出ない。
  ★★**残っている数学は「`p^k` を `p^{(1/(p−1))p^{−k}}` に落とす」ただ 1 点**である。

## ★★`AxDescentStep`(直前の波)との関係

★2 つは**比較不能**である。`AxDescentStep K C` は `g ≡ 1` を要求する代わりに
`c ≡ C` の一様定数でよい。`AxWildDescent K c` は `g = c` を許す代わりに
`c` が**積の意味で総和可能**であることを要求する。
★`AxDescentStep K C → AxWildDescent K (fun _ => C)` は**この形では出せない**
(次数 `[K(x):K]` が下がっても wild 深さ `v_p([K(x):K])` が下がるとは限らない。
★「反例がある」ではなく「導出できない」である ——**測っていない**)。
★逆向きも同様に出せない(`AxWildDescent` は `ε` の保存を言わない)。
★どちらを埋めるかは次の波の選択である。
☆★**本ファイルの見立て**: 円分塔の計算(上)から、
★**`g ≡ 1` は wild な段では成立しない**(`ε` は実際に増える)ので
`AxDescentStep` の方が偽に近い。★`AxWildDescent` を勧める。

## ★★★後の波による訂正 —— 本ファイルの `hdecay` は**指数が 1 つずれており、満たせない**

☆★★`axLemma_of_wildDescent_geometric` の仮説
`c k ≤ p^{(1/(p−1))(1/p)^k}` は、★**`k = 1` で古典的な最良定数 `p^{1/(p−1)}` を下回る**。
`AxWildDescent K c` は深さ 1 の `x` について `d(x,K) ≤ c 1 · ε` を含意し
(`AxEpsilonDecay.axWildDescent_depth_one`)、`K = ℚ_p(ζ_p)`・`x = p^{1/p}` が
`d(x,K) = p^{1/(p−1)}·ε` を実現するので **`c 1 ≥ p^{1/(p−1)}`** が強制される。
⇒ ★★`axLemma_of_wildDescent_geometric` は**真だが空虚**な含意である。
★合図は「仮説の指数の総和 `Σ_{k≥1}(1/(p−1))p^{−k} = 1/(p−1)²` が
**結論の `p/(p−1)²` より小さい**」ことだった(=仮説が強すぎる)。

★★**正しい形は `Found/PGC/AxEpsilonDecay.lean` の `axDecay p k = p^{(1/(p−1))(1/p)^{k−1}}`**
で、指数の総和がちょうど `p/(p−1)²` になる(`axLemma_of_axDecay` / `axSenTate_of_axDecay`)。
★`rpow_geometric_le_axDecay` が「旧仮説 ⇒ 新仮説」を機械で示している。
★★**次の波は `AxWildDescent K (axDecay p)` を掘ること。本ファイルの `hdecay` を掘らないこと。**

## ★残った穴(正直に書く)

★★**`AxLemma K C` の項は作れていない。** 作れたのは

1. 予算つき降下の抽象核(段数によらない勘定の一般形)、
2. tame な基底段を**討ち取った**形の `axLemma_of_wildDescent`、
3. Ax の定数 `p/(p−1)²` が幾何級数 `Σ (1/(p−1))p^{−k}` の値であることの証明、
4. その減衰から `AxLemma K (axConstant p)` が出ること、
5. 非空虚性 `axWildDescent_pow`(定数 `p^k` なら**無条件に真**)、

の 5 つである。残っているのは ★**「wild 深さ `k` の 1 段の損失が
`|p|^{−(1/(p−1))p^{−k}}` で済む」**——これが Sen の補題の本体であり、
★**数学が足りない**(配管ではない)。
★必要な新ノード:「深い段では `ε` が減る」——
`σ^p x − x = Tr_{M/M′}(σx − x)` と、激しく分岐した拡大での跡の評価
(`Tr(O_M) ⊆ 𝔭^{c}`, `c > 0`)。★木の `HerbrandComposition.lean` /
`UpperRamificationGroup.lean` はこの `c` を上付き番号で与える入力になりうる。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. ★原典 (Ax 1970 / Sen 1972) は「次数についての帰納法」を地の文で回し、
   定数の勘定を明示しない。本ファイルは**予算関数 `F`** を導入して勘定を明示した。
   ★これは原典に無い定式化である(`AxDescentStep` と同じ種類の逸脱)。
2. 降下の複雑さを `[K(x):K]` ではなく **`v_p([K(x):K])`(wild 深さ)**で測った。
   ★これは tame な段を基底段に押し込むための読み替えで、
   ★その結果 **tame 側が仮説から消えた**(木の補題で埋まった)。
3. `.src` は `AxLemma.lean` / `SenLemma.lean` と**同じ項目**
   (pGC 物理 p.6 Corollary 3.1)を指す。原典が独立に立てた項目ではない。
4. §1 の抽象核は `M` に群構造しか要求せず、`S` は**ただの集合**でよい。
   原典の「`K` に近づける」という主張より弱い設定である。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued

/-! ## §1 抽象核

★★以下の 7 本には**分岐・付値・Galois・`p` 進の語彙が 1 つも出てこない**。 -/

section TowerCore

/-- ★**抽象核** —— 幾何級数の「**任意の**有限部分和」は `A/(1−r)` を超えない。

★`Finset.range n` ではなく**任意の** `s : Finset ℕ` について言えるのが要点で、
降下で実際に使われる添字の集合が何であっても効く。 -/
theorem finsetSum_geometric_le {A r : ℝ} (hA : 0 ≤ A) (hr0 : 0 ≤ r) (hr1 : r < 1)
    (s : Finset ℕ) : ∑ k ∈ s, A * r ^ k ≤ A / (1 - r) := by
  have hsum : HasSum (fun k : ℕ => A * r ^ k) (A * (1 - r)⁻¹) :=
    (hasSum_geometric_of_lt_one hr0 hr1).mul_left A
  have h := sum_le_hasSum s (fun i _ => mul_nonneg hA (pow_nonneg hr0 i)) hsum
  simpa [div_eq_mul_inv] using h

def finsetSum_geometric_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★**Ax の指数の閉じた形** —— `(1/(P−1))/(1 − 1/P) = P/(P−1)²`。

★`P = p` のとき右辺が原典 (Ax 1970) の指数 `p/(p−1)²` である。
★左辺は「初項 `1/(p−1)`、公比 `1/p` の幾何級数の和」であり、
**1 段あたりの損失 `1/(p−1)` が深さについて `p` 分の 1 ずつ減る**ことを意味する。 -/
theorem axExponent_eq {P : ℝ} (hP : 1 < P) :
    (1 / (P - 1)) / (1 - 1 / P) = P / (P - 1) ^ 2 := by
  have h1 : P ≠ 0 := by linarith
  have h2 : P - 1 ≠ 0 := by
    intro h; rw [sub_eq_zero] at h; linarith [h]
  field_simp

def axExponent_eq.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★`axExponent_eq` の無限和版。 -/
theorem tsum_axExponent {P : ℝ} (hP : 1 < P) :
    ∑' k : ℕ, (1 / (P - 1)) * (1 / P) ^ k = P / (P - 1) ^ 2 := by
  have hP0 : (0:ℝ) < P := by linarith
  have h0 : (0:ℝ) ≤ 1 / P := by positivity
  have h1 : 1 / P < 1 := by rw [div_lt_one hP0]; linarith
  rw [tsum_mul_left, tsum_geometric_of_lt_one h0 h1]
  have h2 : P - 1 ≠ 0 := by intro h; rw [sub_eq_zero] at h; linarith [h]
  field_simp

def tsum_axExponent.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★**幾何級数的に減る係数の積は一様に有界**。

`c k ≤ D^{A r^k}`(`D ≥ 1`, `0 ≤ r < 1`)なら、**どんな有限集合上の積も** `D^{A/(1−r)}` 以下。
★これが「塔に沿った跳びの減衰から一様定数が出る」ことの本体である。 -/
theorem prod_le_rpow_of_geometric {D A r : ℝ} (hD : 1 ≤ D) (hA : 0 ≤ A) (hr0 : 0 ≤ r) (hr1 : r < 1)
    {c : ℕ → ℝ} (hc0 : ∀ k, 0 ≤ c k) (hc : ∀ k, c k ≤ D ^ (A * r ^ k)) (s : Finset ℕ) :
    ∏ k ∈ s, c k ≤ D ^ (A / (1 - r)) := by
  have hD0 : (0:ℝ) < D := lt_of_lt_of_le zero_lt_one hD
  calc ∏ k ∈ s, c k ≤ ∏ k ∈ s, D ^ (A * r ^ k) :=
        Finset.prod_le_prod (fun i _ => hc0 i) (fun i _ => hc i)
    _ = D ^ (∑ k ∈ s, A * r ^ k) := (Real.rpow_sum_of_pos hD0 _ _).symm
    _ ≤ D ^ (A / (1 - r)) :=
        Real.rpow_le_rpow_of_exponent_le hD (finsetSum_geometric_le hA hr0 hr1 s)

def prod_le_rpow_of_geometric.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★**否定的結果** —— 係数が一様に `1` より大きいと、積は**非有界**になる。

★直前の波が地の文で書いた「超距離は**和**を `max` に潰すが**積**は潰さない」を
定理にしたもの。★段ごとに `ε` を `D`(`> 1`)倍する勘定は、段数が伸びると必ず破綻する。
★したがって Ax の一様定数は「1 段の一様な評価」だけからは**絶対に出ない**。 -/
theorem prod_unbounded_of_one_lt {D : ℝ} (hD : 1 < D) {c : ℕ → ℝ} (hc : ∀ k, D ≤ c k) (C : ℝ) :
    ∃ n : ℕ, C < ∏ k ∈ Finset.Icc 1 n, c k := by
  obtain ⟨n, hn⟩ := pow_unbounded_of_one_lt C hD
  refine ⟨n, lt_of_lt_of_le hn ?_⟩
  have hcard : (Finset.Icc 1 n).card = n := by simp
  calc D ^ n = ∏ _k ∈ Finset.Icc 1 n, D := by rw [Finset.prod_const, hcard]
    _ ≤ ∏ k ∈ Finset.Icc 1 n, c k :=
        Finset.prod_le_prod (fun _ _ => le_of_lt (lt_trans zero_lt_one hD))
          (fun i _ => hc i)

def prod_unbounded_of_one_lt.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★**否定的結果(加法版)** —— 各項が正の定数 `δ` 以上なら、部分和は**非有界**。

★★これが「`Σ_j i_j/e_j` は収束しない」の形式的な中身である。
円分塔 `F_n = ℚ_p(μ_{p^n})` では `i_n/e_n = (p^{n−1}−1)/(p^{n−1}(p−1))` で、
`n ≥ 2` なら `δ = 1/(2(p−1))` 以上(手計算。モジュール docstring 参照)。
★したがって**持ち場の目標「`Σ_j i_j/e_j` が `p/(p−1)²` で押さえられる」は
一般の塔については偽である**。 -/
theorem sum_unbounded_of_pos_le {δ : ℝ} (hδ : 0 < δ) {t : ℕ → ℝ} (ht : ∀ k, δ ≤ t k) (C : ℝ) :
    ∃ n : ℕ, C < ∑ k ∈ Finset.range n, t k := by
  obtain ⟨n, hn⟩ := exists_nat_gt (C / δ)
  refine ⟨n, ?_⟩
  have h1 : C < n * δ := by
    rw [div_lt_iff₀ hδ] at hn
    exact hn
  refine lt_of_lt_of_le h1 ?_
  calc (n : ℝ) * δ = ∑ _k ∈ Finset.range n, δ := by
        rw [Finset.sum_const, Finset.card_range, nsmul_eq_mul]
    _ ≤ ∑ k ∈ Finset.range n, t k := Finset.sum_le_sum (fun i _ => ht i)

def sum_unbounded_of_pos_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**抽象核 —— 予算つき降下**(本ファイルの主結果 1)。

`deg : M → ℕ` を「複雑さ」とし、1 段の降下が **2 つの量**を持つとする:

* `c d` —— **距離の損失**(`‖x − x′‖ ≤ c d · ε`)、
* `g d` —— **`ε` の伸び**(`x′` が満たす性質は `P (g d · ε)`)。

「予算関数」`F` が

* `hF1 : 1 ≤ F d`、
* `hFc : d ≠ 0 → c d ≤ F d`(1 段目が予算に収まる)、
* `hFg : d′ < d → g d · F d′ ≤ F d`(予算が合成できる)

を満たせば、`P ε` を満たすすべての `x` が `S` から **`F (deg x) · ε`** 以内にある。

★★**要点**: `g d < 1`(深い段では `ε` が**減る**)なら `F` を有界に取れる。
`g ≡ 1` に固定した `SenLemma.AxDescentStep` も、`g ≡ c` の積の勘定も、
どちらもこの核の特別な場合である。
★超距離性が効くのは 1 箇所(`max` で誤差が吸収され、段数分の増幅が起きない)。
★`M` にも `S` にも構造は要らない(`S` は閉である必要すらない)。 -/
theorem exists_mem_of_descent_budget {M : Type*} [SeminormedAddCommGroup M]
    [IsUltrametricDist M] {S : Set M} (deg : M → ℕ) (P : ℝ → M → Prop) (c g F : ℕ → ℝ)
    (hF1 : ∀ d, 1 ≤ F d) (hg0 : ∀ d, 0 ≤ g d)
    (hFc : ∀ d, d ≠ 0 → c d ≤ F d)
    (hFg : ∀ d d', d' < d → g d * F d' ≤ F d)
    (hbase : ∀ (ε : ℝ) (x : M), 0 ≤ ε → P ε x → deg x = 0 → ∃ y ∈ S, ‖x - y‖ ≤ ε)
    (hstep : ∀ (ε : ℝ) (x : M), 0 ≤ ε → P ε x → deg x ≠ 0 →
      ∃ x', deg x' < deg x ∧ ‖x - x'‖ ≤ c (deg x) * ε ∧ P (g (deg x) * ε) x') :
    ∀ (x : M) (ε : ℝ), 0 ≤ ε → P ε x → ∃ y ∈ S, ‖x - y‖ ≤ F (deg x) * ε := by
  intro x
  generalize hn : deg x = n
  induction n using Nat.strong_induction_on generalizing x with
  | _ n ih =>
    intro ε hε hP
    rcases Nat.eq_zero_or_pos n with h0 | hpos
    · obtain ⟨y, hyS, hy⟩ := hbase ε x hε hP (hn.trans h0)
      exact ⟨y, hyS, hy.trans (by nlinarith [hF1 n])⟩
    · obtain ⟨x', hlt, hd, hP'⟩ := hstep ε x hε hP (by omega)
      rw [hn] at hd hP' hlt
      obtain ⟨y, hyS, hy⟩ :=
        ih (deg x') hlt x' rfl (g n * ε) (mul_nonneg (hg0 n) hε) hP'
      refine ⟨y, hyS, ?_⟩
      have h1 : ‖x - x'‖ ≤ F n * ε :=
        hd.trans (mul_le_mul_of_nonneg_right (hFc n (by omega)) hε)
      have h2 : ‖x' - y‖ ≤ F n * ε := by
        refine hy.trans ?_
        have : F (deg x') * (g n * ε) = (g n * F (deg x')) * ε := by ring
        rw [this]
        exact mul_le_mul_of_nonneg_right (hFg n (deg x') hlt) hε
      have hmax := IsUltrametricDist.norm_add_le_max (x - x') (x' - y)
      have heq : x - x' + (x' - y) = x - y := by abel
      rw [heq] at hmax
      exact hmax.trans (max_le h1 h2)

def exists_mem_of_descent_budget.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★**予算つき降下の系 1** —— `g ≡ c`(`ε` の伸びが距離の損失と同じ)なら
予算は**積** `∏_{k ∈ [1, deg x]} c k` に取れる。 -/
theorem exists_mem_of_descent_prod {M : Type*} [SeminormedAddCommGroup M]
    [IsUltrametricDist M] {S : Set M} (deg : M → ℕ) (P : ℝ → M → Prop) (c : ℕ → ℝ)
    (hc : ∀ k, 1 ≤ c k)
    (hbase : ∀ (ε : ℝ) (x : M), 0 ≤ ε → P ε x → deg x = 0 → ∃ y ∈ S, ‖x - y‖ ≤ ε)
    (hstep : ∀ (ε : ℝ) (x : M), 0 ≤ ε → P ε x → deg x ≠ 0 →
      ∃ x', deg x' < deg x ∧ ‖x - x'‖ ≤ c (deg x) * ε ∧ P (c (deg x) * ε) x') :
    ∀ (x : M) (ε : ℝ), 0 ≤ ε → P ε x →
      ∃ y ∈ S, ‖x - y‖ ≤ (∏ k ∈ Finset.Icc 1 (deg x), c k) * ε := by
  have hc0 : ∀ k, (0:ℝ) ≤ c k := fun k => le_trans zero_le_one (hc k)
  have hprod1 : ∀ j : ℕ, (1:ℝ) ≤ ∏ k ∈ Finset.Icc 1 j, c k :=
    fun j => Finset.one_le_prod (fun i _ => hc i)
  refine exists_mem_of_descent_budget deg P c c (fun j => ∏ k ∈ Finset.Icc 1 j, c k)
    hprod1 hc0 ?_ ?_ hbase hstep
  · intro d hd
    obtain ⟨m, rfl⟩ : ∃ m, d = m + 1 := ⟨d - 1, by omega⟩
    rw [Finset.prod_Icc_succ_top (by omega) c]
    nlinarith [hprod1 m, hc0 (m + 1)]
  · intro d d' hlt
    obtain ⟨m, rfl⟩ : ∃ m, d = m + 1 := ⟨d - 1, by omega⟩
    rw [Finset.prod_Icc_succ_top (by omega) c]
    have hsub : (∏ k ∈ Finset.Icc 1 d', c k) ≤ ∏ k ∈ Finset.Icc 1 m, c k :=
      Finset.prod_le_prod_of_subset_of_one_le
        (Finset.Icc_subset_Icc_right (by omega)) (fun i _ => hc0 i) (fun i _ _ => hc i)
    nlinarith [hc0 (m + 1), hprod1 d']

def exists_mem_of_descent_prod.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★**予算つき降下の系 2 / 整合性の検査** —— `c ≡ g ≡ 1` を入れると
`SenLemma.exists_mem_of_descent`(`ε` を保つ降下)が**再導出**できる。

★`exists_mem_of_descent_budget` が空虚でないことの証拠であり、
同時に「直前の波の核は本核の特別な場合」であることの確認でもある。 -/
theorem exists_mem_of_descent_of_ultrametric {M : Type*} [SeminormedAddCommGroup M]
    [IsUltrametricDist M] {S : Set M} (deg : M → ℕ) (Q : M → Prop) {ε : ℝ} (hε : 0 ≤ ε)
    (hbase : ∀ x, Q x → deg x = 0 → ∃ y ∈ S, ‖x - y‖ ≤ ε)
    (hstep : ∀ x, Q x → deg x ≠ 0 → ∃ x', Q x' ∧ deg x' < deg x ∧ ‖x - x'‖ ≤ ε) :
    ∀ x, Q x → ∃ y ∈ S, ‖x - y‖ ≤ ε := by
  intro x hx
  have h := exists_mem_of_descent_prod (S := S) deg (fun e z => Q z ∧ ε ≤ e) (fun _ => 1)
    (fun _ => le_refl 1) ?_ ?_ x ε hε ⟨hx, le_refl ε⟩
  · simpa using h
  · rintro ε' z _ ⟨hz, hεe⟩ hz0
    obtain ⟨y, hyS, hy⟩ := hbase z hz hz0
    exact ⟨y, hyS, hy.trans hεe⟩
  · rintro ε' z hε' ⟨hz, hεe⟩ hz0
    obtain ⟨z', hz', hlt, hd⟩ := hstep z hz hz0
    exact ⟨z', hlt, by simpa using hd.trans hεe, ⟨hz', by simpa using hεe⟩⟩

def exists_mem_of_descent_of_ultrametric.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end TowerCore

/-! ## §2 具体層 —— wild 深さによる降下 -/

variable {p : ℕ} [Fact p.Prime]

/-- ★**wild 深さ** `v_p([K(x):K])`。

★`wildDepth K x = 0` ⟺ `p ∤ [K(x):K]` ⟺ **tame**。
★これが降下の「複雑さ」であり、`AxDescentStep`(次数そのもの)との違いである。 -/
noncomputable def wildDepth (K : PAdicLocalField p) (x : K.closure) : ℕ :=
  padicValNat p (minpoly K.carrier x).natDegree

def wildDepth.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★`wildDepth K x = 0` ⟺ tame(`p ∤ [K(x):K]`)。 -/
theorem wildDepth_eq_zero_iff (K : PAdicLocalField p) (x : K.closure) :
    wildDepth K x = 0 ↔ ¬ (p ∣ (minpoly K.carrier x).natDegree) := by
  have hint : IsIntegral K.carrier x :=
    (Algebra.IsAlgebraic.isAlgebraic (R := K.carrier) x).isIntegral
  have hdeg : (minpoly K.carrier x).natDegree ≠ 0 := (minpoly.natDegree_pos hint).ne'
  rw [wildDepth, padicValNat.eq_zero_iff]
  constructor
  · rintro (h | h | h)
    · exact absurd h (Fact.out : p.Prime).one_lt.ne'
    · exact absurd h hdeg
    · exact h
  · exact fun h => Or.inr (Or.inr h)

def wildDepth_eq_zero_iff.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★`x ∈ K` なら wild 深さは `0`。 -/
theorem wildDepth_algebraMap (K : PAdicLocalField p) (a : K.carrier) :
    wildDepth K (algebraMap K.carrier K.closure a) = 0 :=
  (wildDepth_eq_zero_iff K _).mpr (natDegree_minpoly_algebraMap_tame K a)

def wildDepth_algebraMap.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★**wild な 1 段の降下**(本ファイルが原典から切り出した形)。

「`p ∣ [K(x):K]`(wild)で `∀σ ‖σx − x‖ ≤ ε` なら、次を満たす `x′` が取れる」:

* `wildDepth K x′ < wildDepth K x`(★**wild 深さ**が真に下がる)、
* `‖x − x′‖ ≤ c (wildDepth K x) · ε`、
* `∀σ ‖σx′ − x′‖ ≤ c (wildDepth K x) · ε`(★`ε` の伸びも同じ `c` で押さえる)。

★★`SenLemma.AxDescentStep` との違いは 2 つ:
1. **複雑さを次数ではなく wild 深さで測る** ⇒ tame な段が仮説から消える
   (木の `exists_norm_sub_algebraMap_le_of_natDegree_tame` で埋まる)。
2. **`ε` が増えることを許す** ⇒ その代償に `c` が総和可能(積が有界)である必要がある。

★★**非空虚である**: `c k = p^k` なら**無条件に真**(`axWildDescent_pow`)。 -/
def AxWildDescent (K : PAdicLocalField p) (c : ℕ → ℝ) : Prop :=
  ∀ ε : ℝ, 0 ≤ ε → ∀ x : K.closure, (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) →
    p ∣ (minpoly K.carrier x).natDegree →
      ∃ x' : K.closure, wildDepth K x' < wildDepth K x ∧
        ‖x - x'‖ ≤ c (wildDepth K x) * ε ∧
        ∀ σ : K.absGal, ‖σ • x' - x'‖ ≤ c (wildDepth K x) * ε

def AxWildDescent.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**wild な降下から Ax の補題が出る**(本ファイルの主結果 2)。

`AxWildDescent K c` と「`c` の**任意の**有限積が `C` 以下」から **`AxLemma K C`**。

★★**tame な基底段は仮説ではない** —— 木の
`exists_norm_sub_algebraMap_le_of_natDegree_tame`(重心。定数 `1`)で埋めてある。
★これが `axLemma_of_descentStep`(`AxDescentStep` から出す版)との実質的な違いで、
**仮説が wild な段だけに絞られている**。 -/
theorem axLemma_of_wildDescent (K : PAdicLocalField p) {c : ℕ → ℝ} {C : ℝ}
    (hc : ∀ k, 1 ≤ c k) (hC : ∀ s : Finset ℕ, ∏ k ∈ s, c k ≤ C)
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
  refine hy.trans (mul_le_mul_of_nonneg_right (hC _) hε)

def axLemma_of_wildDescent.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★**Ax の定数** `|p|^{−p/(p−1)²} = p^{p/(p−1)²}`。

★原典 (Ax 1970) が与える一様定数そのものである。 -/
noncomputable def axConstant (p : ℕ) : ℝ := (p : ℝ) ^ ((p : ℝ) / ((p : ℝ) - 1) ^ 2)

def axConstant.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- `1 ≤ axConstant p`。 -/
theorem one_le_axConstant (p : ℕ) [Fact p.Prime] : 1 ≤ axConstant p := by
  have h1 : (1:ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  refine Real.one_le_rpow (le_of_lt h1) ?_
  positivity

def one_le_axConstant.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**跳びの幾何級数的な減衰から Ax の定数が出る**(本ファイルの主結果 3)。

wild 深さ `k` の 1 段の損失が

  `c k ≤ p^{(1/(p−1))·(1/p)^k}`   (すなわち `|p|^{−(1/(p−1))p^{−k}}`)

で押さえられるなら、**`AxLemma K (axConstant p)`**、
すなわち `d(x, K) ≤ |p|^{−p/(p−1)²} · ε`。

★指数の総和が `Σ_{k≥0}(1/(p−1))p^{−k} = p/(p−1)²`(`axExponent_eq`)であることが
**すべて**である。★★これが持ち場「塔に沿った跳びの減衰」の到達点であり、
残っているのは仮説 `hdecay` を**分岐理論で満たすこと**だけである。 -/
theorem axLemma_of_wildDescent_geometric (K : PAdicLocalField p) {c : ℕ → ℝ}
    (hc : ∀ k, 1 ≤ c k)
    (hdecay : ∀ k, c k ≤ (p : ℝ) ^ ((1 / ((p : ℝ) - 1)) * (1 / (p : ℝ)) ^ k))
    (h : AxWildDescent K c) : AxLemma K (axConstant p) := by
  have h1 : (1:ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  have hp0 : (0:ℝ) < (p : ℝ) := by linarith
  refine axLemma_of_wildDescent K hc (fun s => ?_) h
  have hbound := prod_le_rpow_of_geometric (D := (p : ℝ)) (A := 1 / ((p : ℝ) - 1))
    (r := 1 / (p : ℝ)) (le_of_lt h1) (by positivity) (by positivity)
    (by rw [div_lt_one hp0]; linarith)
    (fun k => le_trans zero_le_one (hc k)) hdecay s
  rwa [axExponent_eq h1] at hbound

def axLemma_of_wildDescent_geometric.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**そこから Ax–Sen–Tate**。

`AxLemma.axSenTate_of_axLemma` と繋いだもの。
★これで pGC §3 の全部が **「wild 深さ `k` の 1 段が `|p|^{−(1/(p−1))p^{−k}}` で済む」
ただ 1 点**に落ちた。 -/
theorem axSenTate_of_wildDescent_geometric (K : PAdicLocalField p) {c : ℕ → ℝ}
    (hc : ∀ k, 1 ≤ c k)
    (hdecay : ∀ k, c k ≤ (p : ℝ) ^ ((1 / ((p : ℝ) - 1)) * (1 / (p : ℝ)) ^ k))
    (h : AxWildDescent K c) : AxSenTate K :=
  axSenTate_of_axLemma (le_trans zero_le_one (one_le_axConstant p))
    (axLemma_of_wildDescent_geometric K hc hdecay h)

def axSenTate_of_wildDescent_geometric.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ### §2.1 非空虚性 —— `c k = p^k` なら無条件に真 -/

/-- `‖([K(x):K] : K)‖ = p^{−wildDepth}`。 -/
theorem norm_natCast_natDegree_eq (K : PAdicLocalField p) (x : K.closure) :
    ‖(((minpoly K.carrier x).natDegree : ℕ) : K.carrier)‖
      = (p : ℝ) ^ (-(wildDepth K x : ℤ)) := by
  have hint : IsIntegral K.carrier x :=
    (Algebra.IsAlgebraic.isAlgebraic (R := K.carrier) x).isIntegral
  have hdeg : (minpoly K.carrier x).natDegree ≠ 0 := (minpoly.natDegree_pos hint).ne'
  have hne : (((minpoly K.carrier x).natDegree : ℕ) : ℚ_[p]) ≠ 0 :=
    Nat.cast_ne_zero.mpr hdeg
  rw [norm_natCast_carrier, Padic.norm_eq_zpow_neg_valuation hne, Padic.valuation_natCast]
  rfl

def norm_natCast_natDegree_eq.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★**非空虚性(無条件)** —— `AxWildDescent K (fun k => p^k)` は**真**である。

証明は「共役の重心を 1 回取る」だけ(木の
`exists_norm_sub_algebraMap_le_div_norm_natDegree`)。得られる `x′` は `K` の元なので
wild 深さは `0` に落ち、`∀σ ‖σx′ − x′‖ = 0` も自動である。

★★もちろん `∏_{k ∈ [1,n]} p^k = p^{n(n+1)/2}` は**非有界**なので、
これだけでは `AxLemma` は出ない(`prod_unbounded_of_one_lt` 参照)。
★★**残っている数学は「`p^k` を `p^{(1/(p−1))p^{−k}}` に落とす」ただ 1 点**である。
★`AxWildDescent` の仮説が空集合について語っているのではないことの証拠。 -/
theorem axWildDescent_pow (K : PAdicLocalField p) :
    AxWildDescent K (fun k => (p : ℝ) ^ k) := by
  intro ε hε x hx hdvd
  have h1 : (1:ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  have hp0 : (0:ℝ) < (p : ℝ) := by linarith
  obtain ⟨a, ha⟩ := exists_norm_sub_algebraMap_le_div_norm_natDegree K hε hx
  refine ⟨algebraMap K.carrier K.closure a, ?_, ?_, ?_⟩
  · rw [wildDepth_algebraMap]
    exact Nat.pos_of_ne_zero (fun h0 => (wildDepth_eq_zero_iff K x).mp h0 hdvd)
  · refine ha.trans (le_of_eq ?_)
    rw [norm_natCast_natDegree_eq, zpow_neg, zpow_natCast, div_eq_mul_inv, inv_inv, mul_comm]
  · intro σ
    rw [smul_closure_def, AlgEquiv.commutes, sub_self, norm_zero]
    positivity

def axWildDescent_pow.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §3 使っている公理の一覧 -/

#print axioms finsetSum_geometric_le
#print axioms axExponent_eq
#print axioms tsum_axExponent
#print axioms prod_le_rpow_of_geometric
#print axioms prod_unbounded_of_one_lt
#print axioms sum_unbounded_of_pos_le
#print axioms exists_mem_of_descent_budget
#print axioms exists_mem_of_descent_prod
#print axioms exists_mem_of_descent_of_ultrametric
#print axioms wildDepth
#print axioms wildDepth_eq_zero_iff
#print axioms wildDepth_algebraMap
#print axioms AxWildDescent
#print axioms axLemma_of_wildDescent
#print axioms axConstant
#print axioms one_le_axConstant
#print axioms axLemma_of_wildDescent_geometric
#print axioms axSenTate_of_wildDescent_geometric
#print axioms norm_natCast_natDegree_eq
#print axioms axWildDescent_pow

end ABC3.Found.PGC
