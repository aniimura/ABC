import ABC3.Found.PGC.AxEpsilonDecay
import ABC3.Found.PGC.WildDepthFieldDescent
import Mathlib.Data.Nat.Choose.Dvd

/-!
# [pGC] 「`x` を含む全分岐巡回 `p` 次の層」は取れない ——「深い段ほど得をする」の出どころ

`Found/PGC/WildDepthFieldDescent.lean` は `AxWildDescent K (fun _ => (p:ℝ))` を無条件に立て、
`Found/PGC/AxEpsilonDecay.lean` は `AxWildDescent K (axDecay p)` さえ来れば
`AxLemma` と `AxSenTate` が出ることを示した。残る 1 点は

  `axDecay p k = p^{(1/(p−1))·p^{1−k}}`

まで 1 段の損失を絞ることである。本ファイルはその 1 点を★測った★。
★★埋まっていない。★以下は「測って分かったこと」であり、見通しではない。

## ★★★測定 1 —— 配られた素朴な形は偽である(★反例を形式化した)

配られた字面はこうだった:

  `wildDepth K x = 1` なら `x` を含む全分岐巡回 `p` 次の層 `K(π)` が在る。

★★これは偽である。理由は分岐でも Galois でもなく★次数★である:

* `x ∈ L`・`[L:K] = p` なら `deg minpoly_K x` は `[L:K] = p` を割る
  (`natDegree_minpoly_dvd_finrank`)。
* ところが `wildDepth K x = 1` は「`p` がちょうど 1 回だけ `deg minpoly_K x` を割る」しか
  言わない。`deg = p·m`(`m > 1`, `p ∤ m`)は `wildDepth = 1` だが `p` を割らない。

⇒ ★素朴な形は `wildDepth K x = 1` より真に強い(`natDegree_minpoly_eq_of_mem_prime_layer`)。

★★★形式化した反例(`exists_wildDepth_one_not_mem_prime_layer`):
`K = ℚ`、`x = ζ₇ = exp(2πi/7)`、`p = 3`。
`deg minpoly_ℚ ζ₇ = φ(7) = 6`、`v₃(6) = 1` なので wild 深さは `1` だが、
`[F:ℚ] = 3` なる中間体 `F ⊆ ℂ` で `ζ₇ ∈ F` なるものは★無い★(`6 ∣ 3` が偽)。
★`sorry` 無しで閉じている。★分岐も付値も出てこない ——
★★素朴な形が壊れるのは次数の段であって、分岐の段ではない。

★同じ壊れ方は `p` 進体でも起きる(`natDegree_minpoly_eq_of_mem_prime_layer` が
`PAdicLocalField` の上でそのまま成り立つ)。★具体的な `p` 進の元
(`ℚ_p` 上 `X^{pq} − p` の根、次数 `pq`)は★構成していない★。
★「反例が在ることは測った、`p` 進での具体元は測っていない」と読むこと。

★★2 つ目の壊れ方(こちらは形式化していない): `deg = p` でも `K(x)/K` は Galois とは限らない。
`K = ℚ_p`(`p` 奇)、`x = p^{1/p}` は `deg = p` だが `ζ_p ∉ K(x)` なので `K(x)/K` は正規でない
(`MinpolyOrbitSplit.lean` 冒頭に手計算がある)。★巡回 `p` 次の層 `L` は
`x ∈ L` なら `K(x) = L` を強制する(次数が `p` を割るから)ので、`K(x)/K` が
正規でない時点で層は取れない。

## ★★★測定 2 —— 正しい形は Sylow 降下が★すでに供給している★(障害ではない)

☆★配られた持ち場は「差し替えの障害は定数ではなく `x` が全分岐巡回 `p` 次の層に入ること」
と書いていた。★★これは外れている。

`WildDepthDescent.exists_pgroup_descent` は `P`(= `Stab(x)` の `p`-Sylow)と
`P ≤ Q`・`[Q:P] = p` を返す。Galois 対応で `M^Q ⊆ M^P` は★次数 `p`★ の拡大であり、
`x ∈ M^H ⊆ M^P` なので

  ★`x` は「`M^Q` の上の次数 `p` の層 `M^P`」に入っている。

`Q` は `p` 群で `[Q:P] = p` だから `P ⊴ Q`、したがって `M^P/M^Q` は★巡回 `p` 次★である。
局所体の `p` 次拡大は `e·f = p` より★不分岐か全分岐のどちらか★なので、
`CyclicJumpNorm` / `RamificationJumpBound` の道はそのまま当たる(全分岐の場合)。
⇒ ★層は取れる。取れないのは「`K` の直上に」という部分だけである。
★底が `K` ではなく `M^Q` になるが、`Δ_{M^Q}(x) ≤ Δ_K(x)` なので損失の評価には効かない。

★★(`P ⊴ Q` は本ファイルでは形式化していない。mathlib に
`Sylow.exists_subgroup_card_pow_succ` は在るが★正規性を返さない★ので、
「`p` 群の指数 `p` の部分群は正規」を自前で書く必要がある。★測ったが書いていない。)

## ★★★測定 3 —— `p^{1−k}` の出どころは「跳びの上界」1 本である

★★これが本ファイルの中心である。抽象核は★群環の恒等式★で、
分岐・付値・Galois・`p` 進の語彙が 1 語も出てこない。

`D := σ − 1` と置くと、hockey-stick 恒等式から

| 恒等式 | 意味 |
|---|---|
| `∑_{j<p} σ^j x − p·x = ∑_{l=1}^{p−1} C(p,l+1)·D^l x` | ★軌道和(=平均の道) |
| `σ^p x − x = ∑_{l=1}^{p} C(p,l)·D^l x` | ★`p` 乗(=`ε` の減衰) |

`1 ≤ l < p` では `p ∣ C(p,l)` なので、超距離で和が `max` に潰れて

  ★`‖∑_{j<p}(σ^j x − x)‖ ≤ max( ‖D^{p−1} x‖ , ‖p‖·‖σx − x‖ )`、
  ★`‖σ^p x − x‖ ≤ max( ‖D^{p} x‖ , ‖p‖·‖σx − x‖ )`。

`D` の収縮率を `θ`(`‖σz − z‖ ≤ θ‖z‖`)とすると

  ★★平均の道の損失 `= max( p·θ^{p−2} , 1 )`(`norm_sub_orbit_average_le_of_contract`)、
  ★★`ε` の縮み `= max( θ^{p−1} , 1/p )`(`norm_smul_pow_prime_sub_self_le_of_contract`)。

★分岐が入るのはここ 1 箇所だけである: `θ = ‖π‖^{i(σ)}` で `i(σ)` は `σ` の跳び。
古典的な上界は 「位数 `p^m` の `σ` について `i(σ) ≤ e/(p^{m−1}(p−1))`」
(`m = 1` が `RamificationJumpBound` の `(p−1)i ≤ e_L`)。これを入れると

  `θ(σ_j) = p^{−1/(p^{k−j}(p−1))}`、  `θ(σ_j)^{p−1} = p^{−p^{−(k−j)}}`。

`σ_1 = τ`(位数 `p^k`)から `σ_k = τ^{p^{k−1}}`(位数 `p`)まで `k−1` 回縮めると
縮みの積は `p^{−∑_{m=1}^{k−1}p^{−m}} = p^{−(1−p^{1−k})/(p−1)}`、
最後の 1 層の損失は `p^{1/(p−1)}` だから、合わせて

  ★★★`p^{1/(p−1)} · p^{−(1−p^{1−k})/(p−1)} = p^{(1/(p−1))·p^{1−k}} = axDecay p k`。

★★ぴったり一致する。★これを機械に検査させたのが
`axDecay_exponent_eq_sub_geomSum` と `axDecay_eq_axDecay_one_mul_gains` である。
⇒ ★`p^{1−k}` は「位数 `p^m` の元の跳びの上界」という 1 本の不等式から来ている。

★★手計算での確認(形式化していない): `K = ℚ_p(ζ_{p²})`、`x = p^{1/p²}`(wild 深さ 2)。
`Δ_K(x) = |ζ_{p²}−1|·|x| = p^{−1/(p(p−1))}|x|`、`d(x,K) = |x|` なので
`d/Δ = p^{1/(p(p−1))} = axDecay p 2`。★降下に使う `σ = τ^p` は
`‖σx − x‖ = p^{−1/(p−1)}|x| = p^{−1/p}·Δ` で、縮み `p^{−1/p} = θ(τ)^{p−1}` が実現している。

## ★★★測定 4 —— 否定的な結果: 平均の道だけでは `k ≥ 2` は出ない

★★`axDecay_one_le_orbit_average_loss`(形式化した)。

跳びの上界 `(p−1)i ≤ e` は `θ = ‖π‖^i ≥ p^{−1/(p−1)}` を★強制する★ので

  `max(p·θ^{p−2}, 1) ≥ p·p^{−(p−2)/(p−1)} = p^{1/(p−1)} = axDecay p 1`。

一方 `axDecay p k < axDecay p 1`(`k ≥ 2`、`axDecay_lt_axDecay_one`)。
⇒ ★★★`exists_natDegree_minpoly_descent_div` の「`p` で割る」箇所を
`(σ−1)^{p−1}` で割り引く、という道は `k ≥ 2` では原理的に `axDecay p k` に届かない。
★持ち場が名指しした「`‖p‖^{−1} → ‖p‖^{−(1/(p−1))p^{1−k}}` に置き換わる」は★起きない★。

★★では何が要るか。上の測定 3 が答えである: 得は★損失の側ではなく `ε` の側★にある。
`Δ(x)` は `Γ_K` 全体の上限だが、降下に使う `σ` は `τ^{p^{k−1}}` という深い元なので
`‖σx − x‖` は `Δ(x)` より `p^{−(1−p^{1−k})/(p−1)}` 倍だけ小さい。
⇒ ★次に要るのは「深さ `k` なら降下生成元 `σ` を `τ^{p^{k−1}}` の形に取れる」(群論)と
「位数 `p^m` の元の跳びの上界」(分岐)の 2 本である。

## ★測定 5 —— 順分岐(prime-to-`p`)の底変換は wild 深さを変えない

★純 `ℕ` の核 `padicValNat_eq_of_dvd_of_dvd_mul` で測った:
`b ∣ a`・`a ∣ b·c`・`p ∤ c` なら `v_p a = v_p b`。
`a = [K(x):K]`、`b = [K'(x):K']`、`c = [K':K]` を入れれば
★`p ∤ [K':K]` のとき `wildDepth_{K'} x = wildDepth_K x` である。
★★2 つの整除(`[K'(x):K'] ∣ [K(x):K]` と `[K(x):K] ∣ [K'(x):K']·[K':K]`)は
体の側の入力であり、本ファイルでは形式化していない(#59 の中間体 2 層に触るため)。
⇒ ★「順分岐に取り替えても深さは変わらない」は算術としては測った。
★代替案(`K` を順分岐拡大に取り替えてから `σ` を選ぶ)は、深さを保つ点では通る。
★しかし測定 2 のとおり 層は最初から取れている ので、この取り替えは要らない。

## 在庫の測定(★自分で測った。★MCP は 1 度も呼んでいない。★引いた索引を残す)

`.cache/mathlib-index.txt` を `IsUltrametricDist.norm_` で引くと
★乗法版 `norm_pow_le` しか出ない。★加法版 `norm_nsmul_le` は索引に無い。
ところが `#check @IsUltrametricDist.norm_nsmul_le` は在る:
`(x : S) (n : ℕ) : ‖n • x‖ ≤ ‖x‖`。
⇒ ★★「索引に無い ⇒ 不在」の反例がまた 1 件。★原因は `to_additive` の生成名を
索引が持たないこと(索引は宣言の字面しか読まない)。

同じ索引を `dvd_choose` で引くと `Nat.Prime.dvd_choose_self`
(`Data/Nat/Choose/Dvd.lean:35`)が出る。★★ABC3 の import 連鎖には入っていないので
`Unknown constant` になる(`lean-idioms.md` #68 そのもの)。
★`import Mathlib.Data.Nat.Choose.Dvd` を足して解決した。

`#check @Sylow.exists_subgroup_card_pow_succ` は
`∃ K, Nat.card K = p^(n+1) ∧ H ≤ K` を返すだけで★正規性を返さない★。
索引を `Normal` で引いて `index` / `coatom` / `maximal` で絞っても空振りする。
⇒ ★「`p` 群の指数 `p` の部分群は正規」は mathlib に見当たらない(測定 2 の穴)。

`Module.End.pow_apply` / `Module.End.natCast_apply` / `Module.End.mul_apply` は全部在る。
★`AddMonoid.End` ではなく `Module.End ℤ M` を使うと `(F^n) x` と `(↑n) x = n • x` が
既にあるので、抽象核が短くなる。
`pow_le_pow_left` は★無い★(`pow_le_pow_left₀` に改名されている)。
`mul_max_of_nonneg` は在る: `0 ≤ a → a * max b c = max (a*b) (a*c)`。
`Complex.isPrimitiveRoot_exp` / `cyclotomic_eq_minpoly_rat` / `natDegree_cyclotomic` も全部在り、
★測定 1 の反例はこの 3 本で 10 行で閉じた。
本ファイルの全 33 宣言名を `.cache/decl-index.txt` で引いて、衝突は 0 件だった。

★★配管の罠を 1 つ踏んだ: `(hp : p.Prime)` と書くと本ファイルの import 環境では
`hp : Irreducible p` に展開されてしまい、`hp.dvd_choose_self` が
`Irreducible.dvd_choose_self` を探して落ちる。逐語:

```
Invalid field `dvd_choose_self`: The environment does not contain
`Irreducible.dvd_choose_self`, so it is not possible to project the field
`dvd_choose_self` from an expression  hp  of type `Irreducible p`
```

★`(hp : Nat.Prime p)` と書き直しても同じで、直し方は
★`Nat.Prime.dvd_choose_self hp` と完全名で書くことである
(そこで初めて `Unknown constant` が出て、import 不足だと分かる)。

## 何が言えたか

### §1–§2 抽象核(★分岐・付値・Galois・`p` 進の語彙が 1 語も出てこない)

| 宣言 | 内容 |
|---|---|
| `add_one_pow_eq_sum_choose` | 任意の環で `(D+1)^n = ∑_{l≤n} C(n,l) D^l` |
| ★`sum_pow_eq_sum_choose` | ★hockey stick `∑_{j<n}(D+1)^j = ∑_{l<n} C(n,l+1) D^l` |
| `iterate_sub_self_eq_sum_choose` | `(D+1)^n x − x = ∑_{l∈[1,n]} C(n,l)•D^l x` |
| ★`sum_iterate_sub_self_eq_sum_choose` | `∑_{j<n}((D+1)^j x − x) = ∑_{l∈[1,n)} C(n,l+1)•D^l x` |
| `norm_pow_apply_le` / `norm_pow_apply_le_self` | `‖D^l z‖ ≤ θ^l‖z‖` |
| ★★`norm_iterate_prime_sub_self_le_max` | ★`ε` の減衰の核 |
| ★★`norm_sum_iterate_sub_self_le_max` | ★平均の道の損失の核 |
| ★`norm_iterate_prime_sub_self_le_of_contract` | 収縮率 `θ` 版 |
| ★`norm_sum_iterate_sub_self_le_of_contract` | 収縮率 `θ` 版 |

### §3 台帳(実数だけ)

| 宣言 | 内容 |
|---|---|
| `sum_inv_pow_succ_eq` | `∑_{m=1}^{k} p^{−m} = (1−p^{−k})/(p−1)` |
| ★★`axDecay_exponent_eq_sub_geomSum` | ★`(1/(p−1))p^{1−k} = 1/(p−1) − ∑_{m=1}^{k−1}p^{−m}` |
| ★★`axDecay_one_le_averaging_loss` | ★`θ ≥ p^{−1/(p−1)}` なら `p·θ^{p−2} ≥ p^{1/(p−1)}` |

### §4 素朴な形の否定

| 宣言 | 内容 |
|---|---|
| ★`natDegree_minpoly_dvd_finrank` | `deg minpoly_K x ∣ [L:K]` |
| ★★★`exists_wildDepth_one_not_mem_prime_layer` | ★反例(`ℚ`, `ζ₇`, `p=3`) |
| ★`natDegree_minpoly_eq_of_mem_prime_layer` | `p` 進でも同じ(層に入るなら次数はちょうど `p`) |
| `axDecay_lt_axDecay_one` | `k ≥ 2` で `axDecay p k < axDecay p 1` |

### §5–§7 底変換・具体層・検算

| 宣言 | 内容 |
|---|---|
| `padicValNat_eq_of_dvd_of_dvd_mul` | 順分岐の底変換は深さを変えない(純 `ℕ`) |
| `smulEnd` / `smulSubOne` / `pow_smulEnd_apply` | `σ` と `σ−1` を `Module.End ℤ K̄` として見る |
| ★★`norm_smul_pow_prime_sub_self_le_of_contract` | ★`ε` の減衰(`AxEpsilonDecay` の同名定理を精密化) |
| ★★`norm_sum_orbit_sub_self_le` | 軌道和の評価 |
| ★★★`norm_sub_orbit_average_le_of_contract` | ★平均の道の 1 段の損失 |
| ★★★`axDecay_one_le_orbit_average_loss` | ★否定的測定(測定 4) |
| ★★★`axDecay_eq_axDecay_one_mul_gains` | ★指数の検算(測定 3) |

## ★残った穴(正直に書く)

★★`AxWildDescent K (axDecay p)` は閉じていない。★残りはちょうど 2 点である:

1. ★群論: 深さ `k` のとき、降下生成元 `σ` を「位数 `p^k` の `τ` の `p^{k−1}` 乗」の形に
   取れること(および `P ⊴ Q`、`Q/P` が巡回 `p` 次であること)。
   ★mathlib に「`p` 群の指数 `p` の部分群は正規」が見当たらないので自前で要る。
2. ★分岐: 位数 `p^m` の `σ` について `i(σ) ≤ e/(p^{m−1}(p−1))`。
   ★`m = 1` は `RamificationJumpBound` に在る。★`m ≥ 2` は無い。

★この 2 本が来れば、本ファイルの `norm_smul_pow_prime_sub_self_le_of_contract` と
`axDecay_eq_axDecay_one_mul_gains` がそのまま `axDecay p k` を与える(測定 3 の勘定)。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. ★配られた素朴な形(「`K` の直上に全分岐巡回 `p` 次の層」)は偽なので捨てた。
   反例を形式化した(測定 1)。★これは一般化ではなく訂正である。
2. ★配られた見込み(「平均の道の `p` で割る箇所が `p^{(1/(p−1))p^{1−k}}` に置き換わる」)も
   起きないことを測った(測定 4)。★得は損失側ではなく `ε` 側にある。
3. §1–§2 の抽象核は `M` に「超距離な半ノルム加法群」しか要求せず、
   `D` はただの `ℤ`-線形自己準同型でよい(★全単射性も等長性も `σ^p = 1` も要らない)。
   原典より弱い設定である。
4. ★反例(測定 1)は `ℚ` と `ℂ` の上で作った。原典は `p` 進体の話なので設定が違うが、
   壊れる理由(次数)は同じである。★`p` 進の具体元は構成していない(上に明記)。
5. ★`.src` は `AxEpsilonDecay.lean` / `WildDepthDescent.lean` と同じ項目
   (pGC 物理 p.6 Corollary 3.1)を指す。原典が独立に立てた項目ではない。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open Polynomial IntermediateField

/-! ## §1 抽象核 A -/

section RingCore

/-- ★★抽象核: 任意の環で `(D+1)^n = ∑_{l ≤ n} C(n,l) D^l`。 -/
theorem add_one_pow_eq_sum_choose {R : Type*} [Ring R] (D : R) (n : ℕ) :
    (D + 1) ^ n = ∑ l ∈ Finset.range (n + 1), (n.choose l : R) * D ^ l := by
  rw [(Commute.one_right D).add_pow n]
  refine Finset.sum_congr rfl (fun l _ => ?_)
  rw [one_pow, mul_one]
  exact ((Nat.cast_commute (n.choose l) (D ^ l)).symm).eq

def add_one_pow_eq_sum_choose.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★抽象核: `∑_{j<n}(D+1)^j = ∑_{l<n} C(n,l+1) D^l`(★hockey stick)。 -/
theorem sum_pow_eq_sum_choose {R : Type*} [Ring R] (D : R) (n : ℕ) :
    ∑ j ∈ Finset.range n, (D + 1) ^ j = ∑ l ∈ Finset.range n, (n.choose (l + 1) : R) * D ^ l := by
  induction n with
  | zero => simp
  | succ m ih =>
      rw [Finset.sum_range_succ, ih, add_one_pow_eq_sum_choose,
        Finset.sum_range_succ (fun l => ((m).choose l : R) * D ^ l) m,
        Finset.sum_range_succ (fun l => (((m + 1)).choose (l + 1) : R) * D ^ l) m,
        ← add_assoc]
      congr 1
      · rw [← Finset.sum_add_distrib]
        refine Finset.sum_congr rfl (fun l _ => ?_)
        rw [Nat.choose_succ_succ m l, Nat.cast_add, add_mul, add_comm]
      · simp

def sum_pow_eq_sum_choose.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end RingCore

/-! ## §2 抽象核 B -/

section ModuleCore

variable {M : Type*} [AddCommGroup M]

/-- `(D+1)^n x − x = ∑_{l ∈ [1,n]} C(n,l) • D^l x`。 -/
theorem iterate_sub_self_eq_sum_choose (D : Module.End ℤ M) (x : M) (n : ℕ) :
    ((D + 1) ^ n) x - x = ∑ l ∈ Finset.Icc 1 n, (n.choose l) • ((D ^ l) x) := by
  have key : ((D + 1) ^ n) x = ∑ l ∈ Finset.range (n + 1), (n.choose l) • ((D ^ l) x) := by
    have h := congrArg (fun F : Module.End ℤ M => F x) (add_one_pow_eq_sum_choose D n)
    simpa [LinearMap.sum_apply, Module.End.mul_apply, Module.End.natCast_apply] using h
  rw [key, Finset.range_eq_Ico, Finset.sum_eq_sum_Ico_succ_bot (Nat.succ_pos n)]
  have hIco : Finset.Ico 1 (n + 1) = Finset.Icc 1 n := by
    ext m; simp only [Finset.mem_Ico, Finset.mem_Icc]; omega
  rw [hIco]
  simp

def iterate_sub_self_eq_sum_choose.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★軌道和(跡)の展開: `∑_{j<n}((D+1)^j x − x) = ∑_{l ∈ [1,n)} C(n,l+1) • D^l x`。 -/
theorem sum_iterate_sub_self_eq_sum_choose (D : Module.End ℤ M) (x : M) (n : ℕ) :
    ∑ j ∈ Finset.range n, (((D + 1) ^ j) x - x)
      = ∑ l ∈ Finset.Ico 1 n, (n.choose (l + 1)) • ((D ^ l) x) := by
  rcases Nat.eq_zero_or_pos n with hn | hn
  · subst hn; simp
  have key : ∑ j ∈ Finset.range n, ((D + 1) ^ j) x
      = ∑ l ∈ Finset.range n, (n.choose (l + 1)) • ((D ^ l) x) := by
    have h := congrArg (fun F : Module.End ℤ M => F x) (sum_pow_eq_sum_choose D n)
    simpa [LinearMap.sum_apply, Module.End.mul_apply, Module.End.natCast_apply] using h
  have hsplit : ∑ j ∈ Finset.range n, (((D + 1) ^ j) x - x)
      = (∑ j ∈ Finset.range n, ((D + 1) ^ j) x) - n • x := by
    rw [Finset.sum_sub_distrib, Finset.sum_const, Finset.card_range]
  rw [hsplit, key, Finset.range_eq_Ico, Finset.sum_eq_sum_Ico_succ_bot hn]
  simp

def sum_iterate_sub_self_eq_sum_choose.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end ModuleCore

section NormCore

variable {M : Type*} [SeminormedAddCommGroup M]

/-- `‖D z‖ ≤ θ‖z‖` なら `‖D^l z‖ ≤ θ^l ‖z‖`。 -/
theorem norm_pow_apply_le (D : Module.End ℤ M) {θ : ℝ} (hθ0 : 0 ≤ θ)
    (hθ : ∀ z : M, ‖D z‖ ≤ θ * ‖z‖) : ∀ (l : ℕ) (z : M), ‖(D ^ l) z‖ ≤ θ ^ l * ‖z‖ := by
  intro l
  induction l with
  | zero => intro z; simp
  | succ m ih =>
      intro z
      rw [pow_succ, Module.End.mul_apply]
      refine le_trans (ih (D z)) ?_
      calc θ ^ m * ‖D z‖ ≤ θ ^ m * (θ * ‖z‖) :=
            mul_le_mul_of_nonneg_left (hθ z) (pow_nonneg hθ0 m)
        _ = θ ^ (m + 1) * ‖z‖ := by ring

def norm_pow_apply_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- `‖D z‖ ≤ ‖z‖` なら `‖D^l z‖ ≤ ‖z‖`。 -/
theorem norm_pow_apply_le_self (D : Module.End ℤ M) (hD : ∀ z : M, ‖D z‖ ≤ ‖z‖) :
    ∀ (l : ℕ) (z : M), ‖(D ^ l) z‖ ≤ ‖z‖ := by
  intro l
  induction l with
  | zero => intro z; simp
  | succ m ih =>
      intro z
      rw [pow_succ, Module.End.mul_apply]
      exact le_trans (ih (D z)) (hD z)

def norm_pow_apply_le_self.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end NormCore

section UltraCore

variable {M : Type*} [SeminormedAddCommGroup M] [IsUltrametricDist M]

/-- ★★★抽象核: `‖(D+1)^p x − x‖ ≤ max ‖D^p x‖ (c‖D x‖)`。 -/
theorem norm_iterate_prime_sub_self_le_max {p : ℕ} (hp : Nat.Prime p) (D : Module.End ℤ M)
    (hD : ∀ z : M, ‖D z‖ ≤ ‖z‖) {c : ℝ} (hc0 : 0 ≤ c)
    (hc : ∀ w : M, ‖(p : ℕ) • w‖ ≤ c * ‖w‖) (x : M) :
    ‖((D + 1) ^ p) x - x‖ ≤ max ‖(D ^ p) x‖ (c * ‖D x‖) := by
  rw [iterate_sub_self_eq_sum_choose]
  refine IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg
    (le_trans (norm_nonneg _) (le_max_left _ _)) (fun l hl => ?_)
  rw [Finset.mem_Icc] at hl
  rcases eq_or_lt_of_le hl.2 with heq | hlt
  · rw [heq, Nat.choose_self, one_smul]
    exact le_max_left _ _
  · obtain ⟨m, hm⟩ := Nat.Prime.dvd_choose_self hp (by omega) hlt
    rw [hm, mul_smul]
    refine le_trans (le_trans (hc _) ?_) (le_max_right _ _)
    refine mul_le_mul_of_nonneg_left ?_ hc0
    refine le_trans (IsUltrametricDist.norm_nsmul_le _ m) ?_
    have hstep : (D ^ l) x = (D ^ (l - 1)) (D x) := by
      rw [← Module.End.mul_apply, ← pow_succ]
      congr 2
      omega
    rw [hstep]
    exact norm_pow_apply_le_self D hD (l - 1) (D x)

def norm_iterate_prime_sub_self_le_max.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★抽象核: `‖∑_{j<p}((D+1)^j x − x)‖ ≤ max ‖D^{p−1} x‖ (c‖D x‖)`。 -/
theorem norm_sum_iterate_sub_self_le_max {p : ℕ} (hp : Nat.Prime p) (D : Module.End ℤ M)
    (hD : ∀ z : M, ‖D z‖ ≤ ‖z‖) {c : ℝ} (hc0 : 0 ≤ c)
    (hc : ∀ w : M, ‖(p : ℕ) • w‖ ≤ c * ‖w‖) (x : M) :
    ‖∑ j ∈ Finset.range p, (((D + 1) ^ j) x - x)‖
      ≤ max ‖(D ^ (p - 1)) x‖ (c * ‖D x‖) := by
  rw [sum_iterate_sub_self_eq_sum_choose]
  refine IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg
    (le_trans (norm_nonneg _) (le_max_left _ _)) (fun l hl => ?_)
  rw [Finset.mem_Ico] at hl
  rcases eq_or_lt_of_le (Nat.le_sub_one_of_lt hl.2) with heq | hlt
  · rw [show l + 1 = p by omega, Nat.choose_self, one_smul, heq]
    exact le_max_left _ _
  · obtain ⟨m, hm⟩ := Nat.Prime.dvd_choose_self hp (by omega) (by omega : l + 1 < p)
    rw [hm, mul_smul]
    refine le_trans (le_trans (hc _) ?_) (le_max_right _ _)
    refine mul_le_mul_of_nonneg_left ?_ hc0
    refine le_trans (IsUltrametricDist.norm_nsmul_le _ m) ?_
    have hstep : (D ^ l) x = (D ^ (l - 1)) (D x) := by
      rw [← Module.End.mul_apply, ← pow_succ]
      congr 2
      omega
    rw [hstep]
    exact norm_pow_apply_le_self D hD (l - 1) (D x)

def norm_sum_iterate_sub_self_le_max.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★収縮率 `θ` を入れた形(★`ε` の減衰)。 -/
theorem norm_iterate_prime_sub_self_le_of_contract {p : ℕ} (hp : Nat.Prime p) (D : Module.End ℤ M)
    {θ : ℝ} (hθ0 : 0 ≤ θ) (hθ1 : θ ≤ 1) (hθ : ∀ z : M, ‖D z‖ ≤ θ * ‖z‖) {c : ℝ} (hc0 : 0 ≤ c)
    (hc : ∀ w : M, ‖(p : ℕ) • w‖ ≤ c * ‖w‖) (x : M) :
    ‖((D + 1) ^ p) x - x‖ ≤ max (θ ^ (p - 1)) c * ‖D x‖ := by
  have h2 : 2 ≤ p := Nat.Prime.two_le hp
  have hD : ∀ z : M, ‖D z‖ ≤ ‖z‖ := fun z =>
    le_trans (hθ z) (by nlinarith [norm_nonneg z])
  refine le_trans (norm_iterate_prime_sub_self_le_max hp D hD hc0 hc x) (max_le ?_ ?_)
  · have hstep : (D ^ p) x = (D ^ (p - 1)) (D x) := by
      rw [← Module.End.mul_apply, ← pow_succ]
      congr 2
      omega
    rw [hstep]
    exact le_trans (norm_pow_apply_le D hθ0 hθ (p - 1) (D x))
      (mul_le_mul_of_nonneg_right (le_max_left _ _) (norm_nonneg _))
  · exact mul_le_mul_of_nonneg_right (le_max_right _ _) (norm_nonneg _)

def norm_iterate_prime_sub_self_le_of_contract.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★収縮率 `θ` を入れた形(★平均の道の損失)。 -/
theorem norm_sum_iterate_sub_self_le_of_contract {p : ℕ} (hp : Nat.Prime p) (D : Module.End ℤ M)
    {θ : ℝ} (hθ0 : 0 ≤ θ) (hθ1 : θ ≤ 1) (hθ : ∀ z : M, ‖D z‖ ≤ θ * ‖z‖) {c : ℝ} (hc0 : 0 ≤ c)
    (hc : ∀ w : M, ‖(p : ℕ) • w‖ ≤ c * ‖w‖) (x : M) :
    ‖∑ j ∈ Finset.range p, (((D + 1) ^ j) x - x)‖ ≤ max (θ ^ (p - 2)) c * ‖D x‖ := by
  have hD : ∀ z : M, ‖D z‖ ≤ ‖z‖ := fun z =>
    le_trans (hθ z) (by nlinarith [norm_nonneg z])
  refine le_trans (norm_sum_iterate_sub_self_le_max hp D hD hc0 hc x) (max_le ?_ ?_)
  · have h2 : 2 ≤ p := Nat.Prime.two_le hp
    have hstep : (D ^ (p - 1)) x = (D ^ (p - 2)) (D x) := by
      rw [← Module.End.mul_apply, ← pow_succ]
      congr 2
      omega
    rw [hstep]
    exact le_trans (norm_pow_apply_le D hθ0 hθ (p - 2) (D x))
      (mul_le_mul_of_nonneg_right (le_max_left _ _) (norm_nonneg _))
  · exact mul_le_mul_of_nonneg_right (le_max_right _ _) (norm_nonneg _)

def norm_sum_iterate_sub_self_le_of_contract.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end UltraCore

/-! ## §3 台帳(実数) -/

section Ledger

/-- 幾何級数(添字を 1 からずらした版)。 -/
theorem sum_inv_pow_succ_eq {p : ℕ} (hp : 1 < (p : ℝ)) (k : ℕ) :
    ∑ m ∈ Finset.range k, ((p : ℝ)⁻¹) ^ (m + 1) = (1 - ((p : ℝ)⁻¹) ^ k) / ((p : ℝ) - 1) := by
  have hp0 : (0:ℝ) < (p : ℝ) := lt_trans zero_lt_one hp
  have hpne : (p : ℝ) ≠ 0 := ne_of_gt hp0
  have hp1 : (p : ℝ) - 1 ≠ 0 := by
    intro h; rw [sub_eq_zero] at h; rw [h] at hp; exact lt_irrefl 1 hp
  have h1p : (1 : ℝ) - (p : ℝ) ≠ 0 := by
    intro h; rw [sub_eq_zero] at h; rw [← h] at hp; exact lt_irrefl 1 hp
  have hne : ((p : ℝ))⁻¹ ≠ 1 := by
    simp only [ne_eq, inv_eq_one]
    intro h; rw [h] at hp; exact lt_irrefl 1 hp
  have h1 : ∑ m ∈ Finset.range k, ((p : ℝ)⁻¹) ^ (m + 1)
      = ((p : ℝ)⁻¹) * ∑ m ∈ Finset.range k, ((p : ℝ)⁻¹) ^ m := by
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl (fun m _ => by rw [pow_succ'])
  have hden : ((p : ℝ)⁻¹ - 1) = (1 - (p : ℝ)) / (p : ℝ) := by field_simp
  rw [h1, geom_sum_eq hne k, hden, div_div_eq_mul_div, ← mul_div_assoc,
    div_eq_div_iff h1p hp1]
  field_simp
  ring

def sum_inv_pow_succ_eq.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★`axDecay` の指数の台帳: `(1/(p−1))·p^{1−k} = 1/(p−1) − ∑_{m=1}^{k−1} p^{−m}`。 -/
theorem axDecay_exponent_eq_sub_geomSum {p : ℕ} (hp : 1 < (p : ℝ)) (k : ℕ) :
    (1 / ((p : ℝ) - 1)) * (1 / (p : ℝ)) ^ (k - 1)
      = 1 / ((p : ℝ) - 1) - ∑ m ∈ Finset.range (k - 1), ((p : ℝ)⁻¹) ^ (m + 1) := by
  have hp1 : (p : ℝ) - 1 ≠ 0 := by
    intro h; rw [sub_eq_zero] at h; rw [h] at hp; exact lt_irrefl 1 hp
  rw [sum_inv_pow_succ_eq hp, one_div ((p : ℝ))]
  field_simp
  ring

def axDecay_exponent_eq_sub_geomSum.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★否定的な測定: 平均の道の 1 段の損失は `axDecay p 1` を下回れない。 -/
theorem axDecay_one_le_averaging_loss {p : ℕ} (hp : 1 < (p : ℝ)) (hp2 : 2 ≤ p) {θ : ℝ}
    (hθ : (p : ℝ) ^ (-(1 / ((p : ℝ) - 1))) ≤ θ) :
    (p : ℝ) ^ (1 / ((p : ℝ) - 1)) ≤ (p : ℝ) * θ ^ (p - 2) := by
  have hp0 : (0:ℝ) < (p : ℝ) := lt_trans zero_lt_one hp
  have hp1 : (p : ℝ) - 1 ≠ 0 := by
    intro h; rw [sub_eq_zero] at h; rw [h] at hp; exact lt_irrefl 1 hp
  have hb0 : (0:ℝ) < (p : ℝ) ^ (-(1 / ((p : ℝ) - 1))) := Real.rpow_pos_of_pos hp0 _
  have hstep : ((p : ℝ) ^ (-(1 / ((p : ℝ) - 1)))) ^ (p - 2) ≤ θ ^ (p - 2) :=
    pow_le_pow_left₀ (le_of_lt hb0) hθ _
  have hrhs : ((p : ℝ) ^ (-(1 / ((p : ℝ) - 1)))) ^ (p - 2)
      = (p : ℝ) ^ (-(1 / ((p : ℝ) - 1)) * ((p : ℝ) - 2)) := by
    rw [← Real.rpow_natCast ((p : ℝ) ^ (-(1 / ((p : ℝ) - 1)))) (p - 2),
      ← Real.rpow_mul (le_of_lt hp0), Nat.cast_sub hp2]
    norm_num
  have hmul : ∀ e : ℝ, (p : ℝ) * (p : ℝ) ^ e = (p : ℝ) ^ (1 + e) := by
    intro e; rw [Real.rpow_add hp0, Real.rpow_one]
  refine le_trans (le_of_eq ?_) (mul_le_mul_of_nonneg_left hstep (le_of_lt hp0))
  rw [hrhs, hmul]
  congr 1
  field_simp
  ring

def axDecay_one_le_averaging_loss.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end Ledger

/-! ## §4 素朴な形の否定 -/

/-- ★抽象核: `deg minpoly_K x` は `[L:K]` を割る。 -/
theorem natDegree_minpoly_dvd_finrank {K L : Type*} [Field K] [Field L] [Algebra K L]
    [FiniteDimensional K L] (x : L) : (minpoly K x).natDegree ∣ Module.finrank K L := by
  have hx : IsIntegral K x := Algebra.IsIntegral.isIntegral x
  refine ⟨Module.finrank K⟮x⟯ L, ?_⟩
  rw [← IntermediateField.adjoin.finrank hx]
  exact (Module.finrank_mul_finrank K K⟮x⟯ L).symm

def natDegree_minpoly_dvd_finrank.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**素朴な形の反例**(完全に形式化した)。 -/
theorem exists_wildDepth_one_not_mem_prime_layer :
    ∃ x : ℂ, padicValNat 3 (minpoly ℚ x).natDegree = 1 ∧
      ∀ F : IntermediateField ℚ ℂ, FiniteDimensional ℚ F → Module.finrank ℚ F = 3 → x ∉ F := by
  set ζ : ℂ := Complex.exp (2 * Real.pi * Complex.I / (7 : ℕ)) with hζdef
  have hζ : IsPrimitiveRoot ζ 7 := Complex.isPrimitiveRoot_exp 7 (by norm_num)
  have hdeg : (minpoly ℚ ζ).natDegree = 6 := by
    rw [← cyclotomic_eq_minpoly_rat hζ (by norm_num), natDegree_cyclotomic]
    decide
  refine ⟨ζ, ?_, ?_⟩
  · rw [hdeg, show (6 : ℕ) = 3 * 2 from rfl, padicValNat.mul (by norm_num) (by norm_num),
      padicValNat.self (by norm_num), padicValNat.eq_zero_of_not_dvd (by norm_num)]
  · intro F hFin hrank hmem
    have hdvd := natDegree_minpoly_dvd_finrank (K := ℚ) (L := F) ⟨ζ, hmem⟩
    have heq : minpoly ℚ (⟨ζ, hmem⟩ : F) = minpoly ℚ ζ :=
      IntermediateField.minpoly_eq (⟨ζ, hmem⟩ : F)
    rw [heq, hdeg, hrank] at hdvd
    exact absurd hdvd (by decide)

def exists_wildDepth_one_not_mem_prime_layer.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★`p` 進の層でも同じ: `p` 次の層に入るなら次数はちょうど `p`。 -/
theorem natDegree_minpoly_eq_of_mem_prime_layer {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x : K.closure) (F : IntermediateField K.carrier K.closure)
    [FiniteDimensional K.carrier F] (hF : Module.finrank K.carrier F = p) (hx : x ∈ F)
    (hdvd : p ∣ (minpoly K.carrier x).natDegree) :
    (minpoly K.carrier x).natDegree = p := by
  have hdvd2 := natDegree_minpoly_dvd_finrank (K := K.carrier) (L := F) ⟨x, hx⟩
  have heq : minpoly K.carrier (⟨x, hx⟩ : F) = minpoly K.carrier x :=
    IntermediateField.minpoly_eq (⟨x, hx⟩ : F)
  rw [heq, hF] at hdvd2
  exact Nat.dvd_antisymm hdvd2 hdvd

def natDegree_minpoly_eq_of_mem_prime_layer.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★`k ≥ 2` では `axDecay p k < axDecay p 1`。 -/
theorem axDecay_lt_axDecay_one (p : ℕ) [Fact p.Prime] {k : ℕ} (hk : 2 ≤ k) :
    axDecay p k < axDecay p 1 := by
  have h1 : (1:ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : Nat.Prime p).one_lt
  have hp0 : (0:ℝ) < (p : ℝ) := lt_trans zero_lt_one h1
  have hA : (0:ℝ) < 1 / ((p : ℝ) - 1) := by apply div_pos one_pos; linarith
  have hr0 : (0:ℝ) ≤ 1 / (p : ℝ) := by positivity
  have hr1 : 1 / (p : ℝ) < 1 := by rw [div_lt_one hp0]; linarith
  have hlt : (1 / (p : ℝ)) ^ (k - 1) < 1 := pow_lt_one₀ hr0 hr1 (by omega)
  rw [axDecay_one, axDecay]
  refine Real.rpow_lt_rpow_left_iff h1 |>.mpr ?_
  calc (1 / ((p : ℝ) - 1)) * (1 / (p : ℝ)) ^ (k - 1) < (1 / ((p : ℝ) - 1)) * 1 :=
        mul_lt_mul_of_pos_left hlt hA
    _ = 1 / ((p : ℝ) - 1) := by ring

def axDecay_lt_axDecay_one.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §5 順分岐(prime-to-p)の底変換は wild 深さを変えない -/

/-- ★抽象核(純 `ℕ`)。 -/
theorem padicValNat_eq_of_dvd_of_dvd_mul {p a b c : ℕ} [Fact p.Prime] (ha : a ≠ 0) (hb : b ≠ 0)
    (hc : c ≠ 0) (hba : b ∣ a) (habc : a ∣ b * c) (hpc : ¬ p ∣ c) :
    padicValNat p a = padicValNat p b := by
  have h1 : padicValNat p b ≤ padicValNat p a := padicValNat_le_of_dvd hb ha hba
  have h2 : padicValNat p a ≤ padicValNat p (b * c) :=
    padicValNat_le_of_dvd ha (Nat.mul_ne_zero hb hc) habc
  rw [padicValNat.mul hb hc, padicValNat.eq_zero_of_not_dvd hpc] at h2
  omega

def padicValNat_eq_of_dvd_of_dvd_mul.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §6 具体層 -/

section Concrete

variable {p : ℕ} [Fact p.Prime]

/-- `σ` の `K̄` への作用を `ℤ`-線形自己準同型として見たもの。 -/
noncomputable def smulEnd (K : PAdicLocalField p) (σ : K.absGal) : Module.End ℤ K.closure :=
  (DistribSMul.toAddMonoidHom K.closure σ).toIntLinearMap

def smulEnd.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

@[simp] theorem smulEnd_apply (K : PAdicLocalField p) (σ : K.absGal) (x : K.closure) :
    smulEnd K σ x = σ • x := rfl

/-- ★`D = σ − 1`。 -/
noncomputable def smulSubOne (K : PAdicLocalField p) (σ : K.absGal) : Module.End ℤ K.closure :=
  smulEnd K σ - 1

def smulSubOne.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

@[simp] theorem smulSubOne_apply (K : PAdicLocalField p) (σ : K.absGal) (x : K.closure) :
    smulSubOne K σ x = σ • x - x := rfl

theorem smulSubOne_add_one (K : PAdicLocalField p) (σ : K.absGal) :
    smulSubOne K σ + 1 = smulEnd K σ := by
  rw [smulSubOne, sub_add_cancel]

def smulSubOne_add_one.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

theorem pow_smulEnd_apply (K : PAdicLocalField p) (σ : K.absGal) (n : ℕ) (x : K.closure) :
    ((smulEnd K σ) ^ n) x = (σ ^ n) • x := by
  induction n generalizing x with
  | zero => simp
  | succ m ih => rw [pow_succ, Module.End.mul_apply, smulEnd_apply, ih, ← mul_smul, ← pow_succ]

def pow_smulEnd_apply.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

theorem norm_smulSubOne_le (K : PAdicLocalField p) (σ : K.absGal) (z : K.closure) :
    ‖smulSubOne K σ z‖ ≤ ‖z‖ := by
  have h := IsUltrametricDist.norm_add_le_max (σ • z) (-z)
  rw [← sub_eq_add_neg, norm_neg, norm_smul_closure] at h
  simpa using h

def norm_smulSubOne_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

theorem norm_nsmul_prime_closure (K : PAdicLocalField p) (w : K.closure) :
    ‖(p : ℕ) • w‖ ≤ ((p : ℝ))⁻¹ * ‖w‖ := by
  rw [nsmul_eq_mul, norm_mul, norm_natCast_p_closure]

def norm_nsmul_prime_closure.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★`‖σ^p x − x‖ ≤ max ‖(σ−1)^p x‖ (‖p‖·‖σ x − x‖)`。 -/
theorem norm_smul_pow_prime_sub_self_le (K : PAdicLocalField p) (σ : K.absGal) (x : K.closure) :
    ‖(σ ^ p) • x - x‖ ≤ max ‖((smulSubOne K σ) ^ p) x‖ (((p : ℝ))⁻¹ * ‖σ • x - x‖) := by
  have h := norm_iterate_prime_sub_self_le_max (M := K.closure) (Fact.out : Nat.Prime p)
    (smulSubOne K σ) (norm_smulSubOne_le K σ) (by positivity) (norm_nsmul_prime_closure K) x
  rwa [smulSubOne_add_one, pow_smulEnd_apply, smulSubOne_apply] at h

def norm_smul_pow_prime_sub_self_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★★**`ε` が減る**(本ファイルのもう 1 つの主結果)。

`σ` の変位を `θ` で押さえると `σ^p` の変位は `max (θ^{p−1}) ‖p‖` 倍に**縮む**。

★★これが `axDecay p k` の `p^{1−k}` の出どころである(モジュール docstring §測定 3)。 -/
theorem norm_smul_pow_prime_sub_self_le_of_contract (K : PAdicLocalField p) (σ : K.absGal)
    (x : K.closure) {θ : ℝ} (hθ0 : 0 ≤ θ) (hθ1 : θ ≤ 1)
    (hθ : ∀ z : K.closure, ‖σ • z - z‖ ≤ θ * ‖z‖) :
    ‖(σ ^ p) • x - x‖ ≤ max (θ ^ (p - 1)) ((p : ℝ))⁻¹ * ‖σ • x - x‖ := by
  have hcon : ∀ z : K.closure, ‖smulSubOne K σ z‖ ≤ θ * ‖z‖ := fun z => by
    rw [smulSubOne_apply]; exact hθ z
  have h := norm_iterate_prime_sub_self_le_of_contract (M := K.closure) (Fact.out : Nat.Prime p)
    (smulSubOne K σ) hθ0 hθ1 hcon (c := ((p : ℝ))⁻¹) (by positivity)
    (norm_nsmul_prime_closure K) x
  rwa [smulSubOne_add_one, pow_smulEnd_apply, smulSubOne_apply] at h

def norm_smul_pow_prime_sub_self_le_of_contract.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★`‖∑_{j<p}(σ^j x − x)‖ ≤ max ‖(σ−1)^{p−1} x‖ (‖p‖·‖σ x − x‖)`。 -/
theorem norm_sum_orbit_sub_self_le (K : PAdicLocalField p) (σ : K.absGal) (x : K.closure) :
    ‖∑ j ∈ Finset.range p, ((σ ^ j) • x - x)‖
      ≤ max ‖((smulSubOne K σ) ^ (p - 1)) x‖ (((p : ℝ))⁻¹ * ‖σ • x - x‖) := by
  have h := norm_sum_iterate_sub_self_le_max (M := K.closure) (Fact.out : Nat.Prime p)
    (smulSubOne K σ) (norm_smulSubOne_le K σ) (by positivity) (norm_nsmul_prime_closure K) x
  rw [smulSubOne_apply] at h
  simp only [smulSubOne_add_one, pow_smulEnd_apply] at h
  exact h

def norm_sum_orbit_sub_self_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★★**平均の道の 1 段の損失**(本ファイルの主結果)。

`y := (1/p)·∑_{j<p} σ^j x` に対し

  `‖x − y‖ ≤ max ( p·‖(σ−1)^{p−1} x‖ , ‖σ x − x‖ )`。

★`WildDepthDescent.exists_natDegree_minpoly_descent_div` の損失 `p·Δ` を
★**`(σ−1)^{p−1}` の分だけ割り引いた**形である。 -/
theorem norm_sub_orbit_average_le (K : PAdicLocalField p) (σ : K.absGal) (x : K.closure) :
    ‖x - ((p : ℕ) : K.closure)⁻¹ * ∑ j ∈ Finset.range p, (σ ^ j) • x‖
      ≤ max ((p : ℝ) * ‖((smulSubOne K σ) ^ (p - 1)) x‖) ‖σ • x - x‖ := by
  have hpp : (0:ℝ) < (p : ℝ) := by
    exact_mod_cast (Fact.out : Nat.Prime p).pos
  have hp0 : ((p : ℕ) : K.closure) ≠ 0 := by
    intro h
    have := norm_natCast_p_closure K
    rw [h, norm_zero] at this
    have : (0:ℝ) < ((p : ℝ))⁻¹ := by positivity
    linarith [this]
  have hsum : (∑ j ∈ Finset.range p, (σ ^ j) • x) - (p : ℕ) • x
      = ∑ j ∈ Finset.range p, ((σ ^ j) • x - x) := by
    rw [Finset.sum_sub_distrib, Finset.sum_const, Finset.card_range]
  have hx : x - ((p : ℕ) : K.closure)⁻¹ * (∑ j ∈ Finset.range p, (σ ^ j) • x)
      = -(((p : ℕ) : K.closure)⁻¹ * ((∑ j ∈ Finset.range p, (σ ^ j) • x) - (p : ℕ) • x)) := by
    rw [mul_sub, nsmul_eq_mul, ← mul_assoc, inv_mul_cancel₀ hp0, one_mul]
    abel
  rw [hx, norm_neg, norm_mul, norm_inv, norm_natCast_p_closure, inv_inv, hsum]
  refine le_trans (mul_le_mul_of_nonneg_left (norm_sum_orbit_sub_self_le K σ x)
    (le_of_lt hpp)) ?_
  rw [mul_max_of_nonneg _ _ (le_of_lt hpp), ← mul_assoc, mul_inv_cancel₀ (ne_of_gt hpp), one_mul]

def norm_sub_orbit_average_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★収縮率 `θ` を入れた形 —— ★**損失は `max (p·θ^{p−2}) 1`**。 -/
theorem norm_sub_orbit_average_le_of_contract (K : PAdicLocalField p) (σ : K.absGal)
    (x : K.closure) {θ : ℝ} (hθ0 : 0 ≤ θ) (hθ1 : θ ≤ 1)
    (hθ : ∀ z : K.closure, ‖σ • z - z‖ ≤ θ * ‖z‖) :
    ‖x - ((p : ℕ) : K.closure)⁻¹ * ∑ j ∈ Finset.range p, (σ ^ j) • x‖
      ≤ max ((p : ℝ) * θ ^ (p - 2)) 1 * ‖σ • x - x‖ := by
  have hpp : (0:ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : Nat.Prime p).pos
  have hcon : ∀ z : K.closure, ‖smulSubOne K σ z‖ ≤ θ * ‖z‖ := fun z => by
    rw [smulSubOne_apply]; exact hθ z
  have h := norm_sum_iterate_sub_self_le_of_contract (M := K.closure) (Fact.out : Nat.Prime p)
    (smulSubOne K σ) hθ0 hθ1 hcon (c := ((p : ℝ))⁻¹) (by positivity)
    (norm_nsmul_prime_closure K) x
  rw [smulSubOne_apply] at h
  simp only [smulSubOne_add_one, pow_smulEnd_apply] at h
  have hp0 : ((p : ℕ) : K.closure) ≠ 0 := by
    intro hz
    have hn := norm_natCast_p_closure K
    rw [hz, norm_zero] at hn
    have hpos : (0:ℝ) < ((p : ℝ))⁻¹ := by positivity
    linarith
  have hsum : (∑ j ∈ Finset.range p, (σ ^ j) • x) - (p : ℕ) • x
      = ∑ j ∈ Finset.range p, ((σ ^ j) • x - x) := by
    rw [Finset.sum_sub_distrib, Finset.sum_const, Finset.card_range]
  have hx : x - ((p : ℕ) : K.closure)⁻¹ * (∑ j ∈ Finset.range p, (σ ^ j) • x)
      = -(((p : ℕ) : K.closure)⁻¹ * ((∑ j ∈ Finset.range p, (σ ^ j) • x) - (p : ℕ) • x)) := by
    rw [mul_sub, nsmul_eq_mul, ← mul_assoc, inv_mul_cancel₀ hp0, one_mul]
    abel
  rw [hx, norm_neg, norm_mul, norm_inv, norm_natCast_p_closure, inv_inv, hsum]
  refine le_trans (mul_le_mul_of_nonneg_left h (le_of_lt hpp)) (le_of_eq ?_)
  rw [← mul_assoc, mul_max_of_nonneg _ _ (le_of_lt hpp), mul_inv_cancel₀ (ne_of_gt hpp)]

def norm_sub_orbit_average_le_of_contract.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**否定的な測定(具体層)** —— 平均の道の 1 段の損失は `axDecay p 1` 以上である。

★跳びの上界 `(p−1)i ≤ e_L`(`RamificationJumpBound`)は
`θ = ‖π‖^i ≥ p^{−1/(p−1)}` を強制するので、`max (p·θ^{p−2}) 1 ≥ axDecay p 1`。
⇒ ★**`k ≥ 2` の段で `axDecay p k` を出すことは、この道では原理的にできない**
(`axDecay_lt_axDecay_one`)。 -/
theorem axDecay_one_le_orbit_average_loss (p : ℕ) [Fact p.Prime] {θ : ℝ}
    (hθ : (p : ℝ) ^ (-(1 / ((p : ℝ) - 1))) ≤ θ) :
    axDecay p 1 ≤ max ((p : ℝ) * θ ^ (p - 2)) 1 := by
  have h1 : (1:ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : Nat.Prime p).one_lt
  rw [axDecay_one]
  exact le_trans (axDecay_one_le_averaging_loss h1 (Fact.out : Nat.Prime p).two_le hθ)
    (le_max_left _ _)

def axDecay_one_le_orbit_average_loss.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end Concrete

/-! ## §7 検算 -/

/-- ★★★★**指数の検算**(本ファイルの主結果)。

`axDecay p k` は「1 層の損失 `axDecay p 1 = p^{1/(p−1)}` に、
各段の `ε` の縮み `p^{−p^{−m}}`(`m = 1..k−1`)を掛けたもの」に**ちょうど一致する**。 -/
theorem axDecay_eq_axDecay_one_mul_gains (p : ℕ) [Fact p.Prime] (k : ℕ) :
    axDecay p k
      = axDecay p 1 * (p : ℝ) ^ (-(∑ m ∈ Finset.range (k - 1), ((p : ℝ)⁻¹) ^ (m + 1))) := by
  have h1 : (1:ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : Nat.Prime p).one_lt
  have hp0 : (0:ℝ) < (p : ℝ) := lt_trans zero_lt_one h1
  rw [axDecay_one, axDecay, ← Real.rpow_add hp0]
  congr 1
  rw [axDecay_exponent_eq_sub_geomSum h1 k]
  ring

def axDecay_eq_axDecay_one_mul_gains.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §8 使っている公理の一覧 -/

#print axioms add_one_pow_eq_sum_choose
#print axioms sum_pow_eq_sum_choose
#print axioms iterate_sub_self_eq_sum_choose
#print axioms sum_iterate_sub_self_eq_sum_choose
#print axioms norm_pow_apply_le
#print axioms norm_pow_apply_le_self
#print axioms norm_iterate_prime_sub_self_le_max
#print axioms norm_sum_iterate_sub_self_le_max
#print axioms norm_iterate_prime_sub_self_le_of_contract
#print axioms norm_sum_iterate_sub_self_le_of_contract
#print axioms sum_inv_pow_succ_eq
#print axioms axDecay_exponent_eq_sub_geomSum
#print axioms axDecay_one_le_averaging_loss
#print axioms natDegree_minpoly_dvd_finrank
#print axioms exists_wildDepth_one_not_mem_prime_layer
#print axioms natDegree_minpoly_eq_of_mem_prime_layer
#print axioms axDecay_lt_axDecay_one
#print axioms padicValNat_eq_of_dvd_of_dvd_mul
#print axioms smulEnd
#print axioms smulEnd_apply
#print axioms smulSubOne
#print axioms smulSubOne_apply
#print axioms smulSubOne_add_one
#print axioms pow_smulEnd_apply
#print axioms norm_smulSubOne_le
#print axioms norm_nsmul_prime_closure
#print axioms norm_smul_pow_prime_sub_self_le
#print axioms norm_smul_pow_prime_sub_self_le_of_contract
#print axioms norm_sum_orbit_sub_self_le
#print axioms norm_sub_orbit_average_le
#print axioms norm_sub_orbit_average_le_of_contract
#print axioms axDecay_one_le_orbit_average_loss
#print axioms axDecay_eq_axDecay_one_mul_gains

end ABC3.Found.PGC
