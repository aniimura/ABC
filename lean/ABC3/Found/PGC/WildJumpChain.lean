import ABC3.Found.PGC.CyclicLayerDescent
import Mathlib.Data.ZMod.Basic
import Mathlib.GroupTheory.PGroup

/-!
# [pGC] 「深さ `k` の降下鎖」—— `CyclicLayerDescent` が残した 2 点を測った

`Found/PGC/CyclicLayerDescent.lean` の docstring「残った穴」は、
`AxWildDescent K (axDecay p)` に残るものを次の 2 点と名指ししていた。

1. 群論: 深さ `k` のとき降下生成元 `σ` を「位数 `p^k` の `τ` の `p^{k−1}` 乗」に取れること。
2. 分岐: 位数 `p^m` の `σ` について `i(σ) ≤ e/(p^{m−1}(p−1))`(`m = 1` は在庫、`m ≥ 2` は無い)。

★★本ファイルの結論は 3 つである。

* ★点 1 の字面は偽である(反例を形式化した)。★正しい形も書いた(§5)。
* ★点 2 は閉じた(§1 + §6 + §7)。★入力は `m = 1` の場合だけである。
* ★★★しかし ★点 1・点 2 が揃っても `axDecay p k` は出ない★。
  `CyclicLayerDescent` の「測定 3」の勘定は★一般には閉じない★——
  不等号の向きが逆で、実在する分岐データが反例になる(§4)。

## ★★★測定 1 —— 点 1 の字面は偽(★一般の `p` 群に位数 `p^k` の元は無い)

`(ℤ/p)²` は位数 `p²` の `p` 群だが、全ての元の位数が `p` を割る
(`addOrderOf_dvd_prime_prod_zmod` / `addOrderOf_ne_sq_prod_zmod` /
`card_prod_zmod`、★`sorry` 無し、★分岐も付値も出てこない)。

★★配られた持ち場は「`exists_pgroup_descent` の `Q` は `H` ではなく `G` から取っている、
だから位数を上げる余地は `G` にある」と見立てていた。★★これも外れている。
★★★逐語で確かめられる反例(体の側。★本ファイルでは形式化していない、手計算である):

  `K = ℚ₂`、`x = ζ₈`。`[ℚ₂(ζ₈):ℚ₂] = φ(8) = 4 = 2²` なので wild 深さは `k = 2`。
  ところが Galois 閉包の群は `Gal(ℚ₂(ζ₈)/ℚ₂) ≅ (ℤ/8)ˣ = {1,3,5,7}` で、
  ★これは Klein 四元群であり、位数 `4` の元を★持たない★。

⇒ ★「深さ `k` なら `G` に位数 `p^k` の `τ` がある」は `p = 2`, `k = 2` で既に偽。
★`G` に取り替えても救われない。

## ★★★測定 2 —— 点 2 は閉じた(★入力は `m = 1` だけ)

★鍵は「跳び」を★ノルムの言葉に翻訳する★ことである。`σ` の収縮率を
`θ`(`∀z ‖σz − z‖ ≤ θ‖z‖`)、`c = ‖p‖` と置くと

| 分岐の言葉 | ノルムの言葉 |
|---|---|
| `θ = ‖π‖^{i(σ)}`、`c = ‖π‖^{e}` | ——(π も e も i も消える) |
| `(p−1)·i(σ) ≤ e`(`RamificationJumpBound`) | ★`c ≤ θ^{p−1}` |
| `p^{m}(p−1)·i(σ) ≤ e`(★点 2) | ★`c ≤ θ^{p^m(p−1)}` |

★★つまり点 2 は「`c ≤ θ^{p^m(p−1)}`」であり、★π も `e` も `i` も要らない。
そして `CyclicLayerDescent.norm_smul_pow_prime_sub_self_le_of_contract` は
「`σ` の収縮率が `θ` なら `σ^p` の収縮率は `max(θ^p, c·θ)`」を★既に与えている★。
⇒ 抽象核は★実数だけ★になる:

  ★`c ≤ b^n`・`b ≤ max(a^p, c·a)`・`0 < a ≤ 1`・`0 < c ≤ 1`・`1 ≤ n` ⟹ `c ≤ a^{p·n}`
  (`le_pow_mul_of_le_max`。★退化枝は `a = 1` を強制するので閉じる)。

これを `m` 回まわすと `le_pow_pow_of_chain`、具体層に落とすと
`norm_natCast_le_pow_of_chain`、付値語に戻すと `mul_le_of_norm_natCast_le_pow` で
★`p^m(p−1)·i ≤ e`、すなわち点 2 の字面そのものである。

★入力は `hbase`(= `σ^{p^m}` は位数 `p` なので `m = 1` の場合)ただ 1 つ。
★`m = 1` は `RamificationJumpBound.sub_one_mul_le_of_norm_natCast_eq_pow` が持っている。

## ★★★★測定 3 —— それでも `axDecay p k` は出ない(★測定 3 の勘定は閉じない)

`CyclicLayerDescent` の「測定 3」は、`σ_j = τ^{p^{j−1}}`(`j = 1..k`)と置き

  総損失の指数 `= s_k/(p−1) − Σ_{j<k} s_j`   (`s_j := (p−1)·i(σ_j)/e`)

が `axDecay p k` の指数 `p^{1−k}/(p−1)` に「ぴったり一致する」と書いていた。
★★ぴったり一致するのは `s_{j+1} = p·s_j` のとき★だけ★である。

★鎖から出るのは `s_{j+1} ≥ p·s_j`(★下から)である
(`norm_smul_pow_prime_sub_self_le_of_contract` は `σ^p` の収縮率の★上界★を与えるので、
跳びの★下界★になる)。ところが台帳が要求するのは `s_{j+1} ≤ p·s_j`(★上から、`ledger_le`)。
★★向きが逆である。

★★★実在する分岐データが反例になる(`not_ledger_of_ge`、★形式化した):

  `p = 3`、`K = ℚ₃(ζ₃)`、`L = ℚ₃(ζ₂₇)`、`G = Gal(L/K) ≅ ℤ/9`、`σ` は生成元。
  `Gal(ℚ₃(ζ₂₇)/ℚ₃)` の下付き分岐群は `G̃_u = Gal(L/ℚ₃(ζ_{3^k}))`(`3^{k−1} ≤ u ≤ 3^k−1`)なので
  `i(σ) = 2`、`i(σ³) = 8`、`e_L = φ(27) = 18`。
  ⇒ `s₀ = 2·2/18 = 2/9`、`s₁ = 2·8/18 = 8/9`。
  `s₁ ≥ 3·s₀`(`8/9 ≥ 6/9`)✓、`s₁ ≤ 1` ✓ ——★鎖の仮定は全部満たす★。
  総損失の指数 `= s₁/2 − s₀ = 4/9 − 2/9 = 2/9`。
  ところが `axDecay 3 2` の指数は `(1/2)·3^{−1} = 1/6`。
  ★`2/9 > 1/6` なので★勘定が閉じない★。

★一般には `i(σ^p) = p·i(σ) + (p−1)`(円分の場合)で、`p ≥ 3` では
総損失の指数が `2/p²` になり `p^{−1}/(p−1)` を超える。
★Kummer 塔(`K = ℚ_p(ζ_{p²})`, `x = p^{1/p²}`)では `s_{j+1} = p·s_j` が★等号で成立する★ので
「測定 3」の手計算は通っていた。★等号が特殊だったのである。

★★これは `axDecay p k` が偽だという意味ではない。上の `ℚ₃(ζ₂₇)` の例でも
実際の降下(`x' = 1` を取る)は損失 `< 1` で済む。★★偽なのは「鎖による勘定」の方である。

## ★★★測定 4 —— 点 1 の正しい形(★鎖を捨てると `p^{1−k}` が出る)

鎖が要らない道がある。★★「最初の跳び」を使う:

* `u₁ :=` `L/K` の最小の跳び。上付き番号付けでは `v₁ = u₁`(`φ` は `u ≤ u₁` で恒等)。
* `G` の指数 `p` の商(`N ⊴ G`, `[G:N] = p`)で `u₁` を見る中間体 `L^N/K` は★次数 `p`★で、
  その跳びは `v₁ = u₁` そのもの。`m = 1` の評価をそこに当てると

  ★★`(p−1)·u₁ ≤ e_{L^N} = p·e_K`  (★底 `K` の分岐指数。`e_L` ではない)

* `e_L = [L:K]·e_K` かつ `p^k ∣ [L:K]` なので、これを掛け合わせて

  ★★★`p^{k}·((p−1)·u₁) ≤ p·e_L`、すなわち `(p−1)·p^{k−1}·u₁ ≤ e_L`

  (`sub_one_mul_pow_le_of_first_jump`、★純 `ℕ`、★形式化した)。

⇒ ★点 2 の結論(`i ≤ e_L/(p^{k−1}(p−1))`)が★位数 `p^k` の元を一切使わずに★出る。
★これが点 1 の「修理」である: 要るのは「位数 `p^k` の `τ`」ではなく
★「最初の跳びを見る指数 `p` の商」と「`p^k ∣ [L:K]`」だけ。
★★`(ℤ/p)^k` でも成立する(元の位数は全部 `p` でよい)。

★★本ファイルが形式化していないのは、この道の
「`(p−1)·u₁ ≤ p·e_K`」(= Herbrand `φ` が `u ≤ u₁` で恒等であること + `m = 1` の評価)
である。★これは新しい節点で、`HerbrandFunction.lean` の言葉で書けるはずである。

## 在庫の測定(★自分で測った。★MCP は 1 度も呼んでいない。★引いたコマンドを残す)

```
grep -c "\tABC3\..*\.<名前>\t" .cache/decl-index.txt   （16 個。全部 0 = 衝突なし）
grep -n "norm_natCast_p_closure|norm_smul_closure" .cache/decl-index.txt
  → norm_natCast_p_closure : ‖((p:ℕ) : K.closure)‖ = ((p:ℕ):ℝ)⁻¹（AbsClosureModules.lean:294）
     ＋ _pos / _lt_one も在る。★c = ‖p‖ の 3 点セットが揃っている。
grep -n "\^ \(p - 1\)" .cache/decl-index.txt
  → ★「σ ∈ K.absGal の収縮率 θ について ‖p‖ ≤ θ^{p−1}」は★木に無い★。
     RamificationJumpBound は π と多項式の言葉で持っている(sub_one_mul_le_of_norm_natCast_eq_pow)。
     ⇒ 本ファイルは `hbase` として仮定に出した。★「測って無かった」である。
```

★配管で 1 つ踏んだ: `le_or_lt` は★もう無い★(`le_or_gt` に改名)。逐語:

```
error(lean.unknownIdentifier): Unknown identifier `le_or_lt`
```

## 逸脱の記録(CLAUDE.md「逸脱」)

1. ★配られた点 1 の字面(位数 `p^k` の `τ`)は偽なので捨てた。反例を形式化した(§3)。
   ★正しい形(最初の跳び + 次数)を §5 に置いた。★これは一般化ではなく訂正である。
2. ★点 2 は配られた向き(`≤`)のまま真で、閉じた。ただし★これだけでは `axDecay p k` は
   出ない★ことを §4 で測った。★配られた持ち場の「この 2 本が来れば axDecay p k を与える」は
   ★起きない★。
3. §1・§2・§5 の抽象核は分岐・付値・Galois・`p` 進の語彙が 1 語も出てこない。
   `σ` も `π` も `e` も現れず、実数と `ℕ` だけである。
4. ★`hbase`(`m = 1` の評価を `K.absGal` の収縮率の言葉で述べたもの)は仮定として置いた。
   ★木に無いことを測ったうえでの仮定である(上の在庫の測定)。
5. `.src` は `CyclicLayerDescent.lean` と同じ項目(pGC 物理 p.6 Corollary 3.1)を指す。
   原典が独立に立てた項目ではない。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC

/-! ## §1 抽象核 A —— 実数だけの鎖(★点 2 の心臓) -/

section Chain

/-- ★★★抽象核(実数だけ): `c ≤ b^n` と `b ≤ max(a^p, c·a)` から `c ≤ a^{p·n}`。

★退化枝(`b ≤ c·a` の側)では `c ≤ c^n a^n ≤ c^n ≤ c` が全部等号になり `a = 1` が
強制されるので、結論は空虚に真になる。★場合分けはこの 2 つだけである。 -/
theorem le_pow_mul_of_le_max {p n : ℕ} (hn : 1 ≤ n) {c a b : ℝ}
    (hc0 : 0 < c) (hc1 : c ≤ 1) (ha0 : 0 < a) (ha1 : a ≤ 1) (hb0 : 0 ≤ b)
    (hstep : b ≤ max (a ^ p) (c * a)) (hprev : c ≤ b ^ n) :
    c ≤ a ^ (p * n) := by
  rcases le_or_gt b (a ^ p) with h | h
  · calc c ≤ b ^ n := hprev
      _ ≤ (a ^ p) ^ n := pow_le_pow_left₀ hb0 h n
      _ = a ^ (p * n) := by rw [← pow_mul]
  · have hb : b ≤ c * a := by
      rcases max_cases (a ^ p) (c * a) with ⟨he, _⟩ | ⟨he, _⟩
      · rw [he] at hstep; linarith
      · rw [he] at hstep; exact hstep
    have h1 : c ≤ (c * a) ^ n := le_trans hprev (pow_le_pow_left₀ hb0 hb n)
    have h2 : (c * a) ^ n = c ^ n * a ^ n := mul_pow c a n
    have h3 : a ^ n ≤ 1 := pow_le_one₀ ha0.le ha1
    have h4 : c ^ n ≤ c := pow_le_of_le_one hc0.le hc1 (by omega)
    have hane : a ^ n = 1 := by nlinarith [pow_pos hc0 n, pow_pos ha0 n]
    have ha : a = 1 := by
      by_contra hne
      have hlt : a < 1 := lt_of_le_of_ne ha1 hne
      have := pow_lt_one₀ ha0.le hlt (by omega : n ≠ 0)
      linarith
    rw [ha, one_pow]
    exact hc1

def le_pow_mul_of_le_max.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★抽象核(実数だけ): 鎖 `θ_{j+1} ≤ max(θ_j^p, c·θ_j)` を `m` 回まわす。

★★これが「点 2 を `m = 1` の場合に帰着させる」全体である。分岐は 1 語も出てこない。 -/
theorem le_pow_pow_of_chain {p : ℕ} (hp : 1 ≤ p) {c : ℝ} (hc0 : 0 < c) (hc1 : c ≤ 1)
    {θ : ℕ → ℝ} (hθ0 : ∀ j, 0 < θ j) (hθ1 : ∀ j, θ j ≤ 1)
    (hstep : ∀ j, θ (j + 1) ≤ max (θ j ^ p) (c * θ j)) :
    ∀ (m n : ℕ), 1 ≤ n → c ≤ θ m ^ n → c ≤ θ 0 ^ (p ^ m * n) := by
  intro m
  induction m with
  | zero => intro n _ h; simpa using h
  | succ k ih =>
      intro n hn h
      have h1 : c ≤ θ k ^ (p * n) :=
        le_pow_mul_of_le_max hn hc0 hc1 (hθ0 k) (hθ1 k) (hθ0 (k + 1)).le (hstep k) h
      have hpn : 1 ≤ p * n := le_trans hp (Nat.le_mul_of_pos_right p hn)
      have h2 := ih (p * n) hpn h1
      have he : p ^ k * (p * n) = p ^ (k + 1) * n := by ring
      rwa [he] at h2

def le_pow_pow_of_chain.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end Chain

/-! ## §2 台帳(実数だけ) —— ★勘定が閉じるための正しい仮定 -/

section Ledger

/-- `s (j+d) ≤ p^d · s j`(`d` 段の反復)。 -/
theorem step_le {p : ℕ} {s : ℕ → ℝ} {m : ℕ}
    (hrec : ∀ j, j < m → s (j + 1) ≤ (p : ℝ) * s j) (hp0 : (0:ℝ) ≤ (p:ℝ)) :
    ∀ (d j : ℕ), j + d ≤ m → s (j + d) ≤ (p : ℝ) ^ d * s j := by
  intro d
  induction d with
  | zero => intro j _; simp
  | succ e ih =>
      intro j hj
      have h1 : s (j + e + 1) ≤ (p:ℝ) * s (j + e) := hrec (j + e) (by omega)
      have h2 : s (j + e) ≤ (p:ℝ) ^ e * s j := ih j (by omega)
      have hje : (j + (e + 1)) = (j + e + 1) := by omega
      rw [hje]
      calc s (j + e + 1) ≤ (p:ℝ) * s (j + e) := h1
        _ ≤ (p:ℝ) * ((p:ℝ) ^ e * s j) := mul_le_mul_of_nonneg_left h2 hp0
        _ = (p:ℝ) ^ (e + 1) * s j := by ring

def step_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★台帳(実数だけ)。`CyclicLayerDescent` 測定 3 の勘定が閉じるための★正しい仮定★は

  `s_{j+1} ≤ p·s_j`(★上から)

であって、鎖から出る `s_{j+1} ≥ p·s_j`(★下から)ではない。

分母を払った形で書いてある: `s_m ≤ (p−1)·Σ_{j<m} s_j + p^{−m}` は
`s_m/(p−1) − Σ_{j<m}s_j ≤ p^{−m}/(p−1)` と同値で、右辺が `axDecay p (m+1)` の指数である。 -/
theorem ledger_le {p : ℕ} (hp : 1 < (p : ℝ)) {s : ℕ → ℝ} (m : ℕ)
    (hrec : ∀ j, j < m → s (j + 1) ≤ (p : ℝ) * s j) (hsm : s m ≤ 1) :
    s m ≤ ((p : ℝ) - 1) * (∑ j ∈ Finset.range m, s j) + ((p : ℝ)⁻¹) ^ m := by
  have hp0 : (0:ℝ) < (p : ℝ) := lt_trans zero_lt_one hp
  have hp1 : (0:ℝ) < (p : ℝ) - 1 := by linarith
  have hkey : ∀ j ∈ Finset.range m, ((p:ℝ)⁻¹) ^ (m - j) * s m ≤ s j := by
    intro j hj
    rw [Finset.mem_range] at hj
    have h := step_le hrec (le_of_lt hp0) (m - j) j (by omega)
    have hmj : j + (m - j) = m := by omega
    rw [hmj] at h
    rw [inv_pow, ← div_eq_inv_mul, div_le_iff₀ (by positivity)]
    calc s m ≤ (p:ℝ) ^ (m - j) * s j := h
      _ = s j * (p:ℝ) ^ (m - j) := by ring
  have hsum : (∑ j ∈ Finset.range m, ((p:ℝ)⁻¹) ^ (m - j)) * s m
      ≤ ∑ j ∈ Finset.range m, s j := by
    rw [Finset.sum_mul]
    exact Finset.sum_le_sum hkey
  have hgeom : (∑ j ∈ Finset.range m, ((p:ℝ)⁻¹) ^ (m - j))
      = (1 - ((p:ℝ)⁻¹) ^ m) / ((p:ℝ) - 1) := by
    rw [← sum_inv_pow_succ_eq hp m,
      ← Finset.sum_range_reflect (fun i => ((p:ℝ)⁻¹) ^ (i + 1)) m]
    refine Finset.sum_congr rfl (fun j hj => ?_)
    rw [Finset.mem_range] at hj
    congr 1
    omega
  rw [hgeom, div_mul_eq_mul_div, div_le_iff₀ hp1] at hsum
  have hq : (0:ℝ) ≤ ((p:ℝ)⁻¹) ^ m := by positivity
  nlinarith [hsum, hsm, hq]

def ledger_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end Ledger

/-! ## §3 ★点 1 の字面の反例(純群論) -/

section PGroupCounterexample

/-- `(ℤ/p)²` の元は全部位数が `p` を割る。 -/
theorem addOrderOf_dvd_prime_prod_zmod (p : ℕ) [Fact p.Prime] (g : ZMod p × ZMod p) :
    addOrderOf g ∣ p := by
  refine addOrderOf_dvd_of_nsmul_eq_zero ?_
  ext <;> simp [nsmul_eq_mul]

def addOrderOf_dvd_prime_prod_zmod.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★点 1 の字面の反例: 位数 `p²` の `p` 群 `(ℤ/p)²` に★位数 `p²` の元は無い★。

★配られた「深さ `k` なら位数 `p^k` の `τ` が取れる」は `k = 2` で既に偽。
★体の側の逐語の証人(手計算): `K = ℚ₂`, `x = ζ₈`。`[ℚ₂(ζ₈):ℚ₂] = 4` で wild 深さ `2`、
`Gal(ℚ₂(ζ₈)/ℚ₂) ≅ (ℤ/8)ˣ` は Klein 四元群なので位数 `4` の元を持たない。 -/
theorem addOrderOf_ne_sq_prod_zmod (p : ℕ) [Fact p.Prime] (g : ZMod p × ZMod p) :
    addOrderOf g ≠ p ^ 2 := by
  have hp : 1 < p := (Fact.out : p.Prime).one_lt
  have h1 : addOrderOf g ≤ p := Nat.le_of_dvd (by omega) (addOrderOf_dvd_prime_prod_zmod p g)
  have h2 : p < p ^ 2 := by nlinarith [sq_nonneg p]
  omega

def addOrderOf_ne_sq_prod_zmod.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- 反例の群の位数は確かに `p²`。 -/
theorem card_prod_zmod (p : ℕ) [Fact p.Prime] : Nat.card (ZMod p × ZMod p) = p ^ 2 := by
  have hp : 0 < p := (Fact.out : p.Prime).pos
  haveI : NeZero p := ⟨by omega⟩
  simp [Nat.card_eq_fintype_card, sq]

def card_prod_zmod.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end PGroupCounterexample

/-! ## §4 ★★測定 3 の勘定の反例(実数だけ) -/

/-- ★★★`CyclicLayerDescent` 測定 3 の勘定は、鎖から出る仮定(`s_{j+1} ≥ p·s_j`)だけでは
★閉じない★。

★数値は実在する分岐データである: `p = 3`、`L = ℚ₃(ζ₂₇)`、`K = ℚ₃(ζ₃)`、
`σ` は `Gal(L/K) ≅ ℤ/9` の生成元。`i(σ) = 2`、`i(σ³) = 8`、`e_L = 18` なので
`s₀ = 2/9`、`s₁ = 8/9`。★鎖の仮定(`3·s₀ ≤ s₁`、`s₁ ≤ 1`)は満たすが、
総損失の指数 `s₁/(3−1) − s₀ = 2/9` が `axDecay 3 2` の指数 `3⁻¹/(3−1) = 1/6` を超える。 -/
theorem not_ledger_of_ge :
    ∃ s₀ s₁ : ℝ, 0 ≤ s₀ ∧ (3 : ℝ) * s₀ ≤ s₁ ∧ s₁ ≤ 1 ∧
      ¬ (s₁ / ((3 : ℝ) - 1) - s₀ ≤ ((3 : ℝ))⁻¹ / ((3 : ℝ) - 1)) :=
  ⟨2 / 9, 8 / 9, by norm_num, by norm_num, by norm_num, by norm_num⟩

def not_ledger_of_ge.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §5 ★点 1 の正しい形(純 `ℕ`) -/

/-- ★★★抽象核(純 `ℕ`): 「最初の跳び」からの評価。

`(p−1)·u ≤ p·e_K`(底 `K` での `m = 1` の評価)と `p^k ∣ d`(`d = [L:K]`)、
`e_L = d·e_K` から `p^k·((p−1)·u) ≤ p·e_L`、すなわち `(p−1)·p^{k−1}·u ≤ e_L`。

★★これが点 1 の修理である: ★位数 `p^k` の元は要らない★。
要るのは「最初の跳びを見る指数 `p` の商」と「`p^k ∣ [L:K]`」だけで、
`(ℤ/p)^k` のように指数 `p` の群でも成立する。 -/
theorem sub_one_mul_pow_le_of_first_jump {p k u eK d eL : ℕ}
    (hd : p ^ k ∣ d) (hd0 : 0 < d) (heL : eL = d * eK)
    (h : (p - 1) * u ≤ p * eK) :
    p ^ k * ((p - 1) * u) ≤ p * eL := by
  have hle : p ^ k ≤ d := Nat.le_of_dvd hd0 hd
  calc p ^ k * ((p - 1) * u) ≤ p ^ k * (p * eK) := Nat.mul_le_mul_left _ h
    _ = p * (p ^ k * eK) := by ring
    _ ≤ p * (d * eK) := Nat.mul_le_mul_left _ (Nat.mul_le_mul_right _ hle)
    _ = p * eL := by rw [heL]

def sub_one_mul_pow_le_of_first_jump.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §6 具体層 —— `σ^{p^j}` の収縮率の列 -/

/-- `σ` の収縮率 `θ` から `σ^{p^j}` の収縮率を作る列
(`CyclicLayerDescent.norm_smul_pow_prime_sub_self_le_of_contract` を `j` 回)。 -/
noncomputable def contractSeq (p : ℕ) (θ : ℝ) : ℕ → ℝ
  | 0 => θ
  | j + 1 => max (contractSeq p θ j ^ p) (((p : ℝ))⁻¹ * contractSeq p θ j)

def contractSeq.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

@[simp] theorem contractSeq_zero (p : ℕ) (θ : ℝ) : contractSeq p θ 0 = θ := rfl

theorem contractSeq_succ (p : ℕ) (θ : ℝ) (j : ℕ) :
    contractSeq p θ (j + 1)
      = max (contractSeq p θ j ^ p) (((p : ℝ))⁻¹ * contractSeq p θ j) := rfl

theorem contractSeq_nonneg (p : ℕ) {θ : ℝ} (hθ : 0 ≤ θ) (j : ℕ) : 0 ≤ contractSeq p θ j := by
  induction j with
  | zero => simpa using hθ
  | succ k ih => rw [contractSeq_succ]; exact le_trans (pow_nonneg ih p) (le_max_left _ _)

theorem contractSeq_pos (p : ℕ) {θ : ℝ} (hθ : 0 < θ) (j : ℕ) : 0 < contractSeq p θ j := by
  induction j with
  | zero => simpa using hθ
  | succ k ih => rw [contractSeq_succ]; exact lt_of_lt_of_le (pow_pos ih p) (le_max_left _ _)

theorem contractSeq_le_one {p : ℕ} (hp : 1 ≤ (p : ℝ)) {θ : ℝ} (hθ0 : 0 ≤ θ) (hθ1 : θ ≤ 1)
    (j : ℕ) : contractSeq p θ j ≤ 1 := by
  have hp0 : (0:ℝ) < (p:ℝ) := lt_of_lt_of_le zero_lt_one hp
  have hinv : ((p:ℝ))⁻¹ ≤ 1 := by rw [inv_le_one₀ hp0]; exact hp
  induction j with
  | zero => simpa using hθ1
  | succ k ih =>
      have hk0 : 0 ≤ contractSeq p θ k := contractSeq_nonneg p hθ0 k
      rw [contractSeq_succ]
      refine max_le (pow_le_one₀ hk0 ih) ?_
      have hinv0 : (0:ℝ) ≤ ((p:ℝ))⁻¹ := le_of_lt (by positivity)
      nlinarith

section Concrete

variable {p : ℕ} [Fact p.Prime]

/-- ★★`contractSeq p θ j` は `σ^{p^j}` の収縮率である。 -/
theorem contractSeq_isRate (K : PAdicLocalField p) (σ : K.absGal) {θ : ℝ}
    (hθ0 : 0 ≤ θ) (hθ1 : θ ≤ 1) (hθ : ∀ z : K.closure, ‖σ • z - z‖ ≤ θ * ‖z‖) :
    ∀ (j : ℕ) (z : K.closure), ‖(σ ^ p ^ j) • z - z‖ ≤ contractSeq p θ j * ‖z‖ := by
  have hp1 : (1:ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  have hple : 1 ≤ p := (Fact.out : p.Prime).one_lt.le
  intro j
  induction j with
  | zero => intro z; simpa using hθ z
  | succ k ih =>
      intro z
      have hk0 : 0 ≤ contractSeq p θ k := contractSeq_nonneg p hθ0 k
      have hk1 : contractSeq p θ k ≤ 1 := contractSeq_le_one (le_of_lt hp1) hθ0 hθ1 k
      have h := norm_smul_pow_prime_sub_self_le_of_contract K (σ ^ p ^ k) z hk0 hk1 ih
      have hpow : (σ ^ p ^ k) ^ p = σ ^ p ^ (k + 1) := by rw [← pow_mul, ← pow_succ]
      rw [hpow] at h
      refine le_trans h ?_
      have h2 : max (contractSeq p θ k ^ (p - 1)) ((p:ℝ))⁻¹ * ‖(σ ^ p ^ k) • z - z‖
          ≤ max (contractSeq p θ k ^ (p - 1)) ((p:ℝ))⁻¹ * (contractSeq p θ k * ‖z‖) := by
        refine mul_le_mul_of_nonneg_left (ih z) ?_
        exact le_trans (by positivity) (le_max_right _ _)
      refine le_trans h2 (le_of_eq ?_)
      rw [contractSeq_succ, ← mul_assoc, max_mul_of_nonneg _ _ hk0]
      congr 2
      · rw [← pow_succ]
        congr 1
        omega

def contractSeq_isRate.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★★点 2(ノルム語)。

`σ` の収縮率が `θ` で、`σ^{p^m}`(位数 `p` の元)について `m = 1` の評価
`‖p‖ ≤ η^{p−1}`(`η` はその収縮率)が使えるなら、

  ★`‖p‖ ≤ θ^{p^m(p−1)}`。

★付値語に直すと `p^m(p−1)·i(σ) ≤ e`、すなわち `i(σ) ≤ e/(p^m(p−1))` ——
配られた点 2 の字面そのものである(`σ` の位数は `p^{m+1}`)。

★仮定 `hbase` は「`m = 1` の場合」だけで、`RamificationJumpBound` が持っている内容である
(★ただし `K.absGal` の収縮率の言葉に翻訳した形は木に無い。モジュール docstring の在庫の測定)。 -/
theorem norm_natCast_le_pow_of_chain (K : PAdicLocalField p) (σ : K.absGal) {θ : ℝ}
    (hθ0 : 0 < θ) (hθ1 : θ ≤ 1) (hθ : ∀ z : K.closure, ‖σ • z - z‖ ≤ θ * ‖z‖) (m : ℕ)
    (hbase : ∀ η : ℝ, 0 ≤ η → (∀ z : K.closure, ‖(σ ^ p ^ m) • z - z‖ ≤ η * ‖z‖) →
      ‖((p : ℕ) : K.closure)‖ ≤ η ^ (p - 1)) :
    ‖((p : ℕ) : K.closure)‖ ≤ θ ^ (p ^ m * (p - 1)) := by
  have hp2 : 2 ≤ p := (Fact.out : p.Prime).two_le
  have hp1 : (1:ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  have hc : ‖((p : ℕ) : K.closure)‖ = ((p:ℝ))⁻¹ := norm_natCast_p_closure K
  have hc0 : 0 < ‖((p : ℕ) : K.closure)‖ := norm_natCast_p_closure_pos K
  have hc1 : ‖((p : ℕ) : K.closure)‖ ≤ 1 := le_of_lt (norm_natCast_p_closure_lt_one K)
  have hstep : ∀ j, contractSeq p θ (j + 1)
      ≤ max (contractSeq p θ j ^ p) (‖((p : ℕ) : K.closure)‖ * contractSeq p θ j) := by
    intro j; rw [hc, contractSeq_succ]
  have hbase' : ‖((p : ℕ) : K.closure)‖ ≤ contractSeq p θ m ^ (p - 1) :=
    hbase _ (contractSeq_nonneg p (le_of_lt hθ0) m)
      (contractSeq_isRate K σ (le_of_lt hθ0) hθ1 hθ m)
  have hmain := le_pow_pow_of_chain (Nat.Prime.one_lt (Fact.out : p.Prime)).le hc0 hc1
    (contractSeq_pos p hθ0) (contractSeq_le_one (le_of_lt hp1) (le_of_lt hθ0) hθ1)
    hstep m (p - 1) (by omega) hbase'
  simpa using hmain

def norm_natCast_le_pow_of_chain.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end Concrete

/-! ## §7 付値語への翻訳(★点 2 の字面) -/

/-- ★★★点 2 の字面(付値語): `‖n‖ = ‖π‖^e` かつ `‖n‖ ≤ (‖π‖^i)^N` なら `N·i ≤ e`。

`N = p^m(p−1)`、`n = p`、`i = i(σ)`、`e = e_L` を入れると

  ★`p^m(p−1)·i(σ) ≤ e_L`、すなわち `i(σ) ≤ e_L/(p^m(p−1))`

——★配られた点 2(位数 `p^{m+1}` の `σ`)そのものである。
★`N = p − 1`(`m = 0`)は `RamificationJumpBound.sub_one_mul_le_of_norm_natCast_eq_pow`。 -/
theorem mul_le_of_norm_natCast_le_pow {L : Type*} [NormedField L] {N i e n : ℕ} {π : L}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (he : ‖(n : L)‖ = ‖π‖ ^ e)
    (h : ‖(n : L)‖ ≤ (‖π‖ ^ i) ^ N) : N * i ≤ e := by
  rw [he, ← pow_mul] at h
  have h2 := (pow_le_pow_iff_right_of_lt_one₀ hπ0 hπ1).mp h
  rw [Nat.mul_comm]
  exact h2

def mul_le_of_norm_natCast_le_pow.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §8 使っている公理の一覧 -/

#print axioms le_pow_mul_of_le_max
#print axioms le_pow_pow_of_chain
#print axioms step_le
#print axioms ledger_le
#print axioms addOrderOf_dvd_prime_prod_zmod
#print axioms addOrderOf_ne_sq_prod_zmod
#print axioms card_prod_zmod
#print axioms not_ledger_of_ge
#print axioms sub_one_mul_pow_le_of_first_jump
#print axioms contractSeq
#print axioms contractSeq_nonneg
#print axioms contractSeq_pos
#print axioms contractSeq_le_one
#print axioms contractSeq_isRate
#print axioms norm_natCast_le_pow_of_chain
#print axioms mul_le_of_norm_natCast_le_pow

end ABC3.Found.PGC
