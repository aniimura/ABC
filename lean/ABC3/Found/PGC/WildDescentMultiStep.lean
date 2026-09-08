import ABC3.Found.PGC.NormalizedTraceDescent

/-!
# [pGC] `AxWildDescent K (axDecay p)` の最後の 1 点 —— ★字面は「1 段の予算」では閉じない

配られた持ち場は

> wild 深さ `k` の `x` に対し、深さ `< k` の `x'` で
> `‖x − x'‖ ≤ p^{j/(m·e)}·ε` かつ `∀σ ‖σx' − x'‖ ≤ 同じ` を満たすものを作る
> (`j` は下の層 `E₁/F` の跳び、`m = e(L/E₁)`、`e = e_{E₁}`)

であった。★★**着手して最初に検算した結果を先に書く。**

## ★★★測定 1 —— 3 つの条件のうち **2 つ目はただで出る**(抽象核 1)

`∀σ ‖σx' − x'‖ ≤ δ` は、`‖x − x'‖ ≤ δ` と `ε ≤ δ` と**超距離性 + 等長性**だけから出る:

```
σx' − x' = σ(x' − x) + ((σx − x) + (x − x'))
‖·‖ ≤ max(‖x − x'‖, max(ε, ‖x − x'‖)) ≤ δ
```

⇒ ★持ち場が「ここが落ちやすい」と名指しした点は**落ちない**。
`UltraCore.norm_smul_sub_self_le_of_norm_sub`(分岐・付値・Galois が 1 語も出ない)。
★したがって残っているのは**距離の 1 本だけ**である。

## ★★★★測定 2 —— `x'` を `E₁` に取る形は **偽**(★反例を数値で確定させた)

持ち場の質問「`x'` は `E₁` に入るのか」に答える。★**入るが、損失は `p^{j/(m·e)}` では効かない。**

反例 `p = 3`、`F = ℚ₃(ζ₃)`、`L = ℚ₃(ζ₈₁)`(`[L:F] = 27` 巡回、`e_L = 54`)、
`E₁ = ℚ₃(ζ₉)`(`[E₁:F] = 3`、`m = e(L/E₁) = 9 = p^{k−1}`、`e = e_{E₁} = 6`、`j = i₁ = 2`)。
`y := π_L³`(`π_L = ζ₈₁ − 1`)は wild 深さ `k = 3` で、★機械で検算した値は

* `min_{σ≠1} v_L(σ y − y) = 9`(`Gal(L/F)` の 26 個すべてで確認) ⇒ `Δ(y) = ‖π_L‖⁹`、
* `v_L(y) = 3`、`v_L(E₁^×) = 9ℤ` ⇒ `d(y, E₁) = ‖π_L‖³`。

配られた字面は `d(y,E₁) ≤ p^{j/(m·e)}·Δ = ‖π_L‖^{−2}·‖π_L‖⁹ = ‖π_L‖⁷` を要求するが、
★★`‖π_L‖³ > ‖π_L‖⁷` なので**偽**(`Zeta81.firstJump_bound_false`、指数で `4` 外す)。

★**なぜ間違っていたか(その 1、`p` 冪の補正)**: `σπ_L/π_L − 1` の付値は `i₁` だが、
`v_L(x) = D` の元では `(1+c)^D − 1` を見るので、`p ∣ D` のとき付値が
★`p^{v_p(D)}` 倍に跳ね上がる★。上の例では `D = 3` で `3·i₁ = 6`、
測った `d/Δ = ‖π_L‖^{3−9}` に**ちょうど一致する**。

## ★★★★測定 2b —— `k = 2` でも偽。★**打ち消し**という第 2 の機構

★前の波が「収まる」と報告した層 `ℚ₃(ζ₂₇)/ℚ₃(ζ₃)`(`k = 2`)自体が反例である。
前の波は `x = π_L` しか見ていない。一般の `x = π_L + b`(`b ∈ E₁`)では
コバウンダリ `σb − b`(`v_L ∈ e(L/E₁)·ℤ`)が `σπ_L − π_L` の先頭項を**打ち消す**。
★総当たりで測った最大は `min_{σ≠1} v_L(σx − x) = 4`(`b = 0` なら `3`)で、
損失は `‖π_L‖^{1−4} = 3^{1/6}`。配られた字面 `p^{j/(m·e)} = 3^{1/9}` を**超える**
(`Zeta27.firstJump_bound_false`)。★`k = 2` では `p` 冪の補正は使えない(`t = 1`)ので、
★★これは測定 2 とは**別の機構**である。

## ★★★★測定 2c —— 2 系列とも上限は **`axDecay p 2 = p^{1/(p(p−1))}` ちょうど**

同じ総当たりを `ℚ₃(ζ₂₇)`(`k = 2`)と `ℚ₃(ζ₈₁)`(`k = 3`、`D = 1..8`)で回すと、
`d(x,E₁)/Δ(x)` の最大値は

* `ℚ₃(ζ₂₇)`: `‖π_L‖^{-3} = 3^{3/18} = 3^{1/6}`、
* `ℚ₃(ζ₈₁)`: `‖π_L‖^{-9} = 3^{9/54} = 3^{1/6}`(`D = 3` で実現)

で ★★**どちらも `axDecay 3 2 = 3^{1/6}` ちょうど**である
(`Zeta27.cancel_cost_eq_axDecay_two` / `zeta81_cancel_cost_eq_axDecay_two`)。
★`e_L` にも `k` にも依らない。⇒ ★★**1 段の損失は `axDecay p 2` で止まる(そして下げられない)**
というのが測定の結論であり、配られた `p^{j/(m·e)}` は**両方向に外れている**
(`k = 1` の層では正しいが、`k ≥ 2` の層では小さすぎる)。

## ★★★測定 3 —— `axDecay p 2` は **1 段の予算に収まらない**。多段の台帳が要る

`axDecay p 2` は `k` に依らないので `k ≥ 3` では `axDecay p k` を超える
(`Zeta81.axDecay_three_three_lt_digitCost`: `axDecay 3 3 = 3^{1/18} < 3^{1/9}`)。
★★しかしこの降下は深さを `k` から **`≤ 1`** へ落とす(着地が `E₁`、`[E₁:F] = p`)ので、
予算は `∏_{i=2}^{k} axDecay p i` であり、★`axDecay p 2` 単独ですでに足りている。

⇒ ★★★**台帳を「1 段」から「多段」に替えれば閉じる。**
`AxTowerDecay.exists_mem_of_descent_budget` は損失を `c (deg x)` としか書けないので、
本ファイルは**対を取る**核 `ProdCore.exists_mem_of_descent_pair`(損失
`∏_{k ∈ [deg x' + 1, deg x]} c k`)を新しく立てた。予算関数は同じ `∏_{[1,deg x]} c k`
なので ★**`AxLemma` / `AxSenTate` の定数は 1 ミリも変わらない**。

## 本ファイルが出したもの(すべて `sorry` 0)

| 宣言 | 内容 |
|---|---|
| ★`UltraCore.norm_smul_sub_self_le_of_norm_sub` | 抽象核 1(`ε` の保存はただ) |
| ★★`ProdCore.exists_mem_of_descent_pair` | 抽象核 2(**多段**降下の予算) |
| `ProdCore.prod_Icc_mul_prod_Icc_succ` | `∏_{[1,d']}·∏_{[d'+1,d]} = ∏_{[1,d]}` |
| ★`AxWildDescentMulti` | 多段の降下(★条件は**距離 1 本だけ**) |
| `axWildDescentMulti_of_axWildDescent` | 1 段 ⇒ 多段(★新仮説は**真に弱い**) |
| ★★`axLemma_of_wildDescentMulti` / `axLemma_of_axDecayMulti` | `AxLemma K (axConstant p)` |
| ★★`axSenTate_of_axDecayMulti` | `AxSenTate K` |
| ★★`JumpArithMulti.rpow_mul_div_le_axDecay_two` | 補正つきの指数の勘定(★ℕ と ℝ だけ) |
| `JumpArithMulti.axDecay_two_le_prod_Icc` | 着地が深さ `≤ 1` なら多段予算に収まる |
| ★★★`DeepDescent.AxDeepDescent` | ★**残っている 1 点の素直な形** |
| ★★★`DeepDescent.axSenTate_of_deepDescent` | ★**その 1 点から `AxSenTate`** |
| `FirstJumpDigit.FirstJumpDigitRoute` / `axDeepDescent_of_firstJumpDigit` | 分岐の言葉で書いた証明書 |
| ★★`Zeta81.*` / `Zeta27.*` | 反例の数値(2 系列、2 つの機構) |

## ★★残っているのはちょうど 1 点(★以前と同じ 1 点だが**形が違う**)

  ★`DeepDescent.AxDeepDescent K`

すなわち「wild な `x` に対し、★**深さ `≤ 1`** の `x'` を
`‖x − x'‖ ≤ axDecay p (min 2 k)·ε` で取ること」(`k = wildDepth K x`)。

* `ε` の保存はもう要らない(測定 1)。
* `k ≥ 2` の損失は `k` に依らない `axDecay p 2 = p^{1/(p(p−1))}` でよい(測定 2c)。
  ★これは配られた `p^{j/(m·e)}` より**大きく**、しかも測定では**等号が実現している**。
* 分岐の言葉で書きたければ `FirstJumpDigit.FirstJumpDigitRoute`
  (`p^{t·j/(m·e)}`、`p^{k−1} ≤ m`、`(p−1)j ≤ e`、`1 ≤ t ≤ p^{k−2}`)を使う。
  ★ただし `(j,m,e,t)` は**存在量化された証明書**であって、
  ★層の分岐データそのものを入れると測定 2b で壊れる。

★★★**本ファイルは `AxDeepDescent` を証明していない。** 埋まっていない。
上の表は「その 1 点さえ来れば `AxSenTate` が出る」ことを `sorry` 0 で示しただけである。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. 原典 (Ax 1970 / pGC Cor. 3.1) は 1 段の勘定を地の文で畳む。「多段の台帳」は
   原典に対応する文が無い**本ファイルの読み替え**である。★消費側(`AxLemma` の定数)は
   変わらないので、CLAUDE.md「逸脱」の条件(後続に影響しない)を満たす。
   `.src` は `AxTowerDecay` / `AxEpsilonDecay` / `NormalizedTraceDescent` と同じ項目を指す。
2. `NormalizedTraceDescent.FirstJumpRoute.axLemma_of_firstJump` は**そのまま正しい**
   (仮説が満たせないだけである)。★本ファイルはそれを消さず、`t` を入れた版を**別に**立てた。
3. 反例の `v_L` は Lean の外(整数多項式の計算)で確かめた。★Lean 側に在るのは
   その数値から出る**実数の不等式**だけで、`ℚ₃(ζ₈₁)` 自体は構成していない
   (§6 の表は「測ったが形式化していない」部分である)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC

/-! ## §1 抽象核 1 —— 「`ε` の保存」は距離から**ただで**出る -/

namespace UltraCore

/-- ★★★**抽象核 1** —— 超距離空間に等長に作用する群では、
`x` の変位が `ε` 以下で `‖x − x'‖ ≤ δ`(`ε ≤ δ`)なら、`x'` の変位も `δ` 以下。 -/
theorem norm_smul_sub_self_le_of_norm_sub {M : Type*} [SeminormedAddCommGroup M]
    [IsUltrametricDist M] {G : Type*} [Monoid G] [DistribMulAction G M]
    (hiso : ∀ (g : G) (y : M), ‖g • y‖ = ‖y‖) {x x' : M} {ε δ : ℝ} (hεδ : ε ≤ δ)
    (hx : ∀ g : G, ‖g • x - x‖ ≤ ε) (hd : ‖x - x'‖ ≤ δ) (g : G) :
    ‖g • x' - x'‖ ≤ δ := by
  have h1 : g • x' - x' = g • (x' - x) + ((g • x - x) + (x - x')) := by
    rw [smul_sub]; abel
  have h2 : ‖g • (x' - x)‖ = ‖x - x'‖ := by rw [hiso, norm_sub_rev]
  rw [h1]
  refine (IsUltrametricDist.norm_add_le_max _ _).trans (max_le ?_ ?_)
  · rw [h2]; exact hd
  · exact (IsUltrametricDist.norm_add_le_max _ _).trans (max_le ((hx g).trans hεδ) hd)

end UltraCore

/-! ## §2 抽象核 2 —— **多段**降下の予算 -/

namespace ProdCore

/-- `∏_{[1,d']} · ∏_{[d'+1,d]} = ∏_{[1,d]}`。 -/
theorem prod_Icc_mul_prod_Icc_succ {M : Type*} [CommMonoid M] (c : ℕ → M) {d' d : ℕ}
    (h : d' ≤ d) :
    (∏ k ∈ Finset.Icc 1 d', c k) * (∏ k ∈ Finset.Icc (d' + 1) d, c k)
      = ∏ k ∈ Finset.Icc 1 d, c k := by
  have h1 : Finset.Icc 1 d' = Finset.Ioc 0 d' := by
    ext k; simp only [Finset.mem_Icc, Finset.mem_Ioc]; omega
  have h2 : Finset.Icc (d' + 1) d = Finset.Ioc d' d := by
    ext k; simp only [Finset.mem_Icc, Finset.mem_Ioc]; omega
  have h3 : Finset.Icc 1 d = Finset.Ioc 0 d := by
    ext k; simp only [Finset.mem_Icc, Finset.mem_Ioc]; omega
  rw [h1, h2, h3]
  exact Finset.prod_Ioc_consecutive c (Nat.zero_le d') h

/-- ★★★**抽象核 2 —— 多段降下の予算**。1 段で `deg` が `d` から `d'` へ落ちるとき、
損失が `∏_{k ∈ [d'+1, d]} c k` までなら、全体の損失は `∏_{k ∈ [1, deg x]} c k`。 -/
theorem exists_mem_of_descent_pair {M : Type*} [SeminormedAddCommGroup M]
    [IsUltrametricDist M] {S : Set M} (deg : M → ℕ) (P : ℝ → M → Prop) (c : ℕ → ℝ)
    (hc : ∀ k, 1 ≤ c k)
    (hbase : ∀ (ε : ℝ) (x : M), 0 ≤ ε → P ε x → deg x = 0 → ∃ y ∈ S, ‖x - y‖ ≤ ε)
    (hstep : ∀ (ε : ℝ) (x : M), 0 ≤ ε → P ε x → deg x ≠ 0 →
      ∃ x', deg x' < deg x ∧
        ‖x - x'‖ ≤ (∏ k ∈ Finset.Icc (deg x' + 1) (deg x), c k) * ε ∧
        P ((∏ k ∈ Finset.Icc (deg x' + 1) (deg x), c k) * ε) x') :
    ∀ (x : M) (ε : ℝ), 0 ≤ ε → P ε x →
      ∃ y ∈ S, ‖x - y‖ ≤ (∏ k ∈ Finset.Icc 1 (deg x), c k) * ε := by
  have hc0 : ∀ k, (0:ℝ) ≤ c k := fun k => le_trans zero_le_one (hc k)
  have hprod1 : ∀ (s : Finset ℕ), (1:ℝ) ≤ ∏ k ∈ s, c k :=
    fun s => Finset.one_le_prod (fun i _ => hc i)
  intro x
  generalize hn : deg x = n
  induction n using Nat.strong_induction_on generalizing x with
  | _ n ih =>
    intro ε hε hP
    rcases Nat.eq_zero_or_pos n with h0 | hpos
    · obtain ⟨y, hyS, hy⟩ := hbase ε x hε hP (hn.trans h0)
      exact ⟨y, hyS, hy.trans (by nlinarith [hprod1 (Finset.Icc 1 n)])⟩
    · obtain ⟨x', hlt, hd, hP'⟩ := hstep ε x hε hP (by omega)
      set A : ℝ := ∏ k ∈ Finset.Icc (deg x' + 1) (deg x), c k with hA
      have hA0 : 0 ≤ A := Finset.prod_nonneg (fun i _ => hc0 i)
      obtain ⟨y, hyS, hy⟩ := ih (deg x') (by omega) x' rfl (A * ε) (mul_nonneg hA0 hε) hP'
      refine ⟨y, hyS, ?_⟩
      have hsplit : (∏ k ∈ Finset.Icc 1 (deg x'), c k) * A = ∏ k ∈ Finset.Icc 1 (deg x), c k :=
        prod_Icc_mul_prod_Icc_succ c (le_of_lt hlt)
      have h1 : ‖x - x'‖ ≤ (∏ k ∈ Finset.Icc 1 (deg x), c k) * ε := by
        refine hd.trans (mul_le_mul_of_nonneg_right ?_ hε)
        rw [← hsplit]
        nlinarith [hprod1 (Finset.Icc 1 (deg x'))]
      have h2 : ‖x' - y‖ ≤ (∏ k ∈ Finset.Icc 1 (deg x), c k) * ε := by
        refine hy.trans ?_
        rw [← hsplit]
        apply le_of_eq
        ring
      have hmax := IsUltrametricDist.norm_add_le_max (x - x') (x' - y)
      have heq : x - x' + (x' - y) = x - y := by abel
      rw [heq] at hmax
      rw [hn] at h1 h2
      exact hmax.trans (max_le h1 h2)

end ProdCore

/-! ## §3 具体層 —— 多段の `AxWildDescent` -/

variable {p : ℕ} [Fact p.Prime]

/-- ★★★**多段の wild 降下**(本ファイルが原典の勘定から切り出した形)。

`AxTowerDecay.AxWildDescent` との違いは 2 つ:

* 1 段の損失に許される予算が `c (wildDepth K x)` ではなく
  ★**`∏_{k ∈ [wildDepth K x' + 1, wildDepth K x]} c k`**(飛ばした段の**積**)である。
* `ε` の伸びの条件を**書かない**(`UltraCore` の抽象核 1 でただで出るため)。 -/
def AxWildDescentMulti (K : PAdicLocalField p) (c : ℕ → ℝ) : Prop :=
  ∀ ε : ℝ, 0 ≤ ε → ∀ x : K.closure, (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) →
    p ∣ (minpoly K.carrier x).natDegree →
      ∃ x' : K.closure, wildDepth K x' < wildDepth K x ∧
        ‖x - x'‖ ≤ (∏ k ∈ Finset.Icc (wildDepth K x' + 1) (wildDepth K x), c k) * ε

def AxWildDescentMulti.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★**1 段の降下は多段の降下である**(`c ≥ 1` なら)。
⇒ `AxWildDescentMulti` は `AxWildDescent` より★真に弱い仮説★である。 -/
theorem axWildDescentMulti_of_axWildDescent (K : PAdicLocalField p) {c : ℕ → ℝ}
    (hc : ∀ k, 1 ≤ c k) (h : AxWildDescent K c) : AxWildDescentMulti K c := by
  intro ε hε x hx hdvd
  obtain ⟨x', hlt, hd, _⟩ := h ε hε x hx hdvd
  refine ⟨x', hlt, hd.trans (mul_le_mul_of_nonneg_right ?_ hε)⟩
  have hmem : wildDepth K x ∈ Finset.Icc (wildDepth K x' + 1) (wildDepth K x) := by
    simp only [Finset.mem_Icc]; omega
  have hsplit := Finset.mul_prod_erase
    (Finset.Icc (wildDepth K x' + 1) (wildDepth K x)) c hmem
  have h1 : (1:ℝ) ≤ ∏ k ∈ (Finset.Icc (wildDepth K x' + 1) (wildDepth K x)).erase
      (wildDepth K x), c k := Finset.one_le_prod (fun i _ => hc i)
  rw [← hsplit]
  nlinarith [hc (wildDepth K x)]

def axWildDescentMulti_of_axWildDescent.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**多段の降下から Ax の補題**(本ファイルの主結果 1)。 -/
theorem axLemma_of_wildDescentMulti (K : PAdicLocalField p) {c : ℕ → ℝ} {C : ℝ}
    (hc : ∀ k, 1 ≤ c k) (hC : ∀ n : ℕ, ∏ k ∈ Finset.Icc 1 n, c k ≤ C)
    (h : AxWildDescentMulti K c) : AxLemma K C := by
  classical
  have hiso : ∀ (g : K.absGal) (y : K.closure), ‖g • y‖ = ‖y‖ := norm_smul_closure K
  have key : ∀ (z : K.closure) (e : ℝ), 0 ≤ e → (∀ σ : K.absGal, ‖σ • z - z‖ ≤ e) →
      ∃ y ∈ Set.range (algebraMap K.carrier K.closure),
        ‖z - y‖ ≤ (∏ k ∈ Finset.Icc 1 (wildDepth K z), c k) * e := by
    refine ProdCore.exists_mem_of_descent_pair (wildDepth K)
      (fun e z => ∀ σ : K.absGal, ‖σ • z - z‖ ≤ e) c hc ?_ ?_
    · intro e z he hz hz0
      obtain ⟨y, hy⟩ := exists_norm_sub_algebraMap_le_of_natDegree_tame K he
        ((wildDepth_eq_zero_iff K z).mp hz0) hz
      exact ⟨algebraMap K.carrier K.closure y, ⟨y, rfl⟩, hy⟩
    · intro e z he hz hz0
      obtain ⟨z', hlt, hd⟩ := h e he z hz
        (by
          by_contra hdvd
          exact hz0 ((wildDepth_eq_zero_iff K z).mpr hdvd))
      have hA1 : (1:ℝ) ≤ ∏ k ∈ Finset.Icc (wildDepth K z' + 1) (wildDepth K z), c k :=
        Finset.one_le_prod (fun i _ => hc i)
      refine ⟨z', hlt, hd, ?_⟩
      exact fun σ => UltraCore.norm_smul_sub_self_le_of_norm_sub hiso
        (by nlinarith) hz hd σ
  intro ε hε x hx
  obtain ⟨y, ⟨a, ha⟩, hy⟩ := key x ε hε hx
  refine ⟨a, ?_⟩
  rw [ha]
  exact hy.trans (mul_le_mul_of_nonneg_right (hC _) hε)

def axLemma_of_wildDescentMulti.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**多段 + `axDecay` ⇒ `AxLemma K (axConstant p)`**(本ファイルの主結果 2)。 -/
theorem axLemma_of_axDecayMulti (K : PAdicLocalField p) {c : ℕ → ℝ}
    (hc : ∀ k, 1 ≤ c k) (hdecay : ∀ k, c k ≤ axDecay p k)
    (h : AxWildDescentMulti K c) : AxLemma K (axConstant p) := by
  have h1 : (1:ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  have hp0 : (0:ℝ) < (p : ℝ) := by linarith
  refine axLemma_of_wildDescentMulti K hc (fun n => ?_) h
  have hbound := prod_le_rpow_of_geometric_shift (D := (p : ℝ)) (A := 1 / ((p : ℝ) - 1))
    (r := 1 / (p : ℝ)) (le_of_lt h1) (by positivity) (by positivity)
    (by rw [div_lt_one hp0]; linarith)
    (fun k => le_trans zero_le_one (hc k)) hdecay n
  rwa [axExponent_eq h1] at hbound

def axLemma_of_axDecayMulti.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**そこから Ax–Sen–Tate**(本ファイルの主結果 3)。 -/
theorem axSenTate_of_axDecayMulti (K : PAdicLocalField p) {c : ℕ → ℝ}
    (hc : ∀ k, 1 ≤ c k) (hdecay : ∀ k, c k ≤ axDecay p k)
    (h : AxWildDescentMulti K c) : AxSenTate K :=
  axSenTate_of_axLemma (le_trans zero_le_one (one_le_axConstant p))
    (axLemma_of_axDecayMulti K hc hdecay h)

def axSenTate_of_axDecayMulti.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §4 抽象核 3 —— 桁の補正 `t` を入れた指数の勘定(★ℕ と ℝ だけ) -/

namespace JumpArithMulti

variable {p : ℕ} [Fact p.Prime]

/-- ★★★**補正つきの閉じる十分条件**(本ファイルの主結果 4)。

`NormalizedTraceDescent.JumpArith.rpow_div_le_axDecay` は損失の指数を `j/(m·e)` と
していたが、★実際の第 1 跳びの降下は分子に **`p` 冪の補正 `t`** が付く(§B)。
`t ≤ p^{k−2}` である限り、損失は `k` に依らず **`axDecay p 2`** で押さえられる:

* `t ≤ p^{k−2}`、`p^{k−1} ≤ m`、`(p−1)j ≤ e` ⇒ `t·j/(m·e) ≤ 1/((p−1)p)`。

★★`t = 1`(補正なし)なら `rpow_div_le_axDecay` の `axDecay p k` に戻る。
★★これが「多段の台帳」を要求する理由である: `axDecay p 2` は
`axDecay p k`(`k ≥ 3`)より**大きい**ので、1 段ぶんの予算には収まらない。 -/
theorem rpow_mul_div_le_axDecay_two {k j m e t : ℕ}
    (hm : p ^ (k - 1) ≤ m) (he : 0 < e) (hjump : (p - 1) * j ≤ e)
    (hk : 2 ≤ k) (ht : t ≤ p ^ (k - 2)) :
    (p : ℝ) ^ (((t : ℝ) * (j : ℝ)) / ((m : ℝ) * (e : ℝ))) ≤ axDecay p 2 := by
  have hp2 : 2 ≤ p := (Fact.out : p.Prime).two_le
  have h1 : (1:ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  have hp0 : (0:ℝ) < (p : ℝ) := by linarith
  have hjumpR : ((p : ℝ) - 1) * (j : ℝ) ≤ (e : ℝ) := by
    have := (Nat.cast_le (α := ℝ)).mpr hjump
    rwa [Nat.cast_mul, Nat.cast_sub (by omega), Nat.cast_one] at this
  have hmR : (p : ℝ) ^ (k - 1) ≤ (m : ℝ) := by
    have := (Nat.cast_le (α := ℝ)).mpr hm
    rwa [Nat.cast_pow] at this
  have htR : (t : ℝ) ≤ (p : ℝ) ^ (k - 2) := by
    have := (Nat.cast_le (α := ℝ)).mpr ht
    rwa [Nat.cast_pow] at this
  have hpk : (0:ℝ) < (p : ℝ) ^ (k - 1) := by positivity
  have hm0 : (0:ℝ) < (m : ℝ) := lt_of_lt_of_le hpk hmR
  have he0 : (0:ℝ) < (e : ℝ) := by exact_mod_cast he
  have hsucc : (p : ℝ) ^ (k - 1) = (p : ℝ) ^ (k - 2) * (p : ℝ) := by
    rw [show k - 1 = (k - 2) + 1 by omega, pow_succ]
  rw [axDecay]
  refine Real.rpow_le_rpow_of_exponent_le (le_of_lt h1) ?_
  have hrhs : (1 / ((p:ℝ) - 1)) * (1 / (p:ℝ)) ^ (2 - 1) = 1 / (((p:ℝ) - 1) * (p:ℝ)) := by
    norm_num
    ring
  rw [hrhs, div_le_div_iff₀ (by positivity) (by positivity)]
  have hj0 : (0:ℝ) ≤ ((p:ℝ) - 1) * (j : ℝ) := by positivity
  calc ((t:ℝ) * (j:ℝ)) * (((p:ℝ) - 1) * (p:ℝ))
      = ((t:ℝ) * (((p:ℝ) - 1) * (j:ℝ))) * (p:ℝ) := by ring
    _ ≤ ((p:ℝ) ^ (k - 2) * (e:ℝ)) * (p:ℝ) := by
        refine mul_le_mul_of_nonneg_right ?_ (le_of_lt hp0)
        exact mul_le_mul htR hjumpR hj0 (by positivity)
    _ = ((p:ℝ) ^ (k - 1)) * (e:ℝ) := by rw [hsucc]; ring
    _ ≤ (m:ℝ) * (e:ℝ) := mul_le_mul_of_nonneg_right hmR (le_of_lt he0)
    _ = 1 * ((m:ℝ) * (e:ℝ)) := (one_mul _).symm

def rpow_mul_div_le_axDecay_two.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★**着地が深さ `≤ 1` なら `axDecay p 2` は多段の予算に収まる**。
(`k' + 1 ≤ 2 ≤ k` すなわち `2 ∈ [k'+1, k]` だけを使う。) -/
theorem axDecay_two_le_prod_Icc {k k' : ℕ} (hk : 2 ≤ k) (hk' : k' ≤ 1) :
    axDecay p 2 ≤ ∏ i ∈ Finset.Icc (k' + 1) k, axDecay p i := by
  have hmem : 2 ∈ Finset.Icc (k' + 1) k := by
    simp only [Finset.mem_Icc]; omega
  have hsplit := Finset.mul_prod_erase (Finset.Icc (k' + 1) k) (axDecay p) hmem
  have h1 : (1:ℝ) ≤ ∏ i ∈ (Finset.Icc (k' + 1) k).erase 2, axDecay p i :=
    Finset.one_le_prod (fun i _ => one_le_axDecay p i)
  rw [← hsplit]
  nlinarith [one_le_axDecay p 2]

def axDecay_two_le_prod_Icc.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★`k = 1` の段は補正が効かない(`t = 1`)ので `axDecay p 1` のまま。 -/
theorem rpow_div_le_prod_Icc_one {j m e : ℕ}
    (hm : p ^ (1 - 1) ≤ m) (he : 0 < e) (hjump : (p - 1) * j ≤ e) :
    (p : ℝ) ^ ((j : ℝ) / ((m : ℝ) * (e : ℝ))) ≤ ∏ i ∈ Finset.Icc 1 1, axDecay p i := by
  simpa using JumpArith.rpow_div_le_axDecay (p := p) (k := 1) hm he hjump

def rpow_div_le_prod_Icc_one.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end JumpArithMulti

/-! ## §5 具体層 —— 「第 1 跳び + 桁の補正」の道 -/

namespace FirstJumpDigit

variable {p : ℕ} [Fact p.Prime]

/-- ★**分岐の言葉で書いた証明書**(本ファイルの主結果 5)。

wild な `x` について、★**深さ `≤ 1` の** `x'` と ★**存在量化された** `(j, m, e, t)` が
`p^{k−1} ≤ m`、`(p−1)j ≤ e`、`1 ≤ t ≤ p^{k−2}`、`‖x − x'‖ ≤ p^{t·j/(m·e)}·ε`
を満たすように取れること。

★★**読み方に注意**(測定 2b):`(j, m, e, t)` に**層の分岐データそのもの**
(`j` = `E₁/F` の跳び、`m = e(L/E₁)`、`e = e_{E₁}`、`t = p^{v_p(v_L(x))}`)を入れた形は
★**偽である**(§6.2 の `ℚ₃(ζ₂₇)`、`k = 2` の打ち消し)。
★★本定義が主張しているのは実質「損失が `axDecay p 2` 以下」だけであり
(`(j,m,e,t)` はその証明書)、素直な形は `DeepDescent.AxDeepDescent`(§5.5)である。
★`NormalizedTraceDescent` の `FirstJumpRoute`(`t = 1` かつ分岐データ固定)は
★**満たせない**(§6.1 と §6.2 の 2 つの反例)。 -/
def FirstJumpDigitRoute (K : PAdicLocalField p) : Prop :=
  ∀ ε : ℝ, 0 ≤ ε → ∀ x : K.closure, (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) →
    p ∣ (minpoly K.carrier x).natDegree →
      ∃ (x' : K.closure) (j m e t : ℕ),
        wildDepth K x' < wildDepth K x ∧ wildDepth K x' ≤ 1 ∧
        p ^ (wildDepth K x - 1) ≤ m ∧ 0 < e ∧ (p - 1) * j ≤ e ∧
        1 ≤ t ∧ t ≤ p ^ (wildDepth K x - 2) ∧
        ‖x - x'‖ ≤ (p : ℝ) ^ (((t : ℝ) * (j : ℝ)) / ((m : ℝ) * (e : ℝ))) * ε

def FirstJumpDigitRoute.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**その道は多段の降下を与える**。 -/
theorem axWildDescentMulti_of_firstJumpDigit (K : PAdicLocalField p)
    (h : FirstJumpDigitRoute K) : AxWildDescentMulti K (axDecay p) := by
  intro ε hε x hx hdvd
  obtain ⟨x', j, m, e, t, hlt, hle1, hm, he, hjump, ht1, ht, hcost⟩ := h ε hε x hx hdvd
  refine ⟨x', hlt, hcost.trans (mul_le_mul_of_nonneg_right ?_ hε)⟩
  have hne : wildDepth K x ≠ 0 := fun h0 => (wildDepth_eq_zero_iff K x).mp h0 hdvd
  rcases Nat.lt_or_ge (wildDepth K x) 2 with hk | hk
  case inr =>
    exact (JumpArithMulti.rpow_mul_div_le_axDecay_two hm he hjump hk ht).trans
      (JumpArithMulti.axDecay_two_le_prod_Icc hk hle1)
  case inl =>
    have hk1 : wildDepth K x = 1 := by omega
    have hk0 : wildDepth K x' = 0 := by omega
    have htone : t = 1 := by
      rw [hk1] at ht
      simp only [show (1:ℕ) - 2 = 0 from rfl, pow_zero] at ht
      omega
    rw [hk1] at hm
    rw [hk0, hk1, htone]
    simpa using JumpArithMulti.rpow_div_le_prod_Icc_one (p := p) hm he hjump


def axWildDescentMulti_of_firstJumpDigit.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**Ax の補題**。 -/
theorem axLemma_of_firstJumpDigit (K : PAdicLocalField p) (h : FirstJumpDigitRoute K) :
    AxLemma K (axConstant p) :=
  axLemma_of_axDecayMulti K (one_le_axDecay p) (fun _ => le_refl _)
    (axWildDescentMulti_of_firstJumpDigit K h)

def axLemma_of_firstJumpDigit.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**Ax–Sen–Tate**。 -/
theorem axSenTate_of_firstJumpDigit (K : PAdicLocalField p) (h : FirstJumpDigitRoute K) :
    AxSenTate K :=
  axSenTate_of_axDecayMulti K (one_le_axDecay p) (fun _ => le_refl _)
    (axWildDescentMulti_of_firstJumpDigit K h)

def axSenTate_of_firstJumpDigit.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end FirstJumpDigit

/-! ## §5.5 ★★残っている 1 点の**素直な形**(測定が支持しているのはこちら) -/

namespace DeepDescent

variable {p : ℕ} [Fact p.Prime]

/-- ★★★**深い降下**(本ファイルが勧める「残っている 1 点」の形)。

wild な `x` に対し、★**深さ `≤ 1` に一気に落ちる** `x'` が
`‖x − x'‖ ≤ axDecay p (min 2 k) · ε`(`k = wildDepth K x`)で取れること。

★★`k ≥ 2` の損失が **`k` に依らない `axDecay p 2 = p^{1/(p(p−1))}`** である点が要。
★§6 で測った 2 系列(`ℚ₃(ζ₂₇)`・`ℚ₃(ζ₈₁)`)では**この値ちょうどが実現する**ので、
★これ以上は絞れない。★1 段の台帳(`axDecay p k`)では `k ≥ 3` で足りず、
★多段の台帳(本ファイル §2)で**ちょうど**足りる。 -/
def AxDeepDescent (K : PAdicLocalField p) : Prop :=
  ∀ ε : ℝ, 0 ≤ ε → ∀ x : K.closure, (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) →
    p ∣ (minpoly K.carrier x).natDegree →
      ∃ x' : K.closure, wildDepth K x' < wildDepth K x ∧ wildDepth K x' ≤ 1 ∧
        ‖x - x'‖ ≤ axDecay p (min 2 (wildDepth K x)) * ε

def AxDeepDescent.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**深い降下 ⇒ 多段の降下**。 -/
theorem axWildDescentMulti_of_deepDescent (K : PAdicLocalField p) (h : AxDeepDescent K) :
    AxWildDescentMulti K (axDecay p) := by
  intro ε hε x hx hdvd
  obtain ⟨x', hlt, hle1, hcost⟩ := h ε hε x hx hdvd
  have hne : wildDepth K x ≠ 0 := fun h0 => (wildDepth_eq_zero_iff K x).mp h0 hdvd
  refine ⟨x', hlt, hcost.trans (mul_le_mul_of_nonneg_right ?_ hε)⟩
  rcases Nat.lt_or_ge (wildDepth K x) 2 with hk | hk
  case inr =>
    rw [show min 2 (wildDepth K x) = 2 by omega]
    exact JumpArithMulti.axDecay_two_le_prod_Icc hk hle1
  case inl =>
    have hk1 : wildDepth K x = 1 := by omega
    have hk0 : wildDepth K x' = 0 := by omega
    rw [hk0, hk1]
    simp

def axWildDescentMulti_of_deepDescent.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**深い降下 ⇒ Ax の補題**。 -/
theorem axLemma_of_deepDescent (K : PAdicLocalField p) (h : AxDeepDescent K) :
    AxLemma K (axConstant p) :=
  axLemma_of_axDecayMulti K (one_le_axDecay p) (fun _ => le_refl _)
    (axWildDescentMulti_of_deepDescent K h)

def axLemma_of_deepDescent.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**深い降下 ⇒ Ax–Sen–Tate**。★これが本ファイルの最終形である。 -/
theorem axSenTate_of_deepDescent (K : PAdicLocalField p) (h : AxDeepDescent K) :
    AxSenTate K :=
  axSenTate_of_axDecayMulti K (one_le_axDecay p) (fun _ => le_refl _)
    (axWildDescentMulti_of_deepDescent K h)

def axSenTate_of_deepDescent.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★`FirstJumpDigit` の証明書は `AxDeepDescent` の十分条件である
(`p^{t·j/(m·e)} ≤ axDecay p 2`、`k = 1` は `axDecay p 1`)。 -/
theorem axDeepDescent_of_firstJumpDigit (K : PAdicLocalField p)
    (h : FirstJumpDigit.FirstJumpDigitRoute K) : AxDeepDescent K := by
  intro ε hε x hx hdvd
  obtain ⟨x', j, m, e, t, hlt, hle1, hm, he, hjump, ht1, ht, hcost⟩ := h ε hε x hx hdvd
  have hne : wildDepth K x ≠ 0 := fun h0 => (wildDepth_eq_zero_iff K x).mp h0 hdvd
  refine ⟨x', hlt, hle1, hcost.trans (mul_le_mul_of_nonneg_right ?_ hε)⟩
  rcases Nat.lt_or_ge (wildDepth K x) 2 with hk | hk
  case inr =>
    rw [show min 2 (wildDepth K x) = 2 by omega]
    exact JumpArithMulti.rpow_mul_div_le_axDecay_two hm he hjump hk ht
  case inl =>
    have hk1 : wildDepth K x = 1 := by omega
    have htone : t = 1 := by
      rw [hk1] at ht
      simp only [show (1:ℕ) - 2 = 0 from rfl, pow_zero] at ht
      omega
    rw [hk1] at hm
    rw [hk1, htone]
    simpa using JumpArith.rpow_div_le_axDecay (p := p) (k := 1) hm he hjump

def axDeepDescent_of_firstJumpDigit.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end DeepDescent

/-! ## §6.1 反例 1(`k = 3`, `ℚ₃(ζ₈₁)/ℚ₃(ζ₃)`)—— ★`p` 冪の補正

★★以下の数値は**機械で検算した**(`v_L` を `ℤ[π]`(`π = ζ₈₁ − 1`)の係数から
`min_i (54·v₃(a_i) + i)` として計算。`Φ₈₁(π+1)` は Eisenstein で `P₀ = 3`)。

| 量 | 値 |
|---|---|
| `e_L = e(ℚ₃(ζ₈₁))` | `54` |
| `min_{σ≠1} v_L(σπ − π)` | `3`  ⇒ 第 1 跳び `i₁ = 2` |
| `y := π³`、`v_L(y)` | `3` |
| `min_{σ≠1} v_L(σy − y)`(σ は 26 個すべて) | ★`9` |
| `d(y, E₁)`(`E₁ = ℚ₃(ζ₉)`、`v_L(E₁^×) = 9ℤ`) | ★`‖π‖³` |
| `v_L(y − π_{ℚ₃(ζ₂₇)})` | `55` ⇒ ★正しい降下は深さ 2 の `π_M` へ |

⇒ `t = 1` の主張「`d(y,E₁) ≤ ‖π_L‖^{−i₁}·Δ(y)`」は `‖π‖³ ≤ ‖π‖⁷` を要求するが**偽**。
★この `y` での補正は `t = p^{v_p(v_L(y))} = 3`(`= p^{k−2}`、`k = 3`)で、
そのとき `p^{t·j/(m·e)} = 3^{6/54} = 3^{1/9}` は `d/Δ = ‖π‖^{3−9}` に★ちょうど一致する★。
★★ただし `b ∈ E₁` を足して打ち消すと同じ層でさらに `‖π‖^{-9} = 3^{1/6}` まで悪くなる
(§6.2 の測定)。★**`p` 冪の補正だけでは足りない** —— 上限は `axDecay p 2` である。 -/

namespace Zeta81

/-- ★★**測った跳びの間隔**(上表): `v_L(y) + i₁ = 3 + 2 = 5` に対し
すべての `σ ≠ 1` で `v_L(σy − y) = 9` 以上。★`t = 1` の道が壊れる幅は `4`。 -/
theorem jump_gap : 3 + 2 < 9 := by norm_num

def jump_gap.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**`t = 1` の第 1 跳びの道は偽**(本ファイルの主結果 6)。

`d(y,E₁) = 3^{−3/54}` は `‖π_L‖^{−i₁}·Δ(y) = 3^{2/54}·3^{−9/54} = 3^{−7/54}` を**超える**。 -/
theorem firstJump_bound_false :
    ¬ ((3:ℝ) ^ (-(3:ℝ)/54) ≤ (3:ℝ) ^ ((2:ℝ)/54) * (3:ℝ) ^ (-(9:ℝ)/54)) := by
  rw [← Real.rpow_add (by norm_num : (0:ℝ) < 3),
    Real.rpow_le_rpow_left_iff (by norm_num : (1:ℝ) < 3)]
  norm_num

def firstJump_bound_false.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

theorem axDecay_three_three_eq : axDecay 3 3 = (3:ℝ) ^ ((1:ℝ)/18) := by
  norm_num [axDecay]

def axDecay_three_three_eq.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★**1 段の予算では足りない**: 実際の損失 `3^{1/9}` は `axDecay 3 3 = 3^{1/18}` を超える。
★これが「多段の台帳」が要る理由の**数値的な証拠**である。 -/
theorem axDecay_three_three_lt_digitCost : axDecay 3 3 < (3:ℝ) ^ ((1:ℝ)/9) := by
  rw [axDecay_three_three_eq]
  exact (Real.rpow_lt_rpow_left_iff (by norm_num)).mpr (by norm_num)

def axDecay_three_three_lt_digitCost.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**多段の予算なら収まる**: `t = 3 = p^{k−2}`、`j = 2`、`m = 9 = p^{k−1}`、`e = 6` を
`rpow_mul_div_le_axDecay_two` に入れると `3^{1/9} ≤ axDecay 3 2`。
★着地は `E₁ = ℚ₃(ζ₉)`(深さ `≤ 1`)なので予算は `axDecay 3 3 · axDecay 3 2` であり、
★`axDecay 3 2` 単独ですでに足りている。 -/
theorem digitCost_le_axDecay_two : (3:ℝ) ^ ((1:ℝ)/9) ≤ axDecay 3 2 := by
  haveI : Fact (Nat.Prime 3) := ⟨Nat.prime_three⟩
  have h := JumpArithMulti.rpow_mul_div_le_axDecay_two (p := 3) (k := 3) (j := 2) (m := 9)
    (e := 6) (t := 3) (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
  have hexp : (((3:ℕ):ℝ) * ((2:ℕ):ℝ)) / (((9:ℕ):ℝ) * ((6:ℕ):ℝ)) = (1:ℝ)/9 := by norm_num
  rw [hexp] at h
  simpa using h

def digitCost_le_axDecay_two.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- `axDecay 3 2 = 3^{1/6}`(`NormalizedTraceDescent` の同名補題の再掲を避けるための別名)。 -/
theorem axDecay_three_two_eq' : axDecay 3 2 = (3:ℝ) ^ ((1:ℝ)/6) := by
  norm_num [axDecay]

def axDecay_three_two_eq'.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end Zeta81

/-! ### §6.2 ★★2 つ目の反例(`k = 2`, `ℚ₃(ζ₂₇)/ℚ₃(ζ₃)`)—— ★**打ち消し**による損失

★★前の波が「収まる」と報告した**まさにその層**で `t = 1` の字面が壊れる。
前の波は `x = π_L`(一様化元)しか見ていない。★一般の `x` は `E₁` 成分を持ち、
その**コバウンダリ** `σb − b`(`b ∈ E₁`、`v_L ∈ (e(L/E₁))ℤ`)が
`σπ_L − π_L` の先頭項を**打ち消す**。

機械で探索した(`b` を `ℤ[ζ₉]/π_{E₁}` の範囲で総当たり、`Gal(L/F)` の 8 元すべてで最小):

| 量 | 値 |
|---|---|
| `x = π_L + b`、`d(x,E₁)` | `‖π_L‖¹` |
| `max_b min_{σ≠1} v_L(σx − x)` | ★`4`(`b = 0` なら `3`) |
| ⇒ 損失 | `‖π_L‖^{1−4} = 3^{3/18} = 3^{1/6}` |
| `p^{j/(m·e)}`(配られた字面、`j=2, m=3, e=6`) | `3^{2/18} = 3^{1/9}` ★**足りない** |
| `axDecay 3 2` | ★`3^{1/6}` ★**ちょうど一致** |

★★`ℚ₃(ζ₈₁)` の `D = 3` でも同じ探索で最大は `‖π_L‖^{-9} = 3^{9/54} = 3^{1/6}`、
★★★**2 系列とも `axDecay p 2` ちょうど**である(`= p^{1/(p(p−1))}`)。
⇒ ★`AxDeepDescent` の定数はこれ以上絞れない。 -/

namespace Zeta27

/-- ★★打ち消しで測った跳び: `d + i₁ = 1 + 2 = 3` に対し実測は `4`。 -/
theorem cancel_jump : 1 + 2 < 4 := by norm_num

def cancel_jump.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**`k = 2` でも `t = 1` の字面は偽**(本ファイルの主結果 7)。

`d(x,E₁) = 3^{−1/18}` は `‖π_L‖^{−i₁}·Δ(x) = 3^{2/18}·3^{−4/18} = 3^{−2/18}` を**超える**。 -/
theorem firstJump_bound_false :
    ¬ ((3:ℝ) ^ (-(1:ℝ)/18) ≤ (3:ℝ) ^ ((2:ℝ)/18) * (3:ℝ) ^ (-(4:ℝ)/18)) := by
  rw [← Real.rpow_add (by norm_num : (0:ℝ) < 3),
    Real.rpow_le_rpow_left_iff (by norm_num : (1:ℝ) < 3)]
  norm_num

def firstJump_bound_false.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★★★**測った損失は `axDecay 3 2` ちょうど**(`3^{3/18} = 3^{1/6}`)。
⇒ `AxDeepDescent` の `axDecay p 2` は**下げられない**(等号が実現している)。 -/
theorem cancel_cost_eq_axDecay_two : (3:ℝ) ^ ((3:ℝ)/18) = axDecay 3 2 := by
  rw [Zeta81.axDecay_three_two_eq']
  norm_num

def cancel_cost_eq_axDecay_two.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-- ★同じ値が `ℚ₃(ζ₈₁)`(`k = 3`)の探索の最大値でもある: `3^{9/54} = 3^{1/6}`。 -/
theorem zeta81_cancel_cost_eq_axDecay_two : (3:ℝ) ^ ((9:ℝ)/54) = axDecay 3 2 := by
  rw [Zeta81.axDecay_three_two_eq']
  norm_num

def zeta81_cancel_cost_eq_axDecay_two.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end Zeta27

end ABC3.Found.PGC

section Audit
open ABC3.Found.PGC
#print axioms UltraCore.norm_smul_sub_self_le_of_norm_sub
#print axioms ProdCore.prod_Icc_mul_prod_Icc_succ
#print axioms ProdCore.exists_mem_of_descent_pair
#print axioms axWildDescentMulti_of_axWildDescent
#print axioms axLemma_of_wildDescentMulti
#print axioms axLemma_of_axDecayMulti
#print axioms axSenTate_of_axDecayMulti
#print axioms JumpArithMulti.rpow_mul_div_le_axDecay_two
#print axioms JumpArithMulti.axDecay_two_le_prod_Icc
#print axioms JumpArithMulti.rpow_div_le_prod_Icc_one
#print axioms FirstJumpDigit.axWildDescentMulti_of_firstJumpDigit
#print axioms FirstJumpDigit.axLemma_of_firstJumpDigit
#print axioms FirstJumpDigit.axSenTate_of_firstJumpDigit
#print axioms DeepDescent.axWildDescentMulti_of_deepDescent
#print axioms DeepDescent.axLemma_of_deepDescent
#print axioms DeepDescent.axSenTate_of_deepDescent
#print axioms DeepDescent.axDeepDescent_of_firstJumpDigit
#print axioms Zeta27.firstJump_bound_false
#print axioms Zeta27.cancel_cost_eq_axDecay_two
#print axioms Zeta27.zeta81_cancel_cost_eq_axDecay_two
#print axioms Zeta81.firstJump_bound_false
#print axioms Zeta81.axDecay_three_three_lt_digitCost
#print axioms Zeta81.digitCost_le_axDecay_two
end Audit
