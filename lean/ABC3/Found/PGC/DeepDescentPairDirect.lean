import ABC3.Found.PGC.DeepDescentRepair

/-!
# [pGC] `AxDeepDescentPair` を「直接」—— ★★**対の予算は 1 ミリも弱くない**(三者は同値)

配られた持ち場は

> `PairDescent.AxDeepDescentPair K` を `AxWildDescent` を経由せずに直接証明する。
> これが来れば `AxWildDescentMulti → AxLemma → AxSenTate` が閉じる。

であった。★★**着手して最初に検算した結果を先に書く。**

## ★★★測定 0 —— 配られた字面は**真**だが、**弱くなっていない**(本ファイルの主結果)

`AxDeepDescentPair K` は Ax–Sen–Tate(古典の定理)から出るので**偽ではない**。
★★しかし「`AxWildDescent` を経由せずに直接」という前提が外れている。本ファイルは

  ★★★`AxDeepDescentPair K` ⟺ `AxWildDescentMulti K (axDecay p)` ⟺ `AxLemmaGraded K`

を `sorry` 0 で示す(`PairDirect.pair_iff_multi` / `pair_iff_graded` / `multi_iff_graded`)。
ここで `AxLemmaGraded K` は「wild 深さ `k` の `x` に対し
`d(x, K) ≤ (∏_{i=1}^{k} axDecay p i)·ε`」、すなわち★**定数を `axConstant p` に潰す前の
Ax の補題そのもの**である。

⇒ ★★**「対の予算に直せば別の道が開ける」という見立ては誤りである。**
`AxDeepDescentPair` を直接証明することは、`AxWildDescentMulti K (axDecay p)` を
証明することと**同じ 1 点**であり、それは(段別に目盛った)Ax の補題そのものである。
★危ういと報告された `AxWildDescent K (axDecay p)`(段ごとに一様な 1 段の予算)は
この 3 つより**真に強い**(`axWildDescentMulti_of_axWildDescent`)ので、
★**「危うい仮説を避ける」という目的自体はこの同値で達成されている**
(避けた先が、残っている 1 点そのものだった、というだけである)。

## ★★★測定 1 —— 反証候補の `x` は対の予算に**収まる**(★厳密整数演算、49,000 件)

`WildDescentDistanceOnly` 測定 3 の反証候補
(`K = ℚ₃(ζ₃)`、`M = ℚ₃(ζ₂₇)`、`e_M = 18`、`H = Gal(M/K) ≅ ℤ/9`、深さ `k = 2`)を
`v_M` の目盛りで測り直した(スクリプト `.../scratchpad/pair/pairsearch.py`。
`ℤ[π]/(g)` 上の厳密な整数演算で、p 進の丸めを一切使っていない):

| 予算 | 着地 | `v_M` での予算 | 測った最悪の必要量 | slack |
|---|---|---|---|---|
| 1 段 `axDecay 3 2` | 深さ `≤ 1` | `3` | ★`4` | ★`−1`(届かない) |
| ★対 `axDecay 3 1·axDecay 3 2` | 深さ `0`(= `K`) | ★`12` | `8` | ★`+4` |

★★**`AxWildDescent` の slack が `−1` になる `x` すべてで、対の予算の slack は `+4` である。**
乱択 49,000 件(`B = 2,3,4`、係数を `mod 3^B`)で ★**対の予算の最小 slack は `+4`**、
反例ゼロ。★これを Lean に落としたのが `Zeta27Pair.*` で、
`measured_fits_pair_budget`(`J = 8` は収まる)と
`pair_budget_sharp`(`J = 13` は収まらない ⇒ 予算はちょうど `12`)。

★他の底でも同じスクリプトを回した(合計 ★**95,298 件**、`.../scratchpad/pair/run2.py`。☆ヒストグラムの度数を足して数えた):

| `M` | 底 `K` | `e_M` | 件数 | 1 段の最小 slack | ★対の最小 slack |
|---|---|---|---|---|---|
| `ℚ₃(ζ₂₇)` | ` ℚ₃(ζ₃)` | 18 | 49,000 | ★`−1` | ★`+4` |
| `ℚ₂(ζ₈)` | `ℚ₂` | 4 | 19,998 | `0` | ★`0`(等号) |
| `ℚ₂(ζ₁₆)` | `ℚ₂` / `ℚ₂(ζ₄)` | 8 | 12,000 | `1` / `0` | `+6` / `+5` |
| `ℚ₂(ζ₃₂)` | `ℚ₂` / `ℚ₂(ζ₄)` | 16 | 2,400 | `1` / `0` | `+22` / `+13` |
| `ℚ₃(ζ₉)` | `ℚ₃(ζ₃)` | 6 | 6,000 | `1` | `+1` |
| `ℚ₃(ζ₈₁)` | `ℚ₃(ζ₃)` / `ℚ₃(ζ₉)` | 54 | 400 | `1` | `+35` / `+23` |
| `ℚ₅(ζ₂₅)` / `ℚ₇(ζ₄₉)` | `ℚ₅` / `ℚ₇` | 20 / 42 | 5,500 | `1` | `+1` |

⇒ ★★**対の予算への反例は 1 件も出ていない**。最小は `0`(`p = 2`, `K = ℚ₂`,
`M = ℚ₂(ζ₈)`, `k = 1`, `vΔ = 10`, `v(d(x,K)) = 6`, 予算 `4`)で ★**等号が実現する**。
★これは `axDecay p 1 = p^{1/(p−1)}` が古典の最良定数であることと整合する。

## ★★★★測定 1b —— 層ごとの sharp な限界を telescope しても**届かない**(★次の波へ)

同じ塔 `ℚ₃(ζ₂₇) ⊃ ℚ₃(ζ₉) ⊃ ℚ₃(ζ₃)` で、`CyclicJumpNorm` の sharp な層の限界
(損失 = その層の跣び)をそのまま合成すると:

| 層 | 跣び | `v_M` での損失 |
|---|---|---|
| `M/E₁` | `i = 8`(`v_M`) | `8` |
| `E₁/K` | `i = 2`(`v_{E₁}`, `e_{E₁} = 6`) | `6` |
| 合計 | | ★`14` |

★★対の予算は `12` なので **`14` は入らない**
(`Zeta27Pair.naive_tower_sum_exceeds_pair_budget`)。
★しかし実測の最悪値は `8` である。
⇒ ★★★**2 つの層の損失は同時に極値を取れない**。
`AxLemmaGraded` を閉じるには、層ごとの限界の単純な telescope では不足で、
**層間の打ち消し**(上の層で `Δ` が小さいと下の層で損失が小さい)を使う議論が要る。
★これが本波で測った中で**次の波にとって最も価値のある負の情報**である。

## ★★★測定 2 —— 対の予算の**閉じた形**(★ℕ と ℝ だけ。本ファイルの主結果 2)

`PairLedger.prod_Icc_axDecay_eq`:

  ★`∏_{i ∈ [k'+1, k'+d]} axDecay p i = p^{ (p^d − 1)·p / ((p−1)²·p^{k'+d}) }`

`PairLedger.rpow_div_le_prod_Icc_axDecay_iff`(★台帳):

  ★`p^{J/E} ≤ ∏_{i ∈ [k'+1, k'+d]} axDecay p i` ⟺ `(p−1)²·p^{k'+d}·J ≤ (p^d − 1)·p·E`

★`d = 1` を入れると `(p−1)·p^{k'}·J ≤ E` に**簡約される**(`pair_ledger_one_step`)ので、
これは `NormalizedTraceDescent.JumpArith.rpow_div_le_axDecay` /
`WildDescentDistanceOnly.rpow_div_le_axDecay_iff` の 1 段の台帳の**真の一般化**である。
★`k' = 0`, `d → ∞` で指数は `p/(p−1)²`(= `axConstant p` の指数)に収束する。

## 本ファイルが出したもの(すべて `sorry` 0、`#print axioms` は下の Audit)

| 宣言 | 内容 |
|---|---|
| ★★`ProdDescendUntil.exists_of_prod_descent_until` | 抽象核(積の予算で `deg ≤ b` まで落とす) |
| ★`PairDirect.AxLemmaGraded` | 段別に目盛った Ax の補題 |
| ★★★`PairDirect.pair_iff_multi` / `pair_iff_graded` / `multi_iff_graded` | ★**三者同値** |
| `PairDirect.axSenTate_of_axLemmaGraded` | graded ⇒ `AxSenTate` |
| ★★`PairLedger.prod_Icc_axDecay_eq` | 対の予算の閉じた形(★ℕ と ℝ だけ) |
| ★★`PairLedger.rpow_div_le_prod_Icc_axDecay_iff` | ★対の台帳(1 段の台帳の一般化) |
| `PairLedger.pair_ledger_one_step` | `d = 1` で 1 段の台帳に戻る |
| ★★`Zeta27Pair.*` | 反証候補が対の予算に収まる(測定を数値で固定) |
| ★★★`Zeta27Pair.naive_tower_sum_exceeds_pair_budget` | 層ごとの telescope では届かない |

### ★抽象核の要点(はまりどころ)

`exists_of_prod_descent_until` の結論を
「`deg x' ≤ b` かつ `‖x − x'‖ ≤ (∏_{[b+1, deg x]} c)·ε`」と書くと**偽**である。
`b = 1`、`deg x = 2`、途中で `deg x' = 0` まで落ちる場合、支払った予算は
`∏_{[1,2]}` で `∏_{[2,2]}` を超える。★正しくは着地の `deg x'` で書く:
`‖x − x'‖ ≤ (∏_{[deg x' + 1, deg x]} c)·ε`。
さらに帰納を回すには結論に ★`deg x' ≤ deg x` を**足しておく**必要がある
(でないと再帰の戻り値が `deg x₁` より上に居る場合を潰せない)。

## ★★★埋まっていないもの(★正確に)

★**`AxDeepDescentPair K` は証明できていない。** 本ファイルが示したのは
「それは `AxWildDescentMulti K (axDecay p)` と同値であり、迂回になっていない」ことである。
残っている 1 点は依然として 1 点で、その最も素直な形は

  ★`PairDirect.AxLemmaGraded K` —— wild 深さ `k` の `x` に対し
  `d(x, K) ≤ (∏_{i=1}^{k} p^{(1/(p−1))p^{1−i}})·ε`

である。★これは Ax(1970) の補題の段別版であって、`axConstant p` に潰す前の形である。

★★**なぜ「直接証明」が迂回にならないか**(1 行): 対の予算で許される着地は
深さ `≤ 1` だが、★**予算が最大になるのは深さ `0`(= 底 `K`)へ落ちる場合**であり、
そのとき予算は `∏_{i=1}^{k} axDecay p i` ちょうど、すなわち Ax の補題の定数そのものになる。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. `AxLemmaGraded` は原典に対応する文が無い**本ファイルの読み替え**である
   (原典は定数を `axConstant p` に潰した形しか書かない)。消費側は
   `axSenTate_of_axLemmaGraded` 1 本で、定数は 1 ミリも変わらない。
2. 測定の `v_M` は Lean の外(整数多項式の計算)で確かめた。Lean 側に在るのは
   その数値から出る**実数の不等式**だけで、`ℚ₃(ζ₂₇)` 自体は構成していない。
   再現スクリプトは `.../scratchpad/pair/pairsearch.py`(`.../scratchpad/awd/cyc.py`,
   `cyc2.py` を import する)。
3. `PairDescent.pair_of_axWildDescent`(危うい仮説からの道)は**消していない**。
   本ファイルはその隣に、より弱い仮説からの道を立てただけである。
-/

namespace ABC3.Found.PGC
open ABC3.Skeleton.PGC

/-! ## §1 抽象核 —— 「積の予算で `deg ≤ b` まで落とす」 -/

namespace ProdDescendUntil

/-- ★★★**抽象核**(分岐・付値・Galois が 1 語も出ない)。

1 段で `deg` が `d` から `d'` へ落ちるとき損失が `∏_{k ∈ [d'+1, d]} c k` までなら、
`deg ≤ b` になるまで回した全体の損失も `∏_{k ∈ [deg x' + 1, deg x]} c k` に収まる。

★`DeepDescentRepair.DescendUntil.exists_of_descent_until`(1 段の予算 `c (deg x)`)と
`WildDescentMultiStep.ProdCore.exists_mem_of_descent_pair`(`b = 0`、着地は集合 `S`)の
**両方の一般化**である。★着地の予算を `∏_{[b+1, deg x]}` ではなく
`∏_{[deg x' + 1, deg x]}` と書くのが要点で、前者にすると
「`b` より下へ落ちすぎた段」で予算を超える(`b = 1`, `deg x = 2`, `deg x' = 0`)。 -/
theorem exists_of_prod_descent_until {M : Type*} [SeminormedAddCommGroup M]
    [IsUltrametricDist M] (deg : M → ℕ) (P : ℝ → M → Prop) (c : ℕ → ℝ) (b : ℕ)
    (hc : ∀ k, 1 ≤ c k)
    (hstep : ∀ (ε : ℝ) (x : M), 0 ≤ ε → P ε x → b < deg x →
      ∃ x', deg x' < deg x ∧
        ‖x - x'‖ ≤ (∏ k ∈ Finset.Icc (deg x' + 1) (deg x), c k) * ε ∧
        P ((∏ k ∈ Finset.Icc (deg x' + 1) (deg x), c k) * ε) x') :
    ∀ (x : M) (ε : ℝ), 0 ≤ ε → P ε x →
      ∃ x', deg x' ≤ b ∧ deg x' ≤ deg x ∧
        ‖x - x'‖ ≤ (∏ k ∈ Finset.Icc (deg x' + 1) (deg x), c k) * ε := by
  have hprod1 : ∀ s : Finset ℕ, (1:ℝ) ≤ ∏ k ∈ s, c k :=
    fun s => Finset.one_le_prod (fun i _ => hc i)
  intro x
  generalize hn : deg x = n
  induction n using Nat.strong_induction_on generalizing x with
  | _ n ih =>
    intro ε hε hP
    subst hn
    rcases Nat.lt_or_ge b (deg x) with hlt | hge
    · obtain ⟨x₁, h1, hd, hP1⟩ := hstep ε x hε hP hlt
      set A : ℝ := ∏ k ∈ Finset.Icc (deg x₁ + 1) (deg x), c k with hA
      have hA0 : (0:ℝ) ≤ A := le_trans zero_le_one (hprod1 _)
      obtain ⟨x', hb', hle, hy⟩ := ih (deg x₁) h1 x₁ rfl (A * ε) (mul_nonneg hA0 hε) hP1
      refine ⟨x', hb', le_trans hle (le_of_lt h1), ?_⟩
      have hsplit : (∏ k ∈ Finset.Icc (deg x' + 1) (deg x₁), c k) * A
          = ∏ k ∈ Finset.Icc (deg x' + 1) (deg x), c k :=
        DescendUntil.prod_Icc_split c hle (le_of_lt h1)
      have hAle : A ≤ ∏ k ∈ Finset.Icc (deg x' + 1) (deg x), c k := by
        rw [← hsplit]
        exact le_mul_of_one_le_left hA0 (hprod1 _)
      have hmax := IsUltrametricDist.norm_add_le_max (x - x₁) (x₁ - x')
      have heq : x - x₁ + (x₁ - x') = x - x' := by abel
      rw [heq] at hmax
      refine hmax.trans (max_le ?_ ?_)
      · exact hd.trans (mul_le_mul_of_nonneg_right hAle hε)
      · have hring : (∏ k ∈ Finset.Icc (deg x' + 1) (deg x₁), c k) * (A * ε)
            = ((∏ k ∈ Finset.Icc (deg x' + 1) (deg x₁), c k) * A) * ε := by ring
        rw [hring, hsplit] at hy
        exact hy
    · exact ⟨x, hge, le_refl _, by
        simp only [sub_self, norm_zero]
        exact mul_nonneg (le_trans zero_le_one (hprod1 _)) hε⟩

end ProdDescendUntil

/-! ## §2 具体層 —— `AxDeepDescentPair` / `AxWildDescentMulti` / graded な Ax の補題は**同値** -/

namespace PairDirect

variable {p : ℕ} [Fact p.Prime]

/-- ★★**段で目盛った Ax の補題**。定数を `axConstant p` に潰す前の形で、
wild 深さ `k` の `x` については `∏_{i=1}^{k} axDecay p i`(★`< axConstant p`)で足りる、と言う。 -/
def AxLemmaGraded (K : PAdicLocalField p) : Prop :=
  ∀ ε : ℝ, 0 ≤ ε → ∀ x : K.closure, (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) →
    ∃ a : K.carrier, ‖x - algebraMap K.carrier K.closure a‖
      ≤ (∏ i ∈ Finset.Icc 1 (wildDepth K x), axDecay p i) * ε

/-- ★★★**多段の 1 歩 ⇒ 対の予算**(抽象核 §1 に代入するだけ)。 -/
theorem pair_of_axWildDescentMulti (K : PAdicLocalField p)
    (h : AxWildDescentMulti K (axDecay p)) : PairDescent.AxDeepDescentPair K := by
  have hiso : ∀ (g : K.absGal) (y : K.closure), ‖g • y‖ = ‖y‖ := norm_smul_closure K
  intro ε hε x hx hdvd
  have hne : wildDepth K x ≠ 0 := fun h0 => (wildDepth_eq_zero_iff K x).mp h0 hdvd
  rcases Nat.lt_or_ge (wildDepth K x) 2 with hk | hk
  case inl =>
    obtain ⟨x', hlt, hd⟩ := h ε hε x hx hdvd
    exact ⟨x', hlt, by omega, hd⟩
  case inr =>
    obtain ⟨x', hb', _, hd⟩ := ProdDescendUntil.exists_of_prod_descent_until (wildDepth K)
      (fun e z => ∀ σ : K.absGal, ‖σ • z - z‖ ≤ e) (axDecay p) 1 (one_le_axDecay p)
      (fun e z he hz hdeg => by
        obtain ⟨z', hlt, hdd⟩ := h e he z hz
          (by by_contra hc; exact absurd ((wildDepth_eq_zero_iff K z).mpr hc) (by omega))
        have hA1 : (1:ℝ) ≤ ∏ k ∈ Finset.Icc (wildDepth K z' + 1) (wildDepth K z), axDecay p k :=
          Finset.one_le_prod (fun i _ => one_le_axDecay p i)
        exact ⟨z', hlt, hdd, fun σ => UltraCore.norm_smul_sub_self_le_of_norm_sub hiso
          (by nlinarith) hz hdd σ⟩)
      x ε hε hx
    exact ⟨x', by omega, hb', hd⟩

/-- ★★**多段の 1 歩 ⇒ 段で目盛った Ax の補題**(`ProdCore` に代入するだけ)。 -/
theorem axLemmaGraded_of_multi (K : PAdicLocalField p)
    (h : AxWildDescentMulti K (axDecay p)) : AxLemmaGraded K := by
  classical
  have hiso : ∀ (g : K.absGal) (y : K.closure), ‖g • y‖ = ‖y‖ := norm_smul_closure K
  have key : ∀ (z : K.closure) (e : ℝ), 0 ≤ e → (∀ σ : K.absGal, ‖σ • z - z‖ ≤ e) →
      ∃ y ∈ Set.range (algebraMap K.carrier K.closure),
        ‖z - y‖ ≤ (∏ k ∈ Finset.Icc 1 (wildDepth K z), axDecay p k) * e := by
    refine ProdCore.exists_mem_of_descent_pair (wildDepth K)
      (fun e z => ∀ σ : K.absGal, ‖σ • z - z‖ ≤ e) (axDecay p) (one_le_axDecay p) ?_ ?_
    · intro e z he hz hz0
      obtain ⟨y, hy⟩ := exists_norm_sub_algebraMap_le_of_natDegree_tame K he
        ((wildDepth_eq_zero_iff K z).mp hz0) hz
      exact ⟨algebraMap K.carrier K.closure y, ⟨y, rfl⟩, hy⟩
    · intro e z he hz hz0
      obtain ⟨z', hlt, hd⟩ := h e he z hz
        (by by_contra hdvd; exact hz0 ((wildDepth_eq_zero_iff K z).mpr hdvd))
      have hA1 : (1:ℝ) ≤ ∏ k ∈ Finset.Icc (wildDepth K z' + 1) (wildDepth K z), axDecay p k :=
        Finset.one_le_prod (fun i _ => one_le_axDecay p i)
      exact ⟨z', hlt, hd, fun σ => UltraCore.norm_smul_sub_self_le_of_norm_sub hiso
        (by nlinarith) hz hd σ⟩
  intro ε hε x hx
  obtain ⟨y, ⟨a, ha⟩, hy⟩ := key x ε hε hx
  refine ⟨a, ?_⟩
  rw [ha]
  exact hy

/-- ★★**段で目盛った Ax の補題 ⇒ 対の予算**(着地を `K` に取るだけ。★予算は最大の `k' = 0`)。 -/
theorem pair_of_axLemmaGraded (K : PAdicLocalField p) (h : AxLemmaGraded K) :
    PairDescent.AxDeepDescentPair K := by
  intro ε hε x hx hdvd
  obtain ⟨a, ha⟩ := h ε hε x hx
  have hne : wildDepth K x ≠ 0 := fun h0 => (wildDepth_eq_zero_iff K x).mp h0 hdvd
  refine ⟨algebraMap K.carrier K.closure a, ?_, ?_, ?_⟩
  · rw [wildDepth_algebraMap]; omega
  · rw [wildDepth_algebraMap]; omega
  · rw [wildDepth_algebraMap]; simpa using ha

/-- ★★★**三者は同値**(本ファイルの主結果 1)。

`AxDeepDescentPair` ⟺ `AxWildDescentMulti K (axDecay p)` ⟺ `AxLemmaGraded`。 -/
theorem pair_iff_multi (K : PAdicLocalField p) :
    PairDescent.AxDeepDescentPair K ↔ AxWildDescentMulti K (axDecay p) :=
  ⟨PairDescent.axWildDescentMulti_of_pair K, pair_of_axWildDescentMulti K⟩

theorem pair_iff_graded (K : PAdicLocalField p) :
    PairDescent.AxDeepDescentPair K ↔ AxLemmaGraded K :=
  ⟨fun hp => axLemmaGraded_of_multi K (PairDescent.axWildDescentMulti_of_pair K hp),
   pair_of_axLemmaGraded K⟩

theorem multi_iff_graded (K : PAdicLocalField p) :
    AxWildDescentMulti K (axDecay p) ↔ AxLemmaGraded K :=
  ⟨axLemmaGraded_of_multi K, fun hg => (pair_iff_multi K).mp (pair_of_axLemmaGraded K hg)⟩

/-- ★`AxLemmaGraded` からも `AxSenTate` が出る(同値なので当然だが、経路を 1 本にしておく)。 -/
theorem axSenTate_of_axLemmaGraded (K : PAdicLocalField p) (h : AxLemmaGraded K) : AxSenTate K :=
  PairDescent.axSenTate_of_pair K (pair_of_axLemmaGraded K h)

end PairDirect

/-! ## §3 抽象核 —— **対の予算の閉じた形**(★ℕ と ℝ だけ) -/

namespace PairLedger

variable {p : ℕ} [Fact p.Prime]

/-- ★★★**対の予算を閉じた形にする**(本ファイルの主結果 2)。

  `∏_{i ∈ [k'+1, k'+d]} axDecay p i = p^{ (p^d − 1)·p / ((p−1)²·p^{k'+d}) }`

★`d = 1` なら `p^{1/((p−1)·p^{k'})} = axDecay p (k'+1)`(1 段の台帳に戻る)。
★`k' = 0`, `d → ∞` なら `p^{p/(p−1)²} = axConstant p`(Ax の定数に収束する)。 -/
theorem prod_Icc_axDecay_eq (p : ℕ) [Fact p.Prime] (k' d : ℕ) :
    ∏ i ∈ Finset.Icc (k' + 1) (k' + d), axDecay p i
      = (p : ℝ) ^ ((((p : ℝ) ^ d - 1) * (p : ℝ)) / (((p : ℝ) - 1) ^ 2 * (p : ℝ) ^ (k' + d))) := by
  have h1 : (1:ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  have hp0 : (p : ℝ) ≠ 0 := by positivity
  have hp1 : ((p : ℝ) - 1) ≠ 0 := ne_of_gt (by linarith)
  induction d with
  | zero =>
    rw [Nat.add_zero, Finset.Icc_eq_empty (by omega), Finset.prod_empty]
    norm_num
  | succ d ih =>
    have hpow : (0:ℝ) < (p : ℝ) ^ (k' + d) := by positivity
    rw [show k' + (d + 1) = (k' + d) + 1 from by omega,
      Finset.prod_Icc_succ_top (by omega), ih, axDecay,
      show (k' + d + 1) - 1 = k' + d from by omega,
      ← Real.rpow_add (by positivity), one_div_pow]
    congr 1
    field_simp
    ring

/-- ★★★**対の予算の台帳**(本ファイルの主結果 3)。

  `p^{J/E} ≤ ∏_{i ∈ [k'+1, k'+d]} axDecay p i` ⟺ `(p−1)²·p^{k'+d}·J ≤ (p^d − 1)·p·E`

★`d = 1` のとき右辺は `(p−1)·p^{k'}·J ≤ E` に**簡約される**ので、
`NormalizedTraceDescent.JumpArith.rpow_div_le_axDecay` /
`WildDescentDistanceOnly.rpow_div_le_axDecay_iff` の 1 段の台帳の**一般化**である。 -/
theorem rpow_div_le_prod_Icc_axDecay_iff {k' d J E : ℕ} (hE : 0 < E) :
    (p : ℝ) ^ ((J : ℝ) / (E : ℝ)) ≤ ∏ i ∈ Finset.Icc (k' + 1) (k' + d), axDecay p i
      ↔ ((p - 1) ^ 2 * p ^ (k' + d)) * J ≤ ((p ^ d - 1) * p) * E := by
  have hp2 : 2 ≤ p := (Fact.out : p.Prime).two_le
  have h1 : (1:ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  have hE0 : (0:ℝ) < (E : ℝ) := by exact_mod_cast hE
  have hd0 : (0:ℝ) < ((p : ℝ) - 1) ^ 2 * (p : ℝ) ^ (k' + d) := by
    have : (0:ℝ) < ((p:ℝ) - 1) := by linarith
    positivity
  rw [prod_Icc_axDecay_eq, Real.rpow_le_rpow_left_iff h1, div_le_div_iff₀ hE0 hd0]
  have hcast : (((p - 1) ^ 2 * p ^ (k' + d)) * J : ℕ) = (J : ℝ) * (((p:ℝ) - 1) ^ 2 * (p:ℝ) ^ (k' + d)) := by
    push_cast [Nat.cast_sub (by omega : 1 ≤ p)]
    ring
  have hcast2 : (((p ^ d - 1) * p) * E : ℕ) = (((p:ℝ) ^ d - 1) * (p:ℝ)) * (E : ℝ) := by
    push_cast [Nat.cast_sub (Nat.one_le_pow d p (by omega))]
    ring
  constructor
  · intro h
    have : ((((p - 1) ^ 2 * p ^ (k' + d)) * J : ℕ) : ℝ) ≤ ((((p ^ d - 1) * p) * E : ℕ) : ℝ) := by
      rw [hcast, hcast2]; exact h
    exact_mod_cast this
  · intro h
    have h' : ((((p - 1) ^ 2 * p ^ (k' + d)) * J : ℕ) : ℝ) ≤ ((((p ^ d - 1) * p) * E : ℕ) : ℝ) := by
      exact_mod_cast h
    rw [hcast, hcast2] at h'
    exact h'

/-- ★`d = 1` のとき対の台帳は 1 段の台帳に戻る(★同値の確認)。 -/
theorem pair_ledger_one_step {k' J E : ℕ} :
    ((p - 1) ^ 2 * p ^ (k' + 1)) * J ≤ ((p ^ 1 - 1) * p) * E ↔ ((p - 1) * p ^ k') * J ≤ E := by
  have hp2 : 2 ≤ p := (Fact.out : p.Prime).two_le
  have hc : 0 < (p - 1) * p := Nat.mul_pos (by omega) (by omega)
  have hL : ((p - 1) ^ 2 * p ^ (k' + 1)) * J = ((p - 1) * p) * (((p - 1) * p ^ k') * J) := by
    ring
  have hR : ((p ^ 1 - 1) * p) * E = ((p - 1) * p) * E := by rw [pow_one]
  constructor
  · intro h; rw [hL, hR] at h; exact Nat.le_of_mul_le_mul_left h hc
  · intro h; rw [hL, hR]; exact Nat.mul_le_mul (le_refl _) h

end PairLedger

/-! ## §4 数値層 —— 反証候補の `x` は**対の予算に収まる**(★台帳から機械的に出る) -/

namespace Zeta27Pair

/-- `K = ℚ₃(ζ₃)`、`M = ℚ₃(ζ₂₇)`(`e_M = 18`)、深さ `k = 2` の反証候補 `x`
(`WildDescentDistanceOnly` 測定 3)の**対の予算はちょうど `3^{12/18}`** である。
★台帳 `(p−1)²p^{k}·J ≤ (p^k−1)·p·E` に `p=3, k'=0, d=2, E=18` を入れると
`36·J ≤ 432` すなわち ★`J ≤ 12`。 -/
theorem pair_budget_eq_twelve :
    ∏ i ∈ Finset.Icc 1 2, axDecay 3 i = (3:ℝ) ^ ((12:ℝ)/18) := by
  haveI : Fact (Nat.Prime 3) := ⟨Nat.prime_three⟩
  have h := PairLedger.prod_Icc_axDecay_eq 3 0 2
  norm_num at h
  rw [h]
  norm_num

/-- ★★**反証候補の `x` は対の予算に収まる**(`J = 8 ≤ 12`)。

測定値: `min_{σ≠1} v_M(σx − x) = 11`、`v_M(d(x,K)) = 3` ⇒ 損失 `3^{8/18}`。 -/
theorem measured_fits_pair_budget :
    (3:ℝ) ^ ((8:ℝ)/18) ≤ ∏ i ∈ Finset.Icc 1 2, axDecay 3 i := by
  haveI : Fact (Nat.Prime 3) := ⟨Nat.prime_three⟩
  have h := (PairLedger.rpow_div_le_prod_Icc_axDecay_iff (p := 3) (k' := 0) (d := 2)
    (J := 8) (E := 18) (by norm_num)).mpr (by norm_num)
  simpa using h

/-- ★**同じ `x` は 1 段の予算には収まらない**(`k' = 1`, `d = 1` ⇒ `J ≤ 3` だが `J = 4`)。
これが `AxWildDescent K (axDecay p)` への反証候補の中身である。 -/
theorem measured_exceeds_one_step_budget :
    ¬ ((3:ℝ) ^ ((4:ℝ)/18) ≤ ∏ i ∈ Finset.Icc 2 2, axDecay 3 i) := by
  haveI : Fact (Nat.Prime 3) := ⟨Nat.prime_three⟩
  intro hle
  have h := (PairLedger.rpow_div_le_prod_Icc_axDecay_iff (p := 3) (k' := 1) (d := 1)
    (J := 4) (E := 18) (by norm_num)).mp (by simpa using hle)
  norm_num at h

/-- ★**予算はここで尽きる** —— `J = 13` はもう収まらない(`12` が sharp)。 -/
theorem pair_budget_sharp :
    ¬ ((3:ℝ) ^ ((13:ℝ)/18) ≤ ∏ i ∈ Finset.Icc 1 2, axDecay 3 i) := by
  haveI : Fact (Nat.Prime 3) := ⟨Nat.prime_three⟩
  intro hle
  have h := (PairLedger.rpow_div_le_prod_Icc_axDecay_iff (p := 3) (k' := 0) (d := 2)
    (J := 13) (E := 18) (by norm_num)).mp (by simpa using hle)
  norm_num at h


/-- ★★★**層ごとの sharp な限界を telescope しても対の予算には入らない**(★次の波への警告)。

`M = ℚ₃(ζ₂₇) ⊃ E₁ = ℚ₃(ζ₉) ⊃ K = ℚ₃(ζ₃)` の 2 層を
`CyclicJumpNorm` の sharp な層の限界(損失 = その層の跣び)で telescope すると

* 上の層 `M/E₁`: 跣び `i = 8`(`v_M`) ⇒ `3^{8/18}`、
* 下の層 `E₁/K`: 跣び `i = 2`(`v_{E₁}`、`e_{E₁} = 6`) ⇒ `3^{2/6} = 3^{6/18}`

で、合計は ★`3^{14/18}`。★★しかし対の予算は `3^{12/18}` しか無いので**足りない**。
★★★実測の最悪値は `3^{8/18}`(slack `+4`)なので、
**2 層の損失は同時に極値を取れない**。
⇒ ★`AxLemmaGraded` を閉じるには、層ごとの限界の単純な合成ではなく
**層間の打ち消しを使う議論**が要る。 -/
theorem naive_tower_sum_exceeds_pair_budget :
    ¬ ((3:ℝ) ^ ((14:ℝ)/18) ≤ ∏ i ∈ Finset.Icc 1 2, axDecay 3 i) := by
  haveI : Fact (Nat.Prime 3) := ⟨Nat.prime_three⟩
  intro hle
  have h := (PairLedger.rpow_div_le_prod_Icc_axDecay_iff (p := 3) (k' := 0) (d := 2)
    (J := 14) (E := 18) (by norm_num)).mp (by simpa using hle)
  norm_num at h

end Zeta27Pair

/-! ## §5 `.src`(原典の対応箇所) -/

namespace ProdDescendUntil
def exists_of_prod_descent_until.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
end ProdDescendUntil

namespace PairDirect
def AxLemmaGraded.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def pair_of_axWildDescentMulti.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def axLemmaGraded_of_multi.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def pair_of_axLemmaGraded.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def pair_iff_multi.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def pair_iff_graded.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def multi_iff_graded.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def axSenTate_of_axLemmaGraded.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
end PairDirect

namespace PairLedger
def prod_Icc_axDecay_eq.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def rpow_div_le_prod_Icc_axDecay_iff.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def pair_ledger_one_step.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
end PairLedger

namespace Zeta27Pair
def pair_budget_eq_twelve.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def measured_fits_pair_budget.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def measured_exceeds_one_step_budget.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def pair_budget_sharp.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
def naive_tower_sum_exceeds_pair_budget.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }
end Zeta27Pair

end ABC3.Found.PGC

section Audit
open ABC3.Found.PGC
#print axioms ProdDescendUntil.exists_of_prod_descent_until
#print axioms PairDirect.pair_of_axWildDescentMulti
#print axioms PairDirect.axLemmaGraded_of_multi
#print axioms PairDirect.pair_of_axLemmaGraded
#print axioms PairDirect.pair_iff_multi
#print axioms PairDirect.pair_iff_graded
#print axioms PairDirect.multi_iff_graded
#print axioms PairDirect.axSenTate_of_axLemmaGraded
#print axioms PairLedger.prod_Icc_axDecay_eq
#print axioms PairLedger.rpow_div_le_prod_Icc_axDecay_iff
#print axioms PairLedger.pair_ledger_one_step
#print axioms Zeta27Pair.pair_budget_eq_twelve
#print axioms Zeta27Pair.measured_fits_pair_budget
#print axioms Zeta27Pair.measured_exceeds_one_step_budget
#print axioms Zeta27Pair.pair_budget_sharp
#print axioms Zeta27Pair.naive_tower_sum_exceeds_pair_budget
end Audit
