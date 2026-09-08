import ABC3.Found.PGC.AxEpsilonDecay
import ABC3.Found.PGC.WildDepthFieldDescent

/-!
# [pGC] `AxWildDescent K (axDecay p)` —— 「`ε` の伸び」は無償／台帳／★字面への反証候補

`Found/PGC/AxEpsilonDecay.lean` 以降、pGC §3 の残りは
**`AxWildDescent K (axDecay p)`(古典の 1 段)ただ 1 点**とされてきた。本ファイルは

1. その仮説の **半分(`ε` が伸びない条件)を無償で落とし**、
2. 残った「距離の条件」を **`ℕ` の不等式 1 本**に翻訳し、
3. ★★その `ℕ` の不等式を **実在する分岐データが破ること**を測った、

の 3 つを出す。★★**`AxWildDescent K (axDecay p)` は閉じていない。**

## ★★★測定 1 —— `AxWildDescent K c` の第 3 条件は `1 ≤ c` から自動で出る

`AxWildDescent K c` の結論は 3 つ組だった:

* `wildDepth K x' < wildDepth K x`、
* `‖x − x'‖ ≤ c k · ε`(距離)、
* ★`∀σ ‖σx' − x'‖ ≤ c k · ε`(`ε` の伸び)。

★★**3 つ目は 1・2 から出る。**超距離と「作用が等長」だけで

  `σx' − x' = σ(x'−x) + (σx − x) + (x−x')`  ⇒  `‖σx'−x'‖ ≤ max(‖x−x'‖, ε) ≤ c k · ε`

(`norm_smul_sub_self_le_max` / `axWildDescent_of_dist`)。
★これは `SenLemma.AxDescentStep`(`ε` を**保つ**ことを要求する)との決定的な違いで、
あちらは `≤ ε` を要求するので同じ論法では出ない。
⇒ ★★**次の波が作るべきものは「距離の評価」1 本だけ**である。

## ★★測定 2 —— 距離の条件は `ℕ` の不等式 1 本に翻訳できる

`rpow_div_le_axDecay_iff`(★実数と `ℕ` だけ。分岐・付値・Galois が 1 語も出ない):

  `p^{j/E} ≤ axDecay p (m+1)` ⟺ `(p−1)·p^m·j ≤ E`。

古典の 1 段はこう書ける。`x ∈ L`、`L/F` が巡回 `p` 次(降下の層)、`τ` をその生成元、
`e := e_L`、`i :=` `L/F` の跳び、`a := v_L(Δ_K(x))`、`b := v_L(τx − x)` と置くと

* 層の最良定数(`CyclicJumpNorm` / `MinpolyOrbitSplit`)から
  `d(x,F) = ‖π_L‖^{b−i}`、要求は `‖π_L‖^{b−i} ≤ axDecay p k · ‖π_L‖^{a}`、
* すなわち **`p^{(i+a−b)/e} ≤ axDecay p k`**、
* 台帳に通すと ★**`(p−1)·p^{k−1}·(i + a − b) ≤ e`**。

`gain := b − a`(★`τ` が深いほど `‖τx−x‖` は `Δ` より小さい)と書けば
★**`(p−1)·p^{k−1}·(i − gain) ≤ e`** が 1 段の条件である。

★★**素元 `x = π_M` ではこれは成り立つ**(`gain = i − u₁`、`u₁` は最初の跳び、
`(p−1)p^{k−1}u₁ ≤ e` は `WildJumpChain.sub_one_mul_pow_le_of_first_jump`)。

## ★★★★測定 3 —— 一般の `x` では `gain` が足りない(★機械計算。★字面への反証候補)

★★**素元でない `x` では `gain` が `i − u₁` より小さくなりうる。**

厳密な整数演算で全部当たった(スクリプトは報告に添付。`p` 進の丸めを一切使っていない):

```
K = ℚ₃(ζ₃)（e_K = 2）,  M = ℚ₃(ζ₂₇)（e_M = 18）,  H = Gal(M/K) ≅ ℤ/9,  深さ k = 2
下付き分岐群:  G_0 = G, G_1 = G_2 = H, G_3..G_8 = Gal(M/ℚ₃(ζ₉)), G_9 = 1
⇒ 位数 9 の σ は i_G(σ) = 3、位数 3 の τ は i_G(τ) = 9、層 M/ℚ₃(ζ₉) の跳び i = 8
x = 素元 のとき     v(σx−x) = 3,  v(τx−x) = 9   ⇒ gain = 6 ≥ i − p^{n−k} = 8 − 3 = 5 ✓
★実在する x では   v(σx−x) = 11, v(τx−x) = 15  ⇒ gain = 4 < 5 ✗
   （この x の d(x, ℚ₃(ζ₉)) は v_M で 7、要求は 8。★ちょうど p^{1/18} 足りない）
   出現率 37/25000, 9/12000, 11/12000（係数を mod 3², 3³, 3⁴ で乱択）。★同じ数値ばかり出る。
```

★この「破れ」を実数の不等式として形式化したのが
`axDecay_three_two_lt_cyclotomic_layer_loss`:  `axDecay 3 2 = 3^{3/18} < 3^{4/18}`。

★★★**ただしこれは `AxWildDescent K (axDecay p)` が偽であることの証明ではない。**
測ったのは「`x' ∈ M` に限れば届かない」ことだけで、
★`M` の外の深さ ≤ 1 の `x'` を排除できていない。Krasner の補題は
`‖x−x'‖ < min_σ‖σx−x‖` のときしか効かず、ここで要求される半径は
`min` より**大きい**ので届かない(`need = 8 ≤ v_max = 15`)。
⇒ ★**「偽」ではなく「反証候補」と読むこと。**★確定させるには
`K = ℚ₃(ζ₃)` の 40 個の巡回 3 次拡大(Kummer)と順分岐拡大の中で `d(x,N)` を測る必要がある。
★測っていない。

★★**合成(`AxLemma`)は無傷である。**同じ `x` で `v(d(x,K)) = 3` であり、
2 段ぶんの予算 `axDecay 3 2 · axDecay 3 1 = 3^{2/3}`(= `v_M` で 12)には
`3 ≥ 11 − 12 = −1` と大きく余っている。
⇒ ★危ういのは **「段ごとに一様な予算 `c k`」という記法**であって、Ax の定数ではない。
★`AxTowerDecay.exists_mem_of_descent_budget` は `c` と `g` を分けた**予算関数 `F`** を
既に持っているので、次の波はそちらに載せ替えるのが自然である(本ファイルは載せ替えていない)。

## 検査した範囲(★反例は出ていない)

同じ厳密計算で `d(x, M^S) ≤ axDecay p k · Δ(x)`(`x' ∈ M` に限った上界)を検査した:

```
p=3: ℚ₃(ζ₉), ℚ₃(ζ₂₇), ℚ₃(ζ₈₁)      底 = ℚ₃ / ℚ₃(ζ₃) / ℚ₃(ζ₉)
p=2: ℚ₂(ζ₈), ℚ₂(ζ₁₆), ℚ₂(ζ₃₂)      底 = ℚ₂ / ℚ₂(ζ₄)
p=5,7: ℚ₅(ζ₂₅), ℚ₅(ζ₁₂₅), ℚ₇(ζ₄₉)
計 約 12 万件。★底が ℚ_p のときは反例ゼロ（最悪でも slack = 0 = 等号)。
★slack < 0 は「底が ℚ₃(ζ₃)」の場合だけに出た（上の測定 3）。
```

## 在庫の測定(★自分で測った。コマンドを残す)

```
grep -rn "IsUltrametricDist" lean/ABC3/Found/PGC/SenLemma.lean
  → IsUltrametricDist.norm_add_le_max は在る（norm_sub_le_max は無い、`to_additive` の穴）
grep -rn "norm_smul_closure" lean/ABC3/Found/PGC/AbsClosureModules.lean
  → norm_smul_closure (AbsClosureModules.lean:270) `‖σ • x‖ = ‖x‖` ★等長性はここ 1 本で足りる
grep -n "div_le_div_iff" .cache/mathlib-index.txt
  → ★`div_le_div_iff` は無い。`div_le_div_iff₀ (hb : 0 < b) (hd : 0 < d) : a/b ≤ c/d ↔ a*d ≤ c*b`
grep -n -i "ramificationGroup|herbrand|lowerIndex" .cache/mathlib-index.txt
  → ★mathlib には **上付き/下付き分岐群も Herbrand 関数も無い**
     (在るのは ValuationSubring.decompositionSubgroup / inertiaSubgroup と differentIdeal だけ)。
     ⇒ 測定 2 の `i`・`u₁` を体の層で扱うには自前の分岐理論が要る。★これが最大の在庫の穴である。
```

## 逸脱の記録(CLAUDE.md「逸脱」)

1. 原典 (Ax 1970) は 1 段の勘定を地の文で畳む。本ファイルの `rpow_div_le_axDecay_iff` は
   原典に対応する文が無い中間結果である。`.src` は上流と同じ項目(pGC 物理 p.6 Corollary 3.1)。
2. ★測定 3 は**反例ではなく反証候補**である(上に明記)。`M` の外を測っていない。
3. `norm_smul_sub_self_le_max` は `M` に「超距離な半ノルム加法群」、`G` に「モノイドの等長作用」
   しか要求しない。★可逆性も忠実性も要らない ——原典より弱い設定である。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-! ## §1 抽象核 —— 「`ε` の伸び」は距離から出る

★分岐・付値・Galois・`p` 進の語彙が 1 語も出てこない。
`M` は超距離な半ノルム加法群、`G` はその上に**等長に**作用するモノイドでよい
（可逆性も忠実性も `g^p = 1` も要らない）。 -/

section AbstractCore

variable {M : Type*} [SeminormedAddCommGroup M] [IsUltrametricDist M]
variable {G : Type*} [Monoid G] [DistribMulAction G M]

/-- ★★**抽象核**: `x` の軌道が `e` 以内に収まり、`x'` が `x` から `d` 以内なら、
`x'` の軌道は `max d e` 以内に収まる。

`σx' − x' = σ(x'−x) + (σx − x) + (x−x')` を超距離で潰すだけである。
★★これが `AxWildDescent K c` の第 3 条件（`ε` の伸び）を無償にする。 -/
theorem norm_smul_sub_self_le_max (hiso : ∀ (g : G) (z : M), ‖g • z‖ = ‖z‖)
    {x x' : M} {d e : ℝ} (hx : ∀ g : G, ‖g • x - x‖ ≤ e) (hd : ‖x - x'‖ ≤ d) (g : G) :
    ‖g • x' - x'‖ ≤ max d e := by
  have hsplit : g • x' - x' = (g • (x' - x) + (g • x - x)) + (x - x') := by
    rw [smul_sub]; abel
  have h1 : ‖g • (x' - x)‖ ≤ d := by
    rw [hiso, norm_sub_rev]; exact hd
  have h2 : ‖(x : M) - x'‖ ≤ d := hd
  calc ‖g • x' - x'‖ = ‖(g • (x' - x) + (g • x - x)) + (x - x')‖ := by rw [hsplit]
    _ ≤ max ‖g • (x' - x) + (g • x - x)‖ ‖x - x'‖ :=
        IsUltrametricDist.norm_add_le_max _ _
    _ ≤ max (max ‖g • (x' - x)‖ ‖g • x - x‖) ‖x - x'‖ := by
        gcongr
        exact IsUltrametricDist.norm_add_le_max _ _
    _ ≤ max (max d e) d := max_le_max (max_le_max h1 (hx g)) h2
    _ = max d e := by
        rcases le_total d e with h | h <;> simp [h]

end AbstractCore

def norm_smul_sub_self_le_max.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §2 具体層 —— `AxWildDescent K c` は「距離の評価」だけで足りる -/

section Concrete

/-- ★★★**`AxWildDescent` の第 3 条件（`ε` の伸び）は無償である。**

`1 ≤ c k` を仮定すると、`AxWildDescent K c` を作るのに必要なのは

  `p ∣ [K(x):K]` かつ `∀σ ‖σx−x‖ ≤ ε` ⇒ `∃ x'`, `wildDepth x' < wildDepth x`,
  `‖x − x'‖ ≤ c(wildDepth x)·ε`

★**距離の評価ただ 1 本**である。

★入力は `AbsClosureModules.norm_smul_closure`（`‖σ • x‖ = ‖x‖`）だけで、
`K.closure` の超距離性と合わせて §1 の抽象核に代入する。
★`SenLemma.AxDescentStep` は `‖σx'−x'‖ ≤ ε`（`c` 倍を許さない）を要求するので、
★**同じ論法はあちらには効かない**。 -/
theorem axWildDescent_of_dist (K : PAdicLocalField p) {c : ℕ → ℝ} (hc : ∀ k, 1 ≤ c k)
    (h : ∀ ε : ℝ, 0 ≤ ε → ∀ x : K.closure, (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) →
      p ∣ (minpoly K.carrier x).natDegree →
        ∃ x' : K.closure, wildDepth K x' < wildDepth K x ∧
          ‖x - x'‖ ≤ c (wildDepth K x) * ε) :
    AxWildDescent K c := by
  intro ε hε x hx hdvd
  obtain ⟨x', hlt, hd⟩ := h ε hε x hx hdvd
  refine ⟨x', hlt, hd, fun σ => ?_⟩
  have hεle : ε ≤ c (wildDepth K x) * ε := by
    nlinarith [hc (wildDepth K x)]
  have := norm_smul_sub_self_le_max (G := K.absGal) (M := K.closure)
    (fun g z => norm_smul_closure K g z) hx hd σ
  exact this.trans (max_le le_rfl hεle)

end Concrete

def axWildDescent_of_dist.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §3 台帳 —— 距離の条件を `ℕ` の不等式に翻訳する -/

section Ledger

/-- ★★**台帳**（★実数と `ℕ` だけ。分岐・付値・Galois が 1 語も出ない）:

  `p^{j/E} ≤ axDecay p (m+1)` ⟺ `(p−1)·p^m·j ≤ E`。

体の層では `E = e_L`、`j = i + a − b`（`i` は層の跳び、`a = v_L(Δ)`、`b = v_L(τx−x)`）
を代入する。★`k = m+1` と書いたのは `axDecay` の `k−1` という `ℕ` の引き算を避けるため。 -/
theorem rpow_div_le_axDecay_iff (p : ℕ) [Fact p.Prime] (m j E : ℕ) (hE : 0 < E) :
    (p : ℝ) ^ ((j : ℝ) / (E : ℝ)) ≤ axDecay p (m + 1) ↔
      ((p : ℝ) - 1) * (p : ℝ) ^ m * j ≤ (E : ℝ) := by
  have hp1 : (1 : ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  have hE0 : (0 : ℝ) < (E : ℝ) := by exact_mod_cast hE
  have hpm : (0 : ℝ) < (p : ℝ) ^ m := by positivity
  have hp10 : (0 : ℝ) < (p : ℝ) - 1 := by linarith
  have hprod : (0 : ℝ) < ((p : ℝ) - 1) * (p : ℝ) ^ m := mul_pos hp10 hpm
  rw [axDecay, Nat.add_sub_cancel, Real.rpow_le_rpow_left_iff hp1]
  rw [show (1 / ((p : ℝ) - 1)) * (1 / (p : ℝ)) ^ m
      = 1 / (((p : ℝ) - 1) * (p : ℝ) ^ m) by rw [one_div_pow, div_mul_div_comm, one_mul]]
  rw [div_le_div_iff₀ hE0 hprod, one_mul]
  constructor <;> intro h <;> nlinarith [h]

/-- 台帳の対偶: `E < (p−1)·p^m·j` なら `axDecay p (m+1) < p^{j/E}`（1 段が届かない）。 -/
theorem axDecay_lt_rpow_div (p : ℕ) [Fact p.Prime] (m j E : ℕ) (hE : 0 < E)
    (h : (E : ℝ) < ((p : ℝ) - 1) * (p : ℝ) ^ m * j) :
    axDecay p (m + 1) < (p : ℝ) ^ ((j : ℝ) / (E : ℝ)) :=
  lt_of_not_ge fun hle => absurd ((rpow_div_le_axDecay_iff p m j E hE).mp hle) (not_le.mpr h)

end Ledger

def rpow_div_le_axDecay_iff.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def axDecay_lt_rpow_div.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §4 測定 —— 実在する分岐データが台帳を破る（★反証候補） -/

section Measurement

/-- ★★★**測定 3 の分岐データ**（★純 `ℕ`）。

`K = ℚ₃(ζ₃)`、`M = ℚ₃(ζ₂₇)`、`e_M = 18`、深さ `k = 2`、層 `M/ℚ₃(ζ₉)` の跳び `i = 8`。

* `2·8 ≤ 18` —— 層の跳びは sharp な `(p−1)i ≤ e` を満たす。
* `2·3·(8−6) ≤ 18` —— ★**素元**の gain `= 6` なら 1 段は通る。
* `¬(2·3·(8−4) ≤ 18)` —— ★★**実在する `x`** の gain `= 4` では通らない。 -/
theorem cyclotomic_three_jump_data :
    2 * 8 ≤ 18 ∧ 2 * 3 * (8 - 6) ≤ 18 ∧ ¬ (2 * 3 * (8 - 4) ≤ 18) := by
  refine ⟨by norm_num, by norm_num, by norm_num⟩

/-- ★★★**測定 3 を実数の不等式にしたもの**: `axDecay 3 2 = 3^{3/18} < 3^{4/18}`。

`gain = 4` の `x`（`K = ℚ₃(ζ₃)`, `M = ℚ₃(ζ₂₇)`, 深さ 2）では、
層 `M/ℚ₃(ζ₉)` を通る古典の 1 段の損失が `3^{(8−4)/18}` になり、
★`axDecay 3 2` を **`3^{1/18}` だけ**超える。

★★これは `AxWildDescent K (axDecay p)` が偽であることの証明ではない
（`M` の外の深さ ≤ 1 の `x'` を排除していない。モジュール docstring の測定 3 を読むこと）。 -/
theorem axDecay_three_two_lt_cyclotomic_layer_loss :
    axDecay 3 2 < (3 : ℝ) ^ (((4 : ℕ) : ℝ) / ((18 : ℕ) : ℝ)) := by
  refine axDecay_lt_rpow_div 3 1 4 18 (by norm_num) ?_
  norm_num

end Measurement

def cyclotomic_three_jump_data.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def axDecay_three_two_lt_cyclotomic_layer_loss.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §5 使っている公理の一覧 -/

#print axioms norm_smul_sub_self_le_max
#print axioms axWildDescent_of_dist
#print axioms rpow_div_le_axDecay_iff
#print axioms axDecay_lt_rpow_div
#print axioms cyclotomic_three_jump_data
#print axioms axDecay_three_two_lt_cyclotomic_layer_loss

end ABC3.Found.PGC
