import ABC3.Found.PGC.LubinTateTowerFIndependent

/-!
# `θ^{(j)}/θ = π′_j/π_j` —— 一様化元の列は 1-コサイクルである(Yoshida 2008 Lemma 4.5)

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) Lemma 4.5(物理 p.7)。
構造化済み原文は `ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
の `section-4.html` の `#lemma-4-5`。

原文 (Yoshida08 p.7):
> Lemma 4.5. If θ ∈ Θ^L_π,π′, then θ^(j)/θ = π′_j/π_j for all j ∈ Z[bb]. Also, π_j ∈ Θ^L_π,π^(j).

原典の証明(p.7–8、`0_Source` の `.txt` の 401–407 行を直読した全文):

> Proof. Using π′_{j+1}/π_{j+1} = (π′_j/π_j)(π′/π)^{(j)} = (π′_j/π_j)(θ^ϕ/θ)^{(j)}
> = (π′_j/π_j)(θ^{(j+1)}/θ^{(j)}), argue by induction in both directions.
> Take π′ = π^ϕ and θ = π for the second claim. □

記号は §4.2(物理 p.7)と Definition 3.9(物理 p.6)による:

> 4.2. Artin map. In this subsection we use the notation ( )^{(i)} := ( )^{ϕ^i} ...
> We extend Definition 3.9 to define π_j ∈ L^× for all j ∈ Z by requiring
> π_{j+j′} = π^{(j)}_{j′} π_j for all j, j′ ∈ Z, i.e. π_j := (π^{−1}_{−j})^{(j)} for j < 0.
> Then v_L(π_j) = j for all j ∈ Z.

> Definition 3.9. ... where we define π_m ∈ O_L by π_m := ∏^{m−1}_{t=0} π^{ϕ^t} and π_0 := 1.

## 何を足したか

### 抽象核(分岐・付値・Lubin-Tate・Frobenius の語彙が 1 つも現れない)

原典の主張から設定を全部抜くと ℤ 上の 1-コサイクルだけが残る。
`Γ = ℤ` が可換群 `M` に `σ : MulAut M` で作用しているとき、
`a : ℤ → M` が `a (j + j′) = σ^j (a j′) · a j` を満たすものを 1-コサイクルと呼ぶ。

| 宣言 | 内容 |
|---|---|
| `IsZCocycle` | `a (j + j′) = (σ ^ j) (a j′) * a j`(★これが §4.2 の `π_{j+j′} = π^{(j)}_{j′} π_j`) |
| `IsZCocycle.eval_zero` | `a 0 = 1`(★仮定ではなく帰結。これが両方向の帰納法の底) |
| `IsZCocycle.ext` | ★★一意性 —— `a 1 = b 1` なら `∀ j : ℤ, a j = b j`。原典の "argue by induction in both directions" そのもの |
| `isZCocycle_of_step` | 逆向き —— `a 0 = 1` と 1 段の漸化式だけからコサイクル条件が出る |
| `IsZCocycle.div` / `.map` | コサイクルの商・`σ` による移送もコサイクル |
| `isZCocycle_conjDiv` | ★`j ↦ (σ^j) θ / θ` はコサイクル(原典左辺) |
| `natProd` / `zpowProd` | ★★存在 —— `π_j`。`j ≥ 0` は `∏_{t<j} σ^t π`、`j < 0` は `σ^j((π_{−j})⁻¹)`。符号の場合分けはこの 1 つの定義と `zpowProd_step` に閉じ込めた |
| `zpowProd_natCast` | `π_n = ∏_{t<n} σ^t π`(Definition 3.9 の `π_m`) |
| `zpowProd_neg_natCast` / `zpowProd_neg` | `π_j = (π_{−j}⁻¹)^{(j)}`(§4.2 の `j < 0` の定義式) |
| `isZCocycle_zpowProd` / `zpowProd_add` | `π_{j+j′} = π^{(j)}_{j′} π_j`(§4.2 の要求) |
| `eq_zpowProd_of_isZCocycle` | コサイクルは `a 1` で決まる |
| `zpowProd_map` | ★`(σπ)_j = σ(π_j)`(第 2 主張の代入に要る) |
| `zpow_div_self_eq_zpowProd_div` | ★★★抽象版 Lemma 4.5 第 1 主張 `σθ/θ = π′/π ⇒ (σ^j)θ/θ = π′_j/π_j` |
| `zpowProd_apply_div_self` | ★★抽象版 第 2 主張 `σ(π_j)/π_j = (σ^j)π/π`(原典どおり `π′ = σπ`, `θ = π` の代入で出る) |
| `zpowProd_one_aut` | ★退化: `σ = 1` なら `π_j = π^j` |
| `map_zpow_apply_of_invariant` / `map_zpowProd` | `σ` 不変な準同型 `v` に対し `v(π_j) = v(π)^j`(§4.2 の `v_L(π_j) = j` の抽象形) |

### 橋(環自己同型 → 単元群の自己同型)

| 宣言 | 内容 |
|---|---|
| `unitsRingAutHom` | `RingAut R →* MulAut Rˣ`。★群準同型として作るので `map_zpow` が使え、`(ϕ^j)` と `(σ^j)` が一致する |
| `coe_unitsRingAutHom_zpow` / `_pow` | `((σ^j) u : R) = (ϕ^j) (u : R)` |

### 具体層(体 `L` と環自己同型 `ϕ`)

| 宣言 | 内容 |
|---|---|
| `uniformizerZ` | ★`π_j ∈ L^×`(`ϕ : RingAut L`、`π : Lˣ`、`j : ℤ`) |
| `coe_uniformizerZ_natCast` | `π_n = ∏_{t<n} ϕ^t(π)`(Definition 3.9 の `π_m`) |
| `coe_uniformizerZ_add` | `π_{j+j′} = ϕ^j(π_{j′}) · π_j`(§4.2 の要求) |
| `coe_uniformizerZ_neg` | `π_j = ϕ^j((π_{−j})⁻¹)`(§4.2 の `j < 0`) |
| `ThetaSet` | `Θ^L_{π,π′}`(Definition 3.3) |
| `zpow_apply_mul_uniformizerZ` | ★★★Lemma 4.5 第 1 主張(積の形。`θ = 0` も込み) |
| `zpow_apply_div_self_eq_uniformizerZ_div` | ★原典の字面どおりの商の形 `θ^{(j)}/θ = π′_j/π_j` |
| `coe_uniformizerZ_mem_thetaSet` | ★★Lemma 4.5 第 2 主張 `π_j ∈ Θ^L_{π,π^{(j)}}`。原典どおり「`π′ = π^ϕ`, `θ = π` を代入」して出した |

## 逸脱の記録

1. Θ^L_{π,π′} を割り算でなく掛け算で書いた。原典 Definition 3.3 は
   `Θ^L_{π,π′} := {θ ∈ O_L | θ^ϕ/θ = π′/π}` だが、同じ行で
   「It is an additive group」と述べている。加法群であるためには `0 ∈ Θ` が要り、
   `0^ϕ/0` は定義できない。したがって原典の意図は交差積
   `θ^ϕ · π = θ · π′` であると読んだ(`ThetaSet`)。★`θ ≠ 0` のときは同値であり、
   商の形も `zpow_apply_div_self_eq_uniformizerZ_div` として別に出している。
2. Θ の「`θ ∈ O_L`」(整性)を落とした。理由は 2 つある。
   (a) 証明のどこにも整性を使わない。
   (b) ★第 2 主張 `π_j ∈ Θ^L_{π,π^{(j)}}` は `j < 0` では `π_j ∉ O_L` なので、
   `Θ ⊂ O_L` のままだと原典の第 2 主張が `j ≥ 0` でしか意味を持たない。
   原典が `∀ j ∈ Z` を意図していると読み、`Θ` を `L` の中で定義した。
3. `π, π′` を「一様化元」でなく「`L^×` の元」にした。証明は付値を一切使わない。
   一様化元であることは `v_L(π_j) = j` を出すときに初めて要る(下記 4)。
4. `v_L(π_j) = j`(§4.2 の最後の文)は具体層では出していない。抽象核の
   `map_zpowProd`(`v(π_j) = v(π)^j`、`v` は `σ` 不変な群準同型)まで作った。
   ★具体層に落とすには `L` の付値を `Lˣ →* Multiplicative ℤ` の形で用意する配管が
   要る。Lemma 4.5 の消費者(Corollary 4.9 後半)は付値を使わないので入れていない。
5. `ϕ` が算術 Frobenius であること・`L` が完備不分岐拡大であることを使っていない。
   `ϕ : RingAut L` は任意の環自己同型でよい。★これは弱めではなく一般化であり、
   `arithFrobenius` を代入すれば原典の設定になる。
6. `uniformizerProd`(`UniformizerExpansion.lean:164`)を ℤ へ延ばすのではなく、
   ℤ 版 `zpowProd` を新設した。理由: `uniformizerProd` は `[Monoid G]`
   `[MulSemiringAction G B]` の上に書かれていて、(i) `σ^j`(`j : ℤ`)が作れない、
   (ii) `B` が環なので `j < 0` の `(π_{−j})⁻¹` が作れない。
   ★`M := Lˣ` と `σ : MulAut M` に移せば両方が同時に解ける。
   `zpowProd_natCast`(`π_n = ∏_{t<n} σ^t π`)が `uniformizerProd` と同じ形なので、
   `j ≥ 0` では両者は一致する。

## ★退化の自己検査

* (D1) `j = 0` で両辺 `1`。`zpowProd_zero`(`π_0 = 1`、Definition 3.9 の `π_0 := 1`)。
  これが両方向の帰納法の底であり、`IsZCocycle.eval_zero` が示すとおり
  仮定ではなくコサイクル条件からの帰結である。
* (D2) 負の `j` を落とすと下流が閉じない。Corollary 4.9 後半は
  `σ ∈ W(K̂^m_f/K) = ϕ^Z` の Frobenius 次数 `j` に対して `σ(α) = [xπ_j](α)` と書くので、
  `j < 0` が本当に現れる。★本ファイルは `j : ℤ` で通してある(`ℕ` に逃げていない)。
* (D3) `θ ∈ Θ^L_{π,π′}` を落とすと偽。`θ` を任意にすると `j = 1` で
  `θ^ϕ/θ = π′/π` を主張することになり、これは仮定そのものである。
* (D4) `M` を可換群にしないと順序が効く。`π_{j+j′} = σ^j(π_{j′}) · π_j` の
  右辺の順序は原典どおりに保ってあるが、`IsZCocycle.div` と
  `zpow_div_self_eq_zpowProd_div` は可換性を使う(商が再びコサイクルになるため)。
  ★可換を落とすなら「商」を「捻れた商」に置き換える必要がある。
* (D5) 除算の分母。`L` は体なので `π_j` が `0` になりうる。これを避けるため
  `π_j` は `Lˣ` の元として定義した(`uniformizerZ : ... → Lˣ`)。
  ★`Units.ne_zero` が分母の非零を無条件に供給するので、
  `π ≠ 0` の仮定を持ち回る必要が無くなった。
* (D6) `σ = 1` のとき `π_j = π^j`(`zpowProd_one_aut`)。
  Frobenius が自明なら望遠鏡積はただの冪に退化する。
-/

namespace ABC3.Found.PGC

def zpow_apply_mul_uniformizerZ.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 7, item := "Lemma 4.5", sectionId := "lemma-4-5" }

/-! ## §1 抽象核 —— ℤ 上の 1-コサイクル

★ここには分岐・付値・Lubin-Tate・Frobenius の語彙が 1 つも現れない。
一般の可換群 `M` と一般の自己同型 `σ : MulAut M` の話である。 -/

section Abstract

variable {M : Type*} [CommGroup M]

/-- `σ^j (σ^k x) = σ^{j+k} x`。`MulAut M` の積が合成 `(f * g) x = f (g x)` であることによる。 -/
theorem zpow_apply_zpow_apply (σ : MulAut M) (j k : ℤ) (x : M) :
    (σ ^ j) ((σ ^ k) x) = (σ ^ (j + k)) x := by
  rw [← MulAut.mul_apply, ← zpow_add]

/-- **ℤ 上の 1-コサイクル** `a (j + j′) = σ^j (a j′) · a j`。

★これは原典 §4.2 の `π_{j+j′} = π^{(j)}_{j′} π_j`(物理 p.7)をそのまま抽象化したものである。
右辺の順序は原典どおり(`σ^j` を掛けたものが左)。 -/
def IsZCocycle (σ : MulAut M) (a : ℤ → M) : Prop :=
  ∀ j j' : ℤ, a (j + j') = (σ ^ j) (a j') * a j

/-- **`a 0 = 1` は仮定ではなく帰結**。`j = j′ = 0` で `a 0 = a 0 * a 0`。

★これが「両方向の帰納法」の底である(退化検査 D1)。 -/
theorem IsZCocycle.eval_zero {σ : MulAut M} {a : ℤ → M} (ha : IsZCocycle σ a) : a 0 = 1 := by
  have h := ha 0 0
  simp at h
  exact h

/-- 1 段の漸化式 `a (j + 1) = σ^j (a 1) · a j`(`j′ = 1` を代入しただけ)。 -/
theorem IsZCocycle.step {σ : MulAut M} {a : ℤ → M} (ha : IsZCocycle σ a) (j : ℤ) :
    a (j + 1) = (σ ^ j) (a 1) * a j := ha j 1

/-- ★★★**コサイクルの一意性** —— `a 1 = b 1` なら `∀ j : ℤ, a j = b j`。

これが原典の "argue by induction in both directions" の中身である。
`j ≥ 0` では `step` を前向きに、`j < 0` では `step` を後ろ向きに(左簡約で)使う。 -/
theorem IsZCocycle.ext {σ : MulAut M} {a b : ℤ → M} (ha : IsZCocycle σ a) (hb : IsZCocycle σ b)
    (h1 : a 1 = b 1) : ∀ j : ℤ, a j = b j := by
  intro j
  induction j using Int.induction_on with
  | zero => rw [ha.eval_zero, hb.eval_zero]
  | succ i ih => rw [ha.step, hb.step, h1, ih]
  | pred i ih =>
      have hA := ha.step (-(i : ℤ) - 1)
      have hB := hb.step (-(i : ℤ) - 1)
      rw [show -(i : ℤ) - 1 + 1 = -(i : ℤ) by ring] at hA hB
      rw [h1] at hA
      rw [ih, hB] at hA
      exact (mul_left_cancel hA).symm

/-- **逆向き** —— `a 0 = 1` と 1 段の漸化式だけからコサイクル条件が出る。

★これも両方向の帰納法。`j` を固定して `j′` について回す。 -/
theorem isZCocycle_of_step {σ : MulAut M} {a : ℤ → M} (h0 : a 0 = 1)
    (hs : ∀ k : ℤ, a (k + 1) = (σ ^ k) (a 1) * a k) : IsZCocycle σ a := by
  intro j j'
  induction j' using Int.induction_on with
  | zero => simp [h0]
  | succ i ih =>
      rw [show j + ((i : ℤ) + 1) = (j + i) + 1 by ring, hs (j + i), ih, hs (i : ℤ), map_mul,
        zpow_apply_zpow_apply, mul_assoc]
  | pred i ih =>
      have hA := hs (j + (-(i : ℤ) - 1))
      rw [show j + (-(i : ℤ) - 1) + 1 = j + -(i : ℤ) by ring] at hA
      have hB := hs (-(i : ℤ) - 1)
      rw [show -(i : ℤ) - 1 + 1 = -(i : ℤ) by ring] at hB
      rw [ih, hB, map_mul, zpow_apply_zpow_apply, mul_assoc] at hA
      exact mul_left_cancel hA.symm

/-- コサイクルの商はコサイクル(★可換性を使う。退化検査 D4)。 -/
theorem IsZCocycle.div {σ : MulAut M} {a b : ℤ → M} (ha : IsZCocycle σ a) (hb : IsZCocycle σ b) :
    IsZCocycle σ (fun j => a j / b j) := by
  intro j j'
  simp only [ha j j', hb j j', map_div, div_mul_div_comm]

/-- コサイクルを `σ` で送ってもコサイクル(`σ` と `σ^j` が可換だから)。 -/
theorem IsZCocycle.map {σ : MulAut M} {a : ℤ → M} (ha : IsZCocycle σ a) :
    IsZCocycle σ (fun j => σ (a j)) := by
  intro j j'
  simp only [ha j j', map_mul, ← MulAut.mul_apply, ← zpow_add_one, ← zpow_one_add,
    show (1 : ℤ) + j = j + 1 from by ring]

/-- ★**原典の左辺** `j ↦ θ^{(j)}/θ` はコサイクルである。

`σ^{j+j′}θ/θ = (σ^j(σ^{j′}θ/θ)) · (σ^jθ/θ)` は `(a/b)(b/c) = a/c` そのもの。 -/
theorem isZCocycle_conjDiv (σ : MulAut M) (θ : M) : IsZCocycle σ (fun j => (σ ^ j) θ / θ) := by
  intro j j'
  simp only [map_div, zpow_apply_zpow_apply]
  rw [div_mul_div_cancel]

/-- **望遠鏡積の `ℕ` 版** `π_n := ∏_{t<n} σ^t π`(原典 Definition 3.9 の `π_m`)。 -/
def natProd (σ : MulAut M) (π : M) (n : ℕ) : M := ∏ t ∈ Finset.range n, (σ ^ t) π

theorem natProd_zero (σ : MulAut M) (π : M) : natProd σ π 0 = 1 := by simp [natProd]

theorem natProd_one (σ : MulAut M) (π : M) : natProd σ π 1 = π := by simp [natProd]

theorem natProd_succ (σ : MulAut M) (π : M) (n : ℕ) :
    natProd σ π (n + 1) = natProd σ π n * (σ ^ n) π := by
  simp [natProd, Finset.prod_range_succ]

/-- 反対側の切り出し `π_{n+1} = π · σ(π_n)`。★`j < 0` の段を作るのに要る。 -/
theorem natProd_succ' (σ : MulAut M) (π : M) (n : ℕ) :
    natProd σ π (n + 1) = π * σ (natProd σ π n) := by
  rw [natProd, natProd, Finset.prod_range_succ' (fun t => (σ ^ t) π) n, map_prod]
  simp [pow_succ', mul_comm]

/-- ★★**望遠鏡積の `ℤ` 版** `π_j`。

`j ≥ 0` では `(-j).toNat = 0` なので `∏_{t<j} σ^t π`、
`j ≤ 0` では `j.toNat = 0` なので `σ^j((π_{−j})⁻¹)` に退化する
(原典 §4.2 の `π_j := (π^{−1}_{−j})^{(j)} for j < 0`)。
★符号の場合分けはこの 1 つの式と `zpowProd_step` の証明だけに閉じ込めてある。 -/
def zpowProd (σ : MulAut M) (π : M) (j : ℤ) : M :=
  natProd σ π j.toNat * (σ ^ j) (natProd σ π (-j).toNat)⁻¹

/-- `j ≥ 0` では Definition 3.9 の `π_m := ∏_{t=0}^{m−1} π^{ϕ^t}` そのもの。 -/
theorem zpowProd_natCast (σ : MulAut M) (π : M) (n : ℕ) :
    zpowProd σ π (n : ℤ) = natProd σ π n := by
  rw [zpowProd]; simp [natProd_zero]

/-- `j ≤ 0` では §4.2 の `π_j := (π^{−1}_{−j})^{(j)}` そのもの。 -/
theorem zpowProd_neg_natCast (σ : MulAut M) (π : M) (n : ℕ) :
    zpowProd σ π (-(n : ℤ)) = (σ ^ (-(n : ℤ))) (natProd σ π n)⁻¹ := by
  have h1 : (-(n : ℤ)).toNat = 0 := by omega
  have h2 : (-(-(n : ℤ))).toNat = n := by omega
  rw [zpowProd, h1, h2, natProd_zero, one_mul]

/-- ★`π_0 = 1`(Definition 3.9 の `π_0 := 1`)。両方向の帰納法の底(退化検査 D1)。 -/
theorem zpowProd_zero (σ : MulAut M) (π : M) : zpowProd σ π 0 = 1 := by
  simpa [natProd_zero] using zpowProd_natCast σ π 0

/-- `π_1 = π`。 -/
theorem zpowProd_one (σ : MulAut M) (π : M) : zpowProd σ π 1 = π := by
  simpa [natProd_one] using zpowProd_natCast σ π 1

/-- ★**1 段** `π_{j+1} = σ^j(π) · π_j`(`j : ℤ`、符号を問わない)。

★符号の場合分けはここだけ。`j ≥ 0` は `Finset.prod_range_succ`、
`j < 0` は `natProd_succ'`(`π_{n+1} = π · σ(π_n)`)で「先頭の `π` を消す」。 -/
theorem zpowProd_step (σ : MulAut M) (π : M) (j : ℤ) :
    zpowProd σ π (j + 1) = (σ ^ j) π * zpowProd σ π j := by
  rcases (by omega : j < 0 ∨ 0 ≤ j) with hj | hj
  · obtain ⟨m, rfl⟩ := Int.eq_negSucc_of_lt_zero hj
    have hneg : (Int.negSucc m) = -((m + 1 : ℕ) : ℤ) := by push_cast [Int.negSucc_eq]; ring
    have key : π * (natProd σ π (m + 1))⁻¹ = (σ (natProd σ π m))⁻¹ := by
      rw [natProd_succ', mul_inv_rev, mul_comm ((σ (natProd σ π m))⁻¹) π⁻¹, mul_inv_cancel_left]
    rw [hneg, show -((m + 1 : ℕ) : ℤ) + 1 = -((m : ℕ) : ℤ) by push_cast; ring,
      zpowProd_neg_natCast, zpowProd_neg_natCast, ← map_mul, key, ← map_inv,
      show -((m + 1 : ℕ) : ℤ) = -((m : ℕ) : ℤ) - 1 by push_cast; ring,
      ← MulAut.mul_apply, ← zpow_add_one,
      show -((m : ℕ) : ℤ) - 1 + 1 = -((m : ℕ) : ℤ) from by ring]
  · obtain ⟨n, rfl⟩ := Int.eq_ofNat_of_zero_le hj
    rw [show (n : ℤ) + 1 = ((n + 1 : ℕ) : ℤ) by push_cast; ring,
      zpowProd_natCast, zpowProd_natCast, natProd_succ, zpow_natCast, mul_comm]

/-- ★★**存在** —— `π_j` は 1-コサイクルである。 -/
theorem isZCocycle_zpowProd (σ : MulAut M) (π : M) : IsZCocycle σ (zpowProd σ π) :=
  isZCocycle_of_step (zpowProd_zero σ π) fun k => by rw [zpowProd_one, zpowProd_step]

/-- 原典 §4.2 が `π_j` に要求している式 `π_{j+j′} = π^{(j)}_{j′} π_j`。 -/
theorem zpowProd_add (σ : MulAut M) (π : M) (j j' : ℤ) :
    zpowProd σ π (j + j') = (σ ^ j) (zpowProd σ π j') * zpowProd σ π j :=
  isZCocycle_zpowProd σ π j j'

/-- `j ≤ 0` の閉じた形 `π_j = (π_{−j})⁻¹ の σ^j 像`(§4.2 の定義式)。 -/
theorem zpowProd_neg (σ : MulAut M) (π : M) {j : ℤ} (hj : j ≤ 0) :
    zpowProd σ π j = (σ ^ j) (zpowProd σ π (-j))⁻¹ := by
  obtain ⟨n, rfl⟩ : ∃ n : ℕ, j = -(n : ℤ) := ⟨(-j).toNat, by omega⟩
  rw [zpowProd_neg_natCast, neg_neg, zpowProd_natCast]

/-- ★**コサイクルは `a 1` で決まる** —— §4.2 が「`π_{j+j′} = π^{(j)}_{j′}π_j` で定める」と
書けるのはこれによる(値 `π_1 = π` を与えれば列が一意に決まる)。 -/
theorem eq_zpowProd_of_isZCocycle {σ : MulAut M} {a : ℤ → M} (ha : IsZCocycle σ a) (j : ℤ) :
    a j = zpowProd σ (a 1) j :=
  ha.ext (isZCocycle_zpowProd σ (a 1)) (zpowProd_one σ (a 1)).symm j

/-- ★**同変性** `(σπ)_j = σ(π_j)`。第 2 主張の「代入」に要る。

★証明は一意性の再利用: 両辺とも `j` のコサイクルで、`j = 1` での値が `σπ` で一致する。 -/
theorem zpowProd_map (σ : MulAut M) (π : M) (j : ℤ) :
    zpowProd σ (σ π) j = σ (zpowProd σ π j) :=
  (isZCocycle_zpowProd σ (σ π)).ext (isZCocycle_zpowProd σ π).map
    (by rw [zpowProd_one]; simp [zpowProd_one]) j

/-- ★★★**抽象版 Lemma 4.5 第 1 主張**。

`σθ/θ = π′/π` ならば `∀ j : ℤ, (σ^j)θ/θ = π′_j/π_j`。
★証明は「両辺がコサイクルで `j = 1` で一致する」だけ —— 原典の
"argue by induction in both directions" は `IsZCocycle.ext` に吸収されている。 -/
theorem zpow_div_self_eq_zpowProd_div {σ : MulAut M} {π π' θ : M} (hθ : σ θ / θ = π' / π)
    (j : ℤ) : (σ ^ j) θ / θ = zpowProd σ π' j / zpowProd σ π j := by
  refine (isZCocycle_conjDiv σ θ).ext
    ((isZCocycle_zpowProd σ π').div (isZCocycle_zpowProd σ π)) ?_ j
  simpa [zpowProd_one, zpow_one] using hθ

/-- ★★**抽象版 Lemma 4.5 第 2 主張** `σ(π_j)/π_j = (σ^j)π/π`。

原典どおり第 1 主張に `π′ = σπ`、`θ = π` を代入して出した
("Take π′ = π^ϕ and θ = π for the second claim")。 -/
theorem zpowProd_apply_div_self (σ : MulAut M) (π : M) (j : ℤ) :
    σ (zpowProd σ π j) / zpowProd σ π j = (σ ^ j) π / π := by
  have h := zpow_div_self_eq_zpowProd_div (σ := σ) (π := π) (π' := σ π) (θ := π) rfl j
  rw [h, zpowProd_map]

/-- ★**退化(D6)** —— `σ = 1` なら望遠鏡積はただの冪 `π_j = π^j`。 -/
theorem zpowProd_one_aut (x : M) (j : ℤ) : zpowProd (1 : MulAut M) x j = x ^ j := by
  have hc : IsZCocycle (1 : MulAut M) (fun j : ℤ => x ^ j) := by
    intro j j'
    simp [zpow_add, mul_comm]
  have h := eq_zpowProd_of_isZCocycle hc j
  simpa using h.symm

/-- `σ` 不変な準同型は `σ^j` でも不変。 -/
theorem map_zpow_apply_of_invariant {N : Type*} [CommGroup N] {σ : MulAut M} (v : M →* N)
    (hv : ∀ x, v (σ x) = v x) (j : ℤ) (x : M) : v ((σ ^ j) x) = v x := by
  have hv' : ∀ x, v (σ⁻¹ x) = v x := by
    intro x
    have h := hv (σ⁻¹ x)
    rw [show σ (σ⁻¹ x) = x from by rw [← MulAut.mul_apply, mul_inv_cancel]; rfl] at h
    exact h.symm
  induction j using Int.induction_on generalizing x with
  | zero => simp
  | succ i ih => rw [zpow_add_one, MulAut.mul_apply, ih, hv]
  | pred i ih => rw [zpow_sub_one, MulAut.mul_apply, ih, hv']

/-- ★`v(π_j) = v(π)^j` —— 原典 §4.2 の「Then `v_L(π_j) = j` for all `j ∈ Z`」の抽象形。

`v` は `σ` 不変な群準同型。具体層への落とし込み(`L` の付値を `Lˣ →* Multiplicative ℤ`
の形で用意する配管)は本ファイルでは行っていない(逸脱 4)。 -/
theorem map_zpowProd {N : Type*} [CommGroup N] {σ : MulAut M} (v : M →* N)
    (hv : ∀ x, v (σ x) = v x) (π : M) (j : ℤ) : v (zpowProd σ π j) = (v π) ^ j := by
  have hc : IsZCocycle (1 : MulAut N) (fun j : ℤ => v (zpowProd σ π j)) := by
    intro j j'
    simp only [zpowProd_add, map_mul, map_zpow_apply_of_invariant v hv, one_zpow]
    simp
  have h := eq_zpowProd_of_isZCocycle hc j
  simpa [zpowProd_one, zpowProd_one_aut] using h

end Abstract

/-! ## §2 橋 —— 環自己同型から単元群の自己同型へ

★`MonoidHom` として作るのが要点。そうすれば `map_zpow` が使えて
`(ϕ^j)` と `(σ^j)` が自動で対応する(場合分けが要らない)。 -/

/-- `RingAut R →* MulAut Rˣ`。★群準同型として作る。 -/
def unitsRingAutHom (R : Type*) [CommRing R] : RingAut R →* MulAut Rˣ where
  toFun ϕ := Units.mapEquiv ϕ.toMulEquiv
  map_one' := by ext u; rfl
  map_mul' f g := by ext u; rfl

theorem coe_unitsRingAutHom_zpow {R : Type*} [CommRing R] (ϕ : RingAut R) (j : ℤ) (u : Rˣ) :
    ((((unitsRingAutHom R ϕ) ^ j) u : Rˣ) : R) = (ϕ ^ j : RingAut R) (u : R) := by
  rw [← map_zpow (unitsRingAutHom R)]; rfl

theorem coe_unitsRingAutHom_pow {R : Type*} [CommRing R] (ϕ : RingAut R) (n : ℕ) (u : Rˣ) :
    ((((unitsRingAutHom R ϕ) ^ n) u : Rˣ) : R) = (ϕ ^ n : RingAut R) (u : R) := by
  rw [← map_pow (unitsRingAutHom R)]; rfl

/-! ## §3 具体層 —— 体 `L` と環自己同型 `ϕ`

原典の `L` は「`K` の完備不分岐拡大」、`ϕ` は「算術 Frobenius」だが、
★本節はそのどちらも使わない(逸脱 5)。`arithFrobenius` を代入すれば原典の設定になる。 -/

section Concrete

variable {L : Type*} [Field L]

/-- ★**`π_j ∈ L^×`**(原典 §4.2、物理 p.7)。

`j ≥ 0` では `∏_{t=0}^{j−1} π^{ϕ^t}`(Definition 3.9)、
`j < 0` では `(π^{−1}_{−j})^{(j)}`。
★値域を `Lˣ` にしたので、割り算の分母の非零が無条件に手に入る(退化検査 D5)。 -/
def uniformizerZ (ϕ : RingAut L) (π : Lˣ) (j : ℤ) : Lˣ :=
  zpowProd (unitsRingAutHom L ϕ) π j

/-- **Definition 3.9 の `π_m := ∏_{t=0}^{m−1} π^{ϕ^t}`**(物理 p.6)。 -/
theorem coe_uniformizerZ_natCast (ϕ : RingAut L) (π : Lˣ) (n : ℕ) :
    ((uniformizerZ ϕ π (n : ℤ) : Lˣ) : L) = ∏ t ∈ Finset.range n, (ϕ ^ t : RingAut L) (π : L) := by
  rw [uniformizerZ, zpowProd_natCast, natProd]
  refine (map_prod (Units.coeHom L) (fun t => ((unitsRingAutHom L ϕ) ^ t) π)
    (Finset.range n)).trans ?_
  exact Finset.prod_congr rfl fun t _ => coe_unitsRingAutHom_pow ϕ t π

/-- **原典 §4.2 が `π_j` に要求する式** `π_{j+j′} = π^{(j)}_{j′} π_j`(∀ `j, j′ ∈ ℤ`)。 -/
theorem coe_uniformizerZ_add (ϕ : RingAut L) (π : Lˣ) (j j' : ℤ) :
    ((uniformizerZ ϕ π (j + j') : Lˣ) : L)
      = (ϕ ^ j : RingAut L) ((uniformizerZ ϕ π j' : Lˣ) : L) * ((uniformizerZ ϕ π j : Lˣ) : L) := by
  rw [uniformizerZ, uniformizerZ, uniformizerZ, zpowProd_add, Units.val_mul,
    coe_unitsRingAutHom_zpow]

/-- **原典 §4.2 の `j < 0` の定義式** `π_j := (π^{−1}_{−j})^{(j)}`。 -/
theorem coe_uniformizerZ_neg (ϕ : RingAut L) (π : Lˣ) {j : ℤ} (hj : j ≤ 0) :
    ((uniformizerZ ϕ π j : Lˣ) : L)
      = (ϕ ^ j : RingAut L) (((uniformizerZ ϕ π (-j) : Lˣ)⁻¹ : Lˣ) : L) := by
  rw [uniformizerZ, uniformizerZ, zpowProd_neg _ _ hj, ← coe_unitsRingAutHom_zpow]

/-- `π_0 = 1`(Definition 3.9)。 -/
theorem uniformizerZ_zero (ϕ : RingAut L) (π : Lˣ) : uniformizerZ ϕ π 0 = 1 :=
  zpowProd_zero _ _

/-- `π_1 = π`。 -/
theorem uniformizerZ_one (ϕ : RingAut L) (π : Lˣ) : uniformizerZ ϕ π 1 = π :=
  zpowProd_one _ _

/-- **`Θ^L_{π,π′}`**(原典 Definition 3.3、物理 p.5)。

原文:
> Definition 3.3. For uniformizers π, π′ of L, set Θ^L_π,π′ := {θ ∈ O_L | θ^ϕ/θ = π′/π}.
> It is an additive group.

★「additive group」であるためには `0 ∈ Θ` が要るので、商ではなく交差積で書いた(逸脱 1)。
★`θ ∈ O_L` は落とした(逸脱 2)——`j < 0` では `π_j ∉ O_L` なので、
第 2 主張 `π_j ∈ Θ^L_{π,π^{(j)}}` を `∀ j ∈ ℤ` で述べるには `Θ` を `L` の中で取るしかない。 -/
def ThetaSet (ϕ : RingAut L) (π π' : L) : Set L := {θ : L | ϕ θ * π = θ * π'}

theorem mem_thetaSet_iff (ϕ : RingAut L) (π π' θ : L) :
    θ ∈ ThetaSet ϕ π π' ↔ ϕ θ * π = θ * π' := Iff.rfl

/-- ★★★**Lemma 4.5 第 1 主張**(Yoshida08 物理 p.7)。

原文 (Yoshida08 p.7):
> Lemma 4.5. If θ ∈ Θ^L_π,π′, then θ^(j)/θ = π′_j/π_j for all j ∈ Z[bb].

★`j : ℤ` である(負の `j` を落としていない)。Corollary 4.9 後半は
`σ ∈ W(K̂^m_f/K) = ϕ^ℤ` の Frobenius 次数 `j` に対してこれを使うので、
`j < 0` が本当に現れる(退化検査 D2)。

★積の形で述べたので `θ = 0` でも正しい(逸脱 1)。商の形は
`zpow_apply_div_self_eq_uniformizerZ_div`。 -/
theorem zpow_apply_mul_uniformizerZ (ϕ : RingAut L) (π π' : Lˣ) {θ : L}
    (hθ : θ ∈ ThetaSet ϕ (π : L) (π' : L)) (j : ℤ) :
    (ϕ ^ j : RingAut L) θ * ((uniformizerZ ϕ π j : Lˣ) : L)
      = θ * ((uniformizerZ ϕ π' j : Lˣ) : L) := by
  rcases eq_or_ne θ 0 with rfl | hθ0
  · simp
  · have hθ' : ϕ θ * (π : L) = θ * (π' : L) := hθ
    have hu : ((Units.mk0 θ hθ0 : Lˣ) : L) = θ := rfl
    have hcond : (unitsRingAutHom L ϕ) (Units.mk0 θ hθ0) / Units.mk0 θ hθ0 = π' / π := by
      rw [div_eq_div_iff_mul_eq_mul]
      ext
      show ϕ θ * (π : L) = (π' : L) * θ
      rw [hθ', mul_comm]
    have h := zpow_div_self_eq_zpowProd_div hcond j
    rw [div_eq_div_iff_mul_eq_mul] at h
    have h' := congrArg (fun u : Lˣ => (u : L)) h
    simp only [Units.val_mul] at h'
    rw [coe_unitsRingAutHom_zpow, hu] at h'
    rw [uniformizerZ, uniformizerZ, h', mul_comm]

/-- ★**原典の字面どおりの商の形** `θ^{(j)}/θ = π′_j/π_j`(`θ ∈ L^×`)。

分母の非零は `Units.ne_zero` から無条件に出る(退化検査 D5)。 -/
theorem zpow_apply_div_self_eq_uniformizerZ_div (ϕ : RingAut L) (π π' θ : Lˣ)
    (hθ : ((θ : L)) ∈ ThetaSet ϕ (π : L) (π' : L)) (j : ℤ) :
    (ϕ ^ j : RingAut L) (θ : L) / (θ : L)
      = ((uniformizerZ ϕ π' j : Lˣ) : L) / ((uniformizerZ ϕ π j : Lˣ) : L) := by
  rw [div_eq_div_iff θ.ne_zero (Units.ne_zero _), zpow_apply_mul_uniformizerZ ϕ π π' hθ j,
    mul_comm]

/-- ★★**Lemma 4.5 第 2 主張**(Yoshida08 物理 p.7)。

原文 (Yoshida08 p.7):
> Also, π_j ∈ Θ^L_π,π^(j).

★原典どおり第 1 主張に `π′ = π^ϕ`、`θ = π` を代入して出した
("Take π′ = π^ϕ and θ = π for the second claim")。代入で出るのは
第 1 主張を `∀ π′, ∀ θ` の一般形で述べたからである。
`(π^ϕ)_j = (π_j)^ϕ` の一致は `zpowProd_map`(同変性)。 -/
theorem coe_uniformizerZ_mem_thetaSet (ϕ : RingAut L) (π : Lˣ) (j : ℤ) :
    ((uniformizerZ ϕ π j : Lˣ) : L) ∈ ThetaSet ϕ (π : L) ((ϕ ^ j : RingAut L) (π : L)) := by
  have hsub := zpow_apply_mul_uniformizerZ ϕ π ((unitsRingAutHom L ϕ) π) (θ := (π : L))
    (show ϕ (π : L) * (π : L) = (π : L) * (((unitsRingAutHom L ϕ) π : Lˣ) : L) from mul_comm _ _) j
  have hmap : uniformizerZ ϕ ((unitsRingAutHom L ϕ) π) j
      = (unitsRingAutHom L ϕ) (uniformizerZ ϕ π j) := zpowProd_map _ _ _
  have hcoe : ((((unitsRingAutHom L ϕ) (uniformizerZ ϕ π j)) : Lˣ) : L)
      = ϕ ((uniformizerZ ϕ π j : Lˣ) : L) := rfl
  rw [hmap, hcoe] at hsub
  show ϕ ((uniformizerZ ϕ π j : Lˣ) : L) * (π : L)
      = ((uniformizerZ ϕ π j : Lˣ) : L) * (ϕ ^ j : RingAut L) (π : L)
  linear_combination -hsub

end Concrete

end ABC3.Found.PGC
