import ABC3.Found.PGC.AxEpsilonDecay

/-!
# [pGC] 幾何減衰 `θ^k` は Hasse–Arf からは**出ない** —— 1 段は出るが、一様な `θ < 1` は取れない

## ★配られた 3 つの問いに答えた

**(1)「`g` が `σ^{p^k}` の形」と「`H`（1 段ごとに `θ` 倍）」は同値か。**

★**同値ではない。そして `H` の方が弱くもない。**
`norm_iterate_pow_sub_self_le`（`Found/PGC/AxEpsilonDecay.lean:294`）の**結論**は
`‖f^[n^k] x − x‖ ≤ θ^k ‖f x − x‖` であって、★**`f^[n^k]`（＝ `σ^{p^k}`）についてのもの**である。
⇒ この補題を降下に使うには、その段で使う `g` が **`σ^{p^k}` そのもの**でなければならない。
★`H` は「`g = σ^{p^k}` を経由せずに済む逃げ道」ではなく、
**`g = σ^{p^k}` を前提にした上でさらに要求される条件**である。

**(2) `JumpStrictMono.lean:237 norm_pow_prime_sub_lt` / `:308 jump_lt_succ` は `H` に化けるか。**

★**1 段だけなら化ける**（§3 `norm_step_le_of_jump`）。
`u m < u (m+1)`（狭義単調）と `p^{m+1} ∣ u(m+1) − u m`（Hasse–Arf の合同）から
`u (m+1) ≥ u m + p^{m+1}` が出る（§1 `add_le_of_dvd_sub_of_lt`、純 `ℤ`）ので
`‖w (m+1)‖ ≤ ‖π‖^{p^{m+1}} · ‖w m‖`、すなわち `θ_m = ‖π‖^{p^{m+1}} < 1`。
★持ち場が疑ったとおり `<` だけでは足りないが、**合同が「どれだけ伸びるか」を与える**ので足りた。

**(3) ★★★しかし `θ^k` は取れない（本波の主結果）。**

`harith` の**上界** `(p−1)·u k ≤ p^{k+1}·e`（`TotallyRamifiedLayer.lean:238` の仮説 (4)）、
`1 ≤ u 0`（同 (1)）、`‖p‖ = ‖π‖^{p^{k+1}e}`（`heM`）**だけ**から

★★**`‖(p : M)‖ · ‖w 0‖^{p−1} ≤ ‖w k‖^{p−1}`**（§5 `norm_pow_ge_of_upper`）

すなわち `‖σ^{p^k}π − π‖` は `‖σπ − π‖` の定数倍**以上**で、その定数 `‖p‖^{1/(p−1)}` は
★**`k` に依らない**。⇒ `‖w k‖ ≤ θ^k ‖w 0‖` を仮定すると `‖p‖ ≤ (θ^k)^{p−1}`
（§5 `theta_pow_ge`）となり、`θ < 1` なら十分大きい `k` で破れる
（§6 `no_uniform_geometric`、§9 `no_uniform_geometric_of_harith`）。
★★**これは出口がすでに持っている仮説だけから出る**ので、仮定を足して逃げられない。

★正確に言うと否定したのは「**`k` に依らない一様な `θ < 1`**」である
（`k` ごとに違う `θ_k` を許す形は否定していない。★測っていないので「無い」とは書かない）。

## ★前波の自分の見込みの訂正

前波の報告で「取れれば `AxEpsilonDecay.lean:294` で `(1/p)^k` が出て `axDecay` に届く」と
書いたが、★**その「取れれば」は取れても届かない**（§9）。
⇒ 前波の `PGroupDescentToAxWild.const_descent_no_go` と合わせて、
**「1 段の損失を深さに応じて絞る」道では `axDecay` に届かない**。
★残るのは `AxWildDescentDecay`（`AxEpsilonDecay.lean:605`、`ε` が増えない形）か、
★**平均化の段で `ε` の伸び `p` を落とす**道である。

## ★木の断定への訂正（他ファイルは書き換えない。名指しでここに書く）

`AxEpsilonDecay.lean:288-292` の docstring は `norm_iterate_pow_sub_self_le` について
「★★これが『`Σ` が収束する』ことの中身であり、`AxTowerDecay.prod_le_rpow_of_geometric` と
対になる」と断定しているが、★**その仮説 `H` が満たされる場面は本波の §9 で塞がれている**
（`harith` を満たす塔では `k` に依らない `θ < 1` の `H` は取れない）。

また同補題の `H : ∀ (m : ℕ) (z : M), …` の ★**`∀ z` は証明で使われていない**
（`AxEpsilonDecay.lean:300` は `H m x` としか呼ばない）。
⇒ §7 に **`x` の 1 点だけで足りる版** `norm_iterate_pow_sub_self_le_of_point` を置いた。

## ★`H` は等長性からは出ない（§8、反例を形式化した）

`f = (z ↦ −z) : ℝ →+ ℝ`、`n = 3`、`x = 1` は**等長な加法自己同型**だが
`‖f^[3^m] x − x‖ = 2`（`m` に依らない）⇒ `H` は `θ ≥ 1` を強いる。
★「反復すれば縮む」は**一般には偽**である（縮ませているのは `p` 進の分岐であって反復ではない）。

## ★①不分岐側との関係

★**本波の否定も不分岐側とは独立である。** 上の否定は `hvalK`（全分岐）を**仮定した上で**
成り立つ（`harith` の (1)(4) と `heM` しか使っていない）。
⇒ `TotallyRamifiedLayer.lean:342 not_valK_of_norm_eq` の穴とは**別の穴**が 3 つ目に増えた
（①不分岐、②一様定数（前波）、③幾何減衰（本波））。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. 本ファイルの主結果（跳びの上界から出る下界）は原典 (Ax 1970 / pGC Cor 3.1) に
   対応する文が無い。`.src` は `AxEpsilonDecay.lean` と同じ項目を指す。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. §1 は純 `ℤ`、§2/§4/§5/§6 は `NormedField` 1 つだけで、
   分岐・付値・Galois の語彙が 1 語も出ない。出るのは docstring だけである。
-/

namespace ABC3.Found.PGC

namespace JumpGeometricDecay

/-! ## §1 抽象核（純 `ℤ`）—— 「真に大きい」＋「差が `d` で割れる」⇒ 「`d` 以上大きい」 -/

section IntCore

/-- ★抽象核: `a < b` と `d ∣ b − a` と `0 < d` から `a + d ≤ b`。

★これが「Hasse–Arf の合同（`p^{m+1} ∣ u(m+1) − u m`）＋ 狭義単調（`u m < u (m+1)`）」
から「跳びが少なくとも `p^{m+1}` 伸びる」を出す唯一の段である。
★分岐・付値・Galois の語彙が 1 語も出ない。 -/
theorem add_le_of_dvd_sub_of_lt {a b d : ℤ} (hd : 0 < d) (hab : a < b) (h : d ∣ b - a) :
    a + d ≤ b := by
  obtain ⟨c, hc⟩ := h
  have hc0 : 1 ≤ c := by
    by_cases h1 : 1 ≤ c
    · exact h1
    · exfalso
      have hle : c ≤ 0 := by omega
      have hdc : d * c ≤ 0 := mul_nonpos_of_nonneg_of_nonpos (le_of_lt hd) hle
      linarith
  nlinarith

end IntCore

/-! ## §2 抽象核（`NormedField` 1 つ）—— 指数の差をノルムの比に直す -/

section NormCore

variable {M : Type*} [NormedField M]

/-- ★抽象核: `‖A‖ = ‖π‖^{a+1}`・`‖B‖ = ‖π‖^{b+1}`・`a + d ≤ b` なら
`‖B‖ ≤ ‖π‖^d · ‖A‖`（`0 < ‖π‖ ≤ 1`）。

★分岐の語彙は 1 語も出ない。 -/
theorem norm_le_zpow_mul {π A B : M} {a b : ℤ} {d : ℤ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≤ 1)
    (hA : ‖A‖ = ‖π‖ ^ (a + 1)) (hB : ‖B‖ = ‖π‖ ^ (b + 1)) (hab : a + d ≤ b) :
    ‖B‖ ≤ ‖π‖ ^ d * ‖A‖ := by
  rw [hA, hB, ← zpow_add₀ (ne_of_gt hπ0)]
  exact zpow_le_zpow_right_of_le_one₀ hπ0 hπ1 (by omega)

end NormCore

/-! ## §3 ★跳びから 1 段の減衰が出る（持ち場の問い (2) への答え） -/

section Jump

variable {M : Type*} [NormedField M] {p : ℕ}

/-- ★★**跳びの 1 段は幾何的に縮む** —— `‖w (m+1)‖ ≤ ‖π‖^{p^{m+1}} · ‖w m‖`。

仮説は `harith` の (2) 狭義単調と (3) Hasse–Arf の合同だけである
（`TotallyRamifiedLayer.lean:238` の `harith`）。
★これが `AxEpsilonDecay.lean:294` の仮説 `H` の形である。
★ただし `θ_m = ‖π‖^{p^{m+1}}` は `m` に依り、一様に取ると `‖π‖^p`
（`k` が大きいと `‖π‖` が 1 に近いので `θ` も 1 に近い）。
★★そのだめなさを定量化したのが §5–§9 である。 -/
theorem norm_step_le_of_jump {π : M} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≤ 1)
    {w : ℕ → M} {u : ℕ → ℤ} {k : ℕ}
    (hw : ∀ j, j ≤ k → ‖w j‖ = ‖π‖ ^ (u j + 1))
    (hult : ∀ m, m < k → u m < u (m + 1))
    (hdvd : ∀ m, m < k → (p : ℤ) ^ (m + 1) ∣ u (m + 1) - u m)
    (hp : 0 < p) {m : ℕ} (hm : m < k) :
    ‖w (m + 1)‖ ≤ ‖π‖ ^ ((p : ℤ) ^ (m + 1)) * ‖w m‖ := by
  have hd0 : (0:ℤ) < (p : ℤ) ^ (m + 1) := by positivity
  have hstep : u m + (p : ℤ) ^ (m + 1) ≤ u (m + 1) :=
    add_le_of_dvd_sub_of_lt hd0 (hult m hm) (hdvd m hm)
  exact norm_le_zpow_mul hπ0 hπ1 (hw m (le_of_lt hm)) (hw (m + 1) hm) hstep

end Jump

/-! ## §4 抽象核 —— 1 段 `θ` 倍を `k` 段で `θ^k` 倍にする（有界版） -/

section Uniform

variable {M : Type*} [NormedField M] {p : ℕ}

/-- 抽象核: 1 段 `θ` 倍なら `j ≤ k` で `‖w j‖ ≤ θ^j · ‖w 0‖`。
★`AxEpsilonDecay.lean:294` の有界版（跳びのデータは `m < k` の範囲にしか無い）。 -/
theorem norm_le_pow_mul_of_step {w : ℕ → M} {θ : ℝ} {k : ℕ} (hθ : 0 ≤ θ)
    (hstep : ∀ m, m < k → ‖w (m + 1)‖ ≤ θ * ‖w m‖) :
    ∀ j, j ≤ k → ‖w j‖ ≤ θ ^ j * ‖w 0‖ := by
  intro j
  induction j with
  | zero => intro _; simp
  | succ m ih =>
      intro hm
      have h1 := ih (by omega)
      have h2 := hstep m (by omega)
      calc ‖w (m + 1)‖ ≤ θ * ‖w m‖ := h2
        _ ≤ θ * (θ ^ m * ‖w 0‖) := mul_le_mul_of_nonneg_left h1 hθ
        _ = θ ^ (m + 1) * ‖w 0‖ := by ring

end Uniform

/-! ## §5 ★★★下界 —— 跳びの**上界**から変位の**下界**が出る -/

section LowerBound

variable {M : Type*} [NormedField M]

/-- ★★★**本ファイルの鍵** —— 跳びの**上界**から変位の**下界**が出る。

`(q : ℤ) * u k ≤ E` と `1 ≤ u 0` と `‖P‖ = ‖π‖^E` から
`‖P‖ · ‖w 0‖^q ≤ ‖w k‖^q`。

★代入すると `q = p−1`、`E = p^{k+1}e`、`P = (p : M)` で、
`harith` の (4) `(p−1)·u k ≤ p^{k+1}·e` と (1) `1 ≤ u 0` と `heM` そのものである。
⇒ ★★`‖w k‖ / ‖w 0‖ ≥ ‖p‖^{1/(p−1)}` で、この下界は **`k` に依らない**。 -/
theorem norm_pow_ge_of_upper {π : M} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≤ 1)
    {w : ℕ → M} {u : ℕ → ℤ} {k q : ℕ} {E : ℤ} {P : M}
    (hw0 : ‖w 0‖ = ‖π‖ ^ (u 0 + 1)) (hwk : ‖w k‖ = ‖π‖ ^ (u k + 1))
    (hu0 : 1 ≤ u 0) (hub : (q : ℤ) * u k ≤ E) (hP : ‖P‖ = ‖π‖ ^ E) :
    ‖P‖ * ‖w 0‖ ^ q ≤ ‖w k‖ ^ q := by
  have hπne : ‖π‖ ≠ 0 := ne_of_gt hπ0
  rw [hw0, hwk, hP, ← zpow_natCast (‖π‖ ^ (u 0 + 1)) q, ← zpow_natCast (‖π‖ ^ (u k + 1)) q,
    ← zpow_mul, ← zpow_mul, ← zpow_add₀ hπne]
  refine zpow_le_zpow_right_of_le_one₀ hπ0 hπ1 ?_
  have hq0 : (0:ℤ) ≤ (q : ℤ) := Int.natCast_nonneg q
  nlinarith

/-- ★下界と幾何減衰をぶつけると `‖P‖ ≤ (θ^k)^q`。 -/
theorem theta_pow_ge {θ : ℝ} {w : ℕ → M} {k q : ℕ} {P : M}
    (hw0 : 0 < ‖w 0‖)
    (hlow : ‖P‖ * ‖w 0‖ ^ q ≤ ‖w k‖ ^ q)
    (hgeo : ‖w k‖ ≤ θ ^ k * ‖w 0‖) :
    ‖P‖ ≤ (θ ^ k) ^ q := by
  have hpow : ‖w k‖ ^ q ≤ (θ ^ k * ‖w 0‖) ^ q :=
    pow_le_pow_left₀ (norm_nonneg _) hgeo q
  have hexp : (θ ^ k * ‖w 0‖) ^ q = (θ ^ k) ^ q * ‖w 0‖ ^ q := mul_pow _ _ q
  have hw0q : (0:ℝ) < ‖w 0‖ ^ q := pow_pos hw0 q
  have := hlow.trans (hpow.trans (le_of_eq hexp))
  exact le_of_mul_le_mul_right (by linarith) hw0q

theorem exists_pow_lt_of_lt_one' {θ c : ℝ} (hθ0 : 0 ≤ θ) (hθ1 : θ < 1) (hc : 0 < c) {q : ℕ}
    (hq : 0 < q) : ∃ k : ℕ, (θ ^ k) ^ q < c := by
  obtain ⟨k, hk⟩ := exists_pow_lt_of_lt_one hc hθ1
  refine ⟨k, ?_⟩
  rcases Nat.eq_zero_or_pos k with hk0 | hk0
  · subst hk0; simpa using hk
  · have h1 : (0:ℝ) ≤ θ ^ k := pow_nonneg hθ0 k
    calc (θ ^ k) ^ q ≤ (θ ^ k) ^ 1 := pow_le_pow_of_le_one h1 (by nlinarith [pow_le_one₀ hθ0 (le_of_lt hθ1) (n := k)]) hq
      _ = θ ^ k := pow_one _
      _ < c := hk

end LowerBound

/-! ## §6 ★★★一様な `θ < 1` は取れない（本波の主結果） -/

section NoGo

variable {M : Type*} [NormedField M]

/-- ★★★**一様な `θ < 1` は取れない**。

`‖P‖ > 0` なら、`‖P‖·‖w 0‖^q ≤ ‖w k‖^q` を満たすどんな `w` も
`‖w k‖ ≤ θ^k · ‖w 0‖` を満たさないような深さ `k` が存在する。 -/
theorem no_uniform_geometric {θ : ℝ} (hθ0 : 0 ≤ θ) (hθ1 : θ < 1) {P : M} (hP : 0 < ‖P‖)
    {q : ℕ} (hq : 0 < q) :
    ∃ k : ℕ, ∀ w : ℕ → M, 0 < ‖w 0‖ → ‖P‖ * ‖w 0‖ ^ q ≤ ‖w k‖ ^ q →
      ¬ (‖w k‖ ≤ θ ^ k * ‖w 0‖) := by
  obtain ⟨k, hk⟩ := exists_pow_lt_of_lt_one' hθ0 hθ1 hP hq
  refine ⟨k, ?_⟩
  intro w hw0 hlow hgeo
  have := theta_pow_ge hw0 hlow hgeo
  linarith

theorem no_uniform_geometric_of_jump {θ : ℝ} (hθ0 : 0 ≤ θ) (hθ1 : θ < 1)
    {π : M} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≤ 1) {P : M} (hP : 0 < ‖P‖) {q : ℕ} (hq : 0 < q) :
    ∃ k : ℕ, ∀ (w : ℕ → M) (u : ℕ → ℤ) (E : ℤ),
      ‖w 0‖ = ‖π‖ ^ (u 0 + 1) → ‖w k‖ = ‖π‖ ^ (u k + 1) → 1 ≤ u 0 →
      (q : ℤ) * u k ≤ E → ‖P‖ = ‖π‖ ^ E → 0 < ‖w 0‖ →
      ¬ (‖w k‖ ≤ θ ^ k * ‖w 0‖) := by
  obtain ⟨k, hk⟩ := no_uniform_geometric (M := M) hθ0 hθ1 hP hq
  refine ⟨k, ?_⟩
  intro w u E hw0 hwk hu0 hub hPE hpos
  exact hk w hpos (norm_pow_ge_of_upper hπ0 hπ1 hw0 hwk hu0 hub hPE)

end NoGo

/-! ## §7 ★訂正 —— `norm_iterate_pow_sub_self_le` の `∀ z` は使われていない -/

section Point

variable {M : Type*} [NormedField M]

/-- ★**訂正版** —— `AxEpsilonDecay.lean:294 norm_iterate_pow_sub_self_le` の
`H : ∀ (m : ℕ) (z : M), …` の **`∀ z` は使われていない**（`:300` は `H m x` だけ）。
★本版は `x` の 1 点だけで仮説を要求する。★既存の宣言は触っていない。 -/
theorem norm_iterate_pow_sub_self_le_of_point (f : M →+ M) (n : ℕ) {θ : ℝ} (hθ : 0 ≤ θ)
    (x : M) (H : ∀ m : ℕ, ‖f^[n ^ (m + 1)] x - x‖ ≤ θ * ‖f^[n ^ m] x - x‖) (k : ℕ) :
    ‖f^[n ^ k] x - x‖ ≤ θ ^ k * ‖f x - x‖ := by
  induction k with
  | zero => simp
  | succ m ih =>
      refine (H m).trans ?_
      calc θ * ‖f^[n ^ m] x - x‖ ≤ θ * (θ ^ m * ‖f x - x‖) :=
            mul_le_mul_of_nonneg_left ih hθ
        _ = θ ^ (m + 1) * ‖f x - x‖ := by ring

end Point

/-! ## §8 ★反例 —— `H` は「等長な加法自己同型」からは出ない -/

section Counter

def negHom : ℝ →+ ℝ where
  toFun := fun z => -z
  map_zero' := neg_zero
  map_add' := fun a b => neg_add a b

theorem iterate_negHom (j : ℕ) (z : ℝ) : negHom^[j] z = (-1) ^ j * z := by
  induction j generalizing z with
  | zero => simp
  | succ m ih =>
      rw [Function.iterate_succ_apply, ih]
      show (-1 : ℝ) ^ m * (-z) = (-1) ^ (m + 1) * z
      ring

theorem negHom_iterate_odd {j : ℕ} (hj : Odd j) (z : ℝ) : negHom^[j] z = -z := by
  rw [iterate_negHom, hj.neg_one_pow]
  ring

/-- ★★**反例** —— 等長な加法自己同型でも `H` は `θ ≥ 1` を強いる。

`f = (z ↦ −z) : ℝ →+ ℝ`、`n = 3`、`x = 1`。`3^m` は奇数なので
`f^[3^m] x − x = −2` で、ノルムは `m` に依らず `2` のままである。
★「反復すれば縮む」は一般には偽である（縮ませているのは分岐であって反復ではない）。 -/
theorem not_forall_theta_lt_one_of_isometry :
    ∃ (f : ℝ →+ ℝ) (x : ℝ), (∀ z : ℝ, ‖f z‖ = ‖z‖) ∧ f x ≠ x ∧
      ∀ θ : ℝ, (∀ m : ℕ, ‖f^[3 ^ (m + 1)] x - x‖ ≤ θ * ‖f^[3 ^ m] x - x‖) → 1 ≤ θ := by
  refine ⟨negHom, 1, ?_, ?_, ?_⟩
  · intro z
    show ‖-z‖ = ‖z‖
    exact norm_neg z
  · show (-1 : ℝ) ≠ 1
    norm_num
  · intro θ H
    have hodd : ∀ m : ℕ, Odd (3 ^ m) := fun m => (by decide : Odd 3).pow
    have hval : ∀ m : ℕ, ‖negHom^[3 ^ m] (1 : ℝ) - 1‖ = 2 := by
      intro m
      rw [negHom_iterate_odd (hodd m)]
      norm_num
    have h0 := H 0
    rw [hval 1, hval 0] at h0
    linarith

end Counter

/-! ## §9 ★★★`harith` の字面への代入（出口が既に持つ仮説だけで塞がる） -/

section Harith

variable {M : Type*} [NormedField M] {p : ℕ}

/-- ★★★★**本波の主結果** —— `harith` の字面への代入。

仮説は `TotallyRamifiedLayer.lean:238`
`exists_norm_sub_algebraMap_le_prod_axDecay_of_totallyRamified` が
★**すでに持っているものだけ**である:
`harith` の (1) `1 ≤ u 0`、(4) `(p−1)·u k ≤ p^{k+1}·e`、そして `heM`。

⇒ ★★★`θ < 1` をどう取っても、その `θ` では `‖w k‖ ≤ θ^k · ‖w 0‖` が
**成り立たない深さ `k`** が存在する。
★否定しているのは「`k` に依らない一様な `θ`」である。 -/
theorem no_uniform_geometric_of_harith (hp : 2 ≤ p)
    {θ : ℝ} (hθ0 : 0 ≤ θ) (hθ1 : θ < 1) {π : M} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≤ 1)
    (hpM : 0 < ‖(p : M)‖) :
    ∃ k : ℕ, ∀ (w : ℕ → M) (u : ℕ → ℤ) (e : ℕ),
      (∀ j, j ≤ k → ‖w j‖ = ‖π‖ ^ (u j + 1)) →
      1 ≤ u 0 →
      ((p : ℤ) - 1) * u k ≤ (p : ℤ) ^ (k + 1) * (e : ℤ) →
      ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e) →
      ¬ (‖w k‖ ≤ θ ^ k * ‖w 0‖) := by
  have hq : 0 < p - 1 := by omega
  obtain ⟨k, hk⟩ := no_uniform_geometric (M := M) (P := (p : M)) hθ0 hθ1 hpM hq
  refine ⟨k, ?_⟩
  intro w u e hw hu0 hub heM hgeo
  have hw0 : ‖w 0‖ = ‖π‖ ^ (u 0 + 1) := hw 0 (Nat.zero_le k)
  have hwk : ‖w k‖ = ‖π‖ ^ (u k + 1) := hw k le_rfl
  have hpos : 0 < ‖w 0‖ := by rw [hw0]; exact zpow_pos hπ0 _
  have hcast : ((p - 1 : ℕ) : ℤ) = (p : ℤ) - 1 := by
    have : (1 : ℕ) ≤ p := by omega
    push_cast [this]
    ring
  have hub' : ((p - 1 : ℕ) : ℤ) * u k ≤ ((p ^ (k + 1) * e : ℕ) : ℤ) := by
    rw [hcast]
    push_cast
    exact hub
  have hP : ‖(p : M)‖ = ‖π‖ ^ ((p ^ (k + 1) * e : ℕ) : ℤ) := by
    rw [heM, zpow_natCast]
  exact hk w hpos (norm_pow_ge_of_upper hπ0 hπ1 hw0 hwk hu0 hub' hP) hgeo

end Harith

/-! ## §10 `.src`（原典の対応箇所）と 使っている公理の一覧 -/

def add_le_of_dvd_sub_of_lt.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_step_le_of_jump.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_pow_ge_of_upper.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def no_uniform_geometric_of_harith.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_iterate_pow_sub_self_le_of_point.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def not_forall_theta_lt_one_of_isometry.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

#print axioms add_le_of_dvd_sub_of_lt
#print axioms norm_le_zpow_mul
#print axioms norm_step_le_of_jump
#print axioms norm_le_pow_mul_of_step
#print axioms norm_pow_ge_of_upper
#print axioms theta_pow_ge
#print axioms exists_pow_lt_of_lt_one'
#print axioms no_uniform_geometric
#print axioms no_uniform_geometric_of_jump
#print axioms norm_iterate_pow_sub_self_le_of_point
#print axioms iterate_negHom
#print axioms not_forall_theta_lt_one_of_isometry
#print axioms no_uniform_geometric_of_harith

end JumpGeometricDecay

end ABC3.Found.PGC
