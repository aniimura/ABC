import ABC3.Found.PGC.PowerSpanTop

/-!
# [pGC] `hvalK` は木に**在った** —— 底 1 本から中間層の全分岐が出る

## 持ち場（前波で「残るのは 1 つだけ」と書いた点）

`hvalK`（`E₁ ⊂ L` の値群が `‖π‖^{pℤ}` に入る）。

## ★★開いたら**木に在った**（★索引ではなく本文を読んで見つけた）

`TotallyRamifiedValueGroup.lean` §4 の

> `exists_zpow_norm_intermediate` —— ★★★★**主定理** —— 底 `K` の値群が `‖π‖^{n·ℤ}`
> (`n = [M:K]`)に入るなら、**どの中間層 `E` の値群も** `‖π‖^{[M:E]·ℤ}` に入る。
> ★★これが `hvalj`(層ごとの全分岐)を**底の 1 本**に落とす中身である。

★**この結論が私の `hvalK` そのもの**である（`q = [L:E₁] = p`）。★docstring の断定を検算した
（5 例目）: 結論の形 `∃ m, ‖algebraMap E M d‖ = ‖π‖^{q·m}` は私の `hvalK` と**一字一句同じ**。
⇒ §1 は 1 行の代入で済んだ。

★**見つけ方**: 索引ではなく `TotallyRamifiedValueGroup.lean` を**上から読んだ**。
`exists_zpow_norm`（§2）を探しに行って、その先の §3・§4 に気づいた。
★#354 の教訓（名前で引いて 0 件なら部品で引く）の延長で、★**ファイルを読むと 2 段先が見える**。

## ★★循環を切った（§3）

素朴に繋ぐと回らない:

* `valK_of_base_layer` は `finrank E M = q` を要求する
* `PowerSpanTop.finrank_eq_of_zeta_pow`（前波）は `hvalK` を要求する

⇒ ★`q` を**未知の層の次数のまま** `valK_of_base_layer` を回し、
★**`q ∣ p`（`E₁` に付値 `p` の元 `μ` があるから）＋ `q ≠ 1`（`L ≠ E₁`）＋ `p` 素数**
で `q = p` を出す（§3 `finrank_eq_of_dvd`）。★`‖μ‖ = ‖π‖^p` は
`ZetaStepRatio.norm_zeta_step`（自分の定理）である。

## 何を埋めたか

| 節 | 宣言 | 内容 |
|---|---|---|
| §1 | `valK_of_base_layer` | ★★★`hvalK` の供給（木の主定理への代入） |
| §2 | `valbase_of_uniformizer` | 底の条件を「素元の冪」に落とす（3 行） |
| §3 | `finrank_eq_of_dvd` | ★循環を切る（`q ∣ p` ＋ `q ≠ 1` ＋ 素数性） |

## ★★★残った入力（★今日いちばん短い表になった）

| 入力 | 出どころ | 状態 |
|---|---|---|
| `‖π‖^n = ‖p‖`（`n = φ(p^{m+2})`） | ★`ZetaSubOnePrime.norm_zeta_sub_one_pow`（自分） | ✓ **在る** |
| `‖μ‖ = ‖π‖^p` | ★`ZetaStepRatio.norm_zeta_step`（自分） | ✓ **在る** |
| `π` の次数 `n` のモニック関係式 | ★円分多項式（`ZetaSubOnePrime` が Eisenstein で使った） | ✓ 材料は在る |
| `Algebra.adjoin ℚ_p {π} = ⊤`、`L ≠ E₁` | 塔の定義 | 定義 |
| ★**`ℚ_p` の値群が `‖p‖^ℤ`** | mathlib の `Padic` | ★**未着手（唯一の外部入力）** |
| `finrank ℚ_p L = n` | ★`PowerSpanTop`（前波）＋ 上の 3 つ | ★導ける（未実装） |

★**`finrank ℚ_p L = n` も `PowerSpanTop.finrank_eq_of_relation` で出る**（`hvalbase` は
`finrank` を要求しないので循環しない）。⇒ ★残る**外部**入力は
★**「`ℚ_p` の値群が `‖p‖^ℤ`」の 1 つだけ**である。

## 逸脱の記録

- §1 は木の定理への代入なので、★**新しい数学は 0**。★価値は「在ることを見つけた」ことにある。
- §3 の `q ≠ 1` は `L ≠ E₁`（`ζ ∉ E₁`）と同値。★塔の定義から明らかだが、
  Lean では `finrank = 1 ↔ 全射` の橋が要る（★未実装）。
-/

namespace ABC3.Found.PGC

namespace ValKFromBase

/-! ## §1 ★`hvalK` は木に**在った** —— 底 1 本から中間層が出る -/

section Layer

/-- ★★★**`hvalK` の供給**。★木の `TotallyRamified.exists_zpow_norm_intermediate`
（`TotallyRamifiedValueGroup.lean:…` §4）にそのまま代入するだけ。

★あちらの docstring は「★★これが `hvalj`（層ごとの全分岐）を**底の 1 本**に落とす中身である」
と言っている。★**本当だった** —— 私の `hvalK`（`E₁ ⊂ L` の値群が `‖π‖^{pℤ}` に入る）は
まさにこの定理の結論そのものである。 -/
theorem valK_of_base_layer {K E M : Type*} [Field K] [Field E] [NormedField M]
    [IsUltrametricDist M] [Algebra K M] [Algebra K E] [Algebra E M] [IsScalarTower K E M]
    [FiniteDimensional K M] [FiniteDimensional E M] [FiniteDimensional K E]
    {π : M} {n p : ℕ} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≠ 1)
    (hn : Module.finrank K M = n) (hq : Module.finrank E M = p)
    (hvalbase : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m)) :
    ∀ a : E, a ≠ 0 → ∃ m : ℤ, ‖algebraMap E M a‖ = ‖π‖ ^ ((p : ℤ) * m) :=
  fun a ha => TotallyRamified.exists_zpow_norm_intermediate hπ0 hπ1 hn hq hvalbase a ha

end Layer

/-! ## §2 底の条件を「素元の冪」に落とす -/

section Base

variable {K M : Type*} [Field K] [NormedField M] [Algebra K M]

/-- ★底 `K` の値群が `‖P‖^ℤ`（`P` は `K` の素元）で、`‖π‖^n = ‖P‖` なら、
`hvalbase` の形になる。★`‖π‖^n = ‖P‖` は `ZetaSubOnePrime.norm_zeta_sub_one_pow`
（`n = φ(p^{m+1})`、`P = p`）そのもの。 -/
theorem valbase_of_uniformizer {π : M} {P : K} {n : ℕ}
    (hpn : ‖π‖ ^ n = ‖algebraMap K M P‖)
    (hK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖algebraMap K M P‖ ^ m) :
    ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m) := by
  intro a ha
  obtain ⟨m, hm⟩ := hK a ha
  exact ⟨m, by rw [hm, ← hpn, ← zpow_natCast, ← zpow_mul]⟩

end Base

/-! ## §3 ★循環を切る —— `hfr` は `hvalK` から**素数性**で出る -/

section NoCycle

/-- ★★`hvalK`（層の次数 `q` を法とする形）と「`E` に付値 `p` の元がある」から `q ∣ p`。

★これが循環を切る鍵である —— `PowerSpanTop.finrank_eq_of_zeta_pow` は `hvalK` を要求し、
`valK_of_base_layer` は `finrank = q` を要求するので、そのままでは回らない。
★`q` を「未知の層の次数」のまま `valK_of_base_layer` を回し、
★**`q ∣ p` と `q ≠ 1` と `p` の素数性**で `q = p` を出す。 -/
theorem finrank_eq_of_dvd {E M : Type*} [Field E] [NormedField M] [Algebra E M]
    {π : M} {q p : ℕ} (hp : p.Prime) (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≠ 1) (hq1 : q ≠ 1)
    (hvalq : ∀ a : E, a ≠ 0 → ∃ m : ℤ, ‖algebraMap E M a‖ = ‖π‖ ^ ((q : ℤ) * m))
    {μ : E} (hμ0 : μ ≠ 0) (hμ : ‖algebraMap E M μ‖ = ‖π‖ ^ p) :
    q = p := by
  obtain ⟨m, hm⟩ := hvalq μ hμ0
  rw [hμ, ← zpow_natCast ‖π‖ p] at hm
  have heq : ((p : ℕ) : ℤ) = (q : ℤ) * m := (zpow_right_inj₀ hπ0 hπ1).mp hm
  have hdvd : (q : ℤ) ∣ ((p : ℕ) : ℤ) := ⟨m, heq⟩
  have hdvdn : q ∣ p := by exact_mod_cast hdvd
  rcases (Nat.Prime.eq_one_or_self_of_dvd hp q hdvdn) with h | h
  · exact absurd h hq1
  · exact h

end NoCycle

/-! ## §4 使っている公理の一覧 -/

#print axioms valK_of_base_layer
#print axioms valbase_of_uniformizer
#print axioms finrank_eq_of_dvd

end ValKFromBase

end ABC3.Found.PGC
