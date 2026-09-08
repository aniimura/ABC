import ABC3.Found.PGC.JumpStrictMono

/-!
# [pGC] `harith` の (3) Hasse–Arf の合同 —— ★**どこで止まるかの測定**と、抽象核

`JumpFromValueGroup.lean:272-276` の `harith` の第 3 連言

    ∀ m, m < k → (p : ℤ) ^ (m + 1) ∣ u (m + 1) - u m

について測った。★**本ファイルは合同そのものは出していない。**
出したのは ①抽象核(Hasse–Arf の整数性 ⇒ 合同)と ②残る 1 点の**正確な位置**である。

## ★★★測定 1 —— 持ち場が指した「橋が無い」は**2 重に古い**(訂正、名指し)

持ち場の証拠 5 は `HasseArfInduction.lean:101-108` の

> 2. ★★`[Algebra A C]` —— **まだ無い**。`C = B^H` は `Subring B` なので
>    `Algebra C B` は付くが、`Algebra A C`(底 `𝒪 ⊆ 𝒪_{K″}`)は木にインスタンスが無い。

を指していた。★**測ったら、どちらも既に埋まっていた**:

1. ★`Found/PGC/FixedRingBaseAlgebra.lean:178` に `fixedRingAlgebra` が在る
   (`Algebra A ↥(fixedRing B H)`、Y17 の `Subring.algebraOfMapsTo` の再利用)。
   同ファイル冒頭は「★★**3 本とも埋まった**」「★**残るのは `hind` だけ**になった」と書いている。
2. ★その `hind`(帰納)も `Found/PGC/HasseArfStrongInduction.lean` が供給しており、
   同 `:447` の

       exists_natCast_herbrandPhiGroup_of_lowerRamificationGroup_one_eq_top
         (habel : ∀ x y : G, x * y = y * x) (h1 : lowerRamificationGroup B G 1 = ⊤)
         (hne : lowerRamificationGroup B G n ≠ lowerRamificationGroup B G (n+1)) :
         ∃ j : ℕ, herbrandPhiGroup G π' (n : ℝ) = (j : ℝ)

   は★**`hind` を持っていない**(同ファイル 106 行「★★`hind` は §3 が供給するので、
   ここに帰納法は無い」)。⇒ ★★**`G` 可換・`G_1 = ⊤` の場合の Hasse–Arf は木で閉じている。**

★★**我々の設定はちょうどその場合である**: `M/K` は全分岐で `G = ⟨g⟩` は位数 `p^{k+1}`。
全分岐なので `G_0 = G`、`G` は `p` 群なので野生慣性 `G_1` も `G` 全体、すなわち `G_1 = ⊤`。
可換も明らか。★`HasseArfStrongInduction.lean:106` が「まだ閉じていない」と言う段 3
(順分岐商 `G_1 ⊊ G_0`)は★**我々には要らない**。

★測ったコマンド:
```
sed -n '95,135p' lean/ABC3/Found/PGC/HasseArfInduction.lean          # 「まだ無い」の字面
grep -n "^theorem\|^def" lean/ABC3/Found/PGC/FixedRingBaseAlgebra.lean # fixedRingAlgebra:178
sed -n '440,470p' lean/ABC3/Found/PGC/HasseArfStrongInduction.lean     # hind を持たない形
```
★元の docstring は他が読んでいるので**直していない**。ここに訂正として書く。

## ★★測定 2 —— 残る 1 点は「**環の言葉と体・ノルムの言葉の橋**」だけ

木の Hasse–Arf は `herbrandPhiGroup G π' (n : ℝ) = (j : ℝ)`(`B` は DVR、`G` は
`MulSemiringAction G B`)で述べられている。我々の `u m` は
`‖(g^{p^m}) π − π‖ = ‖π‖^{u m + 1}`(`M` は `NormedField`)で定義されている。

★**この 2 つを繋ぐには**

    hrec : ∀ m, φ (m+1) = φ m + (u (m+1) − u m) / p^{m+1}

が要る(`φ m := herbrandPhiGroup G π' (u m)`)。★これは Yoshida08 Lemma 6.10 (i)
`φ_G(n) = (1/|G|) Σ_{i=1}^{n} |G_i|`(`HerbrandComposition.lean:449`)と
「`|G_i| = p^{k+1−l}`(`u_{l−1} < i ≤ u_l`)」から**計算で出る**はずだが、
★★**本ファイルはそれを測っていない**(`lowerRamificationGroup` の位数を
ノルムの言葉の `u m` から決める部分が未着手)。★ここが止まった場所である。

★§1 はその `hrec` さえ来れば合同が出ることを示す(＝残りは `hrec` **1 本**)。

★**橋の片側は測れた**: `LowerRamificationGroup.lean:270` の

    σ ∈ lowerRamificationGroup B G n ↔ ∀ x : B, σ • x - x ∈ (maximalIdeal B) ^ (n + 1)

は、ノルムの言葉では `∀ x ∈ 𝒪_M, ‖σ x − x‖ ≤ ‖π‖^{n+1}` である。
★前波の `JumpMono.norm_sub_apply_le_mul`(`‖h z − z‖ ≤ ‖z‖·‖π‖^t`)は
`‖z‖ ≤ 1` で `‖h z − z‖ ≤ ‖π‖^t` を出すので、★**`h ∈ G_{t−1}` までしか出ない**
(真は `h ∈ G_t`)。★1 つずれるのは、単元 `z`(`v(z) = 0`)のとき
`l = 0` の項が消えることを使っていないからである。
⇒ ★**橋の残りは「単元に対する 1 つ分の改良」と、`|G_i|` を `u` で書き下す部分**である。
★本ファイルはどちらも**やっていない**。

## ★★測定 3 —— 前波までの超距離だけの道では (3) は**出ない**(形式化した)

§2 `congruence_not_implied_by_ultrametric` を見ること。
`p = 3, k = 1, e = 1, u₀ = 2, u₁ = 4` は

* `1 ≤ u₀`(前波 (1))、`u₀ < u₁`(前波 `JumpMono.jump_lt_succ` = (2))、
* `(p−1)·u₁ ≤ p^{k+1}·e`(前波 `WildBreakUpper` = (4))

を**すべて満たすのに** `p ∣ u₁ − u₀` は**偽**(`3 ∤ 2`)。
⇒ ★★**(1)(2)(4) から (3) は出ない。**(3) は独立の入力であり、
その入力が Hasse–Arf である。★先行記録(`JumpFromValueGroup.lean` §6 の
`propagated_bound_strictly_weaker` / `chain_upper_false`、
「打ち消しは超距離だけでは見えない」)と同じ結論に別の道で到達した。

## 節の構成

* §1 ★抽象核 —— `φ` が各段で整数 ⇒ `p^{m+1} ∣ u(m+1) − u m`。★分岐・付値・群の語彙が 1 語も出ない。
* §2 ★測定 3 の形式化(独立性)。
-/

namespace ABC3.Found.PGC

namespace HasseArfCongruence

/-! ## §1 抽象核 —— Herbrand の `φ` が整数なら跳びは合同 -/

section Core

/-- ★★★**抽象核** —— `φ` が各段で整数値をとり、`φ(m+1) = φ(m) + (u(m+1) − u m)/p^{m+1}`
という Herbrand の漸化式を満たすなら `p^{m+1} ∣ u(m+1) − u m`。

★これが `harith` の第 3 連言の**中身のすべて**である。
★分岐・付値・Galois・群の語彙が 1 語も出てこない(`p` は単に `0` でない実数)。 -/
theorem dvd_sub_of_phi_intCast {p : ℕ} (hp : p ≠ 0) (u : ℕ → ℤ) (φ : ℕ → ℝ) (j : ℕ → ℤ)
    (hint : ∀ m, φ m = (j m : ℝ))
    (hrec : ∀ m, φ (m + 1) = φ m + ((u (m + 1) - u m : ℤ) : ℝ) / ((p : ℝ) ^ (m + 1)))
    (m : ℕ) : ((p : ℤ) ^ (m + 1)) ∣ (u (m + 1) - u m) := by
  have hp0 : ((p : ℝ)) ≠ 0 := Nat.cast_ne_zero.mpr hp
  have hpow : ((p : ℝ)) ^ (m + 1) ≠ 0 := pow_ne_zero _ hp0
  refine ⟨j (m + 1) - j m, ?_⟩
  have h1 := hrec m
  rw [hint (m + 1), hint m] at h1
  have hd : ((u (m + 1) - u m : ℤ) : ℝ) / ((p : ℝ)) ^ (m + 1)
      = ((j (m + 1) : ℤ) : ℝ) - ((j m : ℤ) : ℝ) := by linarith [h1]
  have hX : ((u (m + 1) - u m : ℤ) : ℝ)
      = (((j (m + 1) : ℤ) : ℝ) - ((j m : ℤ) : ℝ)) * ((p : ℝ)) ^ (m + 1) :=
    (div_eq_iff hpow).mp hd
  have hZ : ((u (m + 1) - u m : ℤ) : ℝ)
      = (((p : ℤ) ^ (m + 1) * (j (m + 1) - j m) : ℤ) : ℝ) := by
    rw [hX]; push_cast; ring
  exact_mod_cast hZ


end Core

/-! ## §2 ★測定 —— (1)(2)(4) から (3) は出ない -/

section Independence

/-- ★★★**独立性(測定 3)** —— `p = 3`, `k = 1`, `e = 1`, `u₀ = 2`, `u₁ = 4` は
前波までに**定理として得た 3 条件**をすべて満たすが、合同 (3) を満たさない。

⇒ ★**(3) は (1)(2)(4) から導けない**。独立の入力(Hasse–Arf)が要る。 -/
theorem congruence_not_implied_by_ultrametric :
    (1 : ℤ) ≤ 2 ∧ (2 : ℤ) < 4 ∧ ((3 : ℤ) - 1) * 4 ≤ (3 : ℤ) ^ (1 + 1) * 1 ∧
      ¬ ((3 : ℤ) ^ (0 + 1) ∣ (4 - 2)) := by
  refine ⟨by norm_num, by norm_num, by norm_num, ?_⟩
  decide

/-- ★同じ組は、前波 §3 の**定量形**(`u₀ + min E u₀ ≤ u₁`, `E = p^{k+1}e = 9`)も満たす。

★すなわち「`h^pπ/π − 1` の 1 次の項の評価」を最大限使っても (3) には届かない。 -/
theorem congruence_not_implied_even_with_quantitative :
    (2 : ℤ) + min 9 2 ≤ 4 ∧ ¬ ((3 : ℤ) ∣ (4 - 2)) := by
  refine ⟨by norm_num, ?_⟩
  decide

end Independence

/-! ## §3 `.src` と 公理 -/

def dvd_sub_of_phi_intCast.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

#print axioms dvd_sub_of_phi_intCast
#print axioms congruence_not_implied_by_ultrametric
#print axioms congruence_not_implied_even_with_quantitative

end HasseArfCongruence

end ABC3.Found.PGC
