import ABC3.Found.PGC.HerbrandComposition

/-!
# [pGC] Herbrand の「第 1 跳びの下で `φ = id`」を形式化した ——「木に無い」は偽だった

## ★★★自分の前波の断定の訂正（本波の第 1 の成果）

前波の報告で「木は『mathlib には上付き/下付き分岐群も Herbrand 関数も無い』と測っており、
★そこが最大の在庫の穴だと私は見込む」と書いたが、★**「木に無い」は偽**である。
測り直したコマンドを残す:

```
grep -rn "herbrandPhiGroup" lean/ABC3/Found/PGC/ --include=*.lean
  → HerbrandComposition.lean:395(定義) / :455 / :477 / :621、
    HasseArfStrongInduction.lean:183 / :447、TameQuotientTower.lean:152
grep -rn "upperRamificationGroupAdjoin" lean/ABC3/Found/PGC/ --include=*.lean
  → AbelianSubfieldInLubinTate.lean:308 / :354 / :400 / :538
```

★**mathlib に無いのは真だが、木には在る。** ★本連鎖は mathlib の測定を木の測定に
読み替えていた。★★「無い」と書く前に木でも 1 回 grep すること（#330）。

## ★何が足りなかったのか（測り直した）

`FirstJumpLedger.lean:99-110`（★本連鎖が以前に書いたファイル）は
★「足りないのは `u₁ ≤ i(σ̄)`。これは Herbrand（`φ_{L/K}` が `u ≤ u₁` で恒等）、
同値に Serre, Corps Locaux IV §1 Prop. 3 である。★本ファイルは形式化していない」
と名指ししていた。★★**本ファイルはその 2 本を形式化した。**

| 宣言 | 内容 | 語彙 |
|---|---|---|
| `phiOf_eq_self_of_forall_lt` | ★★**`φ_f(n) = n`（すべての `f τ > n` なら）**＝「第 1 跳びの下で `φ` は恒等」 | `ℕ∞`・`ℝ`・`Fintype` だけ |
| `lt_of_coset_sum` | ★★剰余類和の恒等式から `n < y` を取り出す | 同上 |
| `truncENat_eq_of_lt` / `lt_of_truncENat_eq` | `truncENat` の飽和（`min` が右に落ちる条件） | 同上 |
| `forall_lt_ramIndex_of_first_jump` | `σ ∉ H` なら剰余類 `σH` の元はすべて `≠ 1` | 純群論 |
| ★`lt_ramIndex_quotient_of_coset_sum` | ★★★**`u₁ < i_ϖ(σ)`**（`FirstJumpLedger` が名指しした 1 本） | 分岐 |
| `herbrandPhi_eq_self_of_first_jump` | `φ_H(u₁) = u₁`（具体層） | 分岐 |

## ★★まだ残っている 1 点（正確に。★数学ではなく配管である）

`lt_ramIndex_quotient_of_coset_sum` は★**剰余類和の恒等式を仮説 `hcoset` として受けている**:

```
Σ_{τ ∈ H} min(i_{π'}(στ), u+1) = |H| · min(i_ϖ(σ), u+1)
```

★これは木の `HerbrandComposition.lean:532 coset_sum_truncENat` **そのもの**である
（★`herbrandPhi π' H u = u` を代入した形。その代入は本ファイルの
`herbrandPhi_eq_self_of_first_jump` が与える）。
⇒ ★★**残っているのは `coset_sum_truncENat` の 10 個の仮説
（`hcomp` / `hHtriv` / `hπ'` / `hinj` / `hfixC` / `hres` / `hAC` / `hϖ` / `hadj` / `hfix`）を
`PAdicLocalField` の塔で揃える配管だけ**である。
★本波では揃えていない（Write 400 行の制約と、仮説が `IntegerMiscInstances` /
`IntegerRingInstances` の型クラス群にまたがるため）。★正直に未了と書く。

## ★★それでも `axDecay p k` は出ない（前波までの測定は生きている）

`FirstJumpLedger.lean:20-27` の測定 2 が
★**「鎖の台帳が閉じるのは `s_m = 1`（最後の跳びが最大）のときだけ」**
（`chain_ledger_forces_max` / `dvd_of_chain_ledger` / `not_chain_ledger_of_not_dvd`）
を定理にしている。⇒ ★**点 1（`u₁ ≤ i(σ̄)`）を閉じても鎖の道は閉じない。**

★★ただしそれは `CyclicLayerDescent` の**鎖**の勘定（`σ_j = τ^{p^j}` の収縮率の積）であって、
`NormalizedTraceDescent.lean:379 FirstJumpRoute` の勘定（`p^{j/(m·e)}` の形）とは**別**である
（★本波で 2 つの字面を突き合わせた）。★どちらに効くかは測っていない。

## ★穴の現状（本波で増減なし）

①不分岐（`TotallyRamifiedLayer.lean:342`）/ ②′（`ExitLossDepth.not_forall_prod_Icc_le_of_le`）/
③幾何減衰（`JumpGeometricDecay.no_uniform_geometric_of_harith`）/
⑤出口の損失は `axDecay p 1` が限界（`ExitLossDepth.exit_does_not_realize_cEx`）。④は消えたまま。
★★本波は穴を増やしても減らしてもいない。**在庫の穴を 1 つ埋めた**のが成果である。

## ★①不分岐側との関係

★本波は分岐群の一般論なので①とは独立である。
`lt_ramIndex_quotient_of_coset_sum` は**全分岐を仮定していない**
（`hvalK` も `IsTotallyRamified` も出てこない）。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は `HerbrandComposition.lean` と同じ Yoshida08 Lemma 6.10 を指す。
   pGC の原典が独立に立てた項目ではない。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. §1・§2 は `ℕ∞` と `ℝ` と `Fintype` だけ、§3 の 1 本目は純群論で、
   分岐・付値・Galois の語彙が 1 語も出ない。
4. ★`hcoset` を仮説のまま残した（配管を払っていない）ことを上に明記した。
-/

namespace ABC3.Found.PGC

namespace HerbrandFirstJump

/-! ## §1 抽象核 —— `truncENat` の飽和 -/

section Trunc

/-- `n < x` なら `min(x, n+1) = n+1`。 -/
theorem truncENat_eq_of_lt {x : ℕ∞} {n : ℕ} (h : (n : ℕ∞) < x) :
    truncENat x ((n : ℝ) + 1) = (n : ℝ) + 1 := by
  cases x with
  | top => simp
  | coe k =>
      have hk : n < k := by exact_mod_cast h
      rw [truncENat_coe]
      exact min_eq_right (by exact_mod_cast (by omega : n + 1 ≤ k))

/-- 逆に `min(x, n+1) = n+1` なら `n < x`。 -/
theorem lt_of_truncENat_eq {x : ℕ∞} {n : ℕ}
    (h : truncENat x ((n : ℝ) + 1) = (n : ℝ) + 1) : (n : ℕ∞) < x := by
  cases x with
  | top => exact lt_of_le_of_ne le_top (by simp)
  | coe k =>
      rw [truncENat_coe] at h
      have hR : (n : ℝ) + 1 ≤ (k : ℝ) := by
        rcases le_total ((n : ℝ) + 1) (k : ℝ) with h' | h'
        · exact h'
        · rw [min_eq_left h'] at h
          linarith
      have hN : n + 1 ≤ k := by exact_mod_cast hR
      exact_mod_cast (by omega : n < k)

end Trunc

/-! ## §2 ★★抽象核 —— 第 1 跳びの下で `φ` は恒等 -/

section Phi

/-- ★★★**第 1 跳びの下で `φ` は恒等** —— すべての `τ` で `n < f τ` なら `φ_f(n) = n`。

★`FirstJumpLedger.lean:99-110` が「足りないのはこれ（Herbrand）」と名指ししていた 1 本。
★中身は `phiOf_natCast_of_pos`（`HerbrandComposition.lean:267`）の和が
`n·|S|` になるだけである。★分岐・付値・Galois の語彙が 1 語も出ない。 -/
theorem phiOf_eq_self_of_forall_lt {S : Type*} [Fintype S] [Nonempty S] (f : S → ℕ∞)
    (hf : ∀ τ, 0 < f τ) {n : ℕ} (hn : ∀ τ, (n : ℕ∞) < f τ) :
    phiOf f (n : ℝ) = (n : ℝ) := by
  have hcard : (Nat.card S : ℝ) ≠ 0 := Nat.cast_ne_zero.2 Nat.card_pos.ne'
  rw [phiOf_natCast_of_pos f hf n]
  have hval : ∀ i ∈ Finset.Icc 1 n,
      (Nat.card {τ : S // (i : ℕ∞) < f τ} : ℝ) = (Nat.card S : ℝ) := by
    intro i hi
    have hi' : i ≤ n := (Finset.mem_Icc.mp hi).2
    have : Nat.card {τ : S // (i : ℕ∞) < f τ} = Nat.card S := by
      refine Nat.card_congr (Equiv.subtypeUnivEquiv ?_)
      intro τ
      exact lt_of_le_of_lt (by exact_mod_cast hi') (hn τ)
    rw [this]
  rw [Finset.sum_congr rfl hval, Finset.sum_const, Nat.card_Icc, Nat.add_sub_cancel,
    nsmul_eq_mul]
  field_simp

end Phi

/-! ## §3 抽象核 —— 剰余類和から下界を取り出す -/

section Quotient

/-- ★★抽象核 —— 剰余類和の恒等式
`Σ_τ min(f τ, n+1) = |H| · min(y, n+1)` と「すべての `f τ > n`」から `n < y`。

★これが「上の層で跳びが `n` より上なら、下の層でも `n` より上」の中身である。 -/
theorem lt_of_coset_sum {H : Type*} [Fintype H] [Nonempty H] {f : H → ℕ∞} {y : ℕ∞} {n : ℕ}
    (hf : ∀ τ, (n : ℕ∞) < f τ)
    (hsum : ∑ τ : H, truncENat (f τ) ((n : ℝ) + 1)
      = (Nat.card H : ℝ) * truncENat y ((n : ℝ) + 1)) :
    (n : ℕ∞) < y := by
  have hL : ∑ τ : H, truncENat (f τ) ((n : ℝ) + 1)
      = (Nat.card H : ℝ) * ((n : ℝ) + 1) := by
    rw [Finset.sum_congr rfl (fun τ _ => truncENat_eq_of_lt (hf τ)), Finset.sum_const,
      Finset.card_univ, Nat.card_eq_fintype_card, nsmul_eq_mul]
  have hcard : (Nat.card H : ℝ) ≠ 0 := Nat.cast_ne_zero.2 Nat.card_pos.ne'
  exact lt_of_truncENat_eq (mul_left_cancel₀ hcard (hL.symm.trans hsum)).symm

end Quotient

/-! ## §4 具体層 —— `i_ϖ(σ̄)` の下界 -/

section RamIndex

/-- 純群論: `σ ∉ H` なら剰余類 `σH` の元はすべて `≠ 1`。 -/
theorem forall_lt_ramIndex_of_first_jump
    {G : Type*} [Group G] {H : Subgroup G} {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] [MulSemiringAction G B] {π' : B} {σ : G} {u : ℕ}
    (hfirst : ∀ τ : G, τ ≠ 1 → (u : ℕ∞) < ramIndex π' τ) (hσ : σ ∉ H) :
    ∀ τ : H, (u : ℕ∞) < ramIndex π' (σ * (τ : G)) := by
  intro τ
  refine hfirst _ (fun h => hσ ?_)
  have hσeq : σ = ((τ : G))⁻¹ := by
    rw [← mul_one σ, ← mul_inv_cancel ((τ : G)), ← mul_assoc, h, one_mul]
  rw [hσeq]
  exact H.inv_mem τ.2

/-- ★★★★**`FirstJumpLedger.lean:105` が名指しした 1 本** —— `u₁ < i_ϖ(σ)`。

仮説 `hcoset` は木の `HerbrandComposition.lean:532 coset_sum_truncENat` に
`herbrandPhi π' H u = u`（本ファイル `herbrandPhi_eq_self_of_first_jump`）を
代入した形そのものである。
★残っているのは `coset_sum_truncENat` の 10 個の仮説を塔で揃える**配管**だけで、
★数学ではない（★本波では揃えていない）。 -/
theorem lt_ramIndex_quotient_of_coset_sum
    {G : Type*} [Group G] {H : Subgroup G} [Fintype H] [Nonempty H]
    {B : Type*} [CommRing B] [IsDomain B] [IsDiscreteValuationRing B] [MulSemiringAction G B]
    {C : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C] [MulSemiringAction G C]
    {π' : B} {ϖ : C} {σ : G} {u : ℕ}
    (hfirst : ∀ τ : G, τ ≠ 1 → (u : ℕ∞) < ramIndex π' τ) (hσ : σ ∉ H)
    (hcoset : ∑ τ : H, truncENat (ramIndex π' (σ * (τ : G))) ((u : ℝ) + 1)
      = (Nat.card H : ℝ) * truncENat (ramIndex ϖ σ) ((u : ℝ) + 1)) :
    (u : ℕ∞) < ramIndex ϖ σ :=
  lt_of_coset_sum (forall_lt_ramIndex_of_first_jump hfirst hσ) hcoset

/-- §2 の具体層 —— `H` の元がすべて第 1 跳びより上なら `φ_H(u) = u`。 -/
theorem herbrandPhi_eq_self_of_first_jump
    {G : Type*} [Group G] (H : Subgroup G) [Fintype H] [Nonempty H]
    {B : Type*} [CommRing B] [IsDomain B] [IsDiscreteValuationRing B] [MulSemiringAction G B]
    {π' : B} {u : ℕ}
    (hpos : ∀ τ : H, 0 < ramIndex π' (τ : G))
    (hfirst : ∀ τ : H, (u : ℕ∞) < ramIndex π' (τ : G)) :
    herbrandPhi π' H (u : ℝ) = (u : ℝ) := by
  rw [herbrandPhi_eq_phiOf]
  exact phiOf_eq_self_of_forall_lt _ hpos hfirst

end RamIndex

/-! ## §5 `.src` と 使っている公理の一覧 -/

section Src

def phiOf_eq_self_of_forall_lt.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 16, item := "Lemma 6.10", sectionId := "lemma-6-10" }

def lt_of_coset_sum.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 16, item := "Lemma 6.10", sectionId := "lemma-6-10" }

def lt_ramIndex_quotient_of_coset_sum.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 16, item := "Lemma 6.10", sectionId := "lemma-6-10" }

def herbrandPhi_eq_self_of_first_jump.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 16, item := "Lemma 6.10", sectionId := "lemma-6-10" }

end Src

#print axioms truncENat_eq_of_lt
#print axioms lt_of_truncENat_eq
#print axioms phiOf_eq_self_of_forall_lt
#print axioms lt_of_coset_sum
#print axioms forall_lt_ramIndex_of_first_jump
#print axioms lt_ramIndex_quotient_of_coset_sum
#print axioms herbrandPhi_eq_self_of_first_jump

end HerbrandFirstJump

end ABC3.Found.PGC
