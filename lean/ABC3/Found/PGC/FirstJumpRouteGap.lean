import ABC3.Found.PGC.ExitLossDepth

/-!
# [pGC] 2 波連続の宿題に答えた —— 点 1 は `FirstJumpRoute` に効く。しかし⑤は消えない

## ★★★問い 1（2 波連続の宿題）: 鎖の勘定と `FirstJumpRoute` のどちらに効くか

★**答え: `FirstJumpRoute` の側に効く。鎖の側には効かない。**

* **鎖の勘定**（`CyclicLayerDescent` の `σ_j = τ^{p^j}` の収縮率の積）: ★効かない。
  `FirstJumpLedger.lean:20-27` の測定 2（`chain_ledger_forces_max` /
  `dvd_of_chain_ledger` / `not_chain_ledger_of_not_dvd`）が
  「閉じるのは `s_m = 1`（最後の跳びが最大）のときだけ」を定理にしており、
  ★これは点 1 とは独立の障害である。
* **`FirstJumpRoute` の勘定**: ★★**効く。前波で閉じた点 1 が `hjump` そのものになる**（問い 2）。

## ★★問い 2: `p` の差（`p·e_K` vs `e k`）は消えるか → ★消える

前波で閉じた `CosetSumFixedRing.sub_one_mul_first_jump_le_of_quotient` は
**`(p−1)·u₁ ≤ p·e_K`**。`FirstJumpRoute` の `hjump` は **`(p−1)·j k ≤ e k`** で、
`NormalizedTraceDescent.lean:214-221` の定義により `e = e_{E₁}`（★**下の層** `E₁` の
絶対分岐指数）である。★`E₁/K` は次数 `p` の完全分岐なので `e_{E₁} = p·e_K`。
⇒ ★★**`p·e_K = e_{E₁} = e k` で差はちょうど消える**（§1 `hjump_of_sub_one_mul_le`）。
★本体が指摘した「右辺に `p` の差がある」は、★**`E₁/K` の分岐指数そのもの**だった。

## ★★★問い 3: `hm` が出ても⑤は消えるか → ★消えない（★自分の 2 定理の突き合わせ）

* `FirstJumpRouteEquiv.exists_firstJump_data_iff`:
  `FirstJumpRoute` の 4 仮説 ⟺ `c k ≤ axDecay p k`。
* `ExitLossDepth.exit_does_not_realize_cEx`（⑤）: 本日の出口の損失は `axDecay p 1` 以上。

⇒ `hm` / `hjump` / `he` の 3 つは★**塔のデータから出る**（§1 `firstJump_three_of_tower`）が、
★★**4 つ目の `hform : c k ≤ p^{u/(m·e)}` は `c k = axDecay p 1` では `k ≥ 2` で満たせない**
（§2 `not_hform_of_axDecay_one`、まとめは `firstJump_gap`）。

★`hm` をいくら供給しても `hform` は変わらない —— `hform` は
★**降下が実際に「第 1 跳びの損失」を達成すること**を要求しており、
本日の出口は**最後の跳び**を使うからである（それが⑤の中身）。
⇒ ★★**本体が「最初に確かめる価値がある」と書いた点は、確かめた結果「消えない」であった。**

## ★★次の 1 点を「仮説 3 本」の形にした（§3）

`axLemma_of_first_jump_loss` / `axSenTate_of_first_jump_loss`:
★塔のデータ（`hm0` / `hdvd` / `heq` / `heK` / `hu`）を渡すと
`NormalizedTraceDescent.lean:379 FirstJumpRoute.axLemma_of_firstJump` の
分岐論的な仮説は**全部落ちて**、残るのは

* ★`hform : c k ≤ p^{u k/(m k·e k)}`（★**第 1 跳びの損失を達成する降下**）、
* `hc : 1 ≤ c k`、`h : AxWildDescent K c`

の 3 つだけになる。★★**これが次のノードの正確な形である。**

## ★①不分岐側との関係 —— ★本日初めて「①に依存する」項目が出た

★本波は分岐指数の等式 `e k = p · eK k`（`e_{E₁} = p·e_K`）を**仮説 `heq` として受けている**。
★これは `E₁/K` が**完全分岐**であることの帰結であり、不分岐なら `e_{E₁} = e_K` で `heq` は偽。
⇒ ★★**`FirstJumpRoute` の道は①（不分岐側）と独立ではない。**
★本日の他の測定（②′・③・⑤）はすべて①と独立だったが、★これだけは違う。

## ★穴の現状（本波で増減なし）

①不分岐（`TotallyRamifiedLayer.lean:342`）/ ②′（`ExitLossDepth.not_forall_prod_Icc_le_of_le`）/
③幾何減衰（`JumpGeometricDecay.no_uniform_geometric_of_harith`）/
⑤出口の限界（`ExitLossDepth.exit_does_not_realize_cEx`。★本波で「`hm` では消えない」と確定）。
④は消えたまま。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は pGC 物理 p.6 Corollary 3.1（上流と同じ）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. §1・§2 は `ℕ` と `ℝ` だけで、分岐・付値・Galois の語彙が 1 語も出ない。
4. ★`heq : e k = p * eK k` は**仮説のまま**で、塔から導いていない。
   ★導くには「`E₁/K` は次数 `p` の完全分岐」を形式化する必要がある（本波では行っていない）。
5. ★`hdvd : p^{k−1} ∣ m k` も**仮説のまま**である（`m = e(L/E₁)` が `p^{k−1}` で割れること）。
-/

namespace ABC3.Found.PGC

namespace FirstJumpGap

/-! ## §1 塔のデータから 3 つの仮説を出す（`ℕ` だけ） -/

section Supply

/-- `hm : p^{k−1} ≤ m` は `p^{k−1} ∣ m` から出る（`Nat.le_of_dvd`）。 -/
theorem hm_of_dvd {p k m : ℕ} (hm0 : 0 < m) (hdvd : p ^ (k - 1) ∣ m) : p ^ (k - 1) ≤ m :=
  Nat.le_of_dvd hm0 hdvd

/-- ★★**`p` の差は消える** —— 前波の `(p−1)·u₁ ≤ p·e_K` が `hjump : (p−1)·u₁ ≤ e_{E₁}`
そのものになる。★`e_{E₁} = p·e_K` は `E₁/K` が次数 `p` の完全分岐だからである。 -/
theorem hjump_of_sub_one_mul_le {p u eK e : ℕ} (he : e = p * eK)
    (h : (p - 1) * u ≤ p * eK) : (p - 1) * u ≤ e := he ▸ h

/-- ★`FirstJumpRoute` の 4 仮説のうち **3 つ**（`hm` / `he` / `hjump`）は
塔のデータから出る。★残るのは `hform` だけである。 -/
theorem firstJump_three_of_tower {p k u eK m e : ℕ}
    (hm0 : 0 < m) (hdvd : p ^ (k - 1) ∣ m) (heq : e = p * eK) (hp : 0 < p) (heK : 0 < eK)
    (hu : (p - 1) * u ≤ p * eK) :
    p ^ (k - 1) ≤ m ∧ 0 < e ∧ (p - 1) * u ≤ e :=
  ⟨hm_of_dvd hm0 hdvd, by subst heq; exact Nat.mul_pos hp heK,
    hjump_of_sub_one_mul_le heq hu⟩

end Supply

/-! ## §2 ★★★しかし 4 つ目（`hform`）は `axDecay p 1` では満たせない -/

section Gap

variable {p : ℕ} [Fact p.Prime]

/-- ★★★**4 つ目は出口が出す `axDecay p 1` では満たせない**（`k ≥ 2`）。

`JumpArith.rpow_div_le_axDecay` が `p^{u/(m·e)} ≤ axDecay p k` を与え、
`PGroupDescentToAxWild.axDecay_succ_lt_axDecay_one` が `axDecay p k < axDecay p 1`
を与えるので、`axDecay p 1 ≤ p^{u/(m·e)}` は矛盾する。
⇒ ★`hm` をいくら供給しても⑤は消えない。 -/
theorem not_hform_of_axDecay_one {k : ℕ} (hk : 2 ≤ k) {j m e : ℕ}
    (hm : p ^ (k - 1) ≤ m) (he : 0 < e) (hjump : (p - 1) * j ≤ e) :
    ¬ (axDecay p 1 ≤ (p : ℝ) ^ ((j : ℝ) / ((m : ℝ) * (e : ℝ)))) := by
  intro h
  have h1 := JumpArith.rpow_div_le_axDecay (p := p) hm he hjump
  obtain ⟨n, rfl⟩ : ∃ n, k = n + 1 := ⟨k - 1, by omega⟩
  have h2 := PGroupToAxWild.axDecay_succ_lt_axDecay_one (p := p) (k := n) (by omega)
  linarith

/-- ★★★★**本ファイルの主結果** —— 3 つはそろうが 4 つ目は満たせない。 -/
theorem firstJump_gap {k : ℕ} (hk : 2 ≤ k) {u eK m e : ℕ}
    (hm0 : 0 < m) (hdvd : p ^ (k - 1) ∣ m) (heq : e = p * eK) (heK : 0 < eK)
    (hu : (p - 1) * u ≤ p * eK) :
    (p ^ (k - 1) ≤ m ∧ 0 < e ∧ (p - 1) * u ≤ e)
      ∧ ¬ (axDecay p 1 ≤ (p : ℝ) ^ ((u : ℝ) / ((m : ℝ) * (e : ℝ)))) := by
  have hp : 0 < p := (Fact.out : p.Prime).pos
  obtain ⟨h1, h2, h3⟩ := firstJump_three_of_tower hm0 hdvd heq hp heK hu
  exact ⟨⟨h1, h2, h3⟩, not_hform_of_axDecay_one hk h1 h2 h3⟩

end Gap

/-! ## §3 ★★次のノードの正確な形（残るのは仮説 3 本） -/

section Target

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-- ★★**次のノードの正確な形** —— 塔のデータを渡すと分岐論的な仮説は全部落ち、
残るのは `hform`（第 1 跳びの損失を達成する降下）と `hc` と `AxWildDescent K c` の 3 本。 -/
theorem axLemma_of_first_jump_loss (K : PAdicLocalField p) {c : ℕ → ℝ}
    {u eK m e : ℕ → ℕ}
    (hc : ∀ k, 1 ≤ c k)
    (hform : ∀ k, c k ≤ (p : ℝ) ^ ((u k : ℝ) / ((m k : ℝ) * (e k : ℝ))))
    (hm0 : ∀ k, 0 < m k) (hdvd : ∀ k, p ^ (k - 1) ∣ m k)
    (heq : ∀ k, e k = p * eK k) (heK : ∀ k, 0 < eK k)
    (hu : ∀ k, (p - 1) * u k ≤ p * eK k)
    (h : AxWildDescent K c) : AxLemma K (axConstant p) := by
  have hp : 0 < p := (Fact.out : p.Prime).pos
  refine FirstJumpRoute.axLemma_of_firstJump K hc hform
    (fun k => hm_of_dvd (hm0 k) (hdvd k)) (fun k => ?_)
    (fun k => hjump_of_sub_one_mul_le (heq k) (hu k)) h
  rw [heq k]
  exact Nat.mul_pos hp (heK k)

/-- 同じものの `AxSenTate K` 版。 -/
theorem axSenTate_of_first_jump_loss (K : PAdicLocalField p) {c : ℕ → ℝ}
    {u eK m e : ℕ → ℕ}
    (hc : ∀ k, 1 ≤ c k)
    (hform : ∀ k, c k ≤ (p : ℝ) ^ ((u k : ℝ) / ((m k : ℝ) * (e k : ℝ))))
    (hm0 : ∀ k, 0 < m k) (hdvd : ∀ k, p ^ (k - 1) ∣ m k)
    (heq : ∀ k, e k = p * eK k) (heK : ∀ k, 0 < eK k)
    (hu : ∀ k, (p - 1) * u k ≤ p * eK k)
    (h : AxWildDescent K c) : AxSenTate K :=
  axSenTate_of_axLemma (le_trans zero_le_one (one_le_axConstant p))
    (axLemma_of_first_jump_loss K hc hform hm0 hdvd heq heK hu h)

end Target

/-! ## §4 `.src` と 使っている公理の一覧 -/

section Src

def firstJump_three_of_tower.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def not_hform_of_axDecay_one.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def firstJump_gap.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def axSenTate_of_first_jump_loss.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end Src

#print axioms hm_of_dvd
#print axioms hjump_of_sub_one_mul_le
#print axioms firstJump_three_of_tower
#print axioms not_hform_of_axDecay_one
#print axioms firstJump_gap
#print axioms axLemma_of_first_jump_loss
#print axioms axSenTate_of_first_jump_loss

end FirstJumpGap

end ABC3.Found.PGC
