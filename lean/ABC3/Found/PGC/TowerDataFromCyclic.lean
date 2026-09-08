import ABC3.Found.PGC.SpectralNormBridge

/-!
# [pGC] ★★★塔のデータ 6 つのうち 4 つは無料だった —— 残るのは**巡回性**

## ★安い順に測った結果（持ち場の問い）

| 項目 | 結果 |
|---|---|
| `τ` / `hτ` | ★★**無料**。`g` の加法部分を取るだけ（`tauOf`、`hτ` は `rfl`） |
| `s` / `hsg` | ★★**無料**。`g^{p^j}` の環準同型部分（`sOf`、`hsg` は `rfl`） |
| `g` / `hg` | ★`IsCyclic` ＋ `Nat.card = p^{k+1}` から出る |
| 巡回性そのもの | ★★★**出ない**（下の否定的結果） |

⇒ ★★**6 つの塔のデータは `IsCyclic` と `Nat.card` の 2 つになった**。

## ★★★否定的な結果（§4、本日の主な成果）

`not_forall_exists_orderOf_eq_card` ——
★**位数が `p^{k+1}` であるだけでは位数 `p^{k+1}` の元は取れない**。
反例は `(ℤ/2)²`（`p = 2`, `k = 1`）。
⇒ ★出口の `hg : orderOf g = p^{k+1}` は
`hnK`（次数）と Galois 性だけからは**得られない**。

★★**体の側でも同じである**（★本波は形式化していない。手での検算として記録する）:
`ℚ₂(√−1, √2)/ℚ₂` は次数 4 の**完全分岐 Galois** 拡大で、Galois 群は `(ℤ/2)²`。
（中間の 2 次体は `√−1` / `√2` / `√−2` の 3 つで、どれも分岐する。
`ℚ₂` の不分岐 2 次拡大は `ℚ₂(√−3)` だけである。）
⇒ ★★**`ht`（完全分岐）＋ `hnK`（次数 `p^{k+1}`）＋ Galois 性 でも巡回性は出ない。**
さらに★`K(x)/K` はそもそも Galois とは限らない。

## ★★どこで止まったか（`file:line`）

★原典の道は「`K(x)` そのものを巡回とする」ではなく、
★**Galois 閉包に移って wild inertia（`p`-群）の正規列を取る**ことである。

* `p`-群性の供給元は `RamificationJumpDivisibility.lean:361
  isPGroup_lowerRamificationGroup_one` であり、その仮説
  （`[SMulCommClass G A B]` / `[FaithfulSMul G B]` / `[CharP (ResidueField B) p]` /
  `hα : maximalIdeal B = span {α}` / `hadj`）は
  ★**本連鎖ですでに全部揃えてある**（`IntegerMiscInstances` / `IntegerRingInstances` /
  `IntegerResidueBase`）。`h1 : G_1 = ⊤` も
  `IntegerMiscInstances.lowerRamificationGroup_one_eq_top` で出る。
  ⇒ ★**`IsPGroup p G` は届く。しかし `IsPGroup` は巡回を意味しない。**
* ★★素朴な部分群降下は既に封じられている:
  `WildDepthDescent.lean:280 not_forall_exists_relIndex_padicValNat_eq`
  （★本人が以前に形式化した `A₄` の反例）。
* ⇒ ★★★**残るのは「`p`-群の中に巡回な商の列を取る」ノード**である。
  本連鎖には `TotallyRamifiedLayer.exists_pgroup_descent_cyclic`（第 1082）が既にあるが、
  ★それを `AxWildDescent` の形に繋ぐのは本波では行っていない。

## ★①不分岐側との関係

★**依然に生きている。本波も触っていない。**
ただし★本波の否定的結果は不分岐側とは**独立**である
（反例 `ℚ₂(√−1,√2)` は完全分岐であり、不分岐層を含まない）。

## 逸脱の記録

1. 体の側の反例（`ℚ₂(√−1,√2)`）は★**形式化していない**。群の側の反例のみ定理である。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
-/

namespace ABC3.Found.PGC

namespace TowerData

open ABC3.Skeleton.PGC
open scoped NormedField Valued

/-! ## §1 ★`τ` と `s` は `g` の詰め替えにすぎない -/

section Free

variable {K M : Type*} [Field K] [Ring M] [Algebra K M]

/-- ★`τ : M →+ M` は `g` の加法部分。 -/
def tauOf (g : M ≃ₐ[K] M) : M →+ M :=
  { toFun := fun z => g z
    map_zero' := map_zero g
    map_add' := fun a b => map_add g a b }

theorem tauOf_apply (g : M ≃ₐ[K] M) (z : M) : tauOf g z = g z := rfl

/-- ★`s : ℕ → (M →+* M)` は `g^{p^j}` の環準同型部分。 -/
def sOf (g : M ≃ₐ[K] M) (p : ℕ) : ℕ → (M →+* M) :=
  fun j => (g ^ p ^ j : M ≃ₐ[K] M).toAlgHom.toRingHom

theorem sOf_apply (g : M ≃ₐ[K] M) (p : ℕ) (j : ℕ) (z : M) :
    sOf g p j z = (g ^ p ^ j) z := rfl

end Free

/-! ## §2 ★`g` と `hg` は「巡回 ＋ 位数」から出る -/

section Cyclic

variable {K M : Type*} [Field K] [Field M] [Algebra K M]

/-- ★`Gal(M/K)` が巡回で位数 `n` なら生成元が取れる。 -/
theorem exists_generator_of_isCyclic [IsCyclic (M ≃ₐ[K] M)] {n : ℕ}
    (hcard : Nat.card (M ≃ₐ[K] M) = n) : ∃ g : M ≃ₐ[K] M, orderOf g = n := by
  obtain ⟨g, hg⟩ := IsCyclic.exists_ofOrder_eq_natCard (α := M ≃ₐ[K] M)
  exact ⟨g, by rw [hg, hcard]⟩

end Cyclic

/-! ## §3 ★★★★塔のデータを「巡回 ＋ 位数」に置き換えた出口 -/

section Exit

variable {p : ℕ} [Fact p.Prime]

/-- ★★★★★**塔のデータ 6 つが 2 つになった出口**。

`SpectralNormBridge.exists_norm_sub_algebraMap_le_prod_axDecay_of_adjoin` は
`g` / `hg` / `τ` / `hτ` / `s` / `hsg` の 6 つを要求していた。
★★**`τ` と `s` は `g` の詰め替えにすぎない**（§1）ので、
残るのは `IsCyclic` と `Nat.card` の 2 つだけである。 -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_adjoin
    (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    [IsCyclic ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))]
    (ht : IsTotallyRamifiedAdjoin K x) {k : ℕ}
    (hnK : Module.finrank K.carrier
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) = p ^ (k + 1))
    (hcard : Nat.card ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
      = p ^ (k + 1))
    (y : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    ∃ g : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)),
      ∃ w : GainedTowerModel.GaloisTower.twr g p k,
        ‖y - algebraMap (GainedTowerModel.GaloisTower.twr g p k)
            (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) w‖
          ≤ (∏ j ∈ Finset.Icc 1 (k + 1), axDecay p j) * ‖g y - y‖ := by
  obtain ⟨g, hg⟩ := exists_generator_of_isCyclic hcard
  refine ⟨g, ?_⟩
  have hmain := SpectralBridge.exists_norm_sub_algebraMap_le_prod_axDecay_of_adjoin
    K x ht hnK g hg (tauOf g) (tauOf_apply g) (sOf g p) (sOf_apply g p) y
  simpa [tauOf_apply] using hmain

end Exit

/-! ## §4 ★★★否定的な結果 —— 巡回性は次数からは**出ない** -/

section Negative

/-- 反例の群: `(ℤ/2)²`。 -/
abbrev V4 : Type := Multiplicative (ZMod 2 × ZMod 2)

theorem card_V4 : Nat.card V4 = 2 ^ 2 := by
  rw [Nat.card_eq_fintype_card]
  decide

theorem sq_eq_one_V4 (g : V4) : g ^ 2 = 1 := by
  revert g
  decide

/-- ★★★**位数 `p^{k+1}` だけからは位数 `p^{k+1}` の元は出ない**。

★すなわち、出口の `hg : orderOf g = p^{k+1}` は
`hnK`（次数）と Galois 性だけからは**得られない**。
反例は `(ℤ/2)²`（`p = 2`, `k = 1`）。 -/
theorem not_forall_exists_orderOf_eq_card :
    ¬ ∀ (G : Type) (_ : Group G), Nat.card G = 2 ^ 2 → ∃ g : G, orderOf g = 2 ^ 2 := by
  intro h
  obtain ⟨g, hg⟩ := h V4 inferInstance card_V4
  have hdvd : orderOf g ∣ 2 := orderOf_dvd_of_pow_eq_one (sq_eq_one_V4 g)
  rw [hg] at hdvd
  omega

end Negative


/-! ## `.src` と 公理 -/

def exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_adjoin.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def not_forall_exists_orderOf_eq_card.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

#print axioms tauOf
#print axioms sOf
#print axioms exists_generator_of_isCyclic
#print axioms exists_norm_sub_algebraMap_le_prod_axDecay_of_cyclic_adjoin
#print axioms not_forall_exists_orderOf_eq_card

end TowerData

end ABC3.Found.PGC
