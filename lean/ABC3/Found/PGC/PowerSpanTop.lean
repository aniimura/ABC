import ABC3.Found.PGC.CyclotomicSubstitution

/-!
# [pGC] `hfr : finrank = p` が定理になった —— ★`IntermediateField` も `minpoly` も使わずに

## 持ち場（残り 2 点のうち `finrank ≤ p`）

前波で私は「`finrank ≤ p` は `IntermediateField.adjoin.finrank` ＋ `minpoly.min` で出るが、
★`IntermediateField` を触るので **#59/#69 の領域**」と書いた。

## ★★道を変えた —— `#296` の助言どおり「作らずに済ませる」

★**`IntermediateField` も `minpoly` も使わなかった。** 代わりに:

1. `π^p = Σ_{j<p} c_j π^j`（モニックな関係式）から ★**すべての冪**が
   `{π^0,…,π^{p−1}}` の張る部分加群に入る（§1、強帰納法）。
2. `Algebra.adjoin_eq_span` ＋ `Submonoid.mem_closure_singleton` で
   `adjoin K {π} = ⊤` から ★**`span = ⊤`**（§1）。
3. `finrank_le_of_span_eq_top`（mathlib、`LinearAlgebra/Dimension/Constructions.lean:491`）で
   ★`finrank ≤ p`（§2）。
4. 前波の `le_finrank_of_valK` と合わせて ★★**`finrank = p`**（§2）。

★`lean-idioms` #296 は「★そもそも `IntermediateField` を作らず …… 触らずに済む」と
言っている。★**本波はその通りになった**（`Submodule.span` だけで閉じた）。

## ★関係式は円分塔から出る（§3）

`π = ζ − 1`、`ζ^p = ξ ∈ E₁` なら `(π+1)^p = ξ` で、二項展開して `π^p` を残せば関係式になる。
★係数は `c_0 = Ξ − 1`、`c_j = −C(p,j)`（`j ≥ 1`）。

★**`TotallyRamifiedLayer.adjoin_eq_top_of_valK` は使えない** —— あちらは `finrank = n` を
**仮定**するので循環する。`hadj`（`L = E₁(ζ)`）は塔の定義として受ける。

## 何を埋めたか

| 節 | 宣言 | 内容 |
|---|---|---|
| §1 | `pow_mem_span_of_relation` | ★抽象核。分岐・付値が 1 語も出ない |
| §1 | `span_eq_top_of_relation` | ★`adjoin = ⊤` ＋ 関係式 ⇒ 冪 `p` 本が張る |
| §2 | `finrank_le_of_relation` / `finrank_eq_of_relation` | ★`hfr` の残り半分と完成 |
| §3 | `relation_of_pow_add_one` | 二項展開（`0 < p` が要る） |
| §3 | `finrank_eq_of_zeta_pow` | ★★★**`hfr` の完成形** |

## ★測って分かったこと（仮定が 1 つ減った）

§1 の 2 本は ★**`0 < p` を要らない**（最初は書いていたが外しても通った）。
`p = 0` のとき関係式は `1 = 0` を意味し、`M` が自明環になるので結論も成り立つ。
★§3 の `relation_of_pow_add_one` では `0 < p` が**要る**（`Ξ` の項を取り出すのに
`0 ∈ range p` が要るため）。

## ★★残るのは `hvalK` 1 つだけになった

前波の棚卸し表の「★具体層」3 行のうち、`hwn` は前波で、`hfr` は本波で定理になった。
★**残るのは `hvalK`（`E₁` の値群が `‖μ‖^ℤ`）だけ**である。
★それには `finrank ℚ_p E₁ = φ(p^{m+1})`（**1 段下の絶対次数**）が要る ——
★本波の §1–§2 がそのまま使える可能性がある（`ℚ_p ⊂ E₁` でも同じ形の関係式
`μ` の最小多項式が Eisenstein）。★開いていないので費用は書かない。

## 逸脱の記録

- `hadj`（`Algebra.adjoin K {π} = ⊤`）は仮説。★塔の定義そのものである。
- ★`p = 2` は本鎖では対象外だが、★本ファイルの定理は `p = 2` でも成り立つ
  （`p` に関する条件は `0 < p` だけ）。
-/

namespace ABC3.Found.PGC

namespace PowerSpanTop

open Finset

/-! ## §1 抽象核 —— モニックな関係式があれば冪は有限個で張れる -/

section Kernel

variable {K M : Type*} [CommRing K] [CommRing M] [Algebra K M]

/-- ★抽象核。`π^p` が低い冪の `K` 係数の和なら、★**すべての冪**が
`{π^0, …, π^{p−1}}` の張る部分加群に入る。★分岐も付値も出てこない。 -/
theorem pow_mem_span_of_relation {π : M} {p : ℕ} (c : ℕ → K)
    (hrel : π ^ p = ∑ j ∈ range p, algebraMap K M (c j) * π ^ j) :
    ∀ k : ℕ, π ^ k ∈ Submodule.span K (Set.range fun l : Fin p => π ^ (l : ℕ)) := by
  intro k
  induction k using Nat.strong_induction_on with
  | _ k ih =>
    by_cases hk : k < p
    · exact Submodule.subset_span ⟨⟨k, hk⟩, rfl⟩
    · have hkp : p ≤ k := not_lt.mp hk
      have hrw : π ^ k = ∑ j ∈ range p, algebraMap K M (c j) * π ^ (k - p + j) := by
        have h1 : π ^ k = π ^ (k - p) * π ^ p := by
          rw [← pow_add]
          congr 1
          omega
        rw [h1, hrel, Finset.mul_sum]
        refine Finset.sum_congr rfl fun j _ => ?_
        rw [← mul_assoc, mul_comm (π ^ (k - p)), mul_assoc, ← pow_add]
      rw [hrw]
      refine Submodule.sum_mem _ fun j hj => ?_
      have hlt : k - p + j < k := by
        have := Finset.mem_range.mp hj
        omega
      have hmem := ih (k - p + j) hlt
      rw [← Algebra.smul_def]
      exact Submodule.smul_mem _ (c j) hmem

/-- ★`Algebra.adjoin K {π} = ⊤` と関係式から、★**冪 `p` 本が張る**。 -/
theorem span_eq_top_of_relation {π : M} {p : ℕ} (c : ℕ → K)
    (hrel : π ^ p = ∑ j ∈ range p, algebraMap K M (c j) * π ^ j)
    (hadj : Algebra.adjoin K ({π} : Set M) = ⊤) :
    Submodule.span K (Set.range fun l : Fin p => π ^ (l : ℕ)) = ⊤ := by
  have hle : Subalgebra.toSubmodule (Algebra.adjoin K ({π} : Set M))
      ≤ Submodule.span K (Set.range fun l : Fin p => π ^ (l : ℕ)) := by
    rw [Algebra.adjoin_eq_span]
    refine Submodule.span_le.mpr ?_
    rintro y hy
    obtain ⟨k, rfl⟩ := Submonoid.mem_closure_singleton.mp hy
    exact pow_mem_span_of_relation c hrel k
  rw [Algebra.toSubmodule_eq_top.mpr hadj] at hle
  exact top_le_iff.mp hle

end Kernel

/-! ## §2 ★`finrank ≤ p`（`hfr` の残り半分） -/

section Finrank

/-- ★★`hfr` の残り半分。★`IntermediateField` も `minpoly` も**使わない**
（#59/#69 を回避、#296 の助言どおり）。`finrank_le_of_span_eq_top` に §1 を渡すだけ。 -/
theorem finrank_le_of_relation {K M : Type*} [Field K] [CommRing M] [Algebra K M]
    {π : M} {p : ℕ} (c : ℕ → K)
    (hrel : π ^ p = ∑ j ∈ range p, algebraMap K M (c j) * π ^ j)
    (hadj : Algebra.adjoin K ({π} : Set M) = ⊤) :
    Module.finrank K M ≤ p := by
  have h := finrank_le_of_span_eq_top (R := K) (span_eq_top_of_relation c hrel hadj)
  simpa using h

/-- ★★★`hfr : finrank K M = p` が**定理になった**。 -/
theorem finrank_eq_of_relation {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M]
    [Algebra K M] [FiniteDimensional K M] {π : M} {p : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≠ 1)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((p : ℤ) * m))
    (c : ℕ → K) (hrel : π ^ p = ∑ j ∈ range p, algebraMap K M (c j) * π ^ j)
    (hadj : Algebra.adjoin K ({π} : Set M) = ⊤) :
    Module.finrank K M = p :=
  CyclotomicSubstitution.finrank_eq_of_le hπ0 hπ1 hvalK (finrank_le_of_relation c hrel hadj)

end Finrank

/-! ## §3 円分塔での関係式 —— `ζ^p = ξ` から作る -/

section Cyclotomic

/-- ★`π = ζ − 1`、`ζ^p = ξ ∈ K` なら、関係式 `π^p = Σ_{j<p} c_j π^j` が**書ける**。

`(π+1)^p = ξ` を二項展開して `π^p` を残すだけ。★係数は
`c_0 = Ξ − 1`、`c_j = −C(p,j)`（`j ≥ 1`）。 -/
theorem relation_of_pow_add_one {K M : Type*} [CommRing K] [CommRing M] [Algebra K M]
    {π : M} {Ξ : K} {p : ℕ} (hp : 0 < p) (h : (π + 1) ^ p = algebraMap K M Ξ) :
    π ^ p = ∑ j ∈ range p,
      algebraMap K M ((if j = 0 then Ξ else 0) - (p.choose j : K)) * π ^ j := by
  have hbin : (π + 1) ^ p = (∑ j ∈ range p, (p.choose j : M) * π ^ j) + π ^ p := by
    rw [add_pow, Finset.sum_range_succ]
    simp [mul_comm]
  have hone : ∑ j ∈ range p, algebraMap K M (if j = 0 then Ξ else 0) * π ^ j
      = algebraMap K M Ξ := by
    rw [Finset.sum_eq_single_of_mem 0 (Finset.mem_range.mpr hp)
      (fun j _ hj => by simp [hj])]
    simp
  have hsum : ∑ j ∈ range p,
      algebraMap K M ((if j = 0 then Ξ else 0) - (p.choose j : K)) * π ^ j
      = algebraMap K M Ξ - ∑ j ∈ range p, (p.choose j : M) * π ^ j := by
    rw [← hone, ← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl fun j _ => by rw [map_sub, map_natCast]; ring
  rw [hsum, ← h, hbin]
  ring

/-- ★★★**`hfr` の完成形** —— `L = E₁(ζ)` と `ζ^p = ξ ∈ E₁` と `hvalK` から
`finrank E₁ L = p` が出る。

★`TotallyRamifiedLayer.adjoin_eq_top_of_valK` は使えない（あちらは `finrank = n` を
**仮定**するので循環する）。★`hadj` は塔の定義（`L` は `ζ` で生成される）として受ける。 -/
theorem finrank_eq_of_zeta_pow {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M]
    [Algebra K M] [FiniteDimensional K M] {π : M} {Ξ : K} {p : ℕ}
    (hp : 0 < p) (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≠ 1)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hadj : Algebra.adjoin K ({π} : Set M) = ⊤)
    (h : (π + 1) ^ p = algebraMap K M Ξ) :
    Module.finrank K M = p :=
  finrank_eq_of_relation hπ0 hπ1 hvalK _ (relation_of_pow_add_one hp h) hadj

end Cyclotomic

/-! ## §4 使っている公理の一覧 -/

#print axioms pow_mem_span_of_relation
#print axioms span_eq_top_of_relation
#print axioms finrank_le_of_relation
#print axioms finrank_eq_of_relation
#print axioms relation_of_pow_add_one
#print axioms finrank_eq_of_zeta_pow

end PowerSpanTop

end ABC3.Found.PGC
