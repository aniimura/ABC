import ABC3.Found.PGC.CoefficientIntegrality

/-!
# [pGC] 最後の穴を**深さでなく剰余**で塞ぐ —— スロットの位の剰余は添字そのもの

## 持ち場（前波で私が (n12) で特定した最後の穴）

前波の測定 (n12) は「`v(σf_j − f_j) ≥ v(f_{j₀}) + 2p`（深さ）は p≥3 で **294/330**、
★破れるときは必ず `j = 0`」と言っていた。⇒ 残る穴は **`E₁` 成分 `f_0` の `σ` による動き**。

★★**深さで押さえるのをやめた。** `f_0` の動きは `E₁` の元だから
★**位の剰余が `0 (mod p)`** で、`j₀ ≢ 0` とは**必ず違う剰余類**にある。
⇒ 打ち消しようがない。深さは要らない。

## ★測定が形を選んだ（3 つの形を同じ母集団で比べた）

| Lean の定理 | 仮定 | p ≥ 3 での充足率 |
|---|---|---|
| `LossExponentMatch.loss_le_of_error_small` | 誤差が `‖f_{j₀}‖‖ρ‖²` 以下（**深さ**） | 294/330 |
| `LossExponentMatch.loss_le_of_error_ne` | 誤差と主部が**同位でない** | 328/330 |
| ★`SlotResidue.loss_le_of_remainder_slots`（本ファイル） | ★`j₀` スロットが**最大でない**（**剰余**） | ★**330/330** |

（`tools/numerology-check.py` の (n12)(n8)(n15)。`p = 2` はどの形でも 12 件破れ、
★その 12 件は真の `loss > 2p−2` になる 12 件と一貫して一致する。）

## 何を埋めたか

| 節 | 宣言 | 内容 |
|---|---|---|
| §1 | `slot_pairwise_ne` | スロットのノルムは相異なる（指数が `n` を法として相異なる） |
| §1 | `exists_slot_of_sum` | ★スロットの和のノルムは**どれか 1 つのスロット**のノルムに等しい |
| §2 | `residue_of_slot` | ★★**`E₁` 係数 × `π^j` の位は `j` (mod n)**。★6 波にわたって非形式に使ってきた「スロットの剰余は全部違う」の形式化 |
| §2 | `residue_of_sum_slots` | 和の位の剰余は「最大を達成するスロットの添字」 |
| §3 | `slot_le_sum` / `exists_slot_ne_of_lt` | ★`j₀` スロットが最大でなければ剰余は `j₀` と違う |
| §4 | `loss_le_of_remainder_slots` | ★★★**到達点**。深さを一切要求せずに `loss ≤ 2p−2` |

## ★木の断定の検算（前波に続き 1 件）

`TotallyRamifiedValueGroup.lean:113` が `IsUltrametricDist.nnnorm_sum_eq_sup_of_pairwise_ne`
を使っている（★索引には出ない `to_additive` 生成名。#354）。★本ファイルはそれを
**そのまま**使った。★「木が実際に使っている箇所を読む」経路が索引より強いことの 2 例目。

## 逸脱の記録

- ★残りを「`𝒪_{E₁}` スロットの和で書ける」ことは**仮説** `hrem` で受けている。
  これは `CoefficientIntegrality.exists_integral_expansion_of_valK`（前波）から出るが、
  ★**展開の一意性**（`A − B·π^{j₀}` の展開が「その」`z` であること）は使っていない
  ——本定理は「そう書ける `z` が 1 つあれば十分」の形にしてある。
- ★`p = 2` はどの形でも 12/200 破れる。`loss ≤ 2p−2` は `p = 2` では**偽**
  （真の最大は `2p−1 = 3`）であり、本ファイルもそれを直していない（直せない）。
-/

namespace ABC3.Found.PGC

namespace SlotResidue

open Finset

section Slots

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

omit [IsUltrametricDist M] in
/-- スロットのノルムは相異なる（指数が `n` を法として相異なるから）。 -/
theorem slot_pairwise_ne {π : M} {n : ℕ} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (z : ℕ → K) (t : Finset ℕ) (ht : ∀ i ∈ t, i < n ∧ z i ≠ 0) :
    (t : Set ℕ).Pairwise fun i i' =>
      ‖algebraMap K M (z i) * π ^ i‖₊ ≠ ‖algebraMap K M (z i') * π ^ i'‖₊ := by
  have hexp : ∀ i ∈ t, ∃ m : ℤ,
      ‖algebraMap K M (z i) * π ^ i‖ = ‖π‖ ^ ((n : ℤ) * m + (i : ℤ)) := by
    intro i hi
    obtain ⟨m, hm⟩ := hvalK (z i) (ht i hi).2
    refine ⟨m, ?_⟩
    simp only [norm_mul, norm_pow, hm]
    rw [zpow_add₀ (ne_of_gt hπ0), zpow_natCast]
  intro i hi i' hi' hii hcon
  have hreal : ‖algebraMap K M (z i) * π ^ i‖ = ‖algebraMap K M (z i') * π ^ i'‖ := by
    simpa using congrArg NNReal.toReal hcon
  obtain ⟨m, hm⟩ := hexp i hi
  obtain ⟨m', hm'⟩ := hexp i' hi'
  rw [hm, hm'] at hreal
  have hzeq := (zpow_right_inj₀ hπ0 (ne_of_lt hπ1)).mp hreal
  refine TotallyRamifiedLayer.not_dvd_sub_of_lt (ht i hi).1 (ht i' hi').1 hii ⟨m' - m, ?_⟩
  linear_combination hzeq

/-- ★★**スロットの和のノルムは、どれか 1 つのスロットのノルムに等しい**。

超距離で「ノルムが相異なる項の和」は最大の項に等しくなる（`p` 進の常識だが、
★ここでは「指数が `n` を法として相異なる」ことから導いている）。 -/
theorem exists_slot_of_sum {π : M} {n : ℕ} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (z : ℕ → K) :
    (∀ j, j < n → z j = 0) ∨
      ∃ j, j < n ∧ z j ≠ 0 ∧
        ‖∑ i ∈ range n, algebraMap K M (z i) * π ^ i‖ = ‖algebraMap K M (z j) * π ^ j‖ := by
  classical
  set y : ℕ → M := fun i => algebraMap K M (z i) * π ^ i with hy
  set t : Finset ℕ := (range n).filter (fun i => z i ≠ 0) with ht
  have hmem : ∀ i, i ∈ t ↔ (i ∈ range n ∧ z i ≠ 0) := fun i => Finset.mem_filter
  by_cases hte : t.Nonempty
  · right
    have hsum : ∑ i ∈ range n, y i = ∑ i ∈ t, y i := by
      refine (Finset.sum_subset (Finset.filter_subset _ _) ?_).symm
      intro i hi hit
      have hzi : z i = 0 := by
        by_contra hne
        exact hit ((hmem i).mpr ⟨hi, hne⟩)
      simp [hy, hzi]
    have hpair := slot_pairwise_ne hπ0 hπ1 hvalK z t
      (fun i hi => ⟨Finset.mem_range.mp ((hmem i).mp hi).1, ((hmem i).mp hi).2⟩)
    have hsup : ‖∑ i ∈ t, y i‖₊ = t.sup fun i => ‖y i‖₊ :=
      IsUltrametricDist.nnnorm_sum_eq_sup_of_pairwise_ne hpair
    obtain ⟨j, hjt, hj⟩ := Finset.exists_mem_eq_sup t hte (fun i => ‖y i‖₊)
    refine ⟨j, Finset.mem_range.mp ((hmem j).mp hjt).1, ((hmem j).mp hjt).2, ?_⟩
    have hnn : ‖∑ i ∈ t, y i‖₊ = ‖y j‖₊ := by rw [hsup, hj]
    have h2 : ‖∑ i ∈ t, y i‖ = ‖y j‖ := by simpa using congrArg NNReal.toReal hnn
    rw [hsum]
    exact h2
  · left
    intro j hj
    by_contra hne
    exact hte ⟨j, (hmem j).mpr ⟨Finset.mem_range.mpr hj, hne⟩⟩

/-! ## §2 ★スロットの位の剰余は添字そのもの -/

omit [IsUltrametricDist M] in
/-- ★★**`E₁` 係数 × `π^j` の位は `j` (mod n)**。★6 波にわたって非形式に使ってきた
「スロットの剰余は全部違う」の、形式化された形。 -/
theorem residue_of_slot {π : M} {n j : ℕ} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : 0 < n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    {a : K} (ha : a ≠ 0) (ha1 : ‖algebraMap K M a‖ ≤ 1) (hj : j < n) :
    ∃ v : ℕ, ‖algebraMap K M a * π ^ j‖ = ‖π‖ ^ v ∧ v % n = j := by
  obtain ⟨m, hm⟩ := hvalK a ha
  have hnm : 0 ≤ (n : ℤ) * m := by
    rw [hm] at ha1
    exact (zpow_le_one_iff_right_of_lt_one₀ hπ0 hπ1).mp ha1
  have hm0 : 0 ≤ m :=
    CoefficientIntegrality.exponent_nonneg_of_lt (n := n) (j := 0) hn (by simpa using hnm)
  refine ⟨n * m.toNat + j, ?_, ?_⟩
  · rw [norm_mul, norm_pow, hm, ← zpow_natCast ‖π‖ j, ← zpow_add₀ (ne_of_gt hπ0),
      ← zpow_natCast ‖π‖ (n * m.toNat + j)]
    congr 1
    push_cast [Int.toNat_of_nonneg hm0]
    ring
  · rw [Nat.mul_add_mod, Nat.mod_eq_of_lt hj]

/-- ★和の位の剰余は「最大を達成するスロットの添字」。 -/
theorem residue_of_sum_slots {π : M} {n : ℕ} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : 0 < n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (z : ℕ → K) (hz1 : ∀ j, ‖algebraMap K M (z j)‖ ≤ 1) :
    (∀ j, j < n → z j = 0) ∨
      ∃ j v : ℕ, j < n ∧ ‖∑ i ∈ range n, algebraMap K M (z i) * π ^ i‖ = ‖π‖ ^ v ∧ v % n = j := by
  rcases exists_slot_of_sum hπ0 hπ1 hvalK z with hzero | ⟨j, hjn, hjz, hj⟩
  · exact Or.inl hzero
  · obtain ⟨v, hv, hvr⟩ := residue_of_slot hπ0 hπ1 hn hvalK hjz (hz1 j) hjn
    exact Or.inr ⟨j, v, hjn, hj.trans hv, hvr⟩

/-! ## §3 ★「`j₀` スロットが最大でない」なら剰余は `j₀` と違う -/

/-- どのスロットも和より大きくならない。 -/
theorem slot_le_sum {π : M} {n i : ℕ} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (z : ℕ → K) (hi : i < n) :
    ‖algebraMap K M (z i) * π ^ i‖ ≤ ‖∑ i' ∈ range n, algebraMap K M (z i') * π ^ i'‖ := by
  classical
  by_cases hz : z i = 0
  · simp [hz]
  set y : ℕ → M := fun i' => algebraMap K M (z i') * π ^ i' with hy
  set t : Finset ℕ := (range n).filter (fun i' => z i' ≠ 0) with ht
  have hmem : ∀ i', i' ∈ t ↔ (i' ∈ range n ∧ z i' ≠ 0) := fun i' => Finset.mem_filter
  have hsum : ∑ i' ∈ range n, y i' = ∑ i' ∈ t, y i' := by
    refine (Finset.sum_subset (Finset.filter_subset _ _) ?_).symm
    intro i' hi' hit
    have hzi : z i' = 0 := by
      by_contra hne
      exact hit ((hmem i').mpr ⟨hi', hne⟩)
    simp [hy, hzi]
  have hpair := slot_pairwise_ne hπ0 hπ1 hvalK z t
    (fun i' hi' => ⟨Finset.mem_range.mp ((hmem i').mp hi').1, ((hmem i').mp hi').2⟩)
  have hit : i ∈ t := (hmem i).mpr ⟨Finset.mem_range.mpr hi, hz⟩
  rw [hsum]
  exact CoefficientIntegrality.norm_term_le_sum_of_pairwise_ne hpair hit

/-- ★★`j₀` スロットより大きいスロットが 1 つでもあれば、和の位の剰余は `j₀` と違う。

★これが最後の穴（`f_0` の動き）を塞ぐ形である。★深さで押さえるのではなく、
★**剰余類が違う**ことで押さえる。 -/
theorem exists_slot_ne_of_lt {π : M} {n j₀ i₀ : ℕ} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : 0 < n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (z : ℕ → K) (hz1 : ∀ j, ‖algebraMap K M (z j)‖ ≤ 1) (hi₀ : i₀ < n)
    (hlt : ‖algebraMap K M (z j₀) * π ^ j₀‖ < ‖algebraMap K M (z i₀) * π ^ i₀‖) :
    ∃ j v : ℕ, j < n ∧ j ≠ j₀ ∧
      ‖∑ i ∈ range n, algebraMap K M (z i) * π ^ i‖ = ‖π‖ ^ v ∧ v % n = j := by
  have hz0 : z i₀ ≠ 0 := by
    intro h
    rw [h] at hlt
    simp only [map_zero, zero_mul, norm_zero] at hlt
    exact absurd hlt (not_lt.mpr (norm_nonneg _))
  rcases exists_slot_of_sum hπ0 hπ1 hvalK z with hzero | ⟨j, hjn, hjz, hj⟩
  · exact absurd (hzero i₀ hi₀) hz0
  · obtain ⟨v, hv, hvr⟩ := residue_of_slot hπ0 hπ1 hn hvalK hjz (hz1 j) hjn
    refine ⟨j, v, hjn, ?_, hj.trans hv, hvr⟩
    intro hcon
    subst hcon
    exact absurd (hj ▸ slot_le_sum hπ0 hπ1 hvalK z hi₀) (not_le.mpr hlt)

/-! ## §4 ★到達点 —— 深さではなく**剰余**で `loss ≤ 2p−2` を出す -/

/-- ★★★**本ファイルの到達点**。残り（`R` も誤差も `f_0` の動きも全部込み）が
`𝒪_{E₁}` スロットの和で書けて、★**`j₀` スロットが最大でない**なら `loss ≤ 2p−2`。

★前波までの形は「誤差が `‖f_{j₀}‖·‖ρ‖²` 以下（深い）」を要求していた。それは
`f_0`（`E₁` 成分）の動きで **36/330 件破れる**（測定 (n12)）。★本定理は深さを要求せず、
**剰余類が違う**ことだけを使う。★測定 (n14) はこの形が **p ≥ 3 で 330/330** と言っている。 -/
theorem loss_le_of_remainder_slots {π A B : M} {p d j₀ jstar i₀ vf0 vfs : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hp : 2 ≤ p)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hj2 : j₀ ≤ p - 1) (hjp : j₀ < p) (hs1 : 1 ≤ jstar)
    (hmin : vf0 ≤ vfs) (hd : d = vfs + jstar) (hdvd : p ∣ vf0)
    (hB : ‖B‖ = ‖π‖ ^ (p + vf0))
    (z : ℕ → K) (hz1 : ∀ j, ‖algebraMap K M (z j)‖ ≤ 1)
    (hrem : A - B * π ^ j₀ = ∑ i ∈ range p, algebraMap K M (z i) * π ^ i)
    (hi₀ : i₀ < p)
    (hlt : ‖algebraMap K M (z j₀) * π ^ j₀‖ < ‖algebraMap K M (z i₀) * π ^ i₀‖) :
    ‖π‖ ^ (d + (2 * p - 2)) ≤ ‖A‖ := by
  obtain ⟨j, v, hjn, hjne, hv, hvr⟩ :=
    exists_slot_ne_of_lt hπ0 hπ1 (by omega) hvalK z hz1 hi₀ hlt
  have hR : ‖A - B * π ^ j₀‖ = ‖π‖ ^ v := by rw [hrem]; exact hv
  have key := LossExponentMatch.loss_le_two_p_sub_two_of_expansion hπ0 hπ1 hp hj2 hjp hs1
    hmin hd hdvd hB hR hvr (Ne.symm hjne)
  simpa using key

end Slots

/-! ## §5 使っている公理の一覧 -/

#print axioms slot_pairwise_ne
#print axioms exists_slot_of_sum
#print axioms residue_of_slot
#print axioms residue_of_sum_slots
#print axioms slot_le_sum
#print axioms exists_slot_ne_of_lt
#print axioms loss_le_of_remainder_slots

end SlotResidue

end ABC3.Found.PGC
