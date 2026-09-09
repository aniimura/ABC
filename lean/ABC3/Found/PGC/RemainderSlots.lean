import ABC3.Found.PGC.SlotResidue

/-!
# [pGC] `hrem` を仮説から外す —— 整数展開を残りに当てる ＋ **展開の一意性**

## 持ち場（前波で私が挙げた残り 3 点の 1 番）

`SlotResidue.loss_le_of_remainder_slots` は「残り `A − B·π^{j₀}` が `𝒪_{E₁}` スロットの和で
書ける」を**仮説 `hrem`** で受けていた。本ファイルがそれを**定理で供給**する。

## ★開いて最初に分かったこと（費用の見立てが 1 つ外れた）

`CoefficientIntegrality.exists_integral_expansion_of_valK`（前波の私のファイル）は
★**`f : ℕ → M` を返す**ので、`SlotResidue` が要求する `z : ℕ → K` に**そのままでは渡せない**。
⇒ §1 で `K` 係数のまま返す版を書いた（中身は同じ、10 行）。

## ★★木の docstring を 2 件検算した（1 件は正しく、1 件は**名前が 1 本ずれている**）

1. `ComponentWitness.lean:50`「（＝展開の**一意性**）。★`linearIndependent_of_ne_mod` が
   その材料である」→ ★**本当だった**。材料はちょうどそれで足り、§2 は 8 行で閉じた。
2. ★`LossTwoPSubTwo.lean:22`「一意性の材料は ★**`TotallyRamifiedLayer.linearIndependent_of_ne_mod`**
   に在る」→ ★**その名前は存在しない**。宣言は
   `TotallyRamifiedValueGroup.lean:133`（名前空間 `TotallyRamified`）にあり、
   `TotallyRamifiedLayer.lean:139` は**使っているだけ**である。
   （測定: `grep -rn "^theorem linearIndependent_of_ne_mod" lean/ABC3/Found/PGC/*.lean` → 1 件、
   `TotallyRamifiedValueGroup.lean:133` のみ）
   ★他ファイルの docstring は直さない規約なので、ここに訂正として書く。
3. 同じ `LossTwoPSubTwo.lean:23`「★`Fin n` 添字の `LinearIndependent` から 2 つの表現の
   一致を出す配管が要る。★費用は書かない」→ ★**測った: 8 行**
   （`Fintype.linearIndependent_iff` → `Algebra.smul_def` → `Fin.sum_univ_eq_sum_range`
   → `Finset.sum_sub_distrib`）。

## 何を埋めたか

| 節 | 宣言 | 内容 |
|---|---|---|
| §1 | `exists_integral_coeffs_of_valK` | `K` 係数のまま返す整数展開 |
| §2 | `coeffs_unique` | ★★**展開の一意性**（`j < n` の範囲で係数は一意） |
| §3 | `norm_sub_le_one` | 超距離の小物 |
| §4 | `loss_le_of_dominated_remainder` | ★★★**到達点**。`hrem` が消え、残る仮定は `‖A‖ ≤ 1` と `hdom` だけ |

## ★残った仮定は 2 つだけ（どちらも測ってある）

- `‖A‖ ≤ 1` —— `A = σx − x` で `x` が整なら自動（超距離）。★具体層で `x` を整に取る段。
- ★`hdom`（`j₀` スロットが最大でない）—— 測定 (n15) が `p ≥ 3` で **330/330**。
  ★数学の中身はここに集まっている。`∀ z` の形に書いてあるが、★§2 の一意性より
  `j < p` の範囲で `z` は一意なので、「どの `z` でも」と「ある `z` で」は**同じこと**である。

## 逸脱の記録

- `hdom` を `∀ z` で書いた（`∃ z` より仮定として強い形）。★一意性があるので実質は同じだが、
  ★**一意性を使って `∃` 版に落とす配管は書いていない**（要るなら `coeffs_unique` で 3 行）。
- `p = 2` は対象外（前波までに機構として確定）。
-/

namespace ABC3.Found.PGC

namespace RemainderSlots

open Finset

section Coeffs

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-- ★`K` 係数のまま返す整数展開（`CoefficientIntegrality` の版は `f : ℕ → M` を返すので
`SlotResidue` の `z : ℕ → K` に渡せない）。 -/
theorem exists_integral_coeffs_of_valK [FiniteDimensional K M] {π : M} {n : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    {x : M} (hx : ‖x‖ ≤ 1) :
    ∃ c : ℕ → K, (∀ j, ‖algebraMap K M (c j)‖ ≤ 1) ∧
      x = ∑ j ∈ range n, algebraMap K M (c j) * π ^ j := by
  obtain ⟨f, hf0, hfa, hfsum⟩ :=
    ComponentBasisExpansion.exists_expansion_of_valK hπ0 (ne_of_lt hπ1) hn hvalK x
  choose c hc using hfa
  have hsum : x = ∑ j ∈ range n, algebraMap K M (c j) * π ^ j := by
    rw [hfsum]
    exact Finset.sum_congr rfl fun i _ => by rw [hc i]
  refine ⟨c, ?_, hsum⟩
  intro j
  by_cases hj : j < n
  · refine CoefficientIntegrality.norm_coeff_le_one_of_valK hπ0 hπ1 hvalK c ?_ hj
    rw [← hsum]
    exact hx
  · have hz : algebraMap K M (c j) = 0 := by
      rw [← hc j]
      exact hf0 j (by omega)
    rw [hz, norm_zero]
    exact zero_le_one

/-! ## §2 ★展開の一意性（木が「材料は在る」と言っていた段） -/

/-- ★★**展開の一意性**。`TotallyRamified.linearIndependent_of_ne_mod` から出る。

★`ComponentWitness.lean:50` の docstring は「（＝展開の**一意性**）。
★`linearIndependent_of_ne_mod` がその材料である」と言っている。★**本当だった** ——
材料はちょうどそれで足り、8 行で閉じた。 -/
theorem coeffs_unique {π : M} {n : ℕ} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≠ 1)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (c c' : ℕ → K)
    (h : ∑ j ∈ range n, algebraMap K M (c j) * π ^ j
       = ∑ j ∈ range n, algebraMap K M (c' j) * π ^ j) :
    ∀ j, j < n → c j = c' j := by
  have hli : LinearIndependent K (fun l : Fin n => π ^ (l : ℕ)) := by
    refine TotallyRamified.linearIndependent_of_ne_mod hπ0 hπ1 hvalK _
      (fun l => ((l : ℕ) : ℤ)) ?_ ?_
    · intro l; simp [zpow_natCast]
    · intro i j hij
      exact TotallyRamifiedLayer.not_dvd_sub_of_lt i.isLt j.isLt (fun hh => hij (Fin.ext hh))
  rw [Fintype.linearIndependent_iff] at hli
  intro j hj
  have hstep : ∀ l : Fin n, (c (l : ℕ) - c' (l : ℕ)) • π ^ (l : ℕ)
      = algebraMap K M (c (l : ℕ)) * π ^ (l : ℕ)
        - algebraMap K M (c' (l : ℕ)) * π ^ (l : ℕ) := by
    intro l
    rw [Algebra.smul_def, map_sub]
    ring
  have hzero : ∑ l : Fin n, (c (l : ℕ) - c' (l : ℕ)) • π ^ (l : ℕ) = 0 := by
    rw [Finset.sum_congr rfl (fun l _ => hstep l), Finset.sum_sub_distrib,
      Fin.sum_univ_eq_sum_range (fun j => algebraMap K M (c j) * π ^ j) n,
      Fin.sum_univ_eq_sum_range (fun j => algebraMap K M (c' j) * π ^ j) n, h, sub_self]
  exact sub_eq_zero.mp (hli (fun l => c (l : ℕ) - c' (l : ℕ)) hzero ⟨j, hj⟩)

/-! ## §3 超距離の小物 -/

omit [Algebra K M] in
theorem norm_sub_le_one {a b : M} (ha : ‖a‖ ≤ 1) (hb : ‖b‖ ≤ 1) : ‖a - b‖ ≤ 1 := by
  rw [sub_eq_add_neg]
  refine le_trans (IsUltrametricDist.norm_add_le_max a (-b)) ?_
  rw [norm_neg]
  exact max_le ha hb

/-! ## §4 ★接続 —— `hrem` を仮説から外す -/

/-- ★★★**本ファイルの到達点** —— `SlotResidue.loss_le_of_remainder_slots` の仮説
`hrem`（残りが `𝒪_{E₁}` スロットの和で書ける）を**定理で供給**する。

★要るのは `‖A‖ ≤ 1`（`A = σx − x` は整）だけ。`B·π^{j₀}` も整なので差も整、
そこに `CoefficientIntegrality` の整数展開を当てる。

★★`hdom`（`j₀` スロットが最大でない）は**残す**。これは測定 (n15) が
`p ≥ 3` で **330/330** と言っている条件で、★数学の中身がそこに集まっている。
`∀ z` の形にしてあるが、★§2 の一意性より `z` は `j < p` の範囲で一意なので
「どの `z` でも」と「ある `z` で」は同じことである。 -/
theorem loss_le_of_dominated_remainder [FiniteDimensional K M] {π A B : M}
    {p d j₀ jstar vf0 vfs : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hp : 2 ≤ p) (hfr : Module.finrank K M = p)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hj2 : j₀ ≤ p - 1) (hjp : j₀ < p) (hs1 : 1 ≤ jstar)
    (hmin : vf0 ≤ vfs) (hd : d = vfs + jstar) (hdvd : p ∣ vf0)
    (hB : ‖B‖ = ‖π‖ ^ (p + vf0)) (hA1 : ‖A‖ ≤ 1)
    (hdom : ∀ z : ℕ → K, (∀ j, ‖algebraMap K M (z j)‖ ≤ 1) →
        A - B * π ^ j₀ = ∑ i ∈ range p, algebraMap K M (z i) * π ^ i →
        ∃ i₀, i₀ < p ∧
          ‖algebraMap K M (z j₀) * π ^ j₀‖ < ‖algebraMap K M (z i₀) * π ^ i₀‖) :
    ‖π‖ ^ (d + (2 * p - 2)) ≤ ‖A‖ := by
  have hB1 : ‖B‖ ≤ 1 := by
    rw [hB]
    exact pow_le_one₀ (norm_nonneg _) hπ1.le
  have hBj : ‖B * π ^ j₀‖ ≤ 1 := by
    rw [norm_mul, norm_pow]
    calc ‖B‖ * ‖π‖ ^ j₀ ≤ 1 * 1 :=
          mul_le_mul hB1 (pow_le_one₀ (norm_nonneg _) hπ1.le) (by positivity) zero_le_one
      _ = 1 := one_mul 1
  obtain ⟨z, hz1, hzsum⟩ :=
    exists_integral_coeffs_of_valK hπ0 hπ1 hfr hvalK (norm_sub_le_one hA1 hBj)
  obtain ⟨i₀, hi₀, hlt⟩ := hdom z hz1 hzsum
  exact SlotResidue.loss_le_of_remainder_slots hπ0 hπ1 hp hvalK hj2 hjp hs1 hmin hd hdvd hB
    z hz1 hzsum hi₀ hlt

end Coeffs

/-! ## §5 使っている公理の一覧 -/

#print axioms exists_integral_coeffs_of_valK
#print axioms coeffs_unique
#print axioms norm_sub_le_one
#print axioms loss_le_of_dominated_remainder

end RemainderSlots

end ABC3.Found.PGC
