import ABC3.Found.PGC.ComponentNormInput

/-!
# [pGC] 全体の誤差評価 → **スロットごと**の誤差評価 → 成分のノルム

## 持ち場（前波で私が挙げた残り 2 点の 1 番）

前波の `ComponentNormInput.loss_le_of_component_norm` は
`hnorm : ‖A_{j₀}‖ = ‖π‖^{p + v(f_{j₀})}` を仮説で受けていた。それを
`ComponentExpansion` の**全体の**誤差評価から導くのが持ち場。

## ★段は 1 本だった —— スロットは全体より大きくならない

私は「成分ごと（スロットごと）の誤差評価が要る」と書いたが、★**新しい段は要らなかった**。
`SlotResidue.slot_le_sum`（前々波）が「どのスロットも和より大きくならない」と言っているので、

`‖e_{j₀}·π^{j₀}‖ ≤ ‖Σ_i e_i π^i‖ ≤ ‖π‖^{v(f_{j₀}) + 2p}`

から `‖e_{j₀}‖ < ‖π‖^{p + v(f_{j₀})}` が出る。★効く不等式は
`p + v(f_{j₀}) + j₀ < v(f_{j₀}) + 2p ⟺ j₀ < p` の**ただ 1 本**で、
★これは前波までと**同じ不等式**である（`exponent_bound` / `error_lt_main` / 本ファイル §1）。
★`j₀ < p` が本鎖で効く 3 度目。

## 何を埋めたか

| 節 | 宣言 | 内容 |
|---|---|---|
| §1 | `slot_coeff_lt_of_global` | ★全体の評価 → `j₀` スロットの係数の評価 |
| §1 | `global_bound_rewrite` | `C·‖ρ‖²` の形を指数の形に直す（`‖ρ‖ = ‖π‖^p`） |
| §2 | `coeffs_add_of_sum_eq` | 一意性で「和の分解」→「係数の分解」 |
| §3 | `component_norm_of_slot_error` | ★`hnorm` を**導く** |
| §4 | `loss_le_of_slot_error` | ★★★到達点 |

## ★2 つの道と、その母集団（★弱い方を主定理にする）

| 入力 | 出典 | p ≥ 3 の充足率 |
|---|---|---|
| ★`‖e_{j₀}‖ < ‖π‖^{p+vf0}`（**スロット**の評価） | §3/§4 が使う | ★**330/330**（(n18)） |
| `‖Σ e_i π^i‖ ≤ ‖π‖^{vf0+2p}`（**全体**の評価） | §1 が使う | 294/330（(n10)） |

★**主定理 §4 はスロットの評価（330/330）を取る。** 全体の評価（294/330）は §1 で
「そこから導ける道」として置いた。★前波と同じ判断基準（弱い仮定を主定理にする）。
★全体の評価が 36/330 で破れるのは、前々波で特定したとおり **`f_0`（`E₁` 成分）の
`σ` による動き**が原因であり、★その破れは `j₀` スロットには波及しない（(n18) は 330/330）。

## ★測って分かったこと

§1 は `hvalK` を使う（`slot_le_sum` 経由）が、★`FiniteDimensional` は要らない。
§3 も要らない。★`FiniteDimensional` が要るのは §4（`coeffs_unique` を経由する鎖）だけである。

## 逸脱の記録

- ★`hsum : ∀ i, i < p → c i = b i + e i`（成分の分解）は**仮説**である。§2 が
  「和のレベルの分解」からそれを出す道を与えるが、★§4 では `hsum` を直接受けている
  （使う側がどちらでも渡せるように）。
- ★`p = 2` は対象外（(n18) が 57% しか成り立たない）。
-/

namespace ABC3.Found.PGC

namespace SlotErrorBound

open Finset

section Slot

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-! ## §1 全体の評価 → スロットの評価（★`j₀ < p` がここでも効く） -/

/-- ★全体の誤差が `‖π‖^{v(f_{j₀}) + 2p}` 以下なら、その `j₀` スロットの係数は
`‖π‖^{p + v(f_{j₀})}` より**真に小さい**。

★理由は `p + vf0 + j₀ < vf0 + 2p ⟺ j₀ < p` のただ 1 本。★前波までと同じ不等式である。 -/
theorem slot_coeff_lt_of_global {π : M} {p j₀ vf0 : ℕ} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((p : ℤ) * m))
    (e : ℕ → K) (hjp : j₀ < p)
    (herr : ‖∑ i ∈ range p, algebraMap K M (e i) * π ^ i‖ ≤ ‖π‖ ^ (vf0 + 2 * p)) :
    ‖algebraMap K M (e j₀)‖ < ‖π‖ ^ (p + vf0) := by
  have hslot : ‖algebraMap K M (e j₀) * π ^ j₀‖ ≤ ‖π‖ ^ (vf0 + 2 * p) :=
    le_trans (SlotResidue.slot_le_sum hπ0 hπ1 hvalK e hjp) herr
  rw [norm_mul, norm_pow] at hslot
  by_contra hcon
  have hge : ‖π‖ ^ (p + vf0) ≤ ‖algebraMap K M (e j₀)‖ := not_lt.mp hcon
  have h1 : ‖π‖ ^ (p + vf0) * ‖π‖ ^ j₀ ≤ ‖algebraMap K M (e j₀)‖ * ‖π‖ ^ j₀ :=
    mul_le_mul_of_nonneg_right hge (by positivity)
  have h2 : ‖π‖ ^ (vf0 + 2 * p) < ‖π‖ ^ (p + vf0) * ‖π‖ ^ j₀ := by
    rw [← pow_add]
    exact pow_lt_pow_right_of_lt_one₀ hπ0 hπ1 (by omega)
  exact lt_irrefl _ (lt_of_lt_of_le (lt_of_lt_of_le h2 h1) hslot)

omit [IsUltrametricDist M] [Algebra K M] in
/-- `C·‖ρ‖²` の形（`LossExponentMatch` の scaled 版が出す形）を指数の形に直す。
★具体層では `‖ρ‖ = ‖π‖^p`（測定 (n1): 5 設定すべてで `v(ρ) = p`）。 -/
theorem global_bound_rewrite {π ρ x : M} {p vf0 : ℕ}
    (hρ : ‖ρ‖ = ‖π‖ ^ p) (h : ‖x‖ ≤ ‖π‖ ^ vf0 * ‖ρ‖ ^ 2) :
    ‖x‖ ≤ ‖π‖ ^ (vf0 + 2 * p) := by
  have hrw : ‖π‖ ^ vf0 * ‖ρ‖ ^ 2 = ‖π‖ ^ (vf0 + 2 * p) := by
    rw [hρ, ← pow_mul, ← pow_add]
    ring_nf
  rw [← hrw]
  exact h

/-! ## §2 一意性で「和の分解」を「係数の分解」に落とす -/

theorem coeffs_add_of_sum_eq {π : M} {p : ℕ} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≠ 1)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((p : ℤ) * m))
    (c b e : ℕ → K)
    (h : ∑ i ∈ range p, algebraMap K M (c i) * π ^ i
       = (∑ i ∈ range p, algebraMap K M (b i) * π ^ i)
         + ∑ i ∈ range p, algebraMap K M (e i) * π ^ i) :
    ∀ i, i < p → c i = b i + e i := by
  refine RemainderSlots.coeffs_unique hπ0 hπ1 hvalK c (fun i => b i + e i) ?_
  rw [h, ← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl fun i _ => by rw [map_add]; ring

/-! ## §3 ★成分のノルムが予測どおりになる -/

/-- ★★前波の `ComponentNormInput.loss_le_of_component_norm` が要求していた
`hnorm : ‖A_{j₀}‖ = ‖π‖^{p+v(f_{j₀})}` を**導く**。

★仮定は「`j₀` スロットの誤差係数が主部より小さい」だけ。★測定 (n18) が `p ≥ 3` で
**330/330** と言っている命題そのもの。 -/
theorem component_norm_of_slot_error {π : M} {p j₀ vf0 : ℕ}
    (c b e : ℕ → K) (hsum : ∀ i, i < p → c i = b i + e i) (hjp : j₀ < p)
    (hb : ‖algebraMap K M (b j₀)‖ = ‖π‖ ^ (p + vf0))
    (hlt : ‖algebraMap K M (e j₀)‖ < ‖π‖ ^ (p + vf0)) :
    ‖algebraMap K M (c j₀)‖ = ‖π‖ ^ (p + vf0) := by
  have hsmall : ‖algebraMap K M (e j₀)‖ < ‖algebraMap K M (b j₀)‖ := by
    rw [hb]; exact hlt
  have hadd := FirstJumpWitness.norm_add_eq_of_norm_lt hsmall
  rw [hsum j₀ hjp, map_add,
    add_comm (algebraMap K M (b j₀)) (algebraMap K M (e j₀)), hadd, hb]

/-! ## §4 ★到達点 -/

/-- ★★★`loss ≤ 2p−2`。★入力は「`j₀` スロットの誤差が主部より小さい」だけになった。 -/
theorem loss_le_of_slot_error [FiniteDimensional K M] {π A : M}
    {p d j₀ jstar i₁ vf0 vfs : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hp : 2 ≤ p) (hfr : Module.finrank K M = p)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hj2 : j₀ ≤ p - 1) (hjp : j₀ < p) (hs1 : 1 ≤ jstar)
    (hmin : vf0 ≤ vfs) (hd : d = vfs + jstar) (hdvd : p ∣ vf0) (hA1 : ‖A‖ ≤ 1)
    (c b e : ℕ → K)
    (hA : A = ∑ i ∈ range p, algebraMap K M (c i) * π ^ i)
    (hsum : ∀ i, i < p → c i = b i + e i)
    (hb : ‖algebraMap K M (b j₀)‖ = ‖π‖ ^ (p + vf0))
    (hlt : ‖algebraMap K M (e j₀)‖ < ‖π‖ ^ (p + vf0))
    (hi₁ : i₁ < p) (hne : i₁ ≠ j₀) (hci : c i₁ ≠ 0) :
    ‖π‖ ^ (d + (2 * p - 2)) ≤ ‖A‖ :=
  ComponentNormInput.loss_le_of_component_norm hπ0 hπ1 hp hfr hvalK hj2 hjp hs1 hmin hd hdvd
    hA1 c hA (component_norm_of_slot_error c b e hsum hjp hb hlt) hi₁ hne hci

end Slot

/-! ## §5 使っている公理の一覧 -/

#print axioms slot_coeff_lt_of_global
#print axioms global_bound_rewrite
#print axioms coeffs_add_of_sum_eq
#print axioms component_norm_of_slot_error
#print axioms loss_le_of_slot_error

end SlotErrorBound

end ABC3.Found.PGC
