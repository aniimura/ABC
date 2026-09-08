import ABC3.Found.PGC.HarithAssembly
import ABC3.Found.PGC.ConcreteNormedModelK1

/-!
# [pGC] ★★★★★`harith` の組み立てを `k = 1` の具体模型で検算した

前のファイル `HarithAssembly.lean` で `harith` を抽象層で組み上げた。
本ファイルはそれを `ℚ₂(i)((1+i)^{1/4}) / ℚ₂(i)`（`p = 2, k = 1, e = 2`）に
★**実際に代入**して確かめる。

## ★★測定したこと

1. ★★**抽象層で唯一残った逸脱 `hiso` は、具体層では 1 行で消える**。

   ```lean
   example (z : M4) : ‖z‖ = spectralNorm K1 M4 z := rfl          -- ★rfl で通る
   theorem hiso_g4 (z : M4) : ‖g4 z‖ = ‖z‖ := (spectralNorm_eq_of_equiv g4 z).symm
   ```

   ★向きだけ違う。`spectralNorm_eq_of_equiv g4 z` の型は
   `spectralNorm K1 M4 z = spectralNorm K1 M4 (g4 z)` である（`.symm` が要る）。
   ★`lean-idioms.md` #329 と同じ手であり、`NormedAlgebra ℚ_[2] M4` は不要。

2. ★★★**`harith_k1` は代入 1 行になる**。
   `ConcreteNormedModelK1.lean:561 harith_k1` は
   `u0_eq_four` / `u1_eq_eight`（具体の跡びの値 `4`, `8`）を使い、
   `interval_cases` と `decide` で 4 条件を 1 つずつ確かめていた。
   ★**本ファイルの `harith_k1_via_norm` は跡びの値を一切使わない**。
   字面は完全に同じである。

3. ★★★★**`k = 1` の出口も `harith` なしで出る**。
   `ConcreteNormedModelK1.lean:587 model_k1_exists` は `harith_k1` と
   `adjoin_root_top` を渡していたが、`model_k1_exists_via_norm` はどちらも渡さない。

## ★検算できなかったもの（`file:line`）

`JumpFromValueGroup.lean:332 harith_zeta81` は
`Zeta81.uz`（具体の数列 `(2,8,26)`）についての**数値の主張**であり、
★ノルムの仮説から出す形ではない（`ℚ₃(ζ₈₁)` の**ノルム付き模型が木にない**）。
⇒ ★**`harith_of_norm` は適用できない**。数値の一致は
`HasseArfCongruenceNorm.lean` §4 で既に確かめている（`φ = (2,4,6)`）。

## 逸脱の記録

1. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
   `harith_k1` / `model_k1_exists` は**残してある**（本ファイルは別の宣言）。
-/

namespace ABC3.Found.PGC

namespace HarithConcreteCheck

open ConcreteNormedModelK1

/-- ★**測定**: `M4` のノルムは `spectralNorm K1 M4` そのもの（`rfl`）。 -/
example (z : M4) : ‖z‖ = spectralNorm K1 M4 z := rfl

/-- ★★★**`hiso` は具体層では 1 行**。

★抽象層（`HarithAssembly.harith_of_norm`）では `hiso` を仮説で受けるしかなかったが
（`PureStepSetup.norm_algEquiv_eq` は `[NontriviallyNormedField K] [CompleteSpace K]` を要求する）、
具体層では `M4` のノルムが `spectralNorm K1 M4` なので
`spectralNorm_eq_of_equiv` が**そのまま**供給する。
★向きだけ違う（`spectralNorm K1 M4 z = spectralNorm K1 M4 (g4 z)`）ので `.symm`。 -/
theorem hiso_g4 (z : M4) : ‖g4 z‖ = ‖z‖ := (spectralNorm_eq_of_equiv g4 z).symm

/-- ★★★★★**検算** —— `ConcreteNormedModelK1.harith_k1`（`:561`）と**字面が同じ**ものが
`HarithAssembly.harith_of_norm` への★**代入 1 行**で出る。

★原型（`harith_k1`）は `u0_eq_four` / `u1_eq_eight`（具体の跡びの値）を使って
`interval_cases` と `decide` で閉じていた。
★本定理は★**跡びの値を一切使わない**（一般論だけ）。 -/
theorem harith_k1_via_norm :
    ∀ u : ℕ → ℤ, (∀ j, j ≤ 1 → ‖(s4 j) pi4 - pi4‖ = ‖pi4‖ ^ (u j + 1)) →
      1 ≤ u 0 ∧ (∀ m, m < 1 → u m < u (m + 1)) ∧
        (∀ m, m < 1 → (2 : ℤ) ^ (m + 1) ∣ u (m + 1) - u m) ∧
        ((2 : ℤ) - 1) * u 1 ≤ (2 : ℤ) ^ (1 + 1) * (2 : ℤ) :=
  HarithAssembly.harith_of_norm (p := 2) (k := 1) (e := 2) (π := pi4)
    Nat.prime_two g4 orderOf_g4 hiso_g4 (by norm_num) finrank_M4 (by simpa using valK4)
    norm_two_M4 heM4 s4 (fun _ _ => rfl)

/-- ★★★★★★**`k = 1` の具体模型の出口も `harith` なしで出る**。

`ConcreteNormedModelK1.model_k1_exists`（`:587`）は `harith_k1` と
`adjoin_root_top` を渡していたが、本定理はどちらも渡さない。 -/
theorem model_k1_exists_via_norm (x : M4) :
    ∃ y : GainedTowerModel.GaloisTower.twr g4 2 1,
      ‖x - algebraMap (GainedTowerModel.GaloisTower.twr g4 2 1) M4 y‖
        ≤ (∏ j ∈ Finset.Icc 1 (1 + 1), axDecay 2 j) * ‖tau4 x - x‖ :=
  HarithAssembly.exists_norm_sub_algebraMap_le_prod_axDecay_of_norm
    (p := 2) (e := 2) (k := 1) (π := pi4)
    g4 orderOf_g4 hiso_g4 tau4 (fun _ => rfl) s4 (fun _ _ => rfl) (by norm_num)
    finrank_M4 (by simpa using valK4) norm_two_M4 heM4 x

/-! ## `.src` と 公理 -/

def harith_k1_via_norm.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

#print axioms hiso_g4
#print axioms harith_k1_via_norm
#print axioms model_k1_exists_via_norm

end HarithConcreteCheck

end ABC3.Found.PGC
