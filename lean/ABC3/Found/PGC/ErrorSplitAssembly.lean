import ABC3.Found.PGC.MainPartCoeffs

/-!
# [pGC] `hsplit` の供給 —— ★これで**抽象側の鎖が閉じた**

## 持ち場（前波で私が挙げた残り 2 点の 1 番）

`MainPartCoeffs.loss_le_of_main_coeffs` は `hsplit`（`A = 主部 + 誤差`）を仮説で受けていた。
★誤差 `e` は「`A − 主部` の整数展開」そのものなので、
`RemainderSlots.exists_integral_coeffs_of_valK` を当てれば出る（§1、4 行）。

## ★★到達点 —— 誤差について要求するのは「全体のノルム評価 1 本」だけ

§2 `loss_le_of_global_error` は、誤差の係数を**外から渡さない**。内部で
(1) 整数展開で `e` を作り、(2) `SlotErrorBound.slot_coeff_lt_of_global` で
`j₀` スロットの評価に落とし、(3) `MainPartCoeffs.loss_le_of_main_coeffs` に渡す。

## ★★★仮説の棚卸し（★本ファイルで抽象側は尽きた）

`loss_le_of_global_error` の仮説を全部並べ、★**どこから来るか**を書く:

| 仮説 | 出どころ | 測定 |
|---|---|---|
| `hπ0` `hπ1` | `π` が素元 | —— |
| `hp` `hple` | `p` は素数、`‖p‖ < 1` | —— |
| `hfr : finrank K M = p` | ★具体層（円分塔） | 未着手 |
| `hvalK` | ★具体層（全分岐） | (n2) 1060/1060 |
| `hj1` `hj2` `hjp` | `MaxMinIndex.exists_max_min_index`（**定理**） | —— |
| `hs1` `hmin` `hd` | `d` の定義（`jstar` は `v(f_j)+j` の最小） | (n3) 530/530 |
| `hdvd : p ∣ vf0` | `ComponentBasisExpansion.dvd_exponent_of_norm_eq`（**定理**） | (n2) |
| `hA1 : ‖A‖ ≤ 1` | `A = σx − x`、`x` 整 | —— |
| `hw` `hwn : ‖w‖ = ‖π‖^p` | ★具体層（`RhoFactorization`） | (n1) 5/5 |
| `hfn` | `vf0` の定義 | —— |
| `hltf` | `j₀` の定義（`MaxMinIndex`、**定理**） | (C) 540/540 |
| `c` `hA` | `ComponentBasisExpansion.exists_expansion_of_valK`（**定理**） | —— |
| `hglob` | ★`ComponentExpansion` ＋ `LossExponentMatch` の scaled 版 | (n10) 294/330 |
| `hi₁` `hne` `hci` | 残りが 0 でない | —— |

★**「★具体層」と書いた 3 行以外は、すべて自分の定理か定義である。**
⇒ ★抽象側（`p ≥ 3` の `loss ≤ 2p−2`）の鎖は**閉じた**。残るのは円分塔への代入だけ。

## ★2 つの入口（★弱い方を選べるようにしてある）

| 定理 | 誤差についての入力 | p ≥ 3 の充足率 |
|---|---|---|
| `MainPartCoeffs.loss_le_of_main_coeffs` | `j₀` **スロット**の評価（係数を外から渡す） | ★**330/330**（(n18)） |
| ★`ErrorSplitAssembly.loss_le_of_global_error`（本ファイル） | **全体**のノルム評価 1 本 | 294/330（(n10)） |

★本ファイルの方が仮定が**強い**。★使い分け: 誤差の係数を持っているなら前波の定理を、
持っていない（全体の評価しかない）なら本ファイルを使う。
★36/330 の差は `f_0`（`E₁` 成分）の `σ` による動きで、★`j₀` スロットには波及しない。

## 逸脱の記録

- ★`hglob` を使う道は 294/330 なので、★**主定理と呼ぶのは前波の方**である。
  本ファイルは「係数を作らずに済ませたい」使い方のための道。
- ★`p = 2` は対象外（機構として確定済み）。
-/

namespace ABC3.Found.PGC

namespace ErrorSplitAssembly

open Finset

section Split

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-! ## §1 `hsplit` の供給 —— 誤差を整数展開する -/

/-- ★`A = 主部 + 誤差` の**誤差の側**を、整数展開で供給する。

`‖A − main‖ ≤ 1` なら `A − main` は `𝒪_{E₁}` スロットの和に書ける
（`RemainderSlots.exists_integral_coeffs_of_valK`）。★それが `hsplit` である。 -/
theorem exists_split_of_integral [FiniteDimensional K M] {π A main : M} {p : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hfr : Module.finrank K M = p)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hsub : ‖A - main‖ ≤ 1) :
    ∃ e : ℕ → K, (∀ j, ‖algebraMap K M (e j)‖ ≤ 1) ∧
      (A - main = ∑ i ∈ range p, algebraMap K M (e i) * π ^ i) ∧
      A = main + ∑ i ∈ range p, algebraMap K M (e i) * π ^ i := by
  obtain ⟨e, he1, hesum⟩ :=
    RemainderSlots.exists_integral_coeffs_of_valK hπ0 hπ1 hfr hvalK hsub
  refine ⟨e, he1, hesum, ?_⟩
  rw [← hesum]
  ring

/-! ## §2 ★到達点 —— 入力は「全体の誤差評価」1 本だけ -/

/-- ★★★`loss ≤ 2p−2`。★誤差については ★**全体のノルム評価 1 本**しか要求しない
（スロットの評価も、誤差の係数も、内部で作る）。

★入力の母集団: `‖A − 主部‖ ≤ ‖π‖^{v(f_{j₀}) + 2p}` は測定 (n10) で `p ≥ 3` の
**294/330**。★これは `MainPartCoeffs.loss_le_of_main_coeffs`（誤差の係数を外から渡す形、
入力は (n18) の **330/330**）より**強い**仮定である。
⇒ ★**弱い方（330/330）を使いたいときは前波の定理を直接使うこと。**
本定理は「誤差の係数を作りたくない」使い方のための道である。 -/
theorem loss_le_of_global_error [FiniteDimensional K M] {π A w : M}
    {p d j₀ jstar i₁ vf0 vfs : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hp : p.Prime) (hple : ‖(p : M)‖ < 1)
    (hfr : Module.finrank K M = p)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hj1 : 1 ≤ j₀) (hj2 : j₀ ≤ p - 1) (hjp : j₀ < p) (hs1 : 1 ≤ jstar)
    (hmin : vf0 ≤ vfs) (hd : d = vfs + jstar) (hdvd : p ∣ vf0) (hA1 : ‖A‖ ≤ 1)
    (W : K) (F : ℕ → K) (hw : w = algebraMap K M W)
    (hwn : ‖w‖ = ‖π‖ ^ p) (hfn : ‖algebraMap K M (F j₀)‖ = ‖π‖ ^ vf0)
    (hltf : ‖algebraMap K M (F (j₀ + 1))‖ < ‖algebraMap K M (F j₀)‖)
    (c : ℕ → K)
    (hA : A = ∑ i ∈ range p, algebraMap K M (c i) * π ^ i)
    (hglob : ‖A - w * ∑ i ∈ range p,
        ((i : M) * algebraMap K M (F i) + ((i + 1 : ℕ) : M) * algebraMap K M (F (i + 1)))
          * π ^ i‖ ≤ ‖π‖ ^ (vf0 + 2 * p))
    (hi₁ : i₁ < p) (hne : i₁ ≠ j₀) (hci : c i₁ ≠ 0) :
    ‖π‖ ^ (d + (2 * p - 2)) ≤ ‖A‖ := by
  have hsub1 : ‖A - w * ∑ i ∈ range p,
      ((i : M) * algebraMap K M (F i) + ((i + 1 : ℕ) : M) * algebraMap K M (F (i + 1)))
        * π ^ i‖ ≤ 1 :=
    le_trans hglob (pow_le_one₀ (norm_nonneg _) hπ1.le)
  obtain ⟨e, _, hesum, hsplit⟩ := exists_split_of_integral hπ0 hπ1 hfr hvalK hsub1
  have hlte : ‖algebraMap K M (e j₀)‖ < ‖π‖ ^ (p + vf0) := by
    refine SlotErrorBound.slot_coeff_lt_of_global hπ0 hπ1 hvalK e hjp ?_
    rw [← hesum]
    exact hglob
  exact MainPartCoeffs.loss_le_of_main_coeffs hπ0 hπ1 hp hple hfr hvalK hj1 hj2 hjp hs1
    hmin hd hdvd hA1 W F hw hwn hfn hltf c e hA hsplit hlte hi₁ hne hci

end Split

/-! ## §3 使っている公理の一覧 -/

#print axioms exists_split_of_integral
#print axioms loss_le_of_global_error

end ErrorSplitAssembly

end ABC3.Found.PGC
