import ABC3.Found.PGC.SlotErrorBound

/-!
# [pGC] `hsum` と `hb` の供給 —— 正体は**層の橋**（`M` の係数 ↔ `K` の係数）だった

## 持ち場（前波で私が挙げた残り 2 点の 1 番）

`SlotErrorBound.loss_le_of_slot_error` は 2 つを仮説で受けていた:
`hsum`（`c i = b i + e i` の成分分解）と `hb`（`‖B_pred‖ = ‖π‖^{p+v(f_{j₀})}`）。

## ★開いて分かったこと —— 数学ではなく**層の橋**だった

`ComponentExpansion` は係数を ★**`M` の元** `f : ℕ → M` として扱い、
`SlotErrorBound` は ★**`K` の元** `c b e : ℕ → K` として扱う。両者は同じものを指しているが、
★**Lean では別の型**である。⇒ 本波の中身は「橋を架けること」だった:

`w · Σ_i ( i·f_i + (i+1)·f_{i+1} ) π^i  =  Σ_i algebraMap( W·( i·F_i + (i+1)·F_{i+1} ) ) π^i`

（`w = algebraMap W`、`f_i = algebraMap F_i`）。★`map_mul` / `map_add` / `map_natCast` と
`Finset.mul_sum` で 4 行（§1）。★新しい数学は**出てこない**。

## 何を埋めたか

| 節 | 宣言 | 内容 |
|---|---|---|
| §1 | `main_coeffs_eq` | ★層の橋。主部を `K` 係数のスロット和に書き直す |
| §2 | `norm_main_coeff_eq` | ★`hb` の供給（`rw` 1 本、`MaxMinIndex.norm_pair_no_cancel` を使う） |
| §3 | `loss_le_of_main_coeffs` | ★★★到達点。`A = 主部 + 誤差` の形から直接 `loss ≤ 2p−2` |

## ★`hb` に要ったもの（すべて測定済み）

1. `‖w‖ = ‖π‖^p` —— 測定 (n1): ★5 設定すべてで `v(w) = p`。
2. `‖f_{j₀}‖ = ‖π‖^{vf0}` —— `vf0` の定義そのもの。
3. ★先頭の対が打ち消さないこと —— `MaxMinIndex.norm_pair_no_cancel`（`p ∤ j₀` と
   `‖f_{j₀+1}‖ < ‖f_{j₀}‖` から）。★測定 (D1): **540/540**（`p = 2` を含む）。

★`p ∤ j₀` は `1 ≤ j₀ ≤ p−1` から出る（`MaxMinIndex.not_dvd_max_min_index`）。

## ★`stmt%` を試した（本体の提案。★測って報告する）

`lean/ABC3/Meta/Stmt.lean` の `stmt%` を scratch で試した:

```lean
import ABC3.Meta.Stmt
import ABC3.Found.PGC.SlotErrorBound
example (h : stmt% ABC3.Found.PGC.SlotErrorBound.component_norm_of_slot_error) : True := trivial
```

★**通った**（9.6 秒）。ただし ★**本波では使わなかった** —— 私の鎖は他ファイルの定理を
「型として仮説に取る」のではなく `exact Foo.bar …` で**直接適用**しているので、
statement を書き写す場面が無い。★道具は動くが、この持ち場には合わなかった。

## 逸脱の記録

- `hsplit`（`A = 主部 + 誤差`）は**仮説**である。★これは誤差 `e` の**定義**でもあるので、
  `RemainderSlots.exists_integral_coeffs_of_valK` を `A − 主部` に当てれば供給できる
  （★その接続は本波では書いていない。次の 1 点）。
- `hlte`（`j₀` スロットの誤差評価）は前波 `SlotErrorBound` の入力そのまま。測定 (n18) で
  `p ≥ 3` は **330/330**。
- ★`p = 2` は対象外。
-/

namespace ABC3.Found.PGC

namespace MainPartCoeffs

open Finset

section Bridge

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-! ## §1 主部を `K` 係数のスロット和に書き直す -/

omit [IsUltrametricDist M] in
/-- ★`ComponentExpansion` の主部 `w·Σ_i (i f_i + (i+1) f_{i+1}) π^i` は、
`w` と `f_i` が `K` の元なら ★**`K` 係数のスロット和**に書ける。

★これが「層の橋」である —— `ComponentExpansion` は係数を `M` の元として扱い、
`SlotErrorBound` は `K` の元として扱う。 -/
theorem main_coeffs_eq {π w : M} {p : ℕ} (W : K) (F : ℕ → K)
    (hw : w = algebraMap K M W) :
    w * ∑ i ∈ range p,
        ((i : M) * algebraMap K M (F i) + ((i + 1 : ℕ) : M) * algebraMap K M (F (i + 1)))
          * π ^ i
      = ∑ i ∈ range p,
          algebraMap K M (W * ((i : K) * F i + ((i + 1 : ℕ) : K) * F (i + 1))) * π ^ i := by
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [hw, map_mul, map_add, map_mul, map_mul, map_natCast, map_natCast]
  ring

/-! ## §2 主部の `j₀` 係数のノルム（`hb` の供給） -/

/-- ★`hb : ‖B_pred‖ = ‖π‖^{p + v(f_{j₀})}` の供給。

要るのは `‖w‖ = ‖π‖^p`（測定 (n1): 5 設定すべてで `v(w) = p`）と
`‖f_{j₀}‖ = ‖π‖^{v(f_{j₀})}`（`vf0` の定義）と、
★先頭の対が打ち消さないこと（`MaxMinIndex.norm_pair_no_cancel`）だけ。 -/
theorem norm_main_coeff_eq {π w : M} {p j₀ vf0 : ℕ} (hp : p.Prime) (hple : ‖(p : M)‖ < 1)
    (hnd : ¬ p ∣ j₀) (W : K) (F : ℕ → K) (hw : w = algebraMap K M W)
    (hwn : ‖w‖ = ‖π‖ ^ p) (hfn : ‖algebraMap K M (F j₀)‖ = ‖π‖ ^ vf0)
    (hlt : ‖algebraMap K M (F (j₀ + 1))‖ < ‖algebraMap K M (F j₀)‖) :
    ‖algebraMap K M (W * ((j₀ : K) * F j₀ + ((j₀ + 1 : ℕ) : K) * F (j₀ + 1)))‖
      = ‖π‖ ^ (p + vf0) := by
  rw [map_mul, map_add, map_mul, map_mul, map_natCast, map_natCast, ← hw, norm_mul,
    MaxMinIndex.norm_pair_no_cancel hp hple hnd hlt, hwn, hfn, ← pow_add]

/-! ## §3 ★到達点 —— `hsum` と `hb` を供給して `loss ≤ 2p−2` -/

/-- ★★★`A = 主部 + 誤差` の形（`ComponentExpansion` が出す形）から直接 `loss ≤ 2p−2`。

★`hsum`（成分の分解）は §1 の橋 ＋ `SlotErrorBound.coeffs_add_of_sum_eq`（一意性）で
供給される。★`hb` は §2 で供給される。 -/
theorem loss_le_of_main_coeffs [FiniteDimensional K M] {π A w : M}
    {p d j₀ jstar i₁ vf0 vfs : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hp : p.Prime) (hple : ‖(p : M)‖ < 1)
    (hfr : Module.finrank K M = p)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hj1 : 1 ≤ j₀) (hj2 : j₀ ≤ p - 1) (hjp : j₀ < p) (hs1 : 1 ≤ jstar)
    (hmin : vf0 ≤ vfs) (hd : d = vfs + jstar) (hdvd : p ∣ vf0) (hA1 : ‖A‖ ≤ 1)
    (W : K) (F : ℕ → K) (hw : w = algebraMap K M W)
    (hwn : ‖w‖ = ‖π‖ ^ p) (hfn : ‖algebraMap K M (F j₀)‖ = ‖π‖ ^ vf0)
    (hltf : ‖algebraMap K M (F (j₀ + 1))‖ < ‖algebraMap K M (F j₀)‖)
    (c e : ℕ → K)
    (hA : A = ∑ i ∈ range p, algebraMap K M (c i) * π ^ i)
    (hsplit : A = (w * ∑ i ∈ range p,
        ((i : M) * algebraMap K M (F i) + ((i + 1 : ℕ) : M) * algebraMap K M (F (i + 1)))
          * π ^ i) + ∑ i ∈ range p, algebraMap K M (e i) * π ^ i)
    (hlte : ‖algebraMap K M (e j₀)‖ < ‖π‖ ^ (p + vf0))
    (hi₁ : i₁ < p) (hne : i₁ ≠ j₀) (hci : c i₁ ≠ 0) :
    ‖π‖ ^ (d + (2 * p - 2)) ≤ ‖A‖ := by
  have hnd : ¬ p ∣ j₀ := MaxMinIndex.not_dvd_max_min_index hp.two_le hj1 hj2
  have hsum : ∀ i, i < p →
      c i = W * ((i : K) * F i + ((i + 1 : ℕ) : K) * F (i + 1)) + e i := by
    refine SlotErrorBound.coeffs_add_of_sum_eq hπ0 (ne_of_lt hπ1) hvalK c _ e ?_
    rw [← main_coeffs_eq W F hw, ← hsplit, hA]
  exact SlotErrorBound.loss_le_of_slot_error hπ0 hπ1 hp.two_le hfr hvalK hj2 hjp hs1 hmin hd
    hdvd hA1 c _ e hA hsum
    (norm_main_coeff_eq hp hple hnd W F hw hwn hfn hltf) hlte hi₁ hne hci

end Bridge

/-! ## §4 使っている公理の一覧 -/

#print axioms main_coeffs_eq
#print axioms norm_main_coeff_eq
#print axioms loss_le_of_main_coeffs

end MainPartCoeffs

end ABC3.Found.PGC
