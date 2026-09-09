import ABC3.Found.PGC.RemainderSlots

/-!
# [pGC] `hdom` を「`i₀` を名指ししない形」に変える ＋ 前波の「3 行」主張の検算

## 持ち場（2 候補から選んだ方）

前波で私は 2 候補を挙げた ——(A) `hdom` を成分の主部から導く、(B) 円分塔への代入。
★**(A) を選んだ。** 理由は前波で自分が「★数学の中身はここに集まっています」と判定した点だから。
（(B) は開いていないので費用は書かない。★選ばなかった理由は「安いと分かったから」ではない。）

## ★何が変わったか

`RemainderSlots.loss_le_of_dominated_remainder` の `hdom` は
「★`i₀` が**存在して** `j₀` スロットより大きい」という形で、★使う側が `i₀` を名指しする必要があった。
本ファイルは ★**`i₀` を知らなくてよい形**に置き換える:

> `‖(残りの j₀ スロット)‖ < ‖残り全体‖`

★これで十分である（§1）。最大を達成するスロットは全体と同じノルムを持つので、
`j₀` スロットが全体より小さければ、達成するスロットは自動的に `j₀` ではない。

## ★測定（`tools/numerology-check.py`、5 設定）

| 測った命題 | p ≥ 3 | p = 2 |
|---|---|---|
| (n17) ★`j₀` スロット `<` 残り全体（★本ファイルの仮定そのもの） | ★**330/330** | 188/200 |
| (n15) `j₀` スロットが最大でない（前波の形） | 330/330 | 188/200 |
| (n16) `‖err‖ < ‖R‖`（§2 の十分条件） | 294/330 | 123/200 |

★(n17) と (n15) は**同じ数**になった（理論どおり、両者は同値）。
★(n16) は真に強い仮定（294/330）なので、★**§2 は「導ける道」としてだけ置き、
主定理 §4 は (n17) の形で書いた**。

## ★★前波の「3 行」主張の検算（本体が「検算していない」と名指しした点）

前波で私はこう書いた ——「`∀ z` と `∃ z` は一意性より同じ、配管は `coeffs_unique` で 3 行」。
★**測った: 一意性の配管は 6 行**（`intro` 1 行 ＋ `coeffs_unique` の適用 2 行 ＋
`refine` 1 行 ＋ `rw` 1 行 ＋ `exact` 1 行）。★さらにその前段の
`dominated_of_slot_lt`（`i₀` を取り出す部分）が **10 行**要る。
⇒ ★**合計 16 行。「3 行」は楽観だった**（配管の見積もりを外したのは本日 3 度目）。

## 何を埋めたか

| 節 | 宣言 | 内容 |
|---|---|---|
| §1 | `dominated_of_slot_lt` | ★`j₀` スロット `<` 全体 ⇒ 上回るスロットが在る |
| §2 | `slot_lt_of_split` | `R` に `j₀` スロットが無いことから導く道（十分条件、294/330） |
| §3 | `hdom_of_slot_lt` | ★一意性で `∀ z` を供給（前波の主張の検算） |
| §4 | `loss_le_of_slot_lt` | ★★★`loss ≤ 2p−2`。★残る仮定は `‖A‖ ≤ 1` と (n17) の 2 つだけ |

## 逸脱の記録

- ★`p = 2` は対象外（機構として確定済み。(n17) も 188/200 で破れる）。
- §2 の `slot_lt_of_split` は**使っていない**（主定理は §1 経由）。★仮定が強すぎる
  （294/330）ため。★残しているのは「`R` に `j₀` スロットが無い」という構造が
  次の波で効く可能性があるからで、★今の鎖には要らない。
-/

namespace ABC3.Found.PGC

namespace DominatedSlot

open Finset

section Dominated

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-- ★`j₀` スロットが**全体より小さい**なら、それを上回るスロットが必ずある。

★`hdom`（「`j₀` スロットが最大でない」）は `i₀` を名指しする形だったが、
本定理により ★**`i₀` を知らなくてよい形**「`j₀` スロット `<` 全体」に置き換わる。 -/
theorem dominated_of_slot_lt {π : M} {p j₀ : ℕ} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((p : ℤ) * m))
    (z : ℕ → K)
    (hlt : ‖algebraMap K M (z j₀) * π ^ j₀‖
        < ‖∑ i ∈ range p, algebraMap K M (z i) * π ^ i‖) :
    ∃ i₀, i₀ < p ∧
      ‖algebraMap K M (z j₀) * π ^ j₀‖ < ‖algebraMap K M (z i₀) * π ^ i₀‖ := by
  rcases SlotResidue.exists_slot_of_sum hπ0 hπ1 hvalK z with hzero | ⟨i₀, hi₀, _, hi⟩
  · exfalso
    have : ∑ i ∈ range p, algebraMap K M (z i) * π ^ i = 0 := by
      refine Finset.sum_eq_zero fun i hi => ?_
      rw [hzero i (Finset.mem_range.mp hi), map_zero, zero_mul]
    rw [this, norm_zero] at hlt
    exact absurd hlt (not_lt.mpr (norm_nonneg _))
  · exact ⟨i₀, hi₀, by rw [← hi]; exact hlt⟩

/-! ## §2 分解から出す —— `R` に `j₀` スロットが無いこと -/

/-- ★`R`（主部の `j₀` 以外の部分）は **`j₀` スロットを持たない**（`r j₀ = 0`）。
誤差 `e` が `R` より深ければ、残り `r + e` の `j₀` スロットは全体より小さい。 -/
theorem slot_lt_of_split {π : M} {p j₀ : ℕ} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((p : ℤ) * m))
    (r e : ℕ → K) (hjp : j₀ < p) (hr0 : r j₀ = 0)
    (hlt : ‖∑ i ∈ range p, algebraMap K M (e i) * π ^ i‖
        < ‖∑ i ∈ range p, algebraMap K M (r i) * π ^ i‖) :
    ‖algebraMap K M (r j₀ + e j₀) * π ^ j₀‖
      < ‖∑ i ∈ range p, algebraMap K M (r i + e i) * π ^ i‖ := by
  have hsplit : ∑ i ∈ range p, algebraMap K M (r i + e i) * π ^ i
      = (∑ i ∈ range p, algebraMap K M (r i) * π ^ i)
        + ∑ i ∈ range p, algebraMap K M (e i) * π ^ i := by
    rw [← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun i _ => by rw [map_add]; ring
  have hnorm : ‖(∑ i ∈ range p, algebraMap K M (r i) * π ^ i)
        + ∑ i ∈ range p, algebraMap K M (e i) * π ^ i‖
      = ‖∑ i ∈ range p, algebraMap K M (r i) * π ^ i‖ := by
    rw [add_comm]
    exact FirstJumpWitness.norm_add_eq_of_norm_lt hlt
  rw [hsplit, hnorm, hr0, zero_add]
  exact lt_of_le_of_lt (SlotResidue.slot_le_sum hπ0 hπ1 hvalK e hjp) hlt

/-! ## §3 ★一意性で `∀ z` を供給する（★前波の「3 行」主張の検算） -/

/-- ★前波で私は「`∀ z` と `∃ z` は一意性より同じ、配管は `coeffs_unique` で 3 行」と書いた。
★**測った: 6 行**（`intro` ＋ `coeffs_unique` 2 行 ＋ `refine` ＋ `rw` ＋ `exact`）。
★さらに §1 の `dominated_of_slot_lt`（10 行）が要るので**合計 16 行**。
★「3 行」は楽観だった。 -/
theorem hdom_of_slot_lt [FiniteDimensional K M] {π A B : M} {p j₀ : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hjp : j₀ < p) (z₀ : ℕ → K)
    (hz₀ : A - B * π ^ j₀ = ∑ i ∈ range p, algebraMap K M (z₀ i) * π ^ i)
    (hlt : ‖algebraMap K M (z₀ j₀) * π ^ j₀‖ < ‖A - B * π ^ j₀‖) :
    ∀ z : ℕ → K, (∀ j, ‖algebraMap K M (z j)‖ ≤ 1) →
      A - B * π ^ j₀ = ∑ i ∈ range p, algebraMap K M (z i) * π ^ i →
      ∃ i₀, i₀ < p ∧
        ‖algebraMap K M (z j₀) * π ^ j₀‖ < ‖algebraMap K M (z i₀) * π ^ i₀‖ := by
  intro z _ hz
  have heq : ∀ j, j < p → z j = z₀ j :=
    RemainderSlots.coeffs_unique hπ0 (ne_of_lt hπ1) hvalK z z₀ (by rw [← hz, hz₀])
  refine dominated_of_slot_lt hπ0 hπ1 hvalK z ?_
  rw [heq j₀ hjp, ← hz]
  exact hlt

/-! ## §4 ★到達点 —— 仮定は「`j₀` スロットが全体より小さい」だけ -/

/-- ★★★`loss ≤ 2p−2`。★残る仮定は `‖A‖ ≤ 1`（`x` が整）と
★**「残りの `j₀` スロットが残り全体より小さい」**の 2 つだけ。 -/
theorem loss_le_of_slot_lt [FiniteDimensional K M] {π A B : M}
    {p d j₀ jstar vf0 vfs : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hp : 2 ≤ p) (hfr : Module.finrank K M = p)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hj2 : j₀ ≤ p - 1) (hjp : j₀ < p) (hs1 : 1 ≤ jstar)
    (hmin : vf0 ≤ vfs) (hd : d = vfs + jstar) (hdvd : p ∣ vf0)
    (hB : ‖B‖ = ‖π‖ ^ (p + vf0)) (hA1 : ‖A‖ ≤ 1)
    (z₀ : ℕ → K)
    (hz₀ : A - B * π ^ j₀ = ∑ i ∈ range p, algebraMap K M (z₀ i) * π ^ i)
    (hlt : ‖algebraMap K M (z₀ j₀) * π ^ j₀‖ < ‖A - B * π ^ j₀‖) :
    ‖π‖ ^ (d + (2 * p - 2)) ≤ ‖A‖ :=
  RemainderSlots.loss_le_of_dominated_remainder hπ0 hπ1 hp hfr hvalK hj2 hjp hs1 hmin hd
    hdvd hB hA1 (hdom_of_slot_lt hπ0 hπ1 hvalK hjp z₀ hz₀ hlt)

end Dominated

/-! ## §5 使っている公理の一覧 -/

#print axioms dominated_of_slot_lt
#print axioms slot_lt_of_split
#print axioms hdom_of_slot_lt
#print axioms loss_le_of_slot_lt

end DominatedSlot

end ABC3.Found.PGC
