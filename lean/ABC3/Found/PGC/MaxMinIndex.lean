import ABC3.Found.PGC.ComponentFormulaScope
import ABC3.Found.PGC.ResidueUnitNorm
import ABC3.Found.PGC.FirstJumpWitness

/-!
# [pGC] `j₀`（尾のうち付値最小をとる**最大の**添字）の存在と、そこで先頭の対が打ち消さないこと

## 持ち場

前波（`ComponentFormulaScope.lean`）で、私が手で導いた成分の公式

> `B_j = w·( j·f_j + (j+1)·f_{j+1} ) + (より深い項)`

が**一般には偽**で、★`j = j₀` でだけ真だと測った。そこで本波の持ち場は
「`j₀` の存在（`Finset.max'` 等）と `v(f_{j₀+1}) > v(f_{j₀})` を型に書く」。

## ★配られた字面の検査 —— **偽ではないが、「余分な事実」ではない**

`v(f_{j₀+1}) > v(f_{j₀})` は**塔について新たに証明すべき事実ではない**。
`j₀` を「`[1,p-1]` で `v(f_j)` を最小にする**最大の**添字」と定義すれば、
これは**定義から出る**（本ファイル §1 の抽象核）。★荷は不等式ではなく
**定義の側**にある。ゆえに型に書くべきものは
「そういう `j₀` が尾の中に存在する」＋「`p ∤ j₀`」＋★**その `j₀` で先頭の対が打ち消さない**
であり、本ファイルはその 3 つを書いた。

★**注意（`Finset.max'` は使っていない）**: `max'` は `Nonempty` を**依存引数**に取るため、
`set T := …` して `rw [hTdef] at hmem` すると

```
Tactic `rewrite` failed: motive is not type correct:
  fun _a ↦ _a.max' hTne ∈ _a
```

で止まる。`Finset.exists_max_image` / `Finset.exists_min_image`（どちらも
`∃ x ∈ s, ∀ x' ∈ s, …` の形で `Nonempty` を**結論に持ち込まない**）で書くと
依存引数が消え、抽象核が 15 行に収まる。

## 何を埋めたか

| 節 | 宣言 | 中身 |
|---|---|---|
| §1 | `exists_max_of_min` | ★抽象核。`α` 線形順序・`β` 線形順序・`s : Finset α` 非空なら、`v` を最小にする最大の元がある。★分岐・付値・素数が 1 語も出ない |
| §2 | `exists_max_min_index` | 具体層。`s = Icc 1 (p-1)` を代入しただけ |
| §2 | `not_dvd_max_min_index` | `1 ≤ j₀ ≤ p-1 ⇒ p ∤ j₀`（`TopIndexSurvives.not_dvd_of_pos_lt` に代入） |
| §3 | `norm_pair_no_cancel` | ★抽象核 2。超距離ノルム体で `p ∤ j` と `‖b‖ < ‖a‖` なら `‖ j·a + (j+1)·b ‖ = ‖a‖` |
| §3 | `no_cancel_at_max_min` | 付値版（`‖f j‖ = ‖π‖^{v j}` を仮定する形） |
| §4 | `exists_no_cancel_index` | ★**付値を持ち出さない形**。仮定は `f p = 0` と「尾に非零成分が 1 つある」だけ |

★§4 が本波の到達点。`v : ℕ → ℕ` を経由せず**ノルムの大小だけ**で `j₀` を選ぶので、
★**成分が 0 になる標本でも使える**（下の測定 (B) を見よ）。

## ★測定（`tools/j0-type-check.py`、厳密整数演算、標本 540）

`j₀` を Lean の定義そのまま（尾で付値最小をとる最大の添字）に取って測った。

| p | n | 標本 | (A) 尾が丸ごと零 | (B) 尾に零成分 | (C) `v(f_{j₀+1}) > v(f_{j₀})` | (D1) 先頭の対が非打消し | (D2) 成分の公式が `j₀` で真 |
|---|---|---|---|---|---|---|---|
| 3 | 3 | 200 | 0 | 0 | 200 | 200 | 200 |
| 3 | 4 | 60 | 0 | 0 | 60 | 60 | 60 |
| 5 | 3 | 40 | 0 | 0 | 40 | 40 | 40 |
| 7 | 2 | 40 | 0 | 0 | 40 | 40 | 40 |
| 2 | 4 | 200 | 0 | 0 | 200 | ★200 | ★114（**86 が不成立**） |

## ★★自分の前波の記録を精密化する（訂正）

`ComponentFormulaScope.lean` に私は「`p=2` では成分の公式が 44% 破れる」と書いた。
★**破れている場所を今回特定した**: `p = 2` でも ★**先頭の対は 200/200 で打ち消さない**（D1）。
破れているのは「より深い項が本当に深いか」の段だけで、反例はすべて
`v(B−P) = v(P)`（**より小さいのではなく、ちょうど同じ**）だった。
⇒ ★`p = 2` の退化は `j₀` の選び方の問題ではない。`TwoIsDegenerate.lean` の
`mul_sub_one_eq_self_iff`（`p(p−1) = p ↔ p = 2`）が言っている退化と整合する。

## ★測定 —— 付値版とノルム版の差は空でない

座標の「成分 `jz` のブロック」を 0 に設計して `y` を作ると（`designed` モード）、
★`f_{jz} = 0` の標本が **50/50 件**作れた（p=3 n=3 で 30、p=5 n=3 で 20）。
この標本では付値版 `no_cancel_at_max_min` の仮定 `‖f j‖ = ‖π‖^{v j}` は
★**そもそも書けない**（`‖0‖ = 0` はどの `‖π‖^v` とも等しくない）が、
§4 のノルム版は使えて、★先頭の対は **50/50 で打ち消さなかった**。
⇒ §4 は §3 の飾りではなく、**真に広い**。

## 逸脱の記録

- `Finset.max'` を使わず `Finset.exists_max_image` で書いた（上記 motive エラーの回避）。
  型が言っていること（最小をとる最大の添字）は同じ。
- §3・§4 は「成分 `f_j` が与えられたとき」の主張である。★**成分そのものの構成
  （`σx − x` の `𝒪_{E₁}` 基底展開）はまだ供給していない**。ここが残る 1 点。
- ★配管の実害を 1 件出した（`tools/lean-idioms.md` #352 に登録）。裸の `python` は
  Store のスタブに解決されて **exit 49 で何もせず**、ヒアドキュメントの編集が
  3 回連続で消えたのに `leanfile.mjs` は `ok` を返した。本ファイルは
  Write/Edit ツールで書き直したうえで測り直してある。
-/

namespace ABC3.Found.PGC

namespace MaxMinIndex

/-! ## §1 抽象核 —— 有限集合の「最小値をとる最大の添字」 -/

section Kernel

/-- ★抽象核。分岐も付値も出てこない、純粋に有限線形順序の事実。 -/
theorem exists_max_of_min {α β : Type*} [LinearOrder α] [LinearOrder β]
    (s : Finset α) (v : α → β) (hs : s.Nonempty) :
    ∃ j₀ ∈ s, (∀ j ∈ s, v j₀ ≤ v j) ∧ (∀ j ∈ s, j₀ < j → v j₀ < v j) := by
  classical
  obtain ⟨b, hb, hbmin⟩ := s.exists_min_image v hs
  have hTne : (s.filter (fun j => v j = v b)).Nonempty :=
    ⟨b, Finset.mem_filter.mpr ⟨hb, rfl⟩⟩
  obtain ⟨j₀, hj₀T, hj₀max⟩ := (s.filter (fun j => v j = v b)).exists_max_image id hTne
  obtain ⟨hj₀s, hj₀v⟩ := Finset.mem_filter.mp hj₀T
  refine ⟨j₀, hj₀s, ?_, ?_⟩
  · intro j hj
    rw [hj₀v]
    exact hbmin j hj
  · intro j hj hlt
    rcases eq_or_lt_of_le (hbmin j hj) with heq | hgt
    · exact absurd (hj₀max j (Finset.mem_filter.mpr ⟨hj, heq.symm⟩)) (not_le.mpr hlt)
    · rw [hj₀v]; exact hgt

end Kernel

/-! ## §2 具体層 —— `j₀ ∈ [1, p-1]` -/

section Concrete

theorem exists_max_min_index {β : Type*} [LinearOrder β] {p : ℕ} (hp : 2 ≤ p) (v : ℕ → β) :
    ∃ j₀, 1 ≤ j₀ ∧ j₀ ≤ p - 1 ∧
      (∀ j, 1 ≤ j → j ≤ p - 1 → v j₀ ≤ v j) ∧
      (∀ j, j₀ < j → j ≤ p - 1 → v j₀ < v j) := by
  have hs : (Finset.Icc 1 (p - 1)).Nonempty := by
    refine ⟨1, ?_⟩
    simp only [Finset.mem_Icc]
    omega
  obtain ⟨j₀, hj₀mem, hmin, hstrict⟩ := exists_max_of_min (Finset.Icc 1 (p - 1)) v hs
  rw [Finset.mem_Icc] at hj₀mem
  refine ⟨j₀, hj₀mem.1, hj₀mem.2, ?_, ?_⟩
  · intro j hj1 hj2
    exact hmin j (Finset.mem_Icc.mpr ⟨hj1, hj2⟩)
  · intro j hlt hj2
    exact hstrict j (Finset.mem_Icc.mpr ⟨by omega, hj2⟩) hlt

theorem not_dvd_max_min_index {p j₀ : ℕ} (hp : 2 ≤ p) (h1 : 1 ≤ j₀) (h2 : j₀ ≤ p - 1) :
    ¬ p ∣ j₀ :=
  TopIndexSurvives.not_dvd_of_pos_lt (by omega) (by omega)

end Concrete

/-! ## §3 `j₀` では先頭の対が打ち消さない -/

section NoCancel

theorem norm_pair_no_cancel {M : Type*} [NormedField M] [IsUltrametricDist M]
    {p j : ℕ} (hp : p.Prime) (hple : ‖(p : M)‖ < 1) (hj : ¬ p ∣ j)
    {a b : M} (hlt : ‖b‖ < ‖a‖) :
    ‖(j : M) * a + ((j + 1 : ℕ) : M) * b‖ = ‖a‖ := by
  have hja : ‖(j : M) * a‖ = ‖a‖ := by
    rw [norm_mul, ResidueUnitNorm.norm_natCast_eq_one_ultra hp hple hj, one_mul]
  have hb1 : ‖((j + 1 : ℕ) : M)‖ ≤ 1 := IsUltrametricDist.norm_natCast_le_one M _
  have hb : ‖((j + 1 : ℕ) : M) * b‖ < ‖(j : M) * a‖ := by
    rw [norm_mul, hja]
    calc ‖((j + 1 : ℕ) : M)‖ * ‖b‖ ≤ 1 * ‖b‖ :=
          mul_le_mul_of_nonneg_right hb1 (norm_nonneg b)
      _ = ‖b‖ := one_mul _
      _ < ‖a‖ := hlt
  rw [add_comm, FirstJumpWitness.norm_add_eq_of_norm_lt hb, hja]

theorem no_cancel_at_max_min {M : Type*} [NormedField M] [IsUltrametricDist M]
    {p : ℕ} (hp : p.Prime) (hple : ‖(p : M)‖ < 1) {π : M}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (f : ℕ → M) (v : ℕ → ℕ)
    (hv : ∀ j, 1 ≤ j → j ≤ p - 1 → ‖f j‖ = ‖π‖ ^ v j) :
    ∃ j₀, 1 ≤ j₀ ∧ j₀ ≤ p - 1 ∧ ¬ p ∣ j₀ ∧
      (∀ j, 1 ≤ j → j ≤ p - 1 → ‖f j‖ ≤ ‖f j₀‖) ∧
      (j₀ + 1 ≤ p - 1 →
        ‖(j₀ : M) * f j₀ + ((j₀ + 1 : ℕ) : M) * f (j₀ + 1)‖ = ‖f j₀‖) := by
  obtain ⟨j₀, h1, h2, hmin, hstrict⟩ := exists_max_min_index hp.two_le v
  have hnd : ¬ p ∣ j₀ := not_dvd_max_min_index hp.two_le h1 h2
  refine ⟨j₀, h1, h2, hnd, ?_, ?_⟩
  · intro j hj1 hj2
    rw [hv j hj1 hj2, hv j₀ h1 h2]
    exact pow_le_pow_of_le_one hπ0.le hπ1.le (hmin j hj1 hj2)
  · intro hsucc
    refine norm_pair_no_cancel hp hple hnd ?_
    rw [hv (j₀ + 1) (by omega) hsucc, hv j₀ h1 h2]
    exact pow_lt_pow_right_of_lt_one₀ hπ0 hπ1 (hstrict (j₀ + 1) (by omega) hsucc)

end NoCancel

/-! ## §4 付値を持ち出さない形（★零成分に強い） -/

section NormForm

theorem exists_no_cancel_index {M : Type*} [NormedField M] [IsUltrametricDist M]
    {p : ℕ} (hp : p.Prime) (hple : ‖(p : M)‖ < 1) (f : ℕ → M)
    (hfp : f p = 0) {a : ℕ} (ha1 : 1 ≤ a) (ha2 : a ≤ p - 1) (hane : f a ≠ 0) :
    ∃ j₀, 1 ≤ j₀ ∧ j₀ ≤ p - 1 ∧ ¬ p ∣ j₀ ∧ f j₀ ≠ 0 ∧
      (∀ j, 1 ≤ j → j ≤ p - 1 → ‖f j‖ ≤ ‖f j₀‖) ∧
      ‖(j₀ : M) * f j₀ + ((j₀ + 1 : ℕ) : M) * f (j₀ + 1)‖ = ‖f j₀‖ := by
  obtain ⟨j₀, h1, h2, hmin, hstrict⟩ :=
    exists_max_min_index hp.two_le (fun j => -‖f j‖)
  have hnd : ¬ p ∣ j₀ := not_dvd_max_min_index hp.two_le h1 h2
  have hmin' : ∀ j, 1 ≤ j → j ≤ p - 1 → ‖f j‖ ≤ ‖f j₀‖ := by
    intro j hj1 hj2
    have h := hmin j hj1 hj2
    simpa using neg_le_neg_iff.mp h
  have hj0ne : f j₀ ≠ 0 := by
    intro hzero
    have hle := hmin' a ha1 ha2
    rw [hzero, norm_zero] at hle
    exact hane (norm_eq_zero.mp (le_antisymm hle (norm_nonneg _)))
  have hpos : 0 < ‖f j₀‖ := norm_pos_iff.mpr hj0ne
  refine ⟨j₀, h1, h2, hnd, hj0ne, hmin', ?_⟩
  refine norm_pair_no_cancel hp hple hnd ?_
  by_cases hcase : j₀ + 1 ≤ p - 1
  · have h := hstrict (j₀ + 1) (by omega) hcase
    simpa using neg_lt_neg_iff.mp h
  · have hj0p : j₀ + 1 = p := by omega
    rw [hj0p, hfp, norm_zero]
    exact hpos

end NormForm

/-! ## §5 使っている公理の一覧 -/

#print axioms exists_max_of_min
#print axioms exists_max_min_index
#print axioms not_dvd_max_min_index
#print axioms norm_pair_no_cancel
#print axioms no_cancel_at_max_min
#print axioms exists_no_cancel_index

end MaxMinIndex

end ABC3.Found.PGC
