import ABC3.Found.PGC.MonicRelation

/-!
# [pGC] ★★出口を測ったら「模型の構成」は最後の 1 歩ではなかった —— `w` は自由、定数は `跳び + p − 2`

## 持ち場（前波で私が数えた「模型の構成」18 本）

★**着手前に出口を開いた**（規律 3）。`AxTowerDecay.lean:473` の `AxWildDescent` を読んだ:

```
def AxWildDescent (K : PAdicLocalField p) (c : ℕ → ℝ) : Prop :=
  ∀ ε : ℝ, 0 ≤ ε → ∀ x : K.closure, (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) →
    p ∣ (minpoly K.carrier x).natDegree → ∃ x' : K.closure, …
```

★★**出口は `K.closure` の任意の `x` についての主張である** —— 円分塔についてではない。
⇒ ★**円分塔の模型を作っても出口は閉じない。** 私の測定（`hform-*.py` 系）も、
本日の鎖も、すべて**円分塔**（`ℚ_p(ζ_{p^n}) ⊃ ℚ_p(ζ_{p^{n−1}})`）で行われている。
★**この差を先に測るべきだった**（前波で「模型を作れば閉じる」と書いたのは**言い過ぎ**）。

## ★では、鎖のどこが円分に依存しているのか（測った）

| 入力 | 円分依存か | 実際 |
|---|---|---|
| `hvalK`（全分岐） | ★**しない** | `TotallyRamified.exists_zpow_norm_intermediate` は任意の中間層 |
| `hfr : finrank = p` | ★**しない** | `PowerSpanTop`（モニック関係式のみ） |
| `hw : ρ = w(1+π)` | ★**しない**（本ファイル §1） | `1+π` は単数なので `w := ρ(1+π)⁻¹`、しかも `‖w‖ = ‖ρ‖` |
| ★`hwn : ‖w‖ = ‖π‖^p` | ★**する** | §1 より `‖w‖ = ‖ρ‖` なので、これは ★**`v(σπ − π) = p`（跳びが `p`）**と同値 |
| 展開・スロット・剰余 | しない | すべて `hvalK` と `finrank` から |

⇒ ★★**円分に依存しているのはただ 1 つ、「跳びが `p` にある」だけ**である。

## ★★定数の一般形（§2）

跳びを `t := v(ρ)` とすると、同じ計算で

`v(B) + j₀ = t + v(f_{j₀}) + j₀ ≤ d + (t + p − 2)`

★**定数は `t + p − 2`** であって `2p − 2` ではない。★円分の層は `t = p`（測定 (n1): 5 設定
すべてで `v(ρ) = p`）なので `2p − 2` に戻る（§2 `jump_p_gives_two_p_sub_two`）。

★これは `LossExponentMatch.exponent_bound` の**一般化**である（あちらは `t = p` に固定）。

## 何を埋めたか

| 節 | 宣言 | 内容 |
|---|---|---|
| §1 | `norm_one_add_of_lt_one` | `‖π‖ < 1 ⇒ ‖1+π‖ = 1` |
| §1 | `exists_w_of_lt_one` | ★`w` は仮説ではない（存在とノルムだけが要る） |
| §1 | `jump_at_p_iff` | ★`hwn ⟺ 跳びが `p`` |
| §2 | `exponent_bound_general` | ★定数の一般形 `t + p − 2` |
| §2 | `jump_p_gives_two_p_sub_two` | `t = p` で `2p−2` |

## ★★★出口までに残っているもの（★前波の「18 本」を訂正する）

1. ★**跳びが `p` でない層**（一般の wild な `x`）。定数は `t + p − 2` になり、
   ★`t` は `1 ≤ t ≤ p·e_K/(p−1)` の範囲で動く（古典的な上界）。
   ⇒ `AxWildDescent` に要る `c k` は `t` に依存する。★ここは**まだ測っていない**。
2. 模型の構成（12 インスタンス＋6 事実）—— ★**円分塔の場合の**具体化。
   出口には直結しないが、`cor_3_1` の主張を具体例で確かめるには要る。

★★**前波の「模型を作れば閉じる」は誤り**だった。★正しくは「模型を作れば**円分塔について**
閉じる。出口（任意の `x`）には跳びの一般化が要る」である。★この訂正が本波の主な成果である。

## 逸脱の記録

- §1 の `exists_w_of_lt_one` は `ρ` を任意に取る（`σπ − π` である必要がない）。
- ★`RhoFactorization`（円分の `w` の具体形）は**不要になった** —— 鎖が要求するのは
  存在とノルムだけだった。★あちらの docstring は誤っていない（具体形も真）が、
  ★**鎖に要る強さではなかった**。
-/

namespace ABC3.Found.PGC

namespace JumpAtP

/-! ## §1 `w` は仮説ではない —— `1 + π` は単数 -/

section W

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

theorem norm_one_add_of_lt_one {π : M} (hπ1 : ‖π‖ < 1) : ‖1 + π‖ = 1 := by
  rw [IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm (by rw [norm_one]; exact ne_of_gt hπ1)]
  rw [norm_one]
  exact max_eq_left hπ1.le

/-- ★★`ComponentExpansion` の `hw : ρ = w·(1+π)` は**仮説ではない**。
`1 + π` は単数（`‖1+π‖ = 1`）なので `w := ρ·(1+π)⁻¹` と取れて、しかも `‖w‖ = ‖ρ‖`。

★これで「`w` が円分塔から来る」という縛りが消える。★`RhoFactorization`（円分）は
`w` の**具体形**を与えていたが、★鎖が要求しているのは**存在とノルム**だけである。 -/
theorem exists_w_of_lt_one {π : M} (hπ1 : ‖π‖ < 1) (ρ : M) :
    ∃ w : M, ρ = w * (1 + π) ∧ ‖w‖ = ‖ρ‖ := by
  have hne : (1 + π) ≠ 0 := by
    intro h
    have hn : ‖(1 : M) + π‖ = 1 := norm_one_add_of_lt_one hπ1
    rw [h, norm_zero] at hn
    exact zero_ne_one hn
  refine ⟨ρ * (1 + π)⁻¹, ?_, ?_⟩
  · field_simp
  · rw [norm_mul, norm_inv, norm_one_add_of_lt_one hπ1, inv_one, mul_one]

/-- ★`hwn : ‖w‖ = ‖π‖^p` は「跳びが `p` にある」ことと同値。 -/
theorem jump_at_p_iff {π ρ w : M} {p : ℕ} (hπ1 : ‖π‖ < 1) (hw : ρ = w * (1 + π)) :
    ‖w‖ = ‖π‖ ^ p ↔ ‖ρ‖ = ‖π‖ ^ p := by
  rw [hw, norm_mul, norm_one_add_of_lt_one hπ1, mul_one]

end W

/-! ## §2 ★定数の一般形 —— `2p−2` は「跳びが `p`」の場合 -/

section Constant

/-- ★★`LossExponentMatch.exponent_bound`（定数 `2p−2`）の**一般形**。

跳びを `t := v(ρ)` とすると、同じ計算で `v(B) + j₀ ≤ d + (t + p − 2)` になる。
★`t = p`（円分の層。測定 (n1) で 5 設定すべて）なら `2p − 2` に戻る。 -/
theorem exponent_bound_general {p d t j₀ jstar vf0 vfs : ℕ}
    (hp : 2 ≤ p) (h2 : j₀ ≤ p - 1) (hs1 : 1 ≤ jstar)
    (hmin : vf0 ≤ vfs) (hd : d = vfs + jstar) :
    (t + vf0) + j₀ ≤ d + (t + p - 2) := by
  omega

/-- ★跳びが `p` なら定数は `2p−2`。 -/
theorem jump_p_gives_two_p_sub_two {p : ℕ} (hp : 2 ≤ p) : p + p - 2 = 2 * p - 2 := by
  omega

end Constant

/-! ## §3 使っている公理の一覧 -/

#print axioms norm_one_add_of_lt_one
#print axioms exists_w_of_lt_one
#print axioms jump_at_p_iff
#print axioms exponent_bound_general
#print axioms jump_p_gives_two_p_sub_two

end JumpAtP

end ABC3.Found.PGC
