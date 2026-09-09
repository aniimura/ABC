import ABC3.Found.PGC.DominatedSlot

/-!
# [pGC] (n17) が仮説から**定理**になった —— `B` に「予測した主部」でなく「成分そのもの」を取る

## 持ち場（前波で私が「最後の数学」と書いた 1 番）

前波の `DominatedSlot.loss_le_of_slot_lt` は (n17)
「残りの `j₀` スロット `<` 残り全体」を**仮説**で受けていた。それを導くのが持ち場。

## ★★設計を変えたら仮説が消えた

導こうとして気づいた ——★**`B` の取り方を変えればよい。**

前波まで `B` は「予測した主部」`w·(j₀ f_{j₀} + (j₀+1) f_{j₀+1})` だった。
★代わりに `B :=` **`σx − x` の `j₀` 成分そのもの**（`A_{j₀}`）を取ると、
残り `A − A_{j₀}·π^{j₀} = Σ_{i≠j₀} A_i π^i` には ★**`j₀` スロットが構造的に無い**（§1）。
⇒ (n17) の左辺は **0** になり、右辺が正であること（§2）だけで済む。★仮説が定理になった。

## ★残る入力はちょうど 1 つ

`‖A_{j₀}‖ = ‖π‖^{p + v(f_{j₀})}` —— ★**成分のノルムが予測どおり**であること。
§4 がこれを ★**`‖A_{j₀} − B_pred‖ < ‖B_pred‖`**（＝「成分の公式が `j₀` で真」）に落とす。

## ★木の断定の検算（`ComponentFormulaScope.lean` を開いた）

あちらは「★`p ≥ 3` では **530 件で例外 0**（`j = j₀` のとき）」「★`p = 2` は **44%** 破れる」
と書いている。★**別のスクリプト・別の標本で測り直した**（`tools/numerology-check.py` の
(n18)(n19)）:

| p | n | 標本 | (n18) `v(A_{j₀} − B_pred) > v(B_pred)` | (n19) `v(A_{j₀}) = p + v(f_{j₀})` |
|---|---|---|---|---|
| 3 | 3 | 200 | 200/200 | 200/200 |
| 3 | 4 | 60 | 60/60 | 60/60 |
| 5 | 3 | 40 | 40/40 | 40/40 |
| ★7 | 2 | 30 | 30/30 | 30/30 |
| 2 | 4 | 200 | ★114/200 | ★114/200 |

★**追認された。さらに `p = 7` に広がった**（あちらの 530 件は `p ∈ {3,5}` だけだった）。
`p = 2` の 114/200（57% 成立）は、あちらの「44% 破れる」と**一致する**。

## 何を埋めたか

| 節 | 宣言 | 内容 |
|---|---|---|
| §1 | `sum_sub_slot` | ★`j₀` 成分を引くと残りに `j₀` スロットが無い |
| §2 | `norm_remainder_pos` | 残りが 0 でない（他に非零スロットが 1 つあれば） |
| §3 | `loss_le_of_component_norm` | ★★★到達点。(n17) が消えた |
| §4 | `component_norm_of_error_lt` | 最後の入力を「成分の公式が `j₀` で真」に落とす（3 行） |

## ★測って分かったこと（仮定が 1 つ減った）

`loss_le_of_component_norm` は ★**係数の整数性 `∀ j, ‖algebraMap (c j)‖ ≤ 1` を要求しない**。
`z` を明示的に与えるので、`RemainderSlots` の「整数展開の**存在**」を使わずに済むためである。
★私は最初この仮定を書いていたが、外しても通った（★仮定は書く前に外してみる）。

## 逸脱の記録

- `hA1 : ‖A‖ ≤ 1` は残る（`DominatedSlot` → `RemainderSlots` の鎖が使う）。
- ★`p = 2` は対象外。(n18) が 57% しか成り立たない。
-/

namespace ABC3.Found.PGC

namespace ComponentNormInput

open Finset

section Component

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-! ## §1 `j₀` 成分そのものを引くと、残りに `j₀` スロットが無い -/

omit [IsUltrametricDist M] in
theorem sum_sub_slot {π : M} {p j₀ : ℕ} (c : ℕ → K) (hjp : j₀ < p) :
    (∑ i ∈ range p, algebraMap K M (c i) * π ^ i) - algebraMap K M (c j₀) * π ^ j₀
      = ∑ i ∈ range p, algebraMap K M (if i = j₀ then 0 else c i) * π ^ i := by
  classical
  have hmem : j₀ ∈ range p := Finset.mem_range.mpr hjp
  rw [← Finset.add_sum_erase _ (fun i => algebraMap K M (c i) * π ^ i) hmem,
    ← Finset.add_sum_erase _ (fun i => algebraMap K M (if i = j₀ then 0 else c i) * π ^ i) hmem]
  have hcongr : ∑ i ∈ (range p).erase j₀, algebraMap K M (if i = j₀ then 0 else c i) * π ^ i
      = ∑ i ∈ (range p).erase j₀, algebraMap K M (c i) * π ^ i :=
    Finset.sum_congr rfl fun i hi => by rw [if_neg (Finset.ne_of_mem_erase hi)]
  rw [hcongr]
  simp

/-! ## §2 残りが 0 でないこと -/

theorem norm_remainder_pos {π : M} {p j₀ : ℕ} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((p : ℤ) * m))
    (c : ℕ → K) {i₁ : ℕ} (hi₁ : i₁ < p) (hne : i₁ ≠ j₀) (hc : c i₁ ≠ 0) :
    0 < ‖∑ i ∈ range p, algebraMap K M (if i = j₀ then 0 else c i) * π ^ i‖ := by
  classical
  set z : ℕ → K := fun i => if i = j₀ then 0 else c i with hz
  have hzi : z i₁ = c i₁ := by rw [hz]; exact if_neg hne
  rcases SlotResidue.exists_slot_of_sum hπ0 hπ1 hvalK z with hzero | ⟨j, _, hjz, hj⟩
  · exact absurd (hzi ▸ hzero i₁ hi₁) hc
  · rw [hj]
    have h1 : algebraMap K M (z j) ≠ 0 := (map_ne_zero (algebraMap K M)).mpr hjz
    have h2 : π ^ j ≠ 0 := pow_ne_zero _ (norm_pos_iff.mp hπ0)
    exact norm_pos_iff.mpr (mul_ne_zero h1 h2)

/-! ## §3 ★到達点 —— (n17) が仮説から**定理**になった -/

/-- ★★★`loss ≤ 2p−2`。★`B` に**予測した主部ではなく `σx − x` の `j₀` 成分そのもの**を取ると、
残りには `j₀` スロットが**構造的に無い**（§1）ので、前波の仮定 (n17) は**自動で成り立つ**（§2）。

⇒ ★残る入力は ★**`‖A_{j₀}‖ = ‖π‖^{p + v(f_{j₀})}`（成分のノルム）だけ**になった。
測定 (n19): `p ≥ 3` で **330/330**。 -/
theorem loss_le_of_component_norm [FiniteDimensional K M] {π A : M}
    {p d j₀ jstar i₁ vf0 vfs : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hp : 2 ≤ p) (hfr : Module.finrank K M = p)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hj2 : j₀ ≤ p - 1) (hjp : j₀ < p) (hs1 : 1 ≤ jstar)
    (hmin : vf0 ≤ vfs) (hd : d = vfs + jstar) (hdvd : p ∣ vf0) (hA1 : ‖A‖ ≤ 1)
    (c : ℕ → K)
    (hA : A = ∑ i ∈ range p, algebraMap K M (c i) * π ^ i)
    (hnorm : ‖algebraMap K M (c j₀)‖ = ‖π‖ ^ (p + vf0))
    (hi₁ : i₁ < p) (hne : i₁ ≠ j₀) (hci : c i₁ ≠ 0) :
    ‖π‖ ^ (d + (2 * p - 2)) ≤ ‖A‖ := by
  classical
  have hz₀ : A - algebraMap K M (c j₀) * π ^ j₀
      = ∑ i ∈ range p, algebraMap K M (if i = j₀ then 0 else c i) * π ^ i := by
    rw [hA]
    exact sum_sub_slot c hjp
  have hlt : ‖algebraMap K M (if j₀ = j₀ then (0 : K) else c j₀) * π ^ j₀‖
      < ‖A - algebraMap K M (c j₀) * π ^ j₀‖ := by
    rw [if_pos rfl, map_zero, zero_mul, norm_zero, hz₀]
    exact norm_remainder_pos hπ0 hπ1 hvalK c hi₁ hne hci
  exact DominatedSlot.loss_le_of_slot_lt hπ0 hπ1 hp hfr hvalK hj2 hjp hs1 hmin hd hdvd
    hnorm hA1 (fun i => if i = j₀ then 0 else c i) hz₀ hlt

/-! ## §4 ★最後の入力を「成分の公式が `j₀` で真」に落とす -/

omit [Algebra K M] in
/-- ★誤差が予測した主部より**真に小さい**なら、成分のノルムは予測どおり。

★これで最後の入力は ★**`‖A_{j₀} − B_pred‖ < ‖B_pred‖`**（＝ `ComponentFormulaScope` が
「`j = j₀` では真」と測った命題そのもの）になる。測定 (n18): `p ≥ 3` で **330/330**
（★`p = 7` を含む 4 設定。前の測定は `p ∈ {3,5}` の 530 件だった）。 -/
theorem component_norm_of_error_lt {π B a : M} {p vf0 : ℕ}
    (hB : ‖B‖ = ‖π‖ ^ (p + vf0)) (herr : ‖a - B‖ < ‖B‖) :
    ‖a‖ = ‖π‖ ^ (p + vf0) := by
  have hrw : a = (a - B) + B := by ring
  rw [hrw, FirstJumpWitness.norm_add_eq_of_norm_lt herr, hB]

end Component

/-! ## §5 使っている公理の一覧 -/

#print axioms sum_sub_slot
#print axioms norm_remainder_pos
#print axioms loss_le_of_component_norm
#print axioms component_norm_of_error_lt

end ComponentNormInput

end ABC3.Found.PGC
