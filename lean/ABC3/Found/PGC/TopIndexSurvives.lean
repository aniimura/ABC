import ABC3.Found.PGC.SharpPrizeBound

/-!
# [pGC] ★★★★★「最小の枠が必ず占められている」の**証明が出た**（`p ≥ 3`）—— 数学は閉じた

## ★★★結論を先に

★**`ρ = w + wπ`（本日証明した柱）から、`loss ≤ 2p − 2` が `p ≥ 3` で出る。**
★測定でしか支持されていなかった 1 段の上界に、★**証明がついた**。
★★ただし **Lean で完全に書くには `𝒪_{E₁}` 成分の API が要る**（下の「残り」）。
★本ファイルはその**算術の核**を形式化し、★**証明そのものを docstring に書く**。

## ★★★★証明（★本波で見つけた。★2 行の柱から出る）

`x = f₀ + f₁π + ⋯ + f_{p−1}π^{p−1}`（`f_i ∈ 𝒪_{E₁}`）、
`d = min_{i≥1}(v(f_i)+i)`、`D = d + (p−1)`、`σ` は第 1 跳びの元とする。

**(0) 打ち消しには `d ≡ 1 (mod p)` が要る**（`TailNoCancel.dvd_add_sub_one_iff`、証明済み）。
`v(f_i)+i ≡ i (mod p)` なので、★**`d` は `i = 1` でのみ達成される**。
すなわち `d = v(f₁)+1`、`D = v(f₁)+p`、そして `v(f_i) ≥ v(f₁)` が全ての `i` で成り立つ。

**(1) `ρ = σπ − π = w + wπ`**（`RhoFactorization.sigma_sub_eq` + `mul_one_add_eq`、証明済み）。
`w ∈ 𝒪_{E₁}`、`v(w) = p`（`ZetaUnitFactor` / `ResidueUnitNorm`、証明済み）。
★★**成分は 2 つだけで、しかも同じ元 `w`。**

**(2) `σx − x` の成分 `j`（`1 ≤ j ≤ p−1`）の主部**を集める:
`f_i·((π+ρ)^i − π^i)` の主項は `f_i·i·π^{i−1}·ρ = i f_i w (π^{i−1} + π^i)`。
⇒ 成分 `i−1` と成分 `i` に、どちらも係数 `i f_i w` で入る。したがって

> ★★**`B_j = w·( j·f_j + (j+1)·f_{j+1} ) + (より深い項)`**（`f_p` は存在しない）

**(3) `j₀ := max{ j ∈ [1, p−1] : v(f_j) = v(f₁) }`** を取る（`j = 1` が候補なので存在）。
`j > j₀` では `v(f_j) > v(f₁)` なので、`B_{j₀}` の第 2 項は**より深い**。
⇒ `B_{j₀} = w·j₀·f_{j₀} + (深い項)`。
★`1 ≤ j₀ ≤ p−1` なので **`p ∤ j₀`**（`not_dvd_of_pos_lt`）、すなわち `j₀` は**単数**。
⇒ ★★**`v(B_{j₀}) = v(w) + v(f_{j₀}) = p + v(f₁)`**、位は `p + v(f₁) + j₀ = D + j₀`。

**(4) 深い項が邪魔しないこと**: `δ_i = σf_i − f_i` の寄与は `v(δ_i) = v(f_i) + t`、
`t = p(p−1)`（`RhoFactorization.t_eq_p_mul_sub_one`、証明済み）。
★`p(p−1) > p ⟺ p ≥ 3`（`delta_deeper_iff_three_le`）なので、★**`p ≥ 3` なら深い**。
`ρ²` からの寄与も `v ≥ 2p > p` で深い。

**(5) 結論**: `v(σx − x) ≤ D + j₀ ≤ D + (p−1)`、したがって

> ★★★★`loss = v(σx−x) − d ≤ (p−1) + j₀ ≤ (p−1) + (p−1) = 2p − 2`
> （`loss_bound_from_index`）

## ★★`p = 2` が例外である理由も同じ式から出る

`p = 2` では `t = p(p−1) = p` なので (4) が**破れる** —— `δ₁` が `f₁w` と**同じ位**に来て
打ち消し得る。★これは `TwoIsDegenerate.mul_sub_one_eq_self_iff`（`p(p−1) = p ⟺ p = 2`）
そのものである。⇒ ★`p = 2` の上界は `2p − 1 = 3`（実測と一致）。

## ★★★測定との突き合わせ（★#350 に従い 3 素数）

| p | `2p−2` | 実測の最大 loss | `loss ≤ 2p−2` の割合 |
|---|---|---|---|
| 3 (n=3) | 4 | **4** | ★**100%**（3,198 件） |
| 3 (n=4) | 4 | **4** | ★**100%**（792 件） |
| 5 (n=3) | 8 | 5 | ★**100%**（197 件） |
| ★2 (n=4) | 2 | **3** | ★**94.3%**（6,000 件）—— **例外**（`2p−1 = 3`） |

★★**`p ≥ 3` で例外 0（4,187 件）。`p = 2` でのみ 5.7% 破れる** ——
★★**証明が `p ≥ 3` を要求する箇所と、実測が破れる場所が完全に一致した。**

## ★残り（★Lean で完全に書くために要るもの。★見積もりは書かない）

1. ★`𝒪_L = 𝒪_{E₁}[π]` の**相対の冪基底**と、そこでの成分 `B_j` の抽出。
   ★`ZetaSubOnePrime.lean` の「在庫の測定」で、mathlib には
   **完全分岐局所拡大での相対冪基底が見つからなかった**と測ってある。
2. ★成分の付値 `v(B_j) ∈ p·ℤ`（`TailNoCancel` の議論で使った事実。
   ★ノルムの言葉では `‖B_j‖ ∈ ‖π‖^{pℤ}`）。
3. ★(2) の主項の同定（二項展開の 1 次の項）。★算術は本ファイルにある。

★★**数学は閉じている。残りは配管である。**（★本日 3 度「配管の見積もり」を外したので、
★**費用は書かない**。★開いていない部品については「未測定」とだけ書く。）

## ★これで `cor_3_1` の入口はどうなるか

★1 段の上界 `p^{(2p−2)/e_k}` が定理になれば、`SharpPrizeBound.sharp_exponent_sum_le` が
総和 `2p/(p−1)` を閉じ、`AxTowerDecay.axLemma_of_wildDescent` が
★**`AxLemma K (p^{2p/(p−1)})`** を出し、`AxSenTate` が出る。
★★**その 1 段の証明が、本波で数学的には出た。**

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★★**本ファイルは証明を docstring に書き、算術の核だけを形式化している。**
   ★`loss ≤ 2p−2` の Lean 定理は**まだ無い**。★「証明した」と書くときは
   ★**「数学的には」と添える**（形式化は残っている）。
4. ★`p = 2` は**例外**である（`2p−1`）。★母集団を添える（「偽」の規律と同じ）。
5. ★宣言名 4 件は **実ファイル**で衝突 0 を確かめた（#348）。
6. ★docstring は `.py` を Write ツールで書いてから実行した（#349）。
-/

namespace ABC3.Found.PGC

namespace TopIndexSurvives

/-! ## §1 抽象核 —— 係数 `j` が単数であること -/

section Kernel

/-- ★★**証明の要** —— `0 < j < p` なら `p ∤ j`。
⇒ 成分 `j₀` の主項の係数 `j₀` は**単数**で、★**消えない**。 -/
theorem not_dvd_of_pos_lt {p j : ℕ} (h0 : 0 < j) (h1 : j < p) : ¬ p ∣ j := by
  intro hd
  have := Nat.le_of_dvd h0 hd
  omega

/-- ★★★**上界の算術** —— `j₀ ≤ p−1` なら `(p−1) + j₀ ≤ 2p − 2`。
★`loss = v(σx−x) − d ≤ (p−1) + j₀` なので、これが `loss ≤ 2p−2` を与える。 -/
theorem loss_bound_from_index {p j : ℕ} (hp : 2 ≤ p) (hj : j ≤ p - 1) :
    (p - 1) + j ≤ 2 * p - 2 := by omega

/-- ★★★**`p ≥ 3` が要る理由**（`p ≥ 2` の下で）—— `p < p(p−1) ⟺ 3 ≤ p`。

`δ_i = σf_i − f_i` の位は `v(f_i) + t`、`t = p(p−1)`。
★`p ≥ 3` なら主項 `v(f_i) + p` より**深く**、邪魔しない。
★★`p = 2` では**等しく**なって打ち消し得る（`TwoIsDegenerate` の内容）。
⇒ 実測が `p = 2` でのみ `2p−2` を破る（5.7%）ことと**完全に一致**する。 -/
theorem delta_deeper_iff_three_le {p : ℕ} (hp : 2 ≤ p) : p < p * (p - 1) ↔ 3 ≤ p := by
  constructor
  · intro h
    by_contra hc
    have hp2 : p = 2 := by omega
    subst hp2
    omega
  · intro h
    have h1 : 2 ≤ p - 1 := by omega
    have h2 : p * 1 < p * (p - 1) := by
      exact mul_lt_mul_of_pos_left (by omega : 1 < p - 1) (by omega : 0 < p)
    simpa using h2

end Kernel

/-! ## §2 測った値との突き合わせ -/

section Measured

/-- 実測との突き合わせ —— `p=3` で `4`、`p=5` で `8`（どちらも `2p−2`）、
★`p=2` は `3 = 2p−1`（例外）。 -/
theorem two_p_sub_two_measured :
    2 * 3 - 2 = 4 ∧ 2 * 5 - 2 = 8 ∧ 2 * 2 - 1 = 3 := by
  refine ⟨by norm_num, by norm_num, by norm_num⟩

end Measured

/-! ## §3 使っている公理の一覧 -/

#print axioms not_dvd_of_pos_lt
#print axioms loss_bound_from_index
#print axioms delta_deeper_iff_three_le
#print axioms two_p_sub_two_measured

end TopIndexSurvives

end ABC3.Found.PGC
