import ABC3.Found.PGC.LossExponentMatch
import ABC3.Found.PGC.TotallyRamifiedLayer

/-!
# [pGC] 展開の存在 —— `x = Σ_{j<p} f_j π^j` は**仮説ではなく定理**

## 持ち場（前波で私が挙げた「残り 2 点」の 2 番目）

`ComponentExpansion.lean` / `LossExponentMatch.lean` は展開 `x = Σ_{j<n} f_j π^j` を
**仮説で受けていた**。本ファイルはそれを**供給する**。

## ★木の断定を測った（`TotallyRamifiedLayer.lean` を開いた）

`TotallyRamifiedLayer.adjoin_eq_top_of_valK`（`:134`）の docstring は
「★★★**`htop` の供給** —— 底が全分岐なら `M = K(π)`」と言っている。★**本当だった**が、
★**出口に要るのはこれではない**。あちらが出すのは `Algebra.adjoin K {π} = ⊤`
（**体として**生成）であり、係数を取り出せない。★要るのは
`Submodule.span K (range fun l : Fin n => π^l) = ⊤`（**加群として**生成）の方で、
それは同じ証明の**途中**にある（`hli.span_eq_top_of_card_eq_finrank' hcard`）。
⇒ 本ファイルはその途中の項をそのまま使って**係数を取り出す**。

★在庫の測定: `grep -rn "mem_span_range_iff_exists_fun\|span_eq_top_of_card_eq_finrank"
lean/ABC3/ --include=*.lean` —— `Found/CorrHyp/FieldLimit.lean` と `Found/Falt1/Section1.lean`
に用例があるが、★**`PGC` には無かった**（本ファイルが最初）。

## 何を埋めたか

| 節 | 宣言 | 内容 |
|---|---|---|
| §1 | `exists_coeffs_of_span` | ★抽象核。冪が張れば `x = Σ algebraMap(c l)·π^l`。★分岐も付値も出ない |
| §1 | `exists_expansion_of_span` | ★`ℕ` 添字に伸ばした形。`n` 以上を `0` にするので **`f n = 0` が自動で付く**（`ComponentExpansion` が要求する形） |
| §2 | `exists_expansion_of_valK` | ★★**展開の存在**。仮定は `0 < ‖π‖`、`‖π‖ ≠ 1`、`finrank = n`、`Γ_K ⊆ ‖π‖^{nℤ}` だけ |
| §3 | `exists_dvd_exponent` | 整な元の指数は非負に取れる ⇒ `‖b‖ = ‖π‖^{n·k}` |
| §3 | `dvd_exponent_of_norm_eq` | ★`LossExponentMatch` の `hdvd : p ∣ vf0` が出る形 |

## ★仮定は塔で満たされるか（測定）

- `Γ_{E₁} ⊆ ‖π‖^{pℤ}`（`hvalK`）—— `tools/numerology-check.py` の (n2) が
  ★**1060/1060 件で `p ∣ v(f_j)`** と言っている（5 設定）。
- `finrank E₁ L = p` —— 円分塔の標準事実（`[ℚ_p(ζ_{p^n}) : ℚ_p(ζ_{p^{n-1}})] = p`）。
  ★本ファイルでは仮説 `hn` のまま受けている（具体層の代入はしていない）。

## 逸脱の記録

- `f j` は `M` の元として返す（`∀ j, ∃ a : K, f j = algebraMap K M a` を併せて返す）。
  ★`ComponentExpansion` が `f : ℕ → M` を要求するのに合わせた。
- ★整数性（`‖f j‖ ≤ 1`）は**出していない**。`x` が整でも係数が整とは限らない
  （一般には `𝒪_M ≠ 𝒪_K[π]`。★全分岐なら等しいが、その段はここでは証明していない）。
  ⇒ 下流の `hf : ∀ j, ‖f j‖ ≤ C` は**まだ仮説**である。★ここが次の 1 点。
-/

namespace ABC3.Found.PGC

namespace ComponentBasisExpansion

open Finset

/-! ## §1 抽象核 —— 冪が張るなら、任意の元は冪の和に書ける -/

section Span

theorem exists_coeffs_of_span {K M : Type*} [Field K] [Field M] [Algebra K M]
    {π : M} {n : ℕ}
    (hspan : Submodule.span K (Set.range (fun l : Fin n => π ^ (l : ℕ))) = ⊤) (x : M) :
    ∃ c : Fin n → K, ∑ l : Fin n, algebraMap K M (c l) * π ^ (l : ℕ) = x := by
  have hx : x ∈ Submodule.span K (Set.range (fun l : Fin n => π ^ (l : ℕ))) := by
    rw [hspan]; exact Submodule.mem_top
  obtain ⟨c, hc⟩ := (Submodule.mem_span_range_iff_exists_fun K).mp hx
  refine ⟨c, ?_⟩
  simpa only [Algebra.smul_def] using hc

/-- ★`ℕ` 添字に伸ばした形（`ComponentExpansion` が要求する形）。
`n` 以上では `0` に伸ばすので `f n = 0` が自動で付く。 -/
theorem exists_expansion_of_span {K M : Type*} [Field K] [Field M] [Algebra K M]
    {π : M} {n : ℕ}
    (hspan : Submodule.span K (Set.range (fun l : Fin n => π ^ (l : ℕ))) = ⊤) (x : M) :
    ∃ f : ℕ → M, (∀ j, n ≤ j → f j = 0) ∧ (∀ j, ∃ a : K, f j = algebraMap K M a) ∧
      x = ∑ j ∈ range n, f j * π ^ j := by
  classical
  obtain ⟨c, hc⟩ := exists_coeffs_of_span hspan x
  refine ⟨fun j => if h : j < n then algebraMap K M (c ⟨j, h⟩) else 0, ?_, ?_, ?_⟩
  · intro j hj
    exact dif_neg (by omega)
  · intro j
    by_cases h : j < n
    · exact ⟨c ⟨j, h⟩, by simp [h]⟩
    · exact ⟨0, by simp [h]⟩
  · rw [← hc, ← Fin.sum_univ_eq_sum_range
      (fun j => (if h : j < n then algebraMap K M (c ⟨j, h⟩) else 0) * π ^ j) n]
    exact Finset.sum_congr rfl fun l _ => by simp [l.isLt]

end Span

/-! ## §2 具体層 —— 全分岐なら展開は**存在する**（仮説ではない） -/

section ValK

/-- ★★**展開の存在**。底が全分岐（`Γ_K ⊆ ‖π‖^{nℤ}`、`n = [M:K]`）なら、
`M` の任意の元は `Σ_{j<n} f_j π^j`（`f_j ∈ K`）に書ける。

★`TotallyRamifiedLayer.adjoin_eq_top_of_valK` と**同じ 1 次独立**を使うが、
あちらは `adjoin = ⊤`（体として生成）を出す。★こちらは `span = ⊤`（加群として生成）から
**係数を取り出す**。★出口に要るのは後者である。 -/
theorem exists_expansion_of_valK {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M]
    [Algebra K M] [FiniteDimensional K M] {π : M} {n : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≠ 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m)) (x : M) :
    ∃ f : ℕ → M, (∀ j, n ≤ j → f j = 0) ∧ (∀ j, ∃ a : K, f j = algebraMap K M a) ∧
      x = ∑ j ∈ range n, f j * π ^ j := by
  have hli : LinearIndependent K (fun l : Fin n => π ^ (l : ℕ)) := by
    refine TotallyRamified.linearIndependent_of_ne_mod hπ0 hπ1 hvalK _
      (fun l => ((l : ℕ) : ℤ)) ?_ ?_
    · intro l; simp [zpow_natCast]
    · intro i j hij
      exact TotallyRamifiedLayer.not_dvd_sub_of_lt i.isLt j.isLt (fun h => hij (Fin.ext h))
  have hcard : Fintype.card (Fin n) = Module.finrank K M := by simp [hn]
  exact exists_expansion_of_span (hli.span_eq_top_of_card_eq_finrank' hcard) x

end ValK

/-! ## §3 ★`p ∣ v(f_j)` を供給する（`LossExponentMatch` の `hdvd`） -/

section Dvd

/-- ★整数の元なら、`Γ_K ⊆ ‖π‖^{nℤ}` の指数は**非負**に取れる。
⇒ `‖b‖ = ‖π‖^{n·k}`（`k : ℕ`）となり、★`n ∣ 指数` が出る。 -/
theorem exists_dvd_exponent {M : Type*} [NormedField M] {π b : M} {n : ℕ} {m : ℤ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : 0 < n)
    (h : ‖b‖ = ‖π‖ ^ ((n : ℤ) * m)) (hb : ‖b‖ ≤ 1) :
    ∃ k : ℕ, ‖b‖ = ‖π‖ ^ (n * k) := by
  have hnm : 0 ≤ (n : ℤ) * m := by
    rw [h] at hb
    exact (zpow_le_one_iff_right_of_lt_one₀ hπ0 hπ1).mp hb
  have hm : 0 ≤ m := by
    by_contra hcon
    have : (n : ℤ) * m < 0 := mul_neg_of_pos_of_neg (by exact_mod_cast hn) (not_le.mp hcon)
    omega
  refine ⟨m.toNat, ?_⟩
  rw [h, ← zpow_natCast ‖π‖ (n * m.toNat)]
  congr 1
  push_cast [Int.toNat_of_nonneg hm]
  ring

/-- ★上の帰結 —— `LossExponentMatch` の `hdvd : p ∣ vf0` が出る形。 -/
theorem dvd_exponent_of_norm_eq {M : Type*} [NormedField M] {π b : M} {n k : ℕ}
    (h : ‖b‖ = ‖π‖ ^ (n * k)) : ∃ e : ℕ, ‖b‖ = ‖π‖ ^ e ∧ n ∣ e :=
  ⟨n * k, h, Dvd.intro k rfl⟩

end Dvd

/-! ## §4 使っている公理の一覧 -/

#print axioms exists_coeffs_of_span
#print axioms exists_expansion_of_span
#print axioms exists_expansion_of_valK
#print axioms exists_dvd_exponent
#print axioms dvd_exponent_of_norm_eq

end ComponentBasisExpansion

end ABC3.Found.PGC
