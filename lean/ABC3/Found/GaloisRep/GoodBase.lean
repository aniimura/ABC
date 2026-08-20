import ABC3.Found.GaloisRep.SepRoots

/-!
# Galois (G1) 第 64 ブロック —— **★★★★★★強めた良い基点 `c`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★ファイバーを数えるのに要る 3 条件

`#[n]⁻¹(Q) = n²` を示すには、`Q = (c, d)` の `c` に 3 つの条件が要る:

| 条件 | 何のため |
|---|---|
| (1) `Ψ₂Sq(c) ≠ 0` | ★`Q ≠ −Q`——`R ↦ x_R` が単射になる |
| (2) `g_c` の根は単根で `ΨSq_k(r) ≠ 0` (`k ≤ n`) | ★★乗法公式が使える(全射) |
| (3) `c ≠ Φ_j(x)/ΨSq_j(x)` (`x` は退化点、`j < n`) | ★★★ファイバーに**位数 ≤ n の点が入らない**(well-defined) |

★(3) が要る理由: `ord(R) = K ≤ n` なら `nR = (n mod K)R` で、その x 座標は
`Φ_j(x_R)/ΨSq_j(x_R)`(`j = n mod K < n`)である。★★これを避ければ
そのような `R` は `[n]⁻¹(Q)` に入らない。

## ★★悪い値は有限個

★`Wr_n`(第 55)の根、`∏_{k≤n}ΨSq_k` の根、`Ψ₂Sq` の根——すべて有限。
★★その像も有限で、代数閉体は無限なので避けられる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_good_c'` | ★★★★★★**3 条件を満たす `c` の存在** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial WeierstrassCurve.Affine

universe u

variable {F : Type u} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

/-- ★★★★★★**強めた良い `c`**——3 種類の悪い値を避ける。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem exists_good_c' [IsAlgClosed F] (hΔ : W.Δ ≠ 0) (n : ℕ) (hn : 2 ≤ n)
    (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0) :
    ∃ c : F, W.Ψ₂Sq.eval c ≠ 0
      ∧ (∀ r : F, (gPoly W (n : ℤ) c).IsRoot r →
          (¬ (derivative (gPoly W (n : ℤ) c)).IsRoot r)
          ∧ (∀ k : ℕ, 1 ≤ k → k ≤ n → (W.ΨSq (k : ℤ)).eval r ≠ 0))
      ∧ (∀ x : F, (∃ k : ℕ, 1 ≤ k ∧ k ≤ n ∧ (W.ΨSq (k : ℤ)).eval x = 0) →
          ∀ j : ℕ, j < n → c ≠ (W.Φ (j : ℤ)).eval x / (W.ΨSq (j : ℤ)).eval x) := by
  classical
  have hnabs : 2 ≤ (n : ℤ).natAbs := by simpa using hn
  have h2 : (2 : F) ≠ 0 := by exact_mod_cast hchar 2 (by omega) (by omega)
  have h4 : (4 : F) ≠ 0 := by
    intro hc
    have hz : (2 : F) * 2 = 0 := by linear_combination hc
    exact h2 ((mul_eq_zero.1 hz).resolve_right h2)
  set Wr : F[X] := derivative (W.Φ (n : ℤ)) * W.ΨSq (n : ℤ)
    - W.Φ (n : ℤ) * derivative (W.ΨSq (n : ℤ)) with hWr
  have hWrne : Wr ≠ 0 := by
    rw [hWr]; exact wronskian_ne_zero W (n : ℤ) hnabs (by exact_mod_cast hchar n (by omega) le_rfl)
  set Pk : F[X] := ∏ k ∈ Finset.Icc 1 n, W.ΨSq (k : ℤ) with hPk
  have hPkne : Pk ≠ 0 := by
    rw [hPk]
    refine Finset.prod_ne_zero_iff.2 (fun k hk => ?_)
    simp only [Finset.mem_Icc] at hk
    exact W.ΨSq_ne_zero (by exact_mod_cast hchar k hk.1 hk.2)
  set badRoots : Finset F := Wr.roots.toFinset ∪ Pk.roots.toFinset with hbr
  set bad : Finset F :=
    badRoots.image (fun r => (W.Φ (n : ℤ)).eval r / (W.ΨSq (n : ℤ)).eval r)
    ∪ badRoots.biUnion (fun x => (Finset.range n).image
        (fun j : ℕ => (W.Φ (j : ℤ)).eval x / (W.ΨSq (j : ℤ)).eval x))
    ∪ W.Ψ₂Sq.roots.toFinset with hbad
  obtain ⟨c, hc⟩ := Infinite.exists_notMem_finset bad
  rw [hbad, Finset.mem_union, Finset.mem_union] at hc
  push_neg at hc
  obtain ⟨⟨hc1, hc2⟩, hc3⟩ := hc
  refine ⟨c, ?_, ?_, ?_⟩
  · intro h0
    exact hc3 (by rw [Multiset.mem_toFinset, Polynomial.mem_roots (W.Ψ₂Sq_ne_zero h4)]; exact h0)
  · intro r hr
    have hh : (W.Φ (n : ℤ)).eval r = c * (W.ΨSq (n : ℤ)).eval r := by
      have hrr := hr
      simp only [gPoly, Polynomial.IsRoot, eval_sub, eval_mul, eval_C, sub_eq_zero] at hrr
      exact hrr
    have hPSqr : (W.ΨSq (n : ℤ)).eval r ≠ 0 := PSq_ne_zero_of_root W hΔ n (by omega) c r hh
    have hcr : (W.Φ (n : ℤ)).eval r / (W.ΨSq (n : ℤ)).eval r = c := by rw [hh]; field_simp
    have hnotin : r ∉ badRoots := by
      intro hmem
      apply hc1
      have hi := Finset.mem_image_of_mem
        (fun r => (W.Φ (n : ℤ)).eval r / (W.ΨSq (n : ℤ)).eval r) hmem
      rwa [hcr] at hi
    rw [hbr, Finset.mem_union] at hnotin
    push_neg at hnotin
    constructor
    · intro hd
      refine hnotin.1 ?_
      rw [Multiset.mem_toFinset, Polynomial.mem_roots hWrne]
      exact wronskian_root_of_double W (n : ℤ) c r hr hd
    · intro k hk1 hk2 hzero
      refine hnotin.2 ?_
      rw [Multiset.mem_toFinset, Polynomial.mem_roots hPkne]
      rw [hPk, Polynomial.IsRoot, eval_prod]
      exact Finset.prod_eq_zero (Finset.mem_Icc.mpr ⟨hk1, hk2⟩) hzero
  · rintro x ⟨k, hk1, hk2, hkz⟩ j hj hceq
    refine hc2 ?_
    refine Finset.mem_biUnion.mpr ⟨x, ?_, ?_⟩
    · rw [hbr, Finset.mem_union]
      right
      rw [Multiset.mem_toFinset, Polynomial.mem_roots hPkne, hPk, Polynomial.IsRoot, eval_prod]
      exact Finset.prod_eq_zero (Finset.mem_Icc.mpr ⟨hk1, hk2⟩) hkz
    · exact hceq ▸ Finset.mem_image_of_mem _ (Finset.mem_range.mpr hj)

/-! ## ★出典の紐付け(`.src`) -/

def exists_good_c'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(ファイバーを数えるための良い基点)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
