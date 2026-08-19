import ABC3.Found.GaloisRep.AlgClosedPoint

/-!
# Galois (G1) 第 63 ブロック —— **★★★★★`Φ_n − c·ΨSq_n` は `n²` 個の相異なる根を持つ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★`#E[n] = n²` の (a)–(e) 段

| 段 | 内容 | 本ブロック |
|---|---|---|
| (a) | `g_c := Φ_n − c·ΨSq_n` は monic で次数 `n²` | ★`gPoly_monic` |
| (b) | 重根は Wronskian の根 | ★`wronskian_root_of_double` |
| (c)(d) | 良い `c` を取る | ★★★`exists_good_c` |
| (e) | `n²` 個の相異なる根 | ★★★★`gPoly_toFinset_card` |

## ★★良い `c` の作り方

★`Wr_n`(第 55、`≠ 0`)の根と、`∏_{k≤n} ΨSq_k` の根は**有限個**である。
★★その像 `Φ_n(r)/ΨSq_n(r)` を避ける `c` を取る(代数閉体は無限)。
★★★そうすると `g_c` の根 `r` は

* 重根でない(`Wr_n(r) ≠ 0`)
* `ΨSq_k(r) ≠ 0` (`1 ≤ k ≤ n`)——**乗法公式が使える**

を満たす。

## ★★★逸脱(記録、§9-398)

`∏_{k≤n} ΨSq_k ≠ 0` に mathlib の `ΨSq_ne_zero` を使うので

    ∀ k, 1 ≤ k ≤ n → (k : F) ≠ 0     (標数 0 または 標数 > n)

を仮定する。★ABC の応用は `ℚ̄` なので消費側に影響は無い。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `gPoly` | ★`Φ_n − c·ΨSq_n` |
| `gPoly_monic` / `gPoly_natDegree` | ★monic・次数 `n²` |
| `wronskian_root_of_double` | ★★重根は Wronskian の根 |
| `exists_good_c` | ★★★★★**良い `c` の存在** |
| `gPoly_toFinset_card` | ★★★★**相異なる根が `n²` 個** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial WeierstrassCurve.Affine

universe u

variable {F : Type u} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

/-- ★`g_c := Φ_n − c·ΨSq_n`。 -/
noncomputable def gPoly (W : WeierstrassCurve F) (n : ℤ) (c : F) : F[X] :=
  W.Φ n - C c * W.ΨSq n

theorem gPoly_coeff (n : ℤ) (hn : 2 ≤ n.natAbs) (c : F) :
    (gPoly W n c).coeff (n.natAbs ^ 2) = 1 := by
  have hpos : 0 < n.natAbs ^ 2 := by positivity
  have h1 : (W.ΨSq n).coeff (n.natAbs ^ 2) = 0 :=
    Polynomial.coeff_eq_zero_of_natDegree_lt
      (lt_of_le_of_lt (W.natDegree_ΨSq_le n) (Nat.sub_lt hpos Nat.one_pos))
  simp only [gPoly, Polynomial.coeff_sub, Polynomial.coeff_C_mul, h1, W.coeff_Φ n]
  ring

theorem gPoly_natDegree_le (n : ℤ) (c : F) :
    (gPoly W n c).natDegree ≤ n.natAbs ^ 2 := by
  refine le_trans (Polynomial.natDegree_sub_le _ _) ?_
  refine max_le (W.natDegree_Φ_le n) ?_
  refine le_trans (Polynomial.natDegree_mul_le) ?_
  simp only [Polynomial.natDegree_C, zero_add]
  exact le_trans (W.natDegree_ΨSq_le n) (by omega)

theorem gPoly_monic (n : ℤ) (hn : 2 ≤ n.natAbs) (c : F) : (gPoly W n c).Monic :=
  Polynomial.monic_of_natDegree_le_of_coeff_eq_one _ (gPoly_natDegree_le W n c)
    (gPoly_coeff W n hn c)

theorem gPoly_natDegree (n : ℤ) (hn : 2 ≤ n.natAbs) (c : F) :
    (gPoly W n c).natDegree = n.natAbs ^ 2 :=
  le_antisymm (gPoly_natDegree_le W n c)
    (Polynomial.le_natDegree_of_ne_zero (by rw [gPoly_coeff W n hn c]; exact one_ne_zero))

/-- ★★**重根は Wronskian の根である**。 -/
theorem wronskian_root_of_double (n : ℤ) (c r : F)
    (h1 : (gPoly W n c).IsRoot r) (h2 : (derivative (gPoly W n c)).IsRoot r) :
    (derivative (W.Φ n) * W.ΨSq n - W.Φ n * derivative (W.ΨSq n)).eval r = 0 := by
  simp only [gPoly, Polynomial.IsRoot, eval_sub, eval_mul, eval_C, derivative_sub,
    derivative_C_mul] at h1 h2
  simp only [eval_sub, eval_mul]
  linear_combination (W.ΨSq n).eval r * h2 - (derivative (W.ΨSq n)).eval r * h1

/-- ★★★★★**良い `c` が取れる**——重根も退化点も現れない。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem exists_good_c [IsAlgClosed F] (hΔ : W.Δ ≠ 0) (n : ℕ) (hn : 2 ≤ n)
    (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0) :
    ∃ c : F, ∀ r : F, (gPoly W (n : ℤ) c).IsRoot r →
      (¬ (derivative (gPoly W (n : ℤ) c)).IsRoot r)
      ∧ (∀ k : ℕ, 1 ≤ k → k ≤ n → (W.ΨSq (k : ℤ)).eval r ≠ 0) := by
  classical
  have hnabs : 2 ≤ (n : ℤ).natAbs := by simpa using hn
  set Wr : F[X] := derivative (W.Φ (n : ℤ)) * W.ΨSq (n : ℤ)
    - W.Φ (n : ℤ) * derivative (W.ΨSq (n : ℤ)) with hWr
  have hWrne : Wr ≠ 0 := by
    rw [hWr]
    refine wronskian_ne_zero W (n : ℤ) hnabs ?_
    exact_mod_cast hchar n (by omega) le_rfl
  set Pk : F[X] := ∏ k ∈ Finset.Icc 1 n, W.ΨSq (k : ℤ) with hPk
  have hPkne : Pk ≠ 0 := by
    rw [hPk]
    refine Finset.prod_ne_zero_iff.2 (fun k hk => ?_)
    simp only [Finset.mem_Icc] at hk
    refine W.ΨSq_ne_zero ?_
    exact_mod_cast hchar k hk.1 hk.2
  set bad : Finset F :=
    (Wr.roots.toFinset ∪ Pk.roots.toFinset).image
      (fun r => (W.Φ (n : ℤ)).eval r / (W.ΨSq (n : ℤ)).eval r) with hbad
  obtain ⟨c, hc⟩ := Infinite.exists_notMem_finset bad
  refine ⟨c, fun r hr => ?_⟩
  have hh : (W.Φ (n : ℤ)).eval r = c * (W.ΨSq (n : ℤ)).eval r := by
    have hrr := hr
    simp only [gPoly, Polynomial.IsRoot, eval_sub, eval_mul, eval_C, sub_eq_zero] at hrr
    exact hrr
  have hPSqr : (W.ΨSq (n : ℤ)).eval r ≠ 0 := PSq_ne_zero_of_root W hΔ n (by omega) c r hh
  have hcr : (W.Φ (n : ℤ)).eval r / (W.ΨSq (n : ℤ)).eval r = c := by rw [hh]; field_simp
  have hnotin : r ∉ Wr.roots.toFinset ∪ Pk.roots.toFinset := by
    intro hmem
    apply hc
    rw [hbad]
    have hi := Finset.mem_image_of_mem
      (fun r => (W.Φ (n : ℤ)).eval r / (W.ΨSq (n : ℤ)).eval r) hmem
    rwa [hcr] at hi
  rw [Finset.mem_union] at hnotin
  constructor
  · intro hd
    refine hnotin (Or.inl ?_)
    rw [Multiset.mem_toFinset, Polynomial.mem_roots hWrne]
    exact wronskian_root_of_double W (n : ℤ) c r hr hd
  · intro k hk1 hk2 hzero
    refine hnotin (Or.inr ?_)
    rw [Multiset.mem_toFinset, Polynomial.mem_roots hPkne]
    rw [hPk, Polynomial.IsRoot, eval_prod]
    exact Finset.prod_eq_zero (Finset.mem_Icc.mpr ⟨hk1, hk2⟩) hzero

/-! ## ★★★★根の個数 -/

theorem gPoly_roots_card [IsAlgClosed F] (n : ℤ) (hn : 2 ≤ n.natAbs) (c : F) :
    (gPoly W n c).roots.card = n.natAbs ^ 2 := by
  have hs : (gPoly W n c).Splits := IsAlgClosed.splits _
  rw [← gPoly_natDegree W n hn c]
  exact Polynomial.splits_iff_card_roots.mp hs

theorem gPoly_roots_nodup [IsAlgClosed F] (n : ℤ) (hn : 2 ≤ n.natAbs) (c : F)
    (hgood : ∀ r : F, (gPoly W n c).IsRoot r → ¬ (derivative (gPoly W n c)).IsRoot r) :
    (gPoly W n c).roots.Nodup := by
  have hne : gPoly W n c ≠ 0 := (gPoly_monic W n hn c).ne_zero
  rw [Multiset.nodup_iff_count_le_one]
  intro a
  by_contra hcon
  have hlt : 1 < (gPoly W n c).roots.count a := Nat.lt_of_not_le hcon
  rw [Polynomial.count_roots] at hlt
  obtain ⟨hr, hd⟩ := (Polynomial.one_lt_rootMultiplicity_iff_isRoot hne).mp hlt
  exact hgood a hr hd

/-- ★★★★**相異なる根がちょうど `n²` 個**。 -/
theorem gPoly_toFinset_card [IsAlgClosed F] (n : ℤ) (hn : 2 ≤ n.natAbs) (c : F)
    (hgood : ∀ r : F, (gPoly W n c).IsRoot r → ¬ (derivative (gPoly W n c)).IsRoot r) :
    (gPoly W n c).roots.toFinset.card = n.natAbs ^ 2 := by
  rw [Multiset.toFinset_card_of_nodup (gPoly_roots_nodup W n hn c hgood),
    gPoly_roots_card W n hn c]

/-! ## ★出典の紐付け(`.src`) -/

def exists_good_c.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(乗法写像の分離性——良い基点の存在)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
