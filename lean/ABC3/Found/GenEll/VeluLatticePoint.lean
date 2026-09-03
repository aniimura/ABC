/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluLatticeSet
import ABC3.Meta.Claim

/-!
# 第 1333 ブロック —— **`ℂ` 上の Vélu の商は楕円（原文の語彙で）**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か——`VeluQuotOK` の楕円性の `ℂ` 版

第 1330 は代表系 `T` の語彙、第 1332 は座標集合の突き合わせ。
★本ブロックはその 2 つを繋ぎ、**原文（と `VeluQuotOK` の消費側）の語彙**

  `veluQuotientFull (latticeCurve P) (((range l).erase 0).image (fun k => pointCoords (k • Q)))`

が**楕円**であることを、`Q` の位数が `l`（`l ≥ 1`）というだけで示す。

☆残るのは `ℂ` へ埋め込んで降ろす段（第 1331 `isElliptic_veluQuotientFull_of_map`）である。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve ABC3.Meta

open scoped Classical

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**`ℂ` 上：`Q` の生成する部分群での Vélu の商は楕円**——★**無条件**（第 1333）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★これが `Skeleton/GenEll/AlphaBridge.lean` の `VeluQuotOK` が要求する
**楕円性そのもの**（の `ℂ` 版）である。 -/
theorem isElliptic_veluQuotientFull_nsmul_lattice (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {Q : (latticeCurve P).toAffine.Point} {l : ℕ} (hl : 0 < l) (hQ : addOrderOf Q = l) :
    (veluQuotientFull (latticeCurve P)
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).IsElliptic := by
  classical
  obtain ⟨z₀, P', A, B, C, D, hz₀, h1, h2, hdet, hP'⟩ :=
    exists_isogeny_periodPair P hΔ hl hQ
  have hord := intCast_mul_mem_lattice_iff P hΔ hQ hz₀
  obtain ⟨T, h0T, hcard, hT, hrep, hTdef⟩ :=
    exists_velu_rep P P' z₀ l hl hP' hord
  have hle : P.lattice ≤ P'.lattice := by rw [hP']; exact le_sup_left
  have hvelu : ∀ z : ℂ, P'.weierstrassP z = veluAnalyticX P T (veluAnalyticC P T) z :=
    fun z => weierstrassP_eq_velu_of_rep P P' hle T h0T hT hrep z
  have heq := latticeCurve_eq_veluQuotientFull P P' T h0T hT hrep
      (g₂_isogeny P P' T h0T hT hrep hvelu) (g₃_isogeny P P' T h0T hT hrep hvelu)
  -- ★ 段 1: `1 ≤ k < l` では `k z₀ ∉ Λ`
  have hnot : ∀ k ∈ (Finset.range l).erase 0, ((k : ℂ)) * z₀ ∉ P.lattice := by
    intro k hk hmem
    have hk0 : k ≠ 0 := Finset.ne_of_mem_erase hk
    have hkl : k < l := Finset.mem_range.1 (Finset.mem_of_mem_erase hk)
    have hmem' : ((k : ℤ) : ℂ) * z₀ ∈ P.lattice := by push_cast; exact hmem
    have hd := (hord (k : ℤ)).1 hmem'
    have hle' : (l : ℤ) ≤ (k : ℤ) :=
      Int.le_of_dvd (by exact_mod_cast Nat.pos_of_ne_zero hk0) hd
    omega
  -- ★ 段 2: 代表系から 0 を抜いたものは `{k z₀ : 1 ≤ k < l}`
  have hTe : T.erase 0 = ((Finset.range l).erase 0).image (fun k : ℕ => (k : ℂ) * z₀) := by
    ext w
    simp only [hTdef, Finset.mem_erase, Finset.mem_image, Finset.mem_range]
    constructor
    · rintro ⟨hw0, k, hkl, rfl⟩
      refine ⟨k, ⟨?_, hkl⟩, rfl⟩
      rintro rfl
      exact hw0 (by simp)
    · rintro ⟨k, ⟨hk0, hkl⟩, rfl⟩
      refine ⟨?_, k, hkl, rfl⟩
      intro hz
      exact hnot k (Finset.mem_erase.2 ⟨hk0, Finset.mem_range.2 hkl⟩)
        (hz ▸ P.lattice.zero_mem)
  -- ★ 段 3: 座標集合の突き合わせ（第 1332）
  have hset : ((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))
      = (T.erase 0).image (fun w => (latticePointX P w, latticePointY P w)) := by
    rw [hTe, Finset.image_image, ← hz₀]
    exact image_pointCoords_nsmul_eq P hΔ z₀ l hnot
  rw [hset, ← heq]
  exact isElliptic_latticeCurve' P'

/-! ## ★出典の紐付け(`.src`) -/

def isElliptic_veluQuotientFull_nsmul_lattice.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(ℂ 上：Q の生成する部分群での Vélu の商は楕円。★無条件)",
    sectionId := "genell-lemma-3-5" }

def isElliptic_veluQuotientFull_nsmul_lattice.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_isogeny_periodPair・exists_velu_rep(在庫)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_isogeny_periodPair") 1,
    .citation "[ABC3]" "image_pointCoords_nsmul_eq(第 1332、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.image_pointCoords_nsmul_eq") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1333）**——`Skeleton/GenEll/AlphaBridge.lean` の " ++
       "`VeluQuotOK` が要求する**楕円性そのもの**（の `ℂ` 版）である。" ++
       "☆残るのは `ℂ` へ埋め込んで降ろす段（第 1331）だけである。") 3 ]

end ABC3.Found.GenEll
