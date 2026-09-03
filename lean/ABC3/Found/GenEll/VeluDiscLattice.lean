/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.LatticeDiscPow
import ABC3.Found.GenEll.VeluLatticePoint
import ABC3.Found.GaloisRep.VeluKernelNorm
import ABC3.Meta.Claim

/-!
# 第 1402 ブロック —— **原文の語彙での判別式の恒等式（格子曲線）**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★これは何か

第 1401 の

    latticeDisc Λ ^ l = latticeDisc Λ′ · ( ∏_{w∈T∖0} ℘′(w) )^4

を、**消費側（`Skeleton/GenEll/VeluDiscIdentity.lean`）の語彙**

    Δ(latticeCurve Λ)^l = Δ(veluQuotientFull (latticeCurve Λ) ⟨Q⟩) · N^4

に直す。☆道具は第 1333（`isElliptic_veluQuotientFull_nsmul_lattice`）と同じ:

| 段 | 内容 | 出どころ |
|---|---|---|
| 1 | `Λ′ = Λ + ℤz₀`・代表系 `T = {kz₀}` | 在庫（第 660-663） |
| 2 | 座標集合の突き合わせ `⟨Q⟩ ↔ T∖{0}` | 第 1332 `image_pointCoords_nsmul_eq` |
| 3 | `Λ′/Λ` に 2-捩れが無いこと | ★本ブロック（`l` が奇であること） |
| 4 | `veluKernelNorm = ∏ ℘′(w)` | ★本ブロック（`a₁ = a₃ = 0`） |

★★★これで `Skeleton/GenEll/VeluDiscIdentity.lean` の `sorry` が埋まる。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve PeriodPair Finset ABC3.Meta ABC3.Found.GaloisRep

open scoped Classical

/-- ★★★★★★★★
**`Λ′ = Λ + ℤz₀` で `l·z₀ ∈ Λ`（`l` 奇）なら `Λ′/Λ` に 2-捩れは無い**（第 1402）。 -/
theorem hodd2_of_span (P P' : PeriodPair) (z₀ : ℂ) {l m : ℕ} (hlm : l = 2 * m + 1)
    (hP' : P'.lattice = P.lattice ⊔ (Submodule.span ℤ {z₀}))
    (hord : ∀ k : ℤ, (k : ℂ) * z₀ ∈ P.lattice ↔ (l : ℤ) ∣ k) :
    ∀ y ∈ P'.lattice, 2 * y ∈ P.lattice → y ∈ P.lattice := by
  intro y hy h2y
  rw [hP', Submodule.mem_sup] at hy
  obtain ⟨a, ha, b, hb, hab⟩ := hy
  obtain ⟨n, hn⟩ := Submodule.mem_span_singleton.1 hb
  have hbz : b = (n : ℂ) * z₀ := by rw [← hn]; simp [zsmul_eq_mul]
  have h2b : (2 : ℂ) * b ∈ P.lattice := by
    have he : (2 : ℂ) * b = 2 * y - 2 * a := by rw [← hab]; ring
    rw [he]
    exact P.lattice.sub_mem h2y (by
      have h2a : (2:ℂ) * a = a + a := by ring
      rw [h2a]; exact P.lattice.add_mem ha ha)
  have h2n : ((n * 2 : ℤ) : ℂ) * z₀ ∈ P.lattice := by
    have he : ((n * 2 : ℤ) : ℂ) * z₀ = 2 * b := by rw [hbz]; push_cast; ring
    rw [he]; exact h2b
  have hdvd : (l : ℤ) ∣ n * 2 := (hord (n * 2)).1 h2n
  have hcop : IsCoprime ((l : ℤ)) 2 := ⟨1, -(m : ℤ), by rw [hlm]; push_cast; ring⟩
  have hln : (l : ℤ) ∣ n := hcop.dvd_of_dvd_mul_right hdvd
  have hnz : (n : ℂ) * z₀ ∈ P.lattice := (hord n).2 hln
  rw [← hab, hbz]
  exact P.lattice.add_mem ha hnz

/-- ★★★★★★**格子曲線では核のノルムは `∏ ℘′(w)`**——`a₁ = a₃ = 0` だから（第 1402）。 -/
theorem veluKernelNorm_lattice (P P' : PeriodPair) (T : Finset ℂ) (h0T : (0:ℂ) ∈ T)
    (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice) :
    veluKernelNorm (latticeCurve P)
        ((T.erase 0).image (fun w => (latticePointX P w, latticePointY P w)))
      = ∏ w ∈ T.erase 0, P.derivWeierstrassP w := by
  have hnot : ∀ w ∈ T.erase 0, w ∉ P.lattice :=
    fun w hw => rep_notMem_lattice P P' T h0T hT hrep hw
  have hinj : ∀ w ∈ T.erase 0, ∀ v ∈ T.erase 0,
      (latticePointX P w, latticePointY P w)
        = (latticePointX P v, latticePointY P v) → w = v := by
    intro w hw v hv he
    refine rep_sub_mem_lattice_imp_eq P P' T hT hrep (Finset.mem_of_mem_erase hw)
      (Finset.mem_of_mem_erase hv) ?_
    exact sub_mem_lattice_of_point_eq P (hnot w hw) (hnot v hv)
      (congrArg Prod.fst he) (congrArg Prod.snd he)
  rw [veluKernelNorm, Finset.prod_image hinj]
  refine Finset.prod_congr rfl fun w _ => ?_
  simp only [latticePointY, latticePointX, latticeCurve]
  ring

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] 判別式の恒等式（格子曲線、原文の語彙）**——★**`l` が奇であることだけ**（第 1402）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

    Δ(latticeCurve Λ)^l
      = Δ(veluQuotientFull (latticeCurve Λ) ⟨Q⟩) · (∏_{P∈⟨Q⟩∖{O}}(2y_P + a₁x_P + a₃))^4

★★★これで `Skeleton/GenEll/VeluDiscIdentity.lean` の
`disc_pow_eq_veluQuot_mul_lattice` の `sorry` が埋まる。 -/
theorem disc_pow_eq_lattice (P : PeriodPair) {l m : ℕ} (hlm : l = 2 * m + 1)
    (Q : (latticeCurve P).toAffine.Point) (hQ : addOrderOf Q = l) :
    (latticeCurve P).Δ ^ l
      = (veluQuotientFull (latticeCurve P)
          (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).Δ
        * (veluKernelNorm (latticeCurve P)
          (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) ^ 4 := by
  classical
  have hΔ : latticeDisc P ≠ 0 := latticeDisc_ne_zero P
  have hl : 0 < l := by omega
  obtain ⟨z₀, P', A, B, C, D, hz₀, h1, h2, hdet, hP'⟩ :=
    exists_isogeny_periodPair P hΔ hl hQ
  have hord := intCast_mul_mem_lattice_iff P hΔ hQ hz₀
  obtain ⟨T, h0T, hcard, hT, hrep, hTdef⟩ := exists_velu_rep P P' z₀ l hl hP' hord
  have hle : P.lattice ≤ P'.lattice := by rw [hP']; exact le_sup_left
  have hvelu : ∀ z : ℂ, P'.weierstrassP z = veluAnalyticX P T (veluAnalyticC P T) z :=
    fun z => weierstrassP_eq_velu_of_rep P P' hle T h0T hT hrep z
  have heq := latticeCurve_eq_veluQuotientFull P P' T h0T hT hrep
      (g₂_isogeny P P' T h0T hT hrep hvelu) (g₃_isogeny P P' T h0T hT hrep hvelu)
  have hnot : ∀ k ∈ (Finset.range l).erase 0, ((k : ℂ)) * z₀ ∉ P.lattice := by
    intro k hk hmem
    have hk0 : k ≠ 0 := Finset.ne_of_mem_erase hk
    have hkl : k < l := Finset.mem_range.1 (Finset.mem_of_mem_erase hk)
    have hmem' : ((k : ℤ) : ℂ) * z₀ ∈ P.lattice := by push_cast; exact hmem
    have hd := (hord (k : ℤ)).1 hmem'
    have hle' : (l : ℤ) ≤ (k : ℤ) :=
      Int.le_of_dvd (by exact_mod_cast Nat.pos_of_ne_zero hk0) hd
    omega
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
  have hset : ((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))
      = (T.erase 0).image (fun w => (latticePointX P w, latticePointY P w)) := by
    rw [hTe, Finset.image_image, ← hz₀]
    exact image_pointCoords_nsmul_eq P hΔ z₀ l hnot
  have hodd2 := hodd2_of_span P P' z₀ hlm hP' hord
  have hmain := latticeDisc_pow_eq P P' hle T h0T hT hrep hodd2 (by rw [hcard, hlm])
  rw [hset, ← heq, latticeCurve_Δ, latticeCurve_Δ, hlm, hmain,
    veluKernelNorm_lattice P P' T h0T hT hrep]

/-! ## ★出典の紐付け(`.src`) -/

def hodd2_of_span.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(l が奇なら Λ′/Λ に 2-捩れは無い。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluKernelNorm_lattice.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(格子曲線では核のノルムは ∏ ℘′(w)。★無条件)",
    sectionId := "genell-lemma-3-5" }

def disc_pow_eq_lattice.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(判別式の恒等式——格子曲線、原文の語彙。★l が奇)",
    sectionId := "genell-lemma-3-5" }

def disc_pow_eq_lattice.needs : List ProofObligation :=
  [ .citation "[ABC3]" "latticeDisc_pow_eq(第 1401、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.latticeDisc_pow_eq") 1,
    .citation "[ABC3]" "image_pointCoords_nsmul_eq(第 1332、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.image_pointCoords_nsmul_eq") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1402）**——第 1401 を消費側の語彙に直しただけである。" ++
       "★★★これで `Skeleton/GenEll/VeluDiscIdentity.lean` の " ++
       "`disc_pow_eq_veluQuot_mul_lattice` の `sorry` が埋まり、" ++
       "**判別式の恒等式が数体の上まで閉じる**。") 17 ]

end ABC3.Found.GenEll
