/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluPointSet
import ABC3.Meta.Claim

/-!
# 第 1213 ブロック —— **`⟨Q⟩∖{O}` は `Q ↦ c·Q` で保たれる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——安定直線が集合を保つことの核

第 1205 は `HasLCyclicJ` から
「`∀ σ, ∃ k, σ Q = k • Q`」を与える。☆この `k` は `l` で割れない
（`σ Q ≠ 0` だから）ので、`σ` は `⟨Q⟩∖{O}` を**それ自身に写す**。

★本ブロックはその**群論の核**を取る——曲線も Galois も出てこない:

    l 素数、`addOrderOf Q = l`、`l ∤ c`  ⟹
      `k • (c • Q) ∈ {j • Q : j ∈ (range l).erase 0}`   (`k ∈ (range l).erase 0`)

☆`k • (c • Q) = (k·c) • Q = ((k·c) mod l) • Q` であり、
`l` が素で `l ∤ k`・`l ∤ c` だから `(k·c) mod l ≠ 0` である。

★★★これが `fixesCoeffs_veluQuotientFull`（第 1153）が受け取る
「`∀ z ∈ S, f z ∈ S`」の形の中身である。
-/

namespace ABC3.Found.GenEll

open ABC3.Meta

/-- ★★★★★★★★★★★★
**`⟨Q⟩∖{O}` は `Q ↦ c·Q`（`l ∤ c`）で保たれる**——★**無条件**（第 1213）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`k • (c • Q) = (k·c) • Q = ((k·c) mod l) • Q` であり、
`l` が素で `l ∤ k`・`l ∤ c` だから `(k·c) mod l ≠ 0`。

★★★これが「安定直線は座標集合を保つ」の**群論の核**である。 -/
theorem nsmul_nsmul_mem_image_of_not_dvd {G : Type*} [AddCommGroup G] [DecidableEq G]
    {Q : G} {l : ℕ} (hl : l.Prime) (hQ : addOrderOf Q = l)
    {c : ℕ} (hc : ¬ (l ∣ c)) {k : ℕ} (hk : k ∈ (Finset.range l).erase 0) :
    (k • (c • Q)) ∈ ((Finset.range l).erase 0).image (fun j : ℕ => j • Q) := by
  have hlQ : l • Q = 0 := by
    rw [← hQ]; exact addOrderOf_nsmul_eq_zero Q
  rw [Finset.mem_erase, Finset.mem_range] at hk
  have hkc : ¬ (l ∣ k * c) := by
    intro h
    rcases (Nat.Prime.dvd_mul hl).1 h with h1 | h2
    · exact absurd (Nat.le_of_dvd (Nat.pos_of_ne_zero hk.1) h1) (by omega)
    · exact hc h2
  have hmlt : (k * c) % l < l := Nat.mod_lt _ hl.pos
  have hmne : (k * c) % l ≠ 0 := fun h0 => hkc (Nat.dvd_of_mod_eq_zero h0)
  refine Finset.mem_image.2
    ⟨(k * c) % l, Finset.mem_erase.2 ⟨hmne, Finset.mem_range.2 hmlt⟩, ?_⟩
  have hsm : k • (c • Q) = (k * c) • Q := (mul_smul k c Q).symm
  rw [hsm]
  have hdiv : k * c = ((k * c) / l) * l + (k * c) % l := (Nat.div_add_mod' (k * c) l).symm
  conv_rhs => rw [hdiv]
  rw [add_smul, mul_smul, hlQ, smul_zero, zero_add]

/-! ## ★出典の紐付け(`.src`) -/

def nsmul_nsmul_mem_image_of_not_dvd.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(⟨Q⟩∖{O} は Q ↦ c·Q（l ∤ c）で保たれる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def nsmul_nsmul_mem_image_of_not_dvd.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1213）**——第 1205 が与える" ++
       "「`∀ σ, ∃ k, σ Q = k • Q`」から「`σ` は `⟨Q⟩∖{O}` を保つ」を出す" ++
       "**群論の核**である。☆`fixesCoeffs_veluQuotientFull`（第 1153）が受け取る" ++
       "「`∀ z ∈ S, f z ∈ S`」の形の中身であり、曲線も Galois も出てこない。") 2 ]

end ABC3.Found.GenEll
