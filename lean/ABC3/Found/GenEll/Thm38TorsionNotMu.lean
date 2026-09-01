/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Thm38PiExists
import ABC3.Meta.Claim

/-!
# 第 1278 ブロック —— **`l`-捩れが `l` 個より多ければ `μ_l` からはみ出す**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——第 1277 の仮説を作る段

第 1277（`exists_pow_eq_of_torsion_not_mu`）は
「`l`-捩れの類が `μ_l` の像に収まらない」ことを仮説に受けていた。

☆本ブロックはそれを**個数の勘定だけ**で出す:

> `ζ^l = 1` なら `{[ζ^a] : a ∈ ℤ}` はたかだか `l` 個。
> したがって `l`-捩れの類が `l` 個より多ければ、はみ出す類がある。

★★★消費側では `E[l]` の個数が `l²`（`torsion_card`、在庫）なので、
`l ≥ 2` ならこの条件は満たされる。
-/

namespace ABC3.Found.GenEll

open ABC3.Meta

/-- ★★★★★★★★**`ζ^l = 1` なら `ζ` の整数冪は `0 ≤ i < l` の冪に落ちる**（第 1278）。 -/
theorem zpow_eq_pow_emod {G : Type*} [CommGroup G] {l : ℕ} (hl : 0 < l) (ζ : G)
    (hζl : ζ ^ l = 1) (a : ℤ) :
    ζ ^ a = ζ ^ ((a % (l : ℤ)).toNat) := by
  have hmod : 0 ≤ a % (l : ℤ) := Int.emod_nonneg a (by exact_mod_cast hl.ne')
  have hsplit : a = (l : ℤ) * (a / (l : ℤ)) + a % (l : ℤ) := (Int.ediv_add_emod a (l : ℤ)).symm
  calc ζ ^ a = ζ ^ ((l : ℤ) * (a / (l : ℤ)) + a % (l : ℤ)) := by rw [← hsplit]
    _ = (ζ ^ (l : ℤ)) ^ (a / (l : ℤ)) * ζ ^ (a % (l : ℤ)) := by rw [zpow_add, zpow_mul]
    _ = ζ ^ (a % (l : ℤ)) := by
        rw [zpow_natCast, hζl, one_zpow, one_mul]
    _ = ζ ^ ((a % (l : ℤ)).toNat) := by
        rw [← zpow_natCast ζ ((a % (l : ℤ)).toNat), Int.toNat_of_nonneg hmod]

/-- ★★★★★★★★★★★★★★★★
**`l`-捩れが `l` 個より多ければ `μ_l` の像からはみ出す類がある**——★**無条件**（第 1278）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`{[ζ^a] : a ∈ ℤ}` は `(range l).image (fun i => [ζ^i])` に含まれるので
たかだか `l` 個である。★個数の勘定だけで出る。 -/
theorem exists_torsion_not_mu {G : Type*} [CommGroup G] {l : ℕ} (hl : 0 < l) (Q ζ : G)
    (hζl : ζ ^ l = 1)
    (T : Finset (G ⧸ Subgroup.zpowers Q)) (hT : ∀ c ∈ T, c ^ l = 1) (hcard : l < T.card) :
    ∃ c : G ⧸ Subgroup.zpowers Q, c ^ l = 1 ∧ ∀ a : ℤ, c ≠ QuotientGroup.mk (ζ ^ a) := by
  classical
  by_contra hcon
  push_neg at hcon
  have hsub : T ⊆ (Finset.range l).image
      (fun i : ℕ => (QuotientGroup.mk (ζ ^ i) : G ⧸ Subgroup.zpowers Q)) := by
    intro c hc
    obtain ⟨a, ha⟩ := hcon c (hT c hc)
    refine Finset.mem_image.2 ⟨(a % (l : ℤ)).toNat, ?_, ?_⟩
    · refine Finset.mem_range.2 ?_
      have hlt : a % (l : ℤ) < (l : ℤ) := Int.emod_lt_of_pos a (by exact_mod_cast hl)
      omega
    · rw [← zpow_eq_pow_emod hl ζ hζl a, ← ha]
  have h1 : T.card ≤ ((Finset.range l).image
      (fun i : ℕ => (QuotientGroup.mk (ζ ^ i) : G ⧸ Subgroup.zpowers Q))).card :=
    Finset.card_le_card hsub
  have h2 : ((Finset.range l).image
      (fun i : ℕ => (QuotientGroup.mk (ζ ^ i) : G ⧸ Subgroup.zpowers Q))).card ≤ l := by
    refine le_trans (Finset.card_image_le) ?_
    simp
  omega

/-! ## ★出典の紐付け(`.src`) -/

def zpow_eq_pow_emod.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(ζ^l = 1 なら整数冪は 0 ≤ i < l の冪に落ちる。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_torsion_not_mu.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l-捩れが l 個より多ければ μ_l の像からはみ出す。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_torsion_not_mu.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1278）**——第 1277 の仮説を**個数の勘定だけ**で作る段である。" ++
       "☆消費側では `E[l]` の個数が `l²`（`torsion_card`、在庫）なので、" ++
       "`l ≥ 2` ならこの条件は満たされる。" ++
       "★これで `π` の存在は「`E[l]` が `K` に載っている」ことだけから出る。") 2 ]

end ABC3.Found.GenEll
