/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateKerLevel
import ABC3.Found.GaloisRep.TateLevelOne
import ABC3.Meta.Claim

/-!
# 第 1270 ブロック —— **`T_l E` 上の条件は `E[l]` 上の条件に落ちる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★これは何か——II 側の最後の配管の第 1 段

`alpha_mem_map_of_galTate`（第 1237）が要求するのは 2 つだけである:

| # | 仮説 | 中身 |
|---|---|---|
| `h2` | `∀ x, ∃ u, σ²x + x = 2σx + l·u` | `σ` は `mod l` で**幂単** |
| `h1` | `∃ x, ∀ u, σx ≠ x + l·u` | `σ` は `mod l` で**非自明** |

★どちらも `T_l E` の言葉だが、**`mod l` の条件なので `E[l]` で確かめれば足りる**
——`T_l E / l·T_l E ≅ E[l]` が第 1203（`tateProj_one_eq_zero_iff`）である。

☆本ブロックはその降下を取る。★★★これで II 側の残りは
「Tate 一意化の `σ` が `E[l]` に `α` として作用する」ことだけになる
——それは第 1174（`tate_sigma_coord_alpha`）と第 1212（`kummer_sigma_coord_alpha`）が持っている。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Interface.GaloisRep ABC3.Meta

variable {K L : Type} [Field K] [DecidableEq K] [Field L] [DecidableEq L] [Algebra K L]

/-- ★★★★★★★★★★★★★★★★
**`E[l]` で幂単なら `T_l E` でも `mod l` 幂単**——★**無条件**（第 1270）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`tateProj W l 1` の核はちょうど `l · T_l E`（第 1203）なので、
差が第 1 層で消えることを見れば `l` 倍の形に書ける。

★★★これが `alpha_mem_map_of_galTate`（第 1237）の `h2` である。 -/
theorem galTate_unipotent_of_galPoint (W : WeierstrassCurve K) (l : ℕ) (σ : L ≃ₐ[K] L)
    (h : ∀ P : (W.baseChange L).toAffine.Point, (l ^ 1) • P = 0 →
      galPoint W σ (galPoint W σ P) + P = galPoint W σ P + galPoint W σ P) :
    ∀ x : tateModule (W.baseChange L) l, ∃ u : tateModule (W.baseChange L) l,
      galTate W l σ (galTate W l σ x) + x
        = galTate W l σ x + galTate W l σ x + l • u := by
  intro x
  have hzero : ((tateProj (W.baseChange L) l 1
      (galTate W l σ (galTate W l σ x) + x - (galTate W l σ x + galTate W l σ x)) :
      (W.baseChange L).toAffine.Point)) = 0 := by
    have hP : (l ^ 1) • ((tateProj (W.baseChange L) l 1 x : (W.baseChange L).toAffine.Point))
        = 0 := (tateProj (W.baseChange L) l 1 x).2
    have hkey := h ((tateProj (W.baseChange L) l 1 x : (W.baseChange L).toAffine.Point)) hP
    simp only [map_sub, map_add, AddSubgroup.coe_sub, AddSubgroup.coe_add, tateProj_galTate]
    rw [hkey]
    abel
  obtain ⟨g, hg⟩ := (tateProj_one_eq_zero_iff (W.baseChange L) l _).1 hzero
  refine ⟨g, ?_⟩
  rw [hg]
  abel

/-- ★★★★★★★★★★★★★★★★
**`E[l]` で非自明なら `T_l E` でも `mod l` 非自明**——★**無条件**（第 1270）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`l`-捩れ点は第 1 層に持ち上がる（第 1202、`exists_tateProj_one_eq`）ので、
その持ち上げが証人になる。★`l·T_l E` の第 1 層への像は `0` である。

★★★これが `alpha_mem_map_of_galTate`（第 1237）の `h1` である。 -/
theorem exists_galTate_ne_of_galPoint [IsAlgClosed L] [CharZero L]
    (W : WeierstrassCurve K) [((W.baseChange L).toAffine).IsElliptic] (l : ℕ) [Fact l.Prime]
    (σ : L ≃ₐ[K] L) (P : (W.baseChange L).toAffine.Point) (hP : l • P = 0)
    (hne : galPoint W σ P ≠ P) :
    ∃ x : tateModule (W.baseChange L) l, ∀ u : tateModule (W.baseChange L) l,
      galTate W l σ x ≠ x + l • u := by
  obtain ⟨f, hf⟩ := exists_tateProj_one_eq W l P hP
  refine ⟨f, ?_⟩
  intro u hu
  refine hne ?_
  have hcoe := congrArg
    (fun z : tateModule (W.baseChange L) l =>
      ((tateProj (W.baseChange L) l 1 z : (W.baseChange L).toAffine.Point))) hu
  have hz : (l : ℕ) • ((tateProj (W.baseChange L) l 1 u :
      (W.baseChange L).toAffine.Point)) = 0 := by
    have h1 : (l ^ 1) • ((tateProj (W.baseChange L) l 1 u :
        (W.baseChange L).toAffine.Point)) = 0 := (tateProj (W.baseChange L) l 1 u).2
    simpa using h1
  have hzc : ((l • (tateProj (W.baseChange L) l 1 u) :
      torsionPoints (W.baseChange L) (l ^ 1)) :
      (W.baseChange L).toAffine.Point) = 0 := by
    simpa using hz
  simp only [map_add, map_nsmul, AddSubgroup.coe_add, tateProj_galTate, hf, hzc,
    add_zero] at hcoe
  exact hcoe

/-! ## ★出典の紐付け(`.src`) -/

def galTate_unipotent_of_galPoint.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(E[l] で幂単なら T_l E でも mod l 幂単。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_galTate_ne_of_galPoint.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(E[l] で非自明なら T_l E でも mod l 非自明。★無条件)",
    sectionId := "genell-thm-3-8" }

def galTate_unipotent_of_galPoint.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tateProj_one_eq_zero_iff(第 1203、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateProj_one_eq_zero_iff") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1270）**——`alpha_mem_map_of_galTate`（第 1237）の" ++
       "`h2`・`h1` を `E[l]` の条件に落とした。☆`T_l E / l·T_l E ≅ E[l]`（第 1203）と" ++
       "`l`-捩れ点の持ち上げ（第 1202）が要るところである。" ++
       "★★★これで II 側の残りは「Tate 一意化の `σ` が `E[l]` に `α` として作用する」だけになった。") 2 ]

end ABC3.Found.GaloisRep
