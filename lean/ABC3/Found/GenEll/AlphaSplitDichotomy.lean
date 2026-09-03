/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.SplitLocalAssembly
import ABC3.Found.GaloisRep.RamifiedQuadExt
import ABC3.Meta.Claim

/-!
# 第 1383 ブロック —— **α 側の分裂性の二者択一**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★★★★★これは何か——α 側の最後の一手

第 1381（`exists_h2_h1_of_split_local`）は**分裂**乗法還元を仮説に置いていた。
★本ブロックはそれを外す:

* 分裂する場合——そのまま第 1381
* 分裂しない場合——2 次式が既約なので `R[X]/(f)` は**非分岐 2 次拡大**であり、
  その剰余体で 2 次式は分裂する（第 1025・第 1033）

☆非分岐なので `hpe` の `e` は変わらない（第 1382）。

★★★これで **`exists_h2_h1_unipotent` の局所側がすべて閉じる**。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine IsDedekindDomain NumberField Polynomial
open ABC3.Interface.GaloisRep ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**乗法還元から `h2`・`h1` が出る**——★**分裂性を仮定しない**（第 1383）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これが α 側の最後の一手である。 -/
theorem exists_h2_h1_of_multRed_local
    (L : Type) [Field L] [NumberField L] {Lv : Type} [Field Lv] [CharZero Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [CharZero R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {e : ℕ} (he : 1 ≤ e) (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (W : WeierstrassCurve L) [W.IsElliptic] (C : WeierstrassCurve.VariableChange L)
    [(C • W).IsElliptic]
    (hC : WeierstrassCurve.IsMinimal (primeSubring p) (C • W))
    (hc4ne : (C • W).c₄ ≠ 0) (hc4 : valAdd p (Units.mk0 ((C • W).c₄) hc4ne) = 0)
    [((C • W).baseChange Lv).IsElliptic]
    [WeierstrassCurve.IsMinimal R ((C • W).baseChange Lv)]
    [WeierstrassCurve.HasMultiplicativeReduction R ((C • W).baseChange Lv)]
    (hj : jExp p W < 0) {l : ℕ} [Fact l.Prime] (hle : ¬ (l ∣ e))
    (hcop : ¬ ((l : ℤ) ∣ jExp p W))
    {ζv : Lvˣ} (hζv : IsPrimitiveRoot ((ζv : Lv)) l) :
    ∃ σbar : AlgebraicClosure L ≃ₐ[L] AlgebraicClosure L,
      (∀ x : tateModule (W.baseChange (AlgebraicClosure L)) l,
          ∃ u : tateModule (W.baseChange (AlgebraicClosure L)) l,
          galTate W l σbar (galTate W l σbar x) + x
            = galTate W l σbar x + galTate W l σbar x + l • u) ∧
        (∃ x : tateModule (W.baseChange (AlgebraicClosure L)) l,
          ∀ u : tateModule (W.baseChange (AlgebraicClosure L)) l,
          galTate W l σbar x ≠ x + l • u) := by
  have hA := integralModel_c4_isUnit (R := R) ((C • W).baseChange Lv)
  by_cases hs : WeierstrassCurve.HasSplitMultiplicativeReduction R ((C • W).baseChange Lv)
  · exact exists_h2_h1_of_split_local L p hpe W C hs hj hle hcop hζv
  · -- ★非分裂——非分岐 2 次拡大へ上げる
    set F := splitQuadPoly (WeierstrassCurve.integralModel R ((C • W).baseChange Lv)) hA
      with hFdef
    have hfm : F.Monic := monic_splitQuadPoly _ hA
    have hfd : F.natDegree = 2 := natDegree_splitQuadPoly _ hA
    have hirr := irreducible_map_residue_of_not_splits hfm hfd
      (fun hsp => hs (hasSplit_of_splits_splitQuadPoly _ hA hsp))
    haveI hdom := isDomain_adjoinRoot hfm hirr
    haveI hlocR := isLocalRing_adjoinRoot hfm hfd hirr
    haveI hdvr := isDiscreteValuationRing_adjoinRoot hfm hfd hirr
    haveI hac := isAdicComplete_maximalIdeal_adjoinRoot hfm hfd hirr
    haveI hcz : CharZero (AdjoinRoot F) :=
      charZero_of_injective_algebraMap (algebraMap_adjoinRoot_injective (R := R) hfm hfd)
    haveI hcz2 : CharZero (FractionRing (AdjoinRoot F)) := by
      refine charZero_of_injective_algebraMap (R := AdjoinRoot F) ?_
      exact IsFractionRing.injective (AdjoinRoot F) (FractionRing (AdjoinRoot F))
    letI : Algebra L (FractionRing (AdjoinRoot F)) :=
      ((quadFieldHom (K := Lv) hfm hfd).comp (algebraMap L Lv)).toAlgebra
    have halg : ∀ x : L, algebraMap L (FractionRing (AdjoinRoot F)) x
        = quadFieldHom hfm hfd (algebraMap L Lv x) := fun _ => rfl
    haveI hmin' := isMinimal_baseChange_ext_ram p hpe hfm hfd hirr halg W C hC hc4ne hc4
    haveI hmult' :=
      hasMultiplicativeReduction_ext_ram he p hpe hfm hfd hirr halg W C hC hc4ne hc4 hj
    have h' := hasSplitMultiplicativeReduction_ext (C • W) hA hFdef halg
    have hζ2 : IsPrimitiveRoot
        ((Units.map (quadFieldHom (K := Lv) hfm hfd).toMonoidHom ζv :
          (FractionRing (AdjoinRoot F))ˣ) : FractionRing (AdjoinRoot F)) l :=
      hζv.map_of_injective (f := (quadFieldHom (K := Lv) hfm hfd))
        (quadFieldHom_injective hfm hfd)
    exact exists_h2_h1_of_split_local L
      (Lv := FractionRing (AdjoinRoot F)) (R := AdjoinRoot F)
      p (valuation_algebraMap_ext_ram p hpe hfm hfd hirr halg) W C h' hj hle hcop hζ2

/-! ## ★出典の紐付け(`.src`) -/

def exists_h2_h1_of_multRed_local.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(乗法還元から h2・h1 が出る。分裂性を仮定しない。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_h2_h1_of_multRed_local.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_h2_h1_of_split_local(第 1381、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_h2_h1_of_split_local") 1,
    .citation "[ABC3]" "hasSplitMultiplicativeReduction_ext(第 1033、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.hasSplitMultiplicativeReduction_ext") 1,
    .citation "[ABC3]" "isAdicComplete_maximalIdeal_adjoinRoot(第 1382、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isAdicComplete_maximalIdeal_adjoinRoot") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1383）**——α 側の最後の一手である。" ++
       "☆非分裂の場合は `R[X]/(f)` へ上げる。非分岐なので `hpe` の `e` は変わらない。") 19 ]

end ABC3.Found.GenEll
