/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.AlphaBridgeFull2
import ABC3.Found.GenEll.LocalInputsBadPrime
import ABC3.Found.GenEll.LocalFixedMoved
import ABC3.Meta.Claim

/-!
# 第 1320 ブロック —— **悪い素点から大域の `h2`・`h1` を出す**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★★★★★★★★★これは何か——II 側の**一本道**

第 1318（悪い素点 → 局所の 2 入力）＋第 1309（局所の固定点・動く点）＋
第 1319（大域へ）を繋ぐ。

☆受け取るのは原文の言葉だけである:

| # | 入力 |
|---|---|
| 1 | `jExp p E < 0`（悪い素点）と `¬ l ∣ jExp p E`（`PrimeToLocalHeights`） |
| 2 | `E ⊗ L_v` は極小かつ**分裂**乗法還元 |
| 3 | `ζ ∈ L̄` が原始 `l` 乗根で、その像が `L_v` に入る |
| 4 | `L̄ ↪ M`（`IsAlgClosed.lift` で取る） |

★★★これで `alpha_mem_map_of_galTate`（第 1237）に渡せる `σ` が出る。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine IsDedekindDomain NumberField
open ABC3.Interface.GaloisRep ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**悪い素点から大域の `h2`・`h1` を出す**——★**無条件**（第 1320）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆第 1318 ＋ 第 1309 ＋ 第 1319 を繋いだだけである。 -/
theorem exists_h2_h1_of_bad_prime
    (L : Type) [Field L] [NumberField L] {Lv M : Type} [Field Lv] [CharZero Lv] [Field M]
    [Algebra L Lv] [Algebra Lv M] [Algebra L M] [IsScalarTower L Lv M]
    [IsAlgClosed M] [IsGalois Lv M]
    {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = (HeightOneSpectrum.valuation L p) x)
    (E : WeierstrassCurve L) [E.IsElliptic]
    [(E.baseChange Lv).IsElliptic] [WeierstrassCurve.IsMinimal R (E.baseChange Lv)]
    [WeierstrassCurve.IsElliptic ((E.baseChange Lv).baseChange M).toAffine]
    (h : WeierstrassCurve.HasSplitMultiplicativeReduction R (E.baseChange Lv))
    (hj : jExp p E < 0) {l : ℕ} [Fact l.Prime] (hcop : ¬ ((l : ℤ) ∣ jExp p E))
    {ζv : Lvˣ} (hζv : IsPrimitiveRoot ((ζv : Lv)) l)
    (e : tateModule (E.baseChange (AlgebraicClosure L)) l ≃+ (Fin 2 → ℤ_[l]))
    (ι : AlgebraicClosure L →ₐ[L] M)
    (ζ : AlgebraicClosure L) (hζ : IsPrimitiveRoot ζ l)
    (z : Lv) (hz : algebraMap Lv M z = ι ζ) :
    ∃ σbar : AlgebraicClosure L ≃ₐ[L] AlgebraicClosure L,
      (∀ x : tateModule (E.baseChange (AlgebraicClosure L)) l,
          ∃ u : tateModule (E.baseChange (AlgebraicClosure L)) l,
          galTate E l σbar (galTate E l σbar x) + x
            = galTate E l σbar x + galTate E l σbar x + l • u) ∧
        (∃ x : tateModule (E.baseChange (AlgebraicClosure L)) l,
          ∀ u : tateModule (E.baseChange (AlgebraicClosure L)) l,
          galTate E l σbar x ≠ x + l • u) := by
  have hl : l.Prime := Fact.out
  -- 局所の 2 入力
  obtain ⟨hord, hcard⟩ :=
    local_inputs_of_bad_prime (R := R) p hp E h hj hl hcop hζv
  obtain ⟨P₀, hP₀⟩ := hord
  -- 局所の固定点・動く点
  obtain ⟨σM, Pfix, Qmov, hord', hfix', hQ', hmov'⟩ :=
    exists_local_fixed_moved (M := M) (E.baseChange Lv) l hl P₀ hP₀ hcard
  -- `L̄ ↪ M` と制限
  letI : Algebra (AlgebraicClosure L) M := ι.toAlgebra
  haveI : IsScalarTower L (AlgebraicClosure L) M :=
    IsScalarTower.of_algebraMap_eq (fun x => (ι.commutes x).symm)
  set σbar := restrictLocalHom L Lv M (AlgebraicClosure L) σM with hσbar
  have hcomm : ∀ x : AlgebraicClosure L, ι (σbar x) = (σM.restrictScalars L) (ι x) := by
    intro x
    exact restrictLocalHom_commutes L Lv M (AlgebraicClosure L) σM x
  refine ⟨σbar, ?_⟩
  exact exists_h2_h1_of_localData' L Lv M E l e ι ζ hζ z hz P₀ hP₀ hcard σbar σM hcomm
    Pfix Qmov hord' hfix' hQ' hmov'

/-! ## ★出典の紐付け(`.src`) -/

def exists_h2_h1_of_bad_prime.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(悪い素点から大域の h2・h1 を出す。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_h2_h1_of_bad_prime.needs : List ProofObligation :=
  [ .citation "[ABC3]" "local_inputs_of_bad_prime(第 1318、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.local_inputs_of_bad_prime") 1,
    .citation "[ABC3]" "exists_local_fixed_moved(第 1309、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_local_fixed_moved") 1,
    .citation "[ABC3]" "exists_h2_h1_of_localData′(第 1319、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_h2_h1_of_localData'") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1320）**——II 側の**一本道**である。" ++
       "☆受け取るのは原文の言葉（悪い素点・`PrimeToLocalHeights`・分裂乗法還元・`ζ ∈ L_v`）だけ。" ++
       "★★★これで `alpha_mem_map_of_galTate`（第 1237）に渡せる `σ` が出る。") 3 ]

end ABC3.Found.GenEll
