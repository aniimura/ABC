/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.AlphaFromBadPrime
import ABC3.Found.GaloisRep.RamifiedBadPrime
import ABC3.Meta.Claim

/-!
# 第 1372 ブロック —— **悪い素点から `h2`・`h1` を出す（分岐版）**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★★★★★これは何か——第 1318／第 1320 の `e` 倍版

第 1370-1371 で「`hp` の分岐版（`hpe`）」の道を通した。
★本ブロックはそれを第 1318（局所の 2 入力）と第 1320（大域の `h2`・`h1`）に流す。

☆受け取るものは第 1320 と同じで、`hp` が `hpe` に、
そこに **`l ∤ e`** が 1 つ足されるだけである:

| # | 入力 |
|---|---|
| 1 | `jExp p E < 0`（悪い素点）と `¬ l ∣ jExp p E` |
| 2 | `hpe`（付値は `e` 倍で両立）と **`l ∤ e`** |
| 3 | `E ⊗ L_v` は極小かつ**分裂**乗法還元 |
| 4 | `ζ ∈ L̄` が原始 `l` 乗根で、その像が `L_v` に入る |
| 5 | `L̄ ↪ M` |

★★★これで `L_v` として **`L_p` の任意の有限拡大**（`ζ_l` を含むもの、
分裂させるもの）が使えるようになった——第 1366-1369 でその完備 DVR 構造を、
第 1369 で `e ≤ [L_v′ : L_v]` を与えてある。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine IsDedekindDomain NumberField
open ABC3.Interface.GaloisRep ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

/-- ★★★★★★★★★★★★★★★★
**悪い素点の言葉で局所の 2 入力を出す**——分岐版（第 1372）。 -/
theorem local_inputs_of_bad_prime_ram {L : Type} [Field L] [NumberField L]
    {Lv : Type} [Field Lv] [CharZero Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {e : ℕ} (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (E : WeierstrassCurve L) [E.IsElliptic]
    [(E.baseChange Lv).IsElliptic] [WeierstrassCurve.IsMinimal R (E.baseChange Lv)]
    (h : WeierstrassCurve.HasSplitMultiplicativeReduction R (E.baseChange Lv))
    (hj : jExp p E < 0) {l : ℕ} (hl : l.Prime) (hle : ¬ (l ∣ e))
    (hcop : ¬ ((l : ℤ) ∣ jExp p E))
    {ζ : Lvˣ} (hζ : IsPrimitiveRoot ((ζ : Lv)) l) :
    (∃ P₀ : (E.baseChange Lv).toAffine.Point, addOrderOf P₀ = l) ∧
      (∀ T : Finset ((E.baseChange Lv).toAffine.Point),
        (∀ q ∈ T, l • q = 0) → T.card ≤ l) := by
  have hnd := not_dvd_vAdd_tateParam_of_not_dvd_jExp_ram p hpe E h hj hl hle hcop
  exact local_inputs_of_split (R := R) (E.baseChange Lv) h hl hζ hnd

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**悪い素点から大域の `h2`・`h1` を出す**——分岐版（第 1372）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★第 1320 の `hp` を `hpe`（`e` 倍）に緩めた形である。
☆`l ∤ e` さえあれば分岐した拡大でも道が通る。 -/
theorem exists_h2_h1_of_bad_prime_ram
    (L : Type) [Field L] [NumberField L] {Lv M : Type} [Field Lv] [CharZero Lv] [Field M]
    [Algebra L Lv] [Algebra Lv M] [Algebra L M] [IsScalarTower L Lv M]
    [IsAlgClosed M] [IsGalois Lv M]
    {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {e : ℕ} (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (E : WeierstrassCurve L) [E.IsElliptic]
    [(E.baseChange Lv).IsElliptic] [WeierstrassCurve.IsMinimal R (E.baseChange Lv)]
    [WeierstrassCurve.IsElliptic ((E.baseChange Lv).baseChange M).toAffine]
    (h : WeierstrassCurve.HasSplitMultiplicativeReduction R (E.baseChange Lv))
    (hj : jExp p E < 0) {l : ℕ} [Fact l.Prime] (hle : ¬ (l ∣ e))
    (hcop : ¬ ((l : ℤ) ∣ jExp p E))
    {ζv : Lvˣ} (hζv : IsPrimitiveRoot ((ζv : Lv)) l)
    (eq : tateModule (E.baseChange (AlgebraicClosure L)) l ≃+ (Fin 2 → ℤ_[l]))
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
    local_inputs_of_bad_prime_ram (R := R) p hpe E h hj hl hle hcop hζv
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
  exact exists_h2_h1_of_localData' L Lv M E l eq ι ζ hζ z hz P₀ hP₀ hcard σbar σM hcomm
    Pfix Qmov hord' hfix' hQ' hmov'

/-! ## ★出典の紐付け(`.src`) -/

def local_inputs_of_bad_prime_ram.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(悪い素点の言葉で局所の 2 入力を出す。分岐版。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_h2_h1_of_bad_prime_ram.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(悪い素点から大域の h2・h1 を出す。分岐版。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_h2_h1_of_bad_prime_ram.needs : List ProofObligation :=
  [ .citation "[ABC3]" "not_dvd_vAdd_tateParam_of_not_dvd_jExp_ram(第 1371、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.not_dvd_vAdd_tateParam_of_not_dvd_jExp_ram") 1,
    .citation "[ABC3]" "exists_h2_h1_of_bad_prime(第 1320、証明済み。e = 1 の場合)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_h2_h1_of_bad_prime") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1372）**——第 1320 の `hp` を `hpe`（`e` 倍）に緩めた形。" ++
       "★★★これで `L_v` として **`L_p` の任意の有限拡大**が使える。" ++
       "☆完備 DVR 構造は第 1366-1368、`e ≤ [L_v′ : L_v]` は第 1369。") 19 ]

end ABC3.Found.GenEll
