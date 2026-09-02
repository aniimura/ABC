/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.GalPointVc
import ABC3.Found.GenEll.UnipTransfer
import ABC3.Found.GenEll.AlphaFromBadPrimeRam
import ABC3.Meta.Claim

/-!
# 第 1380 ブロック —— **局所データを `C • E` から `E` に戻す**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★これは何か——第 1372 の `C • E` 版

第 1372（`exists_h2_h1_of_bad_prime_ram`）は
`[IsMinimal R (E.baseChange Lv)]` を要求するが、
`SemistableAt` が与えるのは **`C • E`** の極小性である。

★★★本ブロックは局所の 2 入力（`P₀`・`hcard`）を
第 1379（`galPoint` は変数変換と可換）で `E` に戻し、
**`h2`・`h1` を `E` について出す**。

☆`exists_local_fixed_moved`（第 1309）以降は `E ⊗ L_v` の側で回るので、
戻すのは `P₀` と `hcard` の 2 つだけでよい。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine IsDedekindDomain NumberField
open ABC3.Interface.GaloisRep ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

section VcTransport

variable {K L : Type} [Field K] [Field L] [Algebra K L]

/-- ★★★★★★★★**変数変換の点の加法同型（底変換つき）**（第 1380）。 -/
noncomputable def vcPointBC (W : WeierstrassCurve K) (C : VariableChange K) :
    (W.baseChange L).toAffine.Point ≃+ ((C • W).baseChange L).toAffine.Point :=
  (ABC3.Found.GaloisRep.vcPointAddEquiv (W.baseChange L) (C.map (algebraMap K L))).trans
    (pointEquivOfEq (baseChange_variableChange W C).symm)

theorem vcPointBC_apply (W : WeierstrassCurve K) (C : VariableChange K)
    (P : (W.baseChange L).toAffine.Point) :
    vcPointBC W C P
      = pointEquivOfEq (baseChange_variableChange W C).symm
          (vcPoint (C.map (algebraMap K L)) (W.baseChange L) P) := rfl

/-- ★★★★★★★★★★★★★★★★
**`galPoint` は `vcPointBC` と可換**——★**無条件**（第 1380）。 -/
theorem galPoint_vcPointBC (W : WeierstrassCurve K) (C : VariableChange K)
    (σ : L ≃ₐ[K] L) (P : (W.baseChange L).toAffine.Point) :
    galPoint (C • W) σ (vcPointBC W C P) = vcPointBC W C (galPoint W σ P) := by
  by_cases hP : P = 0
  · subst hP
    simp
  · have hvc : vcPoint (C.map (algebraMap K L)) (W.baseChange L) P ≠ 0 :=
      vcPoint_ne_zero _ _ hP
    have hvc2 : vcPoint (C.map (algebraMap K L)) (W.baseChange L) (galPoint W σ P) ≠ 0 :=
      vcPoint_ne_zero _ _ (galPoint_ne_zero W σ hP)
    have hL : vcPointBC W C P ≠ 0 := by
      rw [vcPointBC_apply]
      intro hcon
      exact hvc ((pointEquivOfEq (baseChange_variableChange W C).symm).injective
        (by simpa using hcon))
    have hR : vcPointBC W C (galPoint W σ P) ≠ 0 := by
      rw [vcPointBC_apply]
      intro hcon
      exact hvc2 ((pointEquivOfEq (baseChange_variableChange W C).symm).injective
        (by simpa using hcon))
    refine pointCoords_injective (galPoint_ne_zero (C • W) σ hL) hR ?_
    rw [pointCoords_galPoint, vcPointBC_apply, vcPointBC_apply,
      pointCoords_pointEquivOfEq, pointCoords_pointEquivOfEq,
      pointCoords_vcPoint' _ _ P hP,
      pointCoords_vcPoint' _ _ (galPoint W σ P) (galPoint_ne_zero W σ hP),
      pointCoords_galPoint]
    refine Prod.ext ?_ ?_
    · exact map_vcX_fixed (σ : L →+* L) (C.map (algebraMap K L))
        (map_variableChange_algEquiv C σ) _
    · exact map_vcY_fixed (σ : L →+* L) (C.map (algebraMap K L))
        (map_variableChange_algEquiv C σ) _ _

end VcTransport

section BadPrimeVc

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**悪い素点から `h2`・`h1` を出す**——分岐版かつ**変数変換つき**（第 1380）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★第 1372 との違いは、極小性・分裂乗法還元を **`C • W`** の側で受け、
結論は **`W`** について出すことである。
☆`SemistableAt` が与えるのは `C • W` の極小性なので、これが使う形である。 -/
theorem exists_h2_h1_of_bad_prime_ram_vc
    (L : Type) [Field L] [NumberField L] {Lv M : Type} [Field Lv] [CharZero Lv] [Field M]
    [Algebra L Lv] [Algebra Lv M] [Algebra L M] [IsScalarTower L Lv M]
    [IsAlgClosed M] [IsGalois Lv M]
    {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {e : ℕ} (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (W : WeierstrassCurve L) [W.IsElliptic] (C : WeierstrassCurve.VariableChange L)
    [(C • W).IsElliptic]
    [((C • W).baseChange Lv).IsElliptic] [WeierstrassCurve.IsMinimal R ((C • W).baseChange Lv)]
    [WeierstrassCurve.IsElliptic ((W.baseChange Lv).baseChange M).toAffine]
    (h : WeierstrassCurve.HasSplitMultiplicativeReduction R ((C • W).baseChange Lv))
    (hj : jExp p W < 0) {l : ℕ} [Fact l.Prime] (hle : ¬ (l ∣ e))
    (hcop : ¬ ((l : ℤ) ∣ jExp p W))
    {ζv : Lvˣ} (hζv : IsPrimitiveRoot ((ζv : Lv)) l)
    (eq : tateModule (W.baseChange (AlgebraicClosure L)) l ≃+ (Fin 2 → ℤ_[l]))
    (ι : AlgebraicClosure L →ₐ[L] M)
    (ζ : AlgebraicClosure L) (hζ : IsPrimitiveRoot ζ l)
    (z : Lv) (hz : algebraMap Lv M z = ι ζ) :
    ∃ σbar : AlgebraicClosure L ≃ₐ[L] AlgebraicClosure L,
      (∀ x : tateModule (W.baseChange (AlgebraicClosure L)) l,
          ∃ u : tateModule (W.baseChange (AlgebraicClosure L)) l,
          galTate W l σbar (galTate W l σbar x) + x
            = galTate W l σbar x + galTate W l σbar x + l • u) ∧
        (∃ x : tateModule (W.baseChange (AlgebraicClosure L)) l,
          ∀ u : tateModule (W.baseChange (AlgebraicClosure L)) l,
          galTate W l σbar x ≠ x + l • u) := by
  have hl : l.Prime := Fact.out
  have hjC : jExp p (C • W) < 0 := by rw [jExp_variableChange p W C]; exact hj
  have hcopC : ¬ ((l : ℤ) ∣ jExp p (C • W)) := by
    rw [jExp_variableChange p W C]; exact hcop
  -- ★局所の 2 入力（`C • W` の側）
  obtain ⟨hord, hcard⟩ :=
    local_inputs_of_bad_prime_ram (R := R) p hpe (C • W) h hjC hl hle hcopC hζv
  obtain ⟨P₀', hP₀'⟩ := hord
  -- ★変数変換で `W` に戻す
  set φ := vcPointBC (L := Lv) W C with hφ
  refine ?_
  have hP₀ : addOrderOf (φ.symm P₀') = l := by
    have h1 : addOrderOf (φ (φ.symm P₀')) = addOrderOf (φ.symm P₀') :=
      addOrderOf_injective φ.toAddMonoidHom φ.injective (φ.symm P₀')
    rw [φ.apply_symm_apply] at h1
    rw [← h1, hP₀']
  have hcardW : ∀ T : Finset ((W.baseChange Lv).toAffine.Point),
      (∀ q ∈ T, l • q = 0) → T.card ≤ l := by
    intro T hT
    have himg : ∀ q ∈ T.image (φ : _ → _), l • q = 0 := by
      intro q hq
      obtain ⟨q0, hq0, rfl⟩ := Finset.mem_image.1 hq
      rw [← map_nsmul, hT q0 hq0, map_zero]
    have hc := hcard (T.image (φ : _ → _)) himg
    rwa [Finset.card_image_of_injective _ φ.injective] at hc
  -- ★局所の固定点・動く点
  obtain ⟨σM, Pfix, Qmov, hord', hfix', hQ', hmov'⟩ :=
    exists_local_fixed_moved (M := M) (W.baseChange Lv) l hl (φ.symm P₀') hP₀ hcardW
  -- ★`L̄ ↪ M` と制限
  letI : Algebra (AlgebraicClosure L) M := ι.toAlgebra
  haveI : IsScalarTower L (AlgebraicClosure L) M :=
    IsScalarTower.of_algebraMap_eq (fun x => (ι.commutes x).symm)
  set σbar := restrictLocalHom L Lv M (AlgebraicClosure L) σM with hσbar
  have hcomm : ∀ x : AlgebraicClosure L, ι (σbar x) = (σM.restrictScalars L) (ι x) := by
    intro x
    exact restrictLocalHom_commutes L Lv M (AlgebraicClosure L) σM x
  refine ⟨σbar, ?_⟩
  exact exists_h2_h1_of_localData' L Lv M W l eq ι ζ hζ z hz (φ.symm P₀') hP₀ hcardW
    σbar σM hcomm Pfix Qmov hord' hfix' hQ' hmov'

end BadPrimeVc

/-! ## ★出典の紐付け(`.src`) -/

def galPoint_vcPointBC.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(galPoint は vcPointBC と可換。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_h2_h1_of_bad_prime_ram_vc.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(悪い素点から h2・h1 を出す。分岐版かつ変数変換つき。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_h2_h1_of_bad_prime_ram_vc.needs : List ProofObligation :=
  [ .citation "[ABC3]" "galAct_vcPoint_baseChange(第 1379、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.galAct_vcPoint_baseChange") 1,
    .citation "[ABC3]" "exists_h2_h1_of_localData′(第 1319、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_h2_h1_of_localData'") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1380）**——第 1372 との違いは、" ++
       "極小性・分裂乗法還元を `C • W` の側で受け、結論は `W` について出すことである。" ++
       "☆戻すのは `P₀` と `hcard` の 2 つだけでよい。") 19 ]

end ABC3.Found.GenEll
