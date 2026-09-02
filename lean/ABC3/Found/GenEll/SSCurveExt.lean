/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.EllModuliObjects
import ABC3.Found.GaloisRep.JExpCoprime
import ABC3.Found.GenEll.Lemma37Full
import ABC3.Meta.Claim

/-!
# 第 1343 ブロック —— **`SSCurve` を `ℂ` の中の有限拡大へ底変換する**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か——`l`-巡回の生成元は `L` に無い

`HasLCyclicJ` が与えるのは `E[l]` の中の **`Gal`-安定な直線**であり、
その生成元 `Q` が有理になるのは次数が `l−1` を割る拡大 `L(H)` の上である。

★したがって「商の類」を作るには **`E` を `L(H)` へ持ち上げてから** Vélu を回す必要がある。
☆`SSCurve` の定義体は `IntermediateField ℚ ℂ` なので、
`ℂ` の中の有限拡大 `M₀ : IntermediateField E.fld ℂ` を取って `restrictScalars ℚ` する。

★★`j`・`ht^Falt`・`deg∞` は類の函数なので、この持ち上げで**変わらない**
（`degInfJ_eq`・`faltingsHeightJ_eq` は任意の `SSCurve` に当たる）。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve IsDedekindDomain NumberField
open ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

section Ext

variable (E : SSCurve) (M₀ : IntermediateField E.fld ℂ) [FiniteDimensional E.fld M₀]

/-- ★★★★★★**`ℂ` の中の有限拡大を `ℚ` の上の中間体として見る**（第 1343）。 -/
noncomputable def extField : IntermediateField ℚ ℂ := M₀.restrictScalars ℚ

instance finiteDimensional_extField : FiniteDimensional ℚ (extField E M₀) := by
  have h1 : FiniteDimensional ℚ E.fld := inferInstance
  have h2 : FiniteDimensional E.fld M₀ := inferInstance
  exact Module.Finite.trans (R := ℚ) (A := E.fld) (M := M₀)

instance numberField_extField : NumberField (extField E M₀) where
  to_charZero := inferInstance
  to_finiteDimensional := finiteDimensional_extField E M₀

noncomputable instance algebra_extField : Algebra E.fld (extField E M₀) :=
  inferInstanceAs (Algebra E.fld M₀)

instance isScalarTower_extField : IsScalarTower ℚ E.fld (extField E M₀) :=
  IsScalarTower.of_algebraMap_eq fun _ => rfl

/-- ★★★★★★★★★★★★★★★★
**`SSCurve` を `ℂ` の中の有限拡大へ持ち上げる**（第 1343）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆半安定性は `semistableAt_baseChange_all`（第 1220、無条件）が与える。 -/
noncomputable def SSCurve.ext : SSCurve where
  K := extField E M₀
  isNF := numberField_extField E M₀
  W := E.W.baseChange (extField E M₀)
  isEll := by
    show (E.W.map (algebraMap E.fld (extField E M₀))).IsElliptic
    infer_instance
  ss := by
    haveI : (E.W.baseChange (extField E M₀)).IsElliptic := by
      show (E.W.map (algebraMap E.fld (extField E M₀))).IsElliptic
      infer_instance
    exact semistableAt_baseChange_all E.fld (extField E M₀) E.W E.ss

/-- ★★★★★★★★**持ち上げても `j` は変わらない**——★**無条件**（第 1343）。 -/
theorem SSCurve.ext_j : (E.ext M₀).j = E.j := by
  show (extField E M₀).val ((E.W.baseChange (extField E M₀)).j) = E.K.val (E.W.j)
  have h : (E.W.baseChange (extField E M₀)).j
      = algebraMap E.fld (extField E M₀) (E.W.j) := by
    show (E.W.map (algebraMap E.fld (extField E M₀))).j = _
    rw [WeierstrassCurve.map_j]
  rw [h]
  rfl

/-- ★★★★★★★★★★★★★★★★
**`[M₀ : L] < l` なら `PrimeToLocalHeights` は持ち上げても保たれる**——★**無条件**（第 1343）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`jExp` の言葉に直して `not_dvd_jExp_baseChange_all`（第 1220）を当て、戻す。 -/
theorem SSCurve.primeToLocalHeights_ext {l : ℕ} (hl : Nat.Prime l)
    (hdeg : Module.finrank E.fld (extField E M₀) < l)
    (h : E.PrimeToLocalHeights l) :
    (E.ext M₀).PrimeToLocalHeights l := by
  have hj : ∀ p : HeightOneSpectrum (𝓞 E.fld), jExp p E.W < 0 → ¬ ((l : ℤ) ∣ jExp p E.W) :=
    E.not_dvd_jExp_of_primeToLocalHeights hl h
  haveI : (E.W.baseChange (extField E M₀)).IsElliptic := (E.ext M₀).isEll
  have hjb := not_dvd_jExp_baseChange_all E.fld (extField E M₀) E.W hl hdeg hj
  intro P hne
  have hlh : (E.ext M₀).localHeightAt P = minDeltaExp P ((E.ext M₀).W) := rfl
  have hss : SemistableAt P ((E.ext M₀).W) := (E.ext M₀).ss P
  have hmd : minDeltaExp P ((E.ext M₀).W) = max 0 (- jExp P ((E.ext M₀).W)) :=
    minDeltaExp_eq_maxJ_of_semistable P _ hss
  have hnn : 0 ≤ minDeltaExp P ((E.ext M₀).W) := minDeltaExp_nonneg P _
  rw [hlh] at hne
  have hneg : jExp P ((E.ext M₀).W) < 0 := by
    by_cases hc : jExp P ((E.ext M₀).W) < 0
    · exact hc
    · exact absurd (by rw [hmd]; omega) hne
  have hmd2 : minDeltaExp P ((E.ext M₀).W) = - jExp P ((E.ext M₀).W) := by
    rw [hmd]; omega
  rw [hlh]
  refine (Nat.Prime.coprime_iff_not_dvd hl).2 ?_
  intro hdvd
  have hdvdZ : (l : ℤ) ∣ minDeltaExp P ((E.ext M₀).W) := by
    have hz := Int.natCast_dvd_natCast.2 hdvd
    rwa [Int.toNat_of_nonneg hnn] at hz
  rw [hmd2] at hdvdZ
  exact hjb P hneg (dvd_neg.mp hdvdZ)

/-- ★★★★★★**上にある素点は存在する**（整拡大での持ち上げ、第 1347）。 -/
theorem exists_heightOneSpectrum_over (p : HeightOneSpectrum (𝓞 E.fld)) :
    ∃ P : HeightOneSpectrum (𝓞 (extField E M₀)), P.asIdeal.LiesOver p.asIdeal := by
  obtain ⟨P, -, hPp, hPo⟩ :=
    Ideal.exists_ideal_over_prime_of_isIntegral (R := 𝓞 E.fld) (S := 𝓞 (extField E M₀))
      p.asIdeal ⊥ (by simp)
  haveI := hPp
  have hPne : P ≠ ⊥ := by
    intro h
    rw [h] at hPo
    exact p.ne_bot (by simpa using hPo.symm)
  exact ⟨⟨P, hPp, hPne⟩, ⟨hPo.symm⟩⟩

/-- ★★★★★★★★★★★★**乗法還元は持ち上げても保たれる**——★**無条件**（第 1347）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`jExp P (E ⊗ M₀) = e·jExp p E`（`jExp_baseChange`、在庫）で `e > 0` だから、
`jExp p E < 0` なら `jExp P < 0`、したがって `minDeltaExp P ≠ 0` である。 -/
theorem SSCurve.hasMultRed_ext (h : E.HasMultRed) : (E.ext M₀).HasMultRed := by
  obtain ⟨p, hp⟩ := h
  obtain ⟨P, hPo⟩ := exists_heightOneSpectrum_over E M₀ p
  haveI := hPo
  haveI hell : (E.W.baseChange ((extField E M₀ : IntermediateField ℚ ℂ) : Type)).IsElliptic := by
    show (E.W.map (algebraMap E.fld (extField E M₀))).IsElliptic
    infer_instance
  have hmd : minDeltaExp p E.W = max 0 (- jExp p E.W) :=
    minDeltaExp_eq_maxJ_of_semistable p E.W (E.ss p)
  have hneg : jExp p E.W < 0 := by
    by_cases hc : jExp p E.W < 0
    · exact hc
    · exact absurd (by rw [show E.localHeightAt p = minDeltaExp p E.W from rfl, hmd]; omega) hp
  have hjb : jExp P (E.W.baseChange ((extField E M₀ : IntermediateField ℚ ℂ) : Type))
      = (p.asIdeal.ramificationIdx P.asIdeal : ℤ) * jExp p E.W :=
    ABC3.Found.GaloisRep.jExp_baseChange E.fld (extField E M₀) p P E.W
  have hene : p.asIdeal.ramificationIdx P.asIdeal ≠ 0 :=
    Ideal.IsDedekindDomain.ramificationIdx_ne_zero_of_liesOver P.asIdeal p.ne_bot
  have hepos : (0 : ℤ) < (p.asIdeal.ramificationIdx P.asIdeal : ℤ) := by
    exact_mod_cast Nat.pos_of_ne_zero hene
  have hnegP : jExp P (E.W.baseChange ((extField E M₀ : IntermediateField ℚ ℂ) : Type)) < 0 := by
    rw [hjb]
    nlinarith [hepos, hneg]
  refine ⟨P, ?_⟩
  have hmdP : minDeltaExp P (E.W.baseChange ((extField E M₀ : IntermediateField ℚ ℂ) : Type))
      = max 0 (- jExp P (E.W.baseChange ((extField E M₀ : IntermediateField ℚ ℂ) : Type))) :=
    minDeltaExp_eq_maxJ_of_semistable P _ ((E.ext M₀).ss P)
  show minDeltaExp P (E.W.baseChange ((extField E M₀ : IntermediateField ℚ ℂ) : Type)) ≠ 0
  rw [hmdP]
  omega

end Ext

/-! ## ★出典の紐付け(`.src`) -/

def SSCurve.ext.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(SSCurve を ℂ の中の有限拡大へ持ち上げる)",
    sectionId := "genell-lemma-3-5" }

def SSCurve.ext_j.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(持ち上げても j は変わらない。★無条件)",
    sectionId := "genell-lemma-3-5" }

def SSCurve.primeToLocalHeights_ext.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5([M₀ : L] < l なら PrimeToLocalHeights は持ち上げても保たれる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def SSCurve.primeToLocalHeights_ext.needs : List ProofObligation :=
  [ .citation "[ABC3]" "not_dvd_jExp_baseChange_all(第 1220、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.not_dvd_jExp_baseChange_all") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1343）**——`HasLCyclicJ` の生成元は " ++
       "`L(H)`（次数は `l−1` を割る）でないと有理にならないので、" ++
       "商の類を作るには `E` をそこへ持ち上げる必要がある。") 2 ]

def exists_heightOneSpectrum_over.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(上にある素点は存在する——整拡大での持ち上げ。★無条件)",
    sectionId := "genell-lemma-3-5" }

def SSCurve.hasMultRed_ext.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(乗法還元は持ち上げても保たれる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def extField.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(ℂ の中の有限拡大を ℚ の上の中間体として見る)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
