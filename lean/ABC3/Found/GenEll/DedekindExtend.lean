/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ProperExtend

/-!
# [GenEll] Remark 1.5.1 —— **Dedekind 環への延長（局所と一意性）**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.8):
> Then the morphism Spec(OF ) → X determined by x [where we recall that X is proper!] allows one to pull-back the divisor D to Spec(OF ) so as to obtain an eﬀective divisor Dx on Spec(OF ).

## ★★★★★★★★段 2 をさらに 3 つに割った

`ProperExtend.lean` で「付値環なら延びる」（段 1）を閉じた。
残っていた「Dedekind 環 `𝓞_F` へ貼り合わせる」（段 2）は、こう割れる:

| | 内容 | 状態 |
|---|---|---|
| 2a | **各素点 `v` で `(𝓞_F)_v`-点へ延びる** | ★**閉じた**（本ファイル） |
| 2b | **大域延長は一意** | ★**閉じた**（本ファイル） |
| 2c | 局所延長を貼って大域延長を作る | ★★残り |

★2a は `(𝓞_F)_v` が離散付値環（＝付値環）で分数体が `F` であることから、
`ProperExtend.lean` の `exists_unique_extend_of_isProper` に流すだけ。
★★2b は **生成点が稠密**であることと `X′` の分離性から
（mathlib の `ext_of_isDominant_of_isSeparated`）。

## ★★★★2c が素朴に済まない理由（再掲）

`Spec (𝓞_F)_v` は Zariski 開ではないので、局所延長を並べても貼り合わせにならない。
★mathlib の `spread_out_of_isGermInjective`（`SpreadingOut.lean`）で
各 `v` の**開近傍**へ広げてから貼るのが正攻法である
——`Spec 𝓞_F` は整なので `IsGermInjective` は自動で付く。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory IsDedekindDomain

/-! ## ★★★★★★生成点は稠密 -/

/-- ★★★★★★**整域の分数体は `Spec` で稠密な像を持つ**。

★機構は 3 行——像は `⟨⊥⟩`（生成点）を含み、その閉包は
`zeroLocus ⊥ = univ`（`PrimeSpectrum.closure_singleton`）。 -/
theorem isDominant_specMap_fractionRing (R : Type) [CommRing R] [IsDomain R]
    (K : Type) [Field K] [Algebra R K] [IsFractionRing R K] :
    IsDominant (Spec.map (CommRingCat.ofHom (algebraMap R K))) := by
  have hmem : (⟨⊥, Ideal.isPrime_bot⟩ : PrimeSpectrum R) ∈ Set.range (Spec.map (CommRingCat.ofHom
      (algebraMap R K))).base := by
    refine ⟨⟨⊥, Ideal.isPrime_bot⟩, ?_⟩
    apply PrimeSpectrum.ext
    simp only [Spec.map_base]
    show Ideal.comap (algebraMap R K) ⊥ = ⊥
    rw [← RingHom.ker_eq_comap_bot]
    exact (RingHom.injective_iff_ker_eq_bot _).mp (IsFractionRing.injective R K)
  rw [isDominant_iff, denseRange_iff_closure_range]
  apply Set.eq_univ_of_univ_subset
  calc (Set.univ : Set (PrimeSpectrum R))
      = closure {(⟨⊥, Ideal.isPrime_bot⟩ : PrimeSpectrum R)} := by
        rw [PrimeSpectrum.closure_singleton]
        exact (PrimeSpectrum.zeroLocus_bot).symm
    _ ⊆ closure (Set.range (Spec.map (CommRingCat.ofHom (algebraMap R K))).base) :=
        closure_mono (Set.singleton_subset_iff.2 hmem)

/-! ## ★★★★★★★★段 2b —— 大域延長は一意 -/

variable (F : Type) [Field F] [NumberField F]

/-- ★★★★★★★★**大域延長は一意である**。

原文 (GenEll p.8):
> Then the morphism Spec(OF ) → X determined by x [where we recall that X is proper!] allows one to pull-back the divisor D to Spec(OF ) so as to obtain an eﬀective divisor Dx on Spec(OF ).

★★これで `Remark151Sigma.lean` の `remark_1_5_1_bdeq` が受けている点の対応 `ePt` は、
**存在すれば一意**であることが分かる——「与えられたデータ」ではなく
「一意に決まるもの」になった。

★機構は `ext_of_isDominant_of_isSeparated`——`Spec 𝓞_F` は整（したがって被約）、
`X′` は `Spec ℤ` 上分離的、生成点は稠密。 -/
theorem extend_unique {X' : Scheme.{0}} (f' : X' ⟶ Spec (CommRingCat.of ℤ)) [IsSeparated f']
    (xR xR' : specRingOfIntegers F ⟶ X')
    (h : Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) F)) ≫ xR
       = Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) F)) ≫ xR') :
    xR = xR' := by
  haveI := isDominant_specMap_fractionRing (NumberField.RingOfIntegers F) F
  exact ext_of_isDominant_of_isSeparated f' (specZIsTerminal.hom_ext _ _) _ h

/-! ## ★★★★★★★★段 2a —— 各素点で延びる -/

/-- ★`(𝓞_F)_v` から `F` への代数構造（`v` の外の元は分母に許される）。 -/
noncomputable instance algLocAtPrime (v : FinitePlace F) :
    Algebra (Localization.AtPrime v.asIdeal) F :=
  IsLocalization.localizationAlgebraOfSubmonoidLe _ _ v.asIdeal.primeCompl
    (nonZeroDivisors _) v.asIdeal.primeCompl_le_nonZeroDivisors

instance towerLocAtPrime (v : FinitePlace F) :
    IsScalarTower (NumberField.RingOfIntegers F) (Localization.AtPrime v.asIdeal) F :=
  IsLocalization.localization_isScalarTower_of_submonoid_le _ _ v.asIdeal.primeCompl
    (nonZeroDivisors _) v.asIdeal.primeCompl_le_nonZeroDivisors

/-- ★★`F` は `(𝓞_F)_v` の分数体でもある。 -/
instance fracLocAtPrime (v : FinitePlace F) :
    IsFractionRing (Localization.AtPrime v.asIdeal) F :=
  IsFractionRing.isFractionRing_of_isDomain_of_isLocalization v.asIdeal.primeCompl _ _

/-- ★★★★★★★★**素点ごとの延長** —— `F`-点は `(𝓞_F)_v`-点へ**一意に**延びる。

原文 (GenEll p.8):
> Then the morphism Spec(OF ) → X determined by x [where we recall that X is proper!] allows one to pull-back the divisor D to Spec(OF ) so as to obtain an eﬀective divisor Dx on Spec(OF ).

★`(𝓞_F)_v` は Dedekind 環の素点での局所化なので**離散付値環**、
したがって付値環であり、分数体は `F` である。
★★あとは `ProperExtend.lean` の付値環版に流すだけ。 -/
theorem exists_unique_extend_atPrime {X' : Scheme.{0}}
    (f' : X' ⟶ Spec (CommRingCat.of ℤ)) [IsProper f']
    (v : FinitePlace F) (xK : Spec (CommRingCat.of F) ⟶ X') :
    ∃! xv : Spec (CommRingCat.of (Localization.AtPrime v.asIdeal)) ⟶ X',
      Spec.map (CommRingCat.ofHom
        (algebraMap (Localization.AtPrime v.asIdeal) F)) ≫ xv = xK :=
  exists_unique_extend_of_isProper f' _ F xK

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` は置かない。** 残っているのは**段 2c**——
各 `v` の局所延長を開近傍へ広げ（`spread_out_of_isGermInjective`）、貼り合わせる段。 -/

def extend_unique.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (iii)(固有性から延ばした点は一意)",
    sectionId := "genell-def-1-5" }

def exists_unique_extend_atPrime.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (iii)(固有性から点を延ばす——各素点で)",
    sectionId := "genell-def-1-5" }

def exists_unique_extend_atPrime.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_unique_extend_of_isProper(付値環の場合)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_unique_extend_of_isProper") 8,
    .citation "[mathlib]" "Localization.AtPrime の離散付値環性(Dedekind 環)"
      (.inMathlib "IsLocalization.AtPrime.isDiscreteValuationRing_of_dedekind_domain") 8,
    .citation "[mathlib]" "ext_of_isDominant_of_isSeparated(分離的なら生成点で決まる)"
      (.inMathlib "AlgebraicGeometry.ext_of_isDominant_of_isSeparated") 8,
    .implicitStep
      ("★★大域延長の一意性(extend_unique)が付いたので、点の対応 ePt は" ++
       "**存在すれば一意**である——「与えられたデータ」ではなく「一意に決まるもの」になった") 8,
    .implicitStep
      ("★★★残る段 2c: 各 v の局所延長を開近傍へ広げて貼り合わせる。" ++
       "★mathlib の spread_out_of_isGermInjective(SpreadingOut.lean)が使える" ++
       "——Spec 𝓞_F は整なので IsGermInjective は自動で付く。" ++
       "★★貼り合わせの整合は extend_unique と同じ機構(生成点の稠密性)で出る") 8 ]

end ABC3.Found.GenEll
