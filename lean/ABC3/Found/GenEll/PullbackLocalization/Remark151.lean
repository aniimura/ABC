/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.LocalizationBridge
import ABC3.Found.GenEll.ComapAffine
import ABC3.Found.GenEll.PullbackNatural
import ABC3.Found.GenEll.PullbackLocalization.Definition12

/-!
# PullbackLocalization —— `[GenEll] Remark 1.5.1` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open CategoryTheory AlgebraicGeometry

/-! ## ★★★★★★★★★★`LocalizationBridge` が要求する形 -/

variable (F : Type) [Field F] [NumberField F]

/-- ★★★★★★★★**`𝓞_F[1/N]` へ拡大したイデアルは、`Spec 𝓞_F[1/N]` の点での引き戻しである**。

★★これが `LocalizationBridge.lean` の `conductorADiv_fin_eq_of_map_eq` が要求する
仮定 `h` を作るための配線である。 -/
theorem pullbackIdeal_map_algebraMap (A : Type) [CommRing A]
    [Algebra (NumberField.RingOfIntegers F) A]
    {X : Scheme.{0}} (D : X.IdealSheafData) (xF : specRingOfIntegers F ⟶ X) :
    (pullbackIdeal F D xF).map (algebraMap (NumberField.RingOfIntegers F) A)
      = pullbackIdealOf (CommRingCat.of A) D
          (Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) A)) ≫ xF) := by
  rw [← pullbackIdealOf_eq_pullbackIdeal F D xF]
  have h := (pullbackIdealOf_specMap D xF
    (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) A))).symm
  simpa using h

/-- ★★★★★★★★★★**[GenEll] Remark 1.5.1 —— `hagree` を作る形**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

`Spec 𝓞_F[1/N]` の水準で 2 つの因子の引き戻しが一致するなら、
**`N` を割らない素点 `v` で導手の係数が一致する**。

★★★これが `Skeleton/GenEll/Section1.lean` の `remark_1_5_1` が受けている
仮定 `hagree` そのものである。 -/
theorem conductorADiv_fin_eq_of_localized
    (M : Submonoid (NumberField.RingOfIntegers F)) (A : Type) [CommRing A]
    [Algebra (NumberField.RingOfIntegers F) A] [IsLocalization M A]
    (v : FinitePlace F) (hM : ∀ m ∈ M, m ∉ v.asIdeal)
    {X X' : Scheme.{0}} (D : X.IdealSheafData) (D' : X'.IdealSheafData)
    (xF : specRingOfIntegers F ⟶ X) (xF' : specRingOfIntegers F ⟶ X')
    (hI : pullbackIdeal F D xF ≠ 0) (hJ : pullbackIdeal F D' xF' ≠ 0)
    (h : pullbackIdealOf (CommRingCat.of A) D
          (Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) A)) ≫ xF)
        = pullbackIdealOf (CommRingCat.of A) D'
          (Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) A)) ≫ xF')) :
    (conductorADiv F D xF).fin v = (conductorADiv F D' xF').fin v := by
  refine conductorADiv_fin_eq_of_map_eq F M A v hM D D' xF xF' hI hJ ?_
  rw [pullbackIdeal_map_algebraMap, pullbackIdeal_map_algebraMap, h]

def conductorADiv_fin_eq_of_localized.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(hagree を作る形——X_n の点への配線は含まない)",
    sectionId := "genell-rem-1-5-1" }

def conductorADiv_fin_eq_of_localized.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "conductorADiv_fin_eq_of_map_eq(局所化の橋)"
      (.inProject "ABC3" "ABC3.Found.GenEll.conductorADiv_fin_eq_of_map_eq") 9,
    .citation "[ABC3]" "ideal_comap_eq_map_of_isAffine(アフィンでの comap は map)"
      (.inProject "ABC3" "ABC3.Found.GenEll.ideal_comap_eq_map_of_isAffine") 5,
    .citation "[mathlib]" "Scheme.ΓSpecIso_naturality(Spec を取って Γ を取ると戻る)"
      (.inMathlib "AlgebraicGeometry.Scheme.ΓSpecIso_naturality") 5,
    .implicitStep
      ("★既存の pullbackIdeal_specMap は K が数体であることを要求するので、" ++
       "𝓞_F[1/N](数体の整数環ではない)には当たらない。" ++
       "★★そこで任意の可換環の上で定義し直した(pullbackIdealOf)。" ++
       "pullbackIdeal はその特別な場合であり rfl で一致する") 5,
    .implicitStep
      ("★★★残る配線: ConductorDescent.lean の X_n 上の因子の一致から本定理の仮定 h を作る。" ++
       "★機構は Spec 𝓞_F[1/N] ⟶ X_n を作ること" ++
       "(𝓞_F[1/N] は ℤ[1/n!]-代数なので Spec 𝓞_F[1/N] ⟶ Spec ℤ[1/n!] があり、" ++
       "X への射と合わせて引き戻しの普遍性で作れる)") 9 ]

end ABC3.Found.GenEll
