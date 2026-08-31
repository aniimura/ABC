/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ProjectiveModel
import ABC3.Meta.Claim
import Mathlib.RingTheory.Localization.Algebra

/-!
# ★★★★★★★★段 F2c —— 閉包の生成ファイバーは元に戻る（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.8。

原文 (GenEll p.8):
> Remark 1.4.1. Observe that it follows immediately from the definitions, together with Proposition 1.4, (iii), that the theory of

## ★★★★★★★★これは何か —— 段 F の最後

段 F1（閉包の切断が `ℤ`-平坦）と段 F2ab（`ℤ`-射影・`ℤ`-固有）は閉じている。
★残っていたのは「**閉包の生成ファイバーが元の `Y` に戻る**」ことである。

## ★★★★★★機構 —— 局所化は核と可換である

アフィンで読むと、閉包のイデアルは `ker φ`（`φ : A → B`、`A = Γ(U ⊆ ℙᴺ_ℤ)`、
`B = Γ(Y, f⁻¹U)`）である。★生成ファイバーへ移るのは `A` を `ℤ∖{0}` で局所化すること
（`A_ℚ`）であり、`B` は**すでに `ℚ`-代数**なので局所化で変わらない。

★★したがって要るのは

    `ker (A_ℚ → B) = (ker φ) · A_ℚ`   （`IsLocalization.ker_map`、mathlib）

だけである。★★★そして `A_ℚ → B` が全射（`Y ↪ ℙᴺ_ℚ` が閉埋め込み）なら

    `A_ℚ / (ker φ)·A_ℚ ≅ B`

——すなわち**閉包の生成ファイバーは `Y` そのもの**である。

## ★測定の記録

★`B` が局所化で変わらないことは `IsLocalization T B`（`T` が単元のみからなる）で表した。
mathlib に「単元の部分モノイドによる局所化は恒等」の instance が見当たらなかったので、
★★`isLocalization_self_of_isUnit` を**自前で置いた**（3 行）。

★★★`ℤ∖{0}` の像が `B` で単元になるのは「`B` が `ℚ`-代数」だからである
（`isUnit_intCast_of_ratAlgebra`）——`n ≠ 0` なら `(n : B) = algebraMap ℚ B n` が単元。

## ★残っている段（明示）

★★本ファイルは**アフィンでの核**である。`Scheme.IdealSheafData` の水準へ持ち上げるには
「閉包のイデアル層をアフィン開ごとに読む」（`Scheme.Hom.ker_apply`、段 F1 で使った在庫）
を各チャートに当てて貼るだけであり、**新しい数学は要らない**。
-/

namespace ABC3.Found.GenEll

/-! ## ★単元だけからなるモノイドによる局所化は恒等である -/

/-- ★**単元だけからなる部分モノイドによる局所化は恒等である**。

★mathlib に instance が見当たらなかったので自前で置いた（`IsLocalization.atUnits` は
「局所化が与えられているとき同型になる」形で、instance を**作る**形ではない）。 -/
theorem isLocalization_self_of_isUnit {B : Type} [CommRing B] (T : Submonoid B)
    (hT : ∀ t ∈ T, IsUnit t) : IsLocalization T B where
  map_units := fun y => by simpa using hT y.1 y.2
  surj := fun z => ⟨⟨z, 1⟩, by simp⟩
  exists_of_eq := fun {x y} h => ⟨1, by simpa using h⟩

variable {A B : Type} (Aq : Type) [CommRing A] [CommRing B] [CommRing Aq]
  (M : Submonoid A) [Algebra A Aq] [IsLocalization M Aq] (φ : A →+* B)

/-! ## ★★★★★★局所化は核と可換である -/

/-- ★★**局所化した環準同型** `A_M → B`（`M` の像が `B` で単元なら定まる）。 -/
noncomputable def locMap (hM : ∀ m ∈ M, IsUnit (φ m)) : Aq →+* B :=
  letI : IsLocalization (M.map φ) B :=
    isLocalization_self_of_isUnit _ (by rintro _ ⟨m, hm, rfl⟩; exact hM m hm)
  IsLocalization.map (M := M) B φ (Submonoid.le_comap_map M)

/-- ★★★★★★**局所化は核と可換である**。

    `ker (A_M → B) = (ker φ) · A_M`

★これが段 F2c の心臓である——閉包のイデアルを生成ファイバーへ移しても、
**もとのイデアルの生成する分にちょうど等しい**。
★★機構は mathlib の `IsLocalization.ker_map` そのものである。 -/
theorem ker_locMap (hM : ∀ m ∈ M, IsUnit (φ m)) :
    RingHom.ker (locMap Aq M φ hM) = Ideal.map (algebraMap A Aq) (RingHom.ker φ) :=
  letI : IsLocalization (M.map φ) B :=
    isLocalization_self_of_isUnit _ (by rintro _ ⟨m, hm, rfl⟩; exact hM m hm)
  IsLocalization.ker_map (R := A) (S := Aq) (P := B) (Q := B) (M := M) (T := M.map φ) φ rfl

/-- ★★★★★★★★**閉包の生成ファイバーは元に戻る**（アフィンでの核）。

原文 (GenEll p.8):
> Remark 1.4.1. Observe that it follows immediately from the definitions, together with Proposition 1.4, (iii), that the theory of

`A_M → B` が**全射**（＝ `Y` が生成ファイバーの側で閉部分スキーム）なら

    `A_M / (ker φ)·A_M ≅ B`

★すなわち**閉包を生成ファイバーへ落とすと `Y` そのものに戻る**。 -/
noncomputable def quotKerMapEquiv (hM : ∀ m ∈ M, IsUnit (φ m))
    (hsurj : Function.Surjective (locMap Aq M φ hM)) :
    (Aq ⧸ Ideal.map (algebraMap A Aq) (RingHom.ker φ)) ≃+* B :=
  (Ideal.quotEquivOfEq (ker_locMap Aq M φ hM).symm).trans
    (RingHom.quotientKerEquivOfSurjective hsurj)

/-! ## ★★★★★`ℤ → ℚ` の場合 -/

/-- ★**`ℚ`-代数では `0` でない整数は単元である**。 -/
theorem isUnit_intCast_of_ratAlgebra {B : Type} [CommRing B] [Algebra ℚ B] {n : ℤ} (hn : n ≠ 0) :
    IsUnit ((n : ℤ) : B) := by
  have h : ((n : ℤ) : B) = algebraMap ℚ B (n : ℚ) := by simp
  rw [h]
  exact (IsUnit.mk0 (n : ℚ) (Int.cast_ne_zero.2 hn)).map (algebraMap ℚ B)

/-- ★★**したがって `ℤ∖{0}` の像は `ℚ`-代数で単元になる**——`hM` が自動で埋まる。 -/
theorem isUnit_of_mem_intSubmonoid {A B : Type} [CommRing A] [CommRing B] [Algebra ℚ B]
    (φ : A →+* B) {m : A}
    (hm : m ∈ Submonoid.map (Int.castRingHom A) (nonZeroDivisors ℤ)) :
    IsUnit (φ m) := by
  obtain ⟨n, hn, rfl⟩ := hm
  have hn0 : n ≠ 0 := fun h => by simp [h] at hn
  simpa using isUnit_intCast_of_ratAlgebra (B := B) hn0

/-- ★★★★★★★★★**段 F2c の核**（`ℤ → ℚ` の場合）。

原文 (GenEll p.8):
> Remark 1.4.1. Observe that it follows immediately from the definitions, together with Proposition 1.4, (iii), that the theory of

`A` を `ℤ∖{0}` で局所化したもの（＝生成ファイバーの座標環）を `A_ℚ` と書くと、
`A_ℚ → B` が全射なら ★**`A_ℚ / (ker φ)·A_ℚ ≅ B`**。

★★これが「閉包の生成ファイバーが元の `Y` に戻る」ことのアフィンでの中身である。 -/
noncomputable def ratQuotKerEquiv {A B : Type} (Aq : Type) [CommRing A] [CommRing B]
    [Algebra ℚ B] [CommRing Aq] [Algebra A Aq]
    [IsLocalization (Submonoid.map (Int.castRingHom A) (nonZeroDivisors ℤ)) Aq]
    (φ : A →+* B)
    (hsurj : Function.Surjective
      (locMap Aq (Submonoid.map (Int.castRingHom A) (nonZeroDivisors ℤ)) φ
        (fun _ hm => isUnit_of_mem_intSubmonoid φ hm))) :
    (Aq ⧸ Ideal.map (algebraMap A Aq) (RingHom.ker φ)) ≃+* B :=
  quotKerMapEquiv Aq _ φ _ hsurj

/-! ## ★出典の紐付け(`.src`) -/

def isLocalization_self_of_isUnit.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Remark 1.4.1(単元だけからなるモノイドによる局所化は恒等)",
    sectionId := "genell-rem-1-4-1" }

def ker_locMap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Remark 1.4.1(局所化は核と可換である)",
    sectionId := "genell-rem-1-4-1" }

def quotKerMapEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Remark 1.4.1(閉包の生成ファイバーは元に戻る——アフィンでの核)",
    sectionId := "genell-rem-1-4-1" }

def ratQuotKerEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Remark 1.4.1(段 F2c の核——ℤ → ℚ の場合)",
    sectionId := "genell-rem-1-4-1" }

def quotKerMapEquiv.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "IsLocalization.ker_map(局所化は核と可換)"
      (.inMathlib "IsLocalization.ker_map") 8,
    .citation "[mathlib]" "RingHom.quotientKerEquivOfSurjective"
      (.inMathlib "RingHom.quotientKerEquivOfSurjective") 8,
    .citation "[ABC3]" "flat_int_subschemeObj(段 F1——閉包の切断が ℤ-平坦)"
      (.inProject "ABC3" "ABC3.Found.GenEll.flat_int_subschemeObj") 8,
    .implicitStep
      ("★mathlib に「単元の部分モノイドによる局所化は恒等」の instance が見当たらなかったので " ++
       "isLocalization_self_of_isUnit を自前で置いた(3 行)。" ++
       "★IsLocalization.atUnits は『局所化が与えられているとき同型になる』形で、" ++
       "instance を**作る**形ではない(2026-08-28 実測)") 8,
    .implicitStep
      ("★★本ファイルは**アフィンでの核**である。Scheme.IdealSheafData の水準へ持ち上げるには " ++
       "「閉包のイデアル層をアフィン開ごとに読む」(Scheme.Hom.ker_apply、段 F1 で使った在庫)を" ++
       "各チャートに当てて貼るだけであり、**新しい数学は要らない**") 8 ]

end ABC3.Found.GenEll
