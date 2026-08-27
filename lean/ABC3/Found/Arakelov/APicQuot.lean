/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.ArithPic
import ABC3.Found.Arakelov.ArcEvalGlobal
import ABC3.Found.Arakelov.ArcConjInvol
import ABC3.Meta.Claim

/-!
# `APic(X)` —— **同型類**の群（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> of tensor product, thus determine a group APic(X).

## ★★★★★★★★★対の群と同型類の群の差は `Γ(X, 𝒪_X^×)` である

`Found/Arakelov/ArithPic.lean` の

    ArithPic X = Pic(X) × Multiplicative (共役不変な連続関数)

は**対の群**である。★原文は「テンソル積で**同型類**が群 `APic(X)` をなす」と書く。

★★差は単元の作用である——`u ∈ Γ(X, 𝒪_X^×)` による自己同型は計量を `|u|` 倍するので、

    `(L, g)` と `(L, g + log‖u‖)` は**同型な算術直線束**

である。★★★本ファイルはその部分群で商を取り、原文どおりの `APic(X)` を作る。

## ★★★★★★鍵は `evalGlobal` の共役両立

`log‖u‖` が `conjArcCM`（共役不変な連続関数）に入ることを言うには

    `evalGlobal (ι_X p) g = conj (evalGlobal p g)`

が要る。★これは **`Scheme.ΓSpecIso_naturality` を `starRingEnd ℂ` に当てるだけ**で出た
——`continuous_evalGlobal` のようなチャート論法は要らなかった（2026-08-28 実測）。

★★連続性は `continuous_evalGlobal` ＋ `Real.log`（単元は非消失）。

## ★★★引き戻しも降りる

`evalGlobal (p ≫ f) g = evalGlobal p (f^# g)` なので
`arcCMPullback f (log‖u‖) = log‖f^# u‖` であり、単元の像は単元の像へ移る。
★したがって `APicOfPullback` が `QuotientGroup.lift` で作れる（`Definition 1.1, (ii)`）。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite ABC3.Interface.Arakelov ABC3.Found.GenEll

variable {X Y : Scheme.{0}}

/-! ## ★★`evalGlobal` の環としての性質 -/

theorem evalGlobal_one (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    evalGlobal p (1 : ((X.presheaf.obj (op (⊤ : X.Opens))) : Type)) = 1 := by
  show (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom ((Scheme.Hom.appTop p).hom 1) = 1
  rw [map_one, map_one]

theorem evalGlobal_mul (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (g h : ((X.presheaf.obj (op (⊤ : X.Opens))) : Type)) :
    evalGlobal p (g * h) = evalGlobal p g * evalGlobal p h := by
  show (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom ((Scheme.Hom.appTop p).hom (g * h)) = _
  rw [map_mul, map_mul]
  rfl

theorem evalGlobal_unit_ne_zero (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (u : ((X.presheaf.obj (op (⊤ : X.Opens))) : Type)ˣ) :
    evalGlobal p (u : ((X.presheaf.obj (op (⊤ : X.Opens))) : Type)) ≠ 0 := by
  intro h
  have h1 : evalGlobal p ((u : ((X.presheaf.obj (op (⊤ : X.Opens))) : Type)) * ↑u⁻¹) = 1 := by
    rw [u.mul_inv]; exact evalGlobal_one p
  rw [evalGlobal_mul, h, zero_mul] at h1
  exact zero_ne_one h1

theorem evalGlobal_comp (f : X ⟶ Y) (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (g : ((Y.presheaf.obj (op (⊤ : Y.Opens))) : Type)) :
    evalGlobal (p ≫ f) g = evalGlobal p ((Scheme.Hom.appTop f).hom g) := by
  show (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom ((Scheme.Hom.appTop (p ≫ f)).hom g) = _
  have h : Scheme.Hom.appTop (p ≫ f) = Scheme.Hom.appTop f ≫ Scheme.Hom.appTop p := by
    rw [← Scheme.Hom.comp_appTop]
  rw [h]
  rfl

/-- ★★★★**`evalGlobal` は `ι_X` と両立する**——値が複素共役になる。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★機構は `Scheme.ΓSpecIso_naturality` を `starRingEnd ℂ` に当てるだけである
——`continuous_evalGlobal` のようなチャート論法は**要らなかった**。 -/
theorem evalGlobal_conjPoint (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (g : ((X.presheaf.obj (op (⊤ : X.Opens))) : Type)) :
    evalGlobal (conjPoint p) g = starRingEnd ℂ (evalGlobal p g) := by
  have hnat := Scheme.ΓSpecIso_naturality (CommRingCat.ofHom (starRingEnd ℂ))
  show (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom
      ((Scheme.Hom.appTop (conjSpec ≫ p)).hom g) = _
  have hcomp : Scheme.Hom.appTop (conjSpec ≫ p)
      = Scheme.Hom.appTop p ≫ Scheme.Hom.appTop conjSpec := by
    rw [← Scheme.Hom.comp_appTop]
  rw [hcomp]
  show (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom
      ((Scheme.Hom.appTop conjSpec).hom ((Scheme.Hom.appTop p).hom g)) = _
  exact congrArg (fun (m : _ ⟶ _) => (CommRingCat.Hom.hom m) ((Scheme.Hom.appTop p).hom g)) hnat

/-! ## ★★★単元が定める Green 関数 -/

/-- ★★**単元 `u` が定める Green 関数 `log‖u‖`**。 -/
noncomputable def unitLogGreen (X : Scheme.{0})
    (u : ((X.presheaf.obj (op (⊤ : X.Opens))) : Type)ˣ) : arcCM X :=
  arcCM.mk X (fun p => Real.log ‖evalGlobal p (u : ((X.presheaf.obj (op (⊤ : X.Opens))) : Type))‖)
    (by
      letI := arcTopology X
      exact ((continuous_evalGlobal (u : ((X.presheaf.obj (op (⊤ : X.Opens))) : Type))).norm).log
        (fun p => norm_ne_zero_iff.2 (evalGlobal_unit_ne_zero p u)))

/-- ★★★`log‖u‖` は**共役不変**である——`evalGlobal_conjPoint` と `‖conj z‖ = ‖z‖`。 -/
theorem unitLogGreen_mem_conjArcCM (u : ((X.presheaf.obj (op (⊤ : X.Opens))) : Type)ˣ) :
    unitLogGreen X u ∈ conjArcCM X := by
  intro p
  show Real.log ‖evalGlobal (conjPoint p) _‖ = Real.log ‖evalGlobal p _‖
  rw [evalGlobal_conjPoint, RCLike.norm_conj]

theorem unitLogGreen_one : unitLogGreen X 1 = 0 := by
  letI := arcTopology X
  apply ContinuousMap.ext
  intro p
  show Real.log ‖evalGlobal p (1 : ((X.presheaf.obj (op (⊤ : X.Opens))) : Type))‖ = 0
  rw [evalGlobal_one, norm_one, Real.log_one]

theorem unitLogGreen_mul (u v : ((X.presheaf.obj (op (⊤ : X.Opens))) : Type)ˣ) :
    unitLogGreen X (u * v) = unitLogGreen X u + unitLogGreen X v := by
  letI := arcTopology X
  apply ContinuousMap.ext
  intro p
  show Real.log ‖evalGlobal p ((u * v : ((X.presheaf.obj (op (⊤ : X.Opens))) : Type)ˣ)
      : ((X.presheaf.obj (op (⊤ : X.Opens))) : Type))‖ = _
  show _ = Real.log ‖evalGlobal p _‖ + Real.log ‖evalGlobal p _‖
  rw [Units.val_mul, evalGlobal_mul, norm_mul,
    Real.log_mul (norm_ne_zero_iff.2 (evalGlobal_unit_ne_zero p u))
      (norm_ne_zero_iff.2 (evalGlobal_unit_ne_zero p v))]

theorem unitLogGreen_pullback (f : X ⟶ Y)
    (u : ((Y.presheaf.obj (op (⊤ : Y.Opens))) : Type)ˣ) :
    arcCMPullback f (unitLogGreen Y u)
      = unitLogGreen X (Units.map (Scheme.Hom.appTop f).hom u) := by
  letI := arcTopology X
  apply ContinuousMap.ext
  intro p
  show Real.log ‖evalGlobal (p ≫ f) (u : ((Y.presheaf.obj (op (⊤ : Y.Opens))) : Type))‖
    = Real.log ‖evalGlobal p ((Scheme.Hom.appTop f).hom
        (u : ((Y.presheaf.obj (op (⊤ : Y.Opens))) : Type)))‖
  rw [evalGlobal_comp]

/-! ## ★★★★★★★★★同型類の群 `APic(X)` -/

attribute [local instance] arithPicGroup

/-- ★★★**単元から算術直線束の類への写像**——`u ↦ (1, log‖u‖)`。 -/
noncomputable def unitToArithPic (X : Scheme.{0}) :
    ((X.presheaf.obj (op (⊤ : X.Opens))) : Type)ˣ →* ArithPic X :=
  letI := picardDataWitness.group X
  { toFun := fun u => (1, Multiplicative.ofAdd
      (⟨unitLogGreen X u, unitLogGreen_mem_conjArcCM u⟩ : ↥(conjArcCM X)))
    map_one' := by
      refine Prod.ext rfl ?_
      show Multiplicative.ofAdd (⟨unitLogGreen X 1, _⟩ : ↥(conjArcCM X)) = 1
      congr 1
      exact Subtype.ext unitLogGreen_one
    map_mul' := fun u v => by
      refine Prod.ext (by simp) ?_
      show Multiplicative.ofAdd (⟨unitLogGreen X (u * v), _⟩ : ↥(conjArcCM X))
        = Multiplicative.ofAdd (⟨unitLogGreen X u, _⟩ : ↥(conjArcCM X))
          * Multiplicative.ofAdd (⟨unitLogGreen X v, _⟩ : ↥(conjArcCM X))
      congr 1
      exact Subtype.ext (unitLogGreen_mul u v) }

/-- ★★★★**同型による同一視の部分群**——単元の像の全体。 -/
noncomputable def isometrySubgroup (X : Scheme.{0}) : Subgroup (ArithPic X) :=
  (unitToArithPic X).range

/-- ★★★★★★★★★**算術 Picard 群 `APic(X)`** —— 原文どおり**同型類**の群。

原文 (GenEll p.3):
> of tensor product, thus determine a group APic(X). -/
noncomputable def APicOf (X : Scheme.{0}) : Type 1 :=
  ArithPic X ⧸ isometrySubgroup X

noncomputable instance (X : Scheme.{0}) : CommGroup (APicOf X) :=
  inferInstanceAs (CommGroup (ArithPic X ⧸ isometrySubgroup X))

/-- ★類を取る写像（全射）。 -/
noncomputable def APicOf.mk (X : Scheme.{0}) : ArithPic X →* APicOf X :=
  QuotientGroup.mk' (isometrySubgroup X)

theorem APicOf.mk_surjective (X : Scheme.{0}) : Function.Surjective (APicOf.mk X) :=
  QuotientGroup.mk'_surjective _

/-! ## ★★★★★★引き戻しは `APic` に降りる -/

/-- ★引き戻しは群準同型（`pullback_mul` から `MonoidHom.mk'`）。 -/
noncomputable def picPullbackHom (f : X ⟶ Y) :
    letI := picardDataWitness.group X
    letI := picardDataWitness.group Y
    picardDataWitness.Pic Y →* picardDataWitness.Pic X :=
  letI := picardDataWitness.group X
  letI := picardDataWitness.group Y
  MonoidHom.mk' (picardDataWitness.pullback f) (picardDataWitness.pullback_mul f)

theorem picPullback_one (f : X ⟶ Y) :
    letI := picardDataWitness.group X
    letI := picardDataWitness.group Y
    picardDataWitness.pullback f (1 : picardDataWitness.Pic Y) = 1 := by
  letI := picardDataWitness.group X
  letI := picardDataWitness.group Y
  exact (picPullbackHom f).map_one

/-- ★★★★**引き戻しは同型による同一視を保つ**。 -/
theorem arithPicPullback_unitToArithPic (f : X ⟶ Y)
    (u : ((Y.presheaf.obj (op (⊤ : Y.Opens))) : Type)ˣ) :
    arithPicPullback f (unitToArithPic Y u)
      = unitToArithPic X (Units.map (Scheme.Hom.appTop f).hom u) := by
  refine Prod.ext ?_ ?_
  · exact picPullback_one f
  · exact congrArg Multiplicative.ofAdd (Subtype.ext (unitLogGreen_pullback f u))

noncomputable def arithPicPullbackHom (f : X ⟶ Y) : ArithPic Y →* ArithPic X :=
  MonoidHom.mk' (arithPicPullback f) (arithPicPullback_mul f)

/-- ★★★★★★**[GenEll] Definition 1.1, (ii)** —— 引き戻しは `APic` に降りる。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z -proper, Z -flat schemes. Then -/
noncomputable def APicOfPullback (f : X ⟶ Y) : APicOf Y →* APicOf X :=
  QuotientGroup.lift (isometrySubgroup Y)
    ((APicOf.mk X).comp (arithPicPullbackHom f)) (by
      rintro _ ⟨u, rfl⟩
      show (APicOf.mk X) (arithPicPullback f (unitToArithPic Y u)) = 1
      rw [arithPicPullback_unitToArithPic]
      exact (QuotientGroup.eq_one_iff _).2 ⟨_, rfl⟩)

/-- ★★★★**`ι_X` 両立な計量から `APic(X)` の元を作る**。 -/
noncomputable def APicOf.ofMetric (L : picardDataWitness.Pic X)
    (m : TorsorMetric X (picardDataWitness.sheafOf X L)) (hm : m.IsConjCompatible) :
    APicOf X :=
  APicOf.mk X (arithPicOfMetric L m hm)

/-- ★★★★★**すべての類は `ι_X` 両立な計量から来る**。 -/
theorem APicOf.ofMetric_surjective (x : APicOf X) :
    ∃ (L : picardDataWitness.Pic X) (m : TorsorMetric X (picardDataWitness.sheafOf X L))
      (hm : m.IsConjCompatible), APicOf.ofMetric L m hm = x := by
  obtain ⟨y, rfl⟩ := APicOf.mk_surjective X x
  obtain ⟨L, m, hm, hy⟩ := arithPicOfMetric_surjective y
  exact ⟨L, m, hm, congrArg (APicOf.mk X) hy⟩

/-! ### ★出典の紐付け(`.src`) -/

def evalGlobal_conjPoint.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(evalGlobal は ι_X と両立する)",
    sectionId := "genell-def-1-1-i" }

def unitLogGreen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(単元が定める Green 関数 log‖u‖)",
    sectionId := "genell-def-1-1-i" }

def APicOf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(APic(X)——同型類の群)",
    sectionId := "genell-def-1-1-i" }

def APicOfPullback.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(引き戻しは同型類の群に降りる)",
    sectionId := "genell-def-1-1-i" }

def APicOf.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "ArithPic(対の群)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.ArithPic") 3,
    .citation "[ABC3]" "continuous_evalGlobal(大域の正則関数は X^arc 上連続)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.continuous_evalGlobal") 3,
    .citation "[mathlib]" "Scheme.ΓSpecIso_naturality(共役両立の機構)"
      (.inMathlib "AlgebraicGeometry.Scheme.ΓSpecIso_naturality") 3,
    .implicitStep
      ("★原文は「テンソル積で**同型類**が群 APic(X) をなす」と書く。" ++
       "対の群 ArithPic X との差は Γ(X, 𝒪_X^×) の作用であり、" ++
       "本ファイルはその部分群で商を取る") 3,
    .implicitStep
      ("★★deg_F はこの作用で不変(単数の積公式)なので、" ++
       "商を取っても下流(ht)は変わらない") 3 ]

end ABC3.Found.Arakelov
