import ABC3.Interface.Arakelov.APic
import ABC3.Found.Arakelov.ArcMapCont
import ABC3.Found.Arakelov.ArcMetricWitness

/-!
# Arakelov (D1) 第 299 ブロック —— **★★★★★★★`APicData` の非空虚 witness**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★算術直線束の群 `APic(X)`

原文の `L̄ = (L, |−|_L)` そのものである。★捻れ集合の表示(第 294)により

    APic(X) = Pic(X) × C(X^arc, ℝ)

——**可逆層と Green 関数の対**になる。★★群構造は
`Pic` の積と Green 関数の**和**である(mathlib の `Multiplicative` と `Prod.instCommGroup`)。

## ★★引き戻しが Green 関数を運べる理由

`pullback f (L, c) = (f^* L, c ∘ (· ≫ f))` であり、
★★★第 298 の**弧空間の関手性**が `c ∘ (· ≫ f)` の連続性を保証する。

## ★★退化していないこと

`forgetMetric_mk` が「対を忘れると元の層に戻る」を要求するので、
★`APic := PUnit` は通らない((B1) の `equivPicRing` が `Pic` の非自明を強制している)。
-/

namespace ABC3.Interface.Arakelov

open AlgebraicGeometry CategoryTheory ABC3.Found.Arakelov ABC3.Found.GenEll

/-- ★**`X^arc` 上の連続関数のなす加法群**。 -/
@[reducible] noncomputable def arcCM (X : Scheme.{0}) : Type :=
  @ContinuousMap (Spec (CommRingCat.of ℂ) ⟶ X) ℝ (arcTopology X) inferInstance

noncomputable instance instArcCMAddCommGroup (X : Scheme.{0}) : AddCommGroup (arcCM X) :=
  letI := arcTopology X
  inferInstanceAs (AddCommGroup C(Spec (CommRingCat.of ℂ) ⟶ X, ℝ))

/-- ★連続関数から `arcCM` の元を作る。 -/
noncomputable def arcCM.mk (X : Scheme.{0}) (g : (Spec (CommRingCat.of ℂ) ⟶ X) → ℝ)
    (hg : @Continuous _ ℝ (arcTopology X) _ g) : arcCM X :=
  @ContinuousMap.mk _ ℝ (arcTopology X) _ g hg

/-- ★`arcCM` の元を関数と見る。 -/
def arcCM.fn {X : Scheme.{0}} (c : arcCM X) : (Spec (CommRingCat.of ℂ) ⟶ X) → ℝ :=
  @ContinuousMap.toFun _ ℝ (arcTopology X) _ c

theorem arcCM.fn_cont {X : Scheme.{0}} (c : arcCM X) :
    @Continuous _ ℝ (arcTopology X) _ c.fn :=
  @ContinuousMap.continuous_toFun _ ℝ (arcTopology X) _ c

@[simp] theorem arcCM.mk_fn (X : Scheme.{0}) (c : arcCM X) :
    arcCM.mk X c.fn (arcCM.fn_cont c) = c := rfl

/-- ★**射に沿った Green 関数の引き戻し**——第 298 の関手性による。 -/
noncomputable def arcCMPullback {X Y : Scheme.{0}} (f : X ⟶ Y) (c : arcCM Y) : arcCM X :=
  letI := arcTopology X
  letI := arcTopology Y
  arcCM.mk X (fun p => c.fn (p ≫ f)) ((arcCM.fn_cont c).comp (continuous_comp_scheme f))

/-- ★`APic X` の群構造——`Pic` の積と Green 関数の和。 -/
@[reducible] noncomputable def aPicGroup (X : Scheme.{0}) :
    CommGroup (picardDataWitness.Pic X × Multiplicative (arcCM X)) :=
  letI := picardDataWitness.group X
  inferInstance

/-- ★★★★★★★**`APicData` の実装**——可逆層と Green 関数の対。 -/
noncomputable def aPicDataWitness : APicData where
  toHermitianMetricData := hermitianMetricDataWitness
  APic := fun X => picardDataWitness.Pic X × Multiplicative (arcCM X)
  group := aPicGroup
  forgetMetric := fun _ x => x.1
  forgetMetric_mul := fun _ _ _ => rfl
  ofMetric := fun X L m => (L, arcCM.mk X m.green m.green_cont)
  forgetMetric_mk := fun _ _ _ => rfl
  mk_surjective := fun X x =>
    ⟨x.1, ⟨(x.2 : arcCM X).fn, arcCM.fn_cont _, picSheaf_locallyTrivial X x.1⟩, rfl⟩
  pullback := fun f x => (picardDataWitness.pullback f x.1, arcCMPullback f x.2)
  pullback_mul := fun f L M =>
    Prod.ext (picardDataWitness.pullback_mul f L.1 M.1) rfl

/-- ★★★★★★★**`APicData` は非空虚である**——D1 達成。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが Arakelov 理論の 9 件のうち **D1** である。 -/
theorem APicData.nonvacuous : Nonempty APicData :=
  ⟨aPicDataWitness⟩

/-! ## ★出典の紐付け(`.src`) -/

def aPicDataWitness.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(算術直線束の群 APic(X))",
    sectionId := "genell-def-1-1-i" }

end ABC3.Interface.Arakelov
