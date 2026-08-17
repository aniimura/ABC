import ABC3.Found.Arakelov.ArcBasicOpen
import Mathlib.RingTheory.Localization.Away.Basic

/-!
# Arakelov (C1) の配管 —— **基本開集合への持ち上げ(環の水準)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★ここは「新たな数学」ではなく配管である

`ArcBasicOpen.lean` で**解析的核心**(`φ(b)/φ(s)ⁿ` の連続性)を取った。
★本ファイルはそれを **mathlib の局所化 API で包む**段である:

    Localization.awayLift (evalHom A p).hom s hs : Localization.Away s →+* ℂ

★★これが `Arc (D(s))` への持ち上げの環の水準の姿である
(`Arc (Spec B) = Hom_Ring(B, ℂ)` は `ArcEval.lean` の全単射)。

## ★★取るもの

| 定理 | 内容 |
|---|---|
| `awayLiftHom` | 持ち上げの環準同型 |
| `awayLiftHom_algebraMap` | ★もとの `φ` の延長であること |
| `awayLiftHom_mk` | ★★`b/sⁿ ↦ φ(b)/φ(s)ⁿ` という**明示式** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory

/-! ## ★★持ち上げ(環の水準) -/

/-- ★★**基本開集合への持ち上げ**——`φ : A → ℂ`(`φ(s) ≠ 0`)を
`Localization.Away s → ℂ` へ延ばす。 -/
noncomputable def awayLiftHom (A : CommRingCat.{0}) (s : A)
    (p : Spec (CommRingCat.of ℂ) ⟶ Spec A) (h : evalAffine A p s ≠ 0) :
    Localization.Away s →+* ℂ :=
  Localization.awayLift (evalHom A p).hom s (isUnit_iff_ne_zero.2 h)

/-- ★**持ち上げはもとの `φ` の延長である**。 -/
@[simp] theorem awayLiftHom_algebraMap (A : CommRingCat.{0}) (s : A)
    (p : Spec (CommRingCat.of ℂ) ⟶ Spec A) (h : evalAffine A p s ≠ 0) (a : A) :
    awayLiftHom A s p h (algebraMap A (Localization.Away s) a) = evalAffine A p a :=
  IsLocalization.lift_eq _ a

/-- ★★★**明示式** `b/sⁿ ↦ φ(b)/φ(s)ⁿ`。

★★★これと `ArcBasicOpen.lean` の `continuousOn_div_pow_evalAffine` が対になって、
**持ち上げの連続性**を与える。 -/
theorem awayLiftHom_mk (A : CommRingCat.{0}) (s : A)
    (p : Spec (CommRingCat.of ℂ) ⟶ Spec A) (h : evalAffine A p s ≠ 0) (b : A) (n : ℕ) :
    awayLiftHom A s p h (Localization.mk b ⟨s ^ n, n, rfl⟩)
      = evalAffine A p b / (evalAffine A p s) ^ n := by
  rw [awayLiftHom, Localization.mk_eq_mk', Localization.awayLift,
    IsLocalization.Away.lift, IsLocalization.lift_mk'_spec]
  show (evalHom A p).hom b
      = (evalHom A p).hom ((⟨s ^ n, n, rfl⟩ : Submonoid.powers s) : A)
        * (evalAffine A p b / (evalAffine A p s) ^ n)
  have hs : (evalAffine A p s) ^ n ≠ 0 := pow_ne_zero n h
  simp only [map_pow, evalAffine] at hs ⊢
  field_simp

/-! ## ★★★スキームの射への包みと連続性 -/

/-- ★★**持ち上げをスキームの射として見る**。 -/
noncomputable def awayLift (A : CommRingCat.{0}) (s : A)
    (p : arcBasicOpen A s) :
    Spec (CommRingCat.of ℂ) ⟶ Spec (CommRingCat.of (Localization.Away s)) :=
  Spec.map (CommRingCat.ofHom (awayLiftHom A s p.1 p.2))

/-- ★**包んでも値は変わらない**。 -/
@[simp] theorem evalAffine_awayLift (A : CommRingCat.{0}) (s : A)
    (p : arcBasicOpen A s) (x : Localization.Away s) :
    evalAffine (CommRingCat.of (Localization.Away s)) (awayLift A s p) x
      = awayLiftHom A s p.1 p.2 x := by
  rw [awayLift, evalAffine, evalHom_Spec_map]
  rfl

/-- ★★★**持ち上げは連続である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★**これが C1 の残り 1 点の本体である。**

★機構は 2 つの対:
- 代数側 `awayLiftHom_mk`(値は `φ(b)/φ(s)ⁿ`)
- 解析側 `continuousOn_div_pow_evalAffine`(その式が連続)

★★局所化の任意の元は `b/sⁿ` の形なので(`Localization.induction_on`)、
成分ごとにこの対を当てればよい。 -/
theorem continuous_awayLift (A : CommRingCat.{0}) (s : A) :
    @Continuous _ _ (@instTopologicalSpaceSubtype _ _ (arcTopologyAffine A))
      (arcTopologyAffine (CommRingCat.of (Localization.Away s)))
      (awayLift A s) := by
  letI := arcTopologyAffine A
  letI := arcTopologyAffine (CommRingCat.of (Localization.Away s))
  refine continuous_induced_rng.2 (continuous_pi fun x => ?_)
  induction x using Localization.induction_on with
  | H x =>
    obtain ⟨b, y⟩ := x
    obtain ⟨n, hn⟩ := y.2
    have hval : ∀ p : arcBasicOpen A s,
        evalAffine (CommRingCat.of (Localization.Away s)) (awayLift A s p)
            (Localization.mk b y)
          = evalAffine A p.1 b / (evalAffine A p.1 s) ^ n := by
      intro p
      rw [evalAffine_awayLift]
      have : y = (⟨s ^ n, n, rfl⟩ : Submonoid.powers s) := Subtype.ext hn.symm
      rw [this, awayLiftHom_mk]
    simp only [Function.comp_apply, hval]
    exact ((continuousOn_div_pow_evalAffine A b s n).restrict).congr (fun p => rfl)

/-! ## ★出典の紐付け(`.src`) -/

def awayLiftHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——基本開集合への持ち上げ、環の水準)",
    sectionId := "genell-def-1-1-i" }

def awayLiftHom_mk.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——持ち上げの明示式)",
    sectionId := "genell-def-1-1-i" }

def continuous_awayLift.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——基本開集合への持ち上げの連続性)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
