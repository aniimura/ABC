import ABC3.Found.Arakelov.ArcBasicOpen
import Mathlib.RingTheory.Localization.Away.Basic
import ABC3.Found.Arakelov.ArcFunctorial
import ABC3.Found.Arakelov.ArcLift

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

/-! ## ★★★合成の逆であること -/

/-- ★**基本開集合への構造射**。 -/
noncomputable def awayι (A : CommRingCat.{0}) (s : A) :
    Spec (CommRingCat.of (Localization.Away s)) ⟶ Spec A :=
  Spec.map (CommRingCat.ofHom (algebraMap A (Localization.Away s)))

/-- ★★★**`awayLift` は構造射との合成の切断である**。

    `awayLift A s p ≫ awayι A s = p.1`

★★これと `continuous_awayLift` で、**アフィン間の部分空間位相**が出る
——それが `topology_openImmersion` の基底事例である。 -/
theorem awayLift_comp_awayι (A : CommRingCat.{0}) (s : A) (p : arcBasicOpen A s) :
    awayLift A s p ≫ awayι A s = p.1 := by
  refine evalHom_injective A ?_
  rw [evalHom_comp, awayι, Spec.preimage_map, awayLift, evalHom_Spec_map]
  ext a
  exact awayLiftHom_algebraMap A s p.1 p.2 a

/-- ★★**逆に、構造射との合成は基本開集合に入る**。 -/
theorem comp_awayι_mem_arcBasicOpen (A : CommRingCat.{0}) (s : A)
    (q : Spec (CommRingCat.of ℂ) ⟶ Spec (CommRingCat.of (Localization.Away s))) :
    (q ≫ awayι A s) ∈ arcBasicOpen A s := by
  have h : evalAffine A (q ≫ awayι A s) s
      = evalAffine (CommRingCat.of (Localization.Away s)) q
          (algebraMap A (Localization.Away s) s) := by
    rw [evalAffine_comp, awayι, Spec.preimage_map]
    rfl
  rw [arcBasicOpen, Set.mem_setOf_eq, h]
  have hunit : IsUnit (algebraMap (↥A) (Localization.Away s) s) :=
    IsLocalization.Away.algebraMap_isUnit s
  intro hzero
  exact (isUnit_iff_ne_zero.1 (hunit.map (evalHom _ q).hom)) hzero

/-! ## ★★★逆向き —— 持ち上げの一意性 -/

/-- ★★★**構造射との合成を持ち上げるともとに戻る**。

    `awayLift A s ⟨q ≫ awayι A s, _⟩ = q`

★機構は**局所化からの準同型は基底での値で決まる**こと
(`IsLocalization.ringHom_ext`)。 -/
theorem awayLift_comp_awayι_self (A : CommRingCat.{0}) (s : A)
    (q : Spec (CommRingCat.of ℂ) ⟶ Spec (CommRingCat.of (Localization.Away s))) :
    awayLift A s ⟨q ≫ awayι A s, comp_awayι_mem_arcBasicOpen A s q⟩ = q := by
  refine evalHom_injective (CommRingCat.of (Localization.Away s)) ?_
  rw [awayLift, evalHom_Spec_map]
  ext x
  refine congrFun (congrArg (fun (r : Localization.Away s →+* ℂ) => (r : _ → _)) ?_) x
  refine IsLocalization.ringHom_ext (Submonoid.powers s) ?_
  ext a
  have h1 : awayLiftHom A s (q ≫ awayι A s)
        (comp_awayι_mem_arcBasicOpen A s q) (algebraMap _ _ a)
      = evalAffine A (q ≫ awayι A s) a :=
    awayLiftHom_algebraMap A s _ _ a
  have h2 : evalAffine A (q ≫ awayι A s) a
      = (evalHom (CommRingCat.of (Localization.Away s)) q).hom (algebraMap _ _ a) := by
    rw [evalAffine_comp, awayι, Spec.preimage_map]
    rfl
  exact h1.trans h2

/-- ★★**したがって `(· ≫ awayι)` は単射である**。 -/
theorem comp_awayι_injective (A : CommRingCat.{0}) (s : A) :
    Function.Injective
      (fun q : Spec (CommRingCat.of ℂ) ⟶ Spec (CommRingCat.of (Localization.Away s)) =>
        q ≫ awayι A s) := by
  intro q₁ q₂ h
  have h1 := awayLift_comp_awayι_self A s q₁
  have h2 := awayLift_comp_awayι_self A s q₂
  rw [← h1, ← h2]
  congr 1
  exact Subtype.ext h

/-! ## ★★★★基底事例 —— アフィン間の `topology_openImmersion` -/

/-- ★★★★**`topology_openImmersion` の基底事例**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

    arcTopologyAffine (Away s) = induced (· ≫ awayι A s) (arcTopologyAffine A)

★★★**両向き連続な全単射なので部分空間位相になる**——
4 条件(`awayLift_comp_awayι` / `awayLift_comp_awayι_self` /
`continuous_awayLift` / `continuous_comp_affine`)がすべて揃った帰結である。

★★これが `Interface/Arakelov/ArcSpace.lean` の
`topology_openImmersion` の**アフィンの場合**である。 -/
theorem arcTopologyAffine_away (A : CommRingCat.{0}) (s : A) :
    arcTopologyAffine (CommRingCat.of (Localization.Away s))
      = TopologicalSpace.induced
          (fun q : Spec (CommRingCat.of ℂ) ⟶ Spec (CommRingCat.of (Localization.Away s)) =>
            q ≫ awayι A s)
          (arcTopologyAffine A) := by
  letI := arcTopologyAffine A
  letI := arcTopologyAffine (CommRingCat.of (Localization.Away s))
  refine le_antisymm (continuous_iff_le_induced.1 (continuous_comp_affine (awayι A s))) ?_
  intro V hV
  refine isOpen_induced_iff.2 ⟨Subtype.val '' (awayLift A s ⁻¹' V), ?_, ?_⟩
  · exact (isOpen_arcBasicOpen A s).isOpenMap_subtype_val _
      ((continuous_awayLift A s).isOpen_preimage _ hV)
  · ext q
    simp only [Set.mem_preimage, Set.mem_image]
    constructor
    · rintro ⟨p, hp, heq⟩
      have : awayLift A s p = q := by
        have h2 : awayLift A s ⟨q ≫ awayι A s,
            comp_awayι_mem_arcBasicOpen A s q⟩ = q :=
          awayLift_comp_awayι_self A s q
        rw [← h2]
        congr 1
        exact Subtype.ext heq
      rwa [← this]
    · intro hq
      exact ⟨⟨q ≫ awayι A s, comp_awayι_mem_arcBasicOpen A s q⟩,
        by rwa [awayLift_comp_awayι_self], rfl⟩

/-! ## ★★★`awayι` は開埋め込みである -/

/-- ★★**`awayι A s` は開埋め込み**(mathlib、2026-08-17 実測で判明)。

★★★これにより `awayLift` を**一般の枠組み(`openLift`)に載せられる**。 -/
instance isOpenImmersion_awayι (A : CommRingCat.{0}) (s : A) :
    IsOpenImmersion (awayι A s) := by
  rw [awayι]
  exact Scheme.isOpenImmersion_SpecMap_localizationAway s

/-- ★★★**`awayLift` は `openLift` と一致する**。

★★これで一般の持ち上げの連続性が、基本開集合の上では
`continuous_awayLift`(既取得)に帰着する。 -/
theorem awayLift_eq_openLift (A : CommRingCat.{0}) (s : A) (p : arcBasicOpen A s)
    (h : (p.1).base default ∈ Set.range (awayι A s).base) :
    awayLift A s p = openLift (awayι A s) p.1 h := by
  refine comp_openImmersion_injective (awayι A s) ?_
  show awayLift A s p ≫ awayι A s = openLift (awayι A s) p.1 h ≫ awayι A s
  rw [awayLift_comp_awayι, openLift_comp]

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

def awayLift_comp_awayι.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——持ち上げが構造射の切断であること)",
    sectionId := "genell-def-1-1-i" }

def awayLift_eq_openLift.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——基本開集合の持ち上げが一般の持ち上げと一致すること)",
    sectionId := "genell-def-1-1-i" }

def arcTopologyAffine_away.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——開埋め込みで位相が誘導されること、基本開集合の場合)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
