import ABC3.Found.Arakelov.ArcFiberDim
import ABC3.Found.Arakelov.ArcScaleNorm

/-!
# Arakelov (C3) 第 289 ブロック —— **★★★★★計量の比はベクトルに依らない**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★Green 関数の土台

`logMetric`(Green 関数)は**基準計量に対する比**でしか定まらない。
その比が well-defined であるためには、比がどのベクトルで測っても同じでなければならない。

★★第 288 で「ファイバーの 0 でない元は生成元」を得たので、これが出る:

    w' = c • w  ⟹  m(w')/m'(w') = (‖c‖ m(w))/(‖c‖ m'(w)) = m(w)/m'(w)

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_ne_zero_unit` | ★単位加群のファイバーには 0 でない元がある |
| `exists_ne_zero_of_iso` | ★同型で運ぶ |
| `exists_ne_zero_arcFiber` | ★★★局所自明なら任意の複素点で |
| `fiberRatio` | ★比 |
| `fiberRatio_indep` | ★★★★★**比はベクトルに依らない** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite

/-- ★**単位加群のファイバーには 0 でない元がある**。 -/
theorem exists_ne_zero_unit : ∃ x : ((moduleSpecΓFunctor (R := CommRingCat.of ℂ)).obj
    (unitModules (Spec (CommRingCat.of ℂ))) : Type), x ≠ 0 := by
  obtain ⟨x, hx⟩ := unitCoord_bijective.2 1
  refine ⟨x, fun hz => ?_⟩
  rw [hz] at hx
  exact one_ne_zero ((((unitCoord_eq_zero_iff 0).2 rfl).symm.trans hx).symm)

/-- ★**0 でない元の存在は同型で運べる**。 -/
theorem exists_ne_zero_of_iso {R : CommRingCat.{0}} {M N : ModuleCat.{0} R} (φ : M ≅ N)
    (h : ∃ y : (N : Type), y ≠ 0) : ∃ x : (M : Type), x ≠ 0 := by
  obtain ⟨y, hy⟩ := h
  have hb : Function.Bijective φ.hom.hom := ConcreteCategory.bijective_of_isIso _
  obtain ⟨x, hx⟩ := hb.2 y
  exact ⟨x, fun hz => hy (by rw [← hx, hz, map_zero])⟩

variable {X : Scheme.{0}}

/-- ★★★**局所自明な層のファイバーには、任意の複素点で 0 でない元がある**。 -/
theorem exists_ne_zero_arcFiber (F : X.Modules) (hF : IsLocallyTrivial X F.val)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) : ∃ x : ↥(arcFiber p F), x ≠ 0 := by
  classical
  obtain ⟨S, hSmem, htriv⟩ := hF ⊤
  obtain ⟨W, g, hg, hW⟩ := hSmem (p.base default) trivial
  have hgt : S.arrows (homOfLE (le_top : W ≤ ⊤)) :=
    (Subsingleton.elim g (homOfLE (le_top : W ≤ ⊤))) ▸ hg
  have e : (restrictPresheafFunctor X W).obj F.val ≅ 𝟙_ (PresheafModulesOn X W) :=
    Classical.choice (htriv (homOfLE le_top) hgt)
  have htop : p ⁻¹ᵁ W = ⊤ := by
    refine preimage_eq_top_of_mem p W (fun z => ?_)
    rw [Subsingleton.elim z default]
    exact hW
  exact exists_ne_zero_of_iso
    (arcFiberAt W F p htop ≪≫ trivFiberIso (bridgeIso W F e) (liftToOpenOfTop W p htop))
    exists_ne_zero_unit

variable {F : X.Modules}

/-- ★**2 つの計量の比**(1 本のベクトルで測ったもの)。 -/
noncomputable def fiberRatio (m m' : ContArcMetric X F)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (w : ↥(arcFiber p F)) : ℝ :=
  m.nrm p w / m'.nrm p w

/-- ★★★★★**比はベクトルの選び方に依らない**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが Green 関数 `logMetric` を well-defined にする。 -/
theorem fiberRatio_indep (hF : IsLocallyTrivial X F.val) (m m' : ContArcMetric X F)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (w w' : ↥(arcFiber p F))
    (hw : w ≠ 0) (hw' : w' ≠ 0) :
    fiberRatio m m' p w = fiberRatio m m' p w' := by
  obtain ⟨c, hc⟩ := exists_smul_eq_arcFiber F hF p w w' hw
  have hcne : ‖c‖ ≠ 0 := by
    intro hz
    refine hw' ?_
    rw [← hc, show c = 0 from norm_eq_zero.1 hz, zero_smul]
  show m.nrm p w / m'.nrm p w = m.nrm p w' / m'.nrm p w'
  rw [← hc, m.smul, m'.smul, mul_div_mul_left _ _ hcne]

/-! ## ★出典の紐付け(`.src`) -/

def fiberRatio_indep.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C3——計量の比がベクトルに依らないこと)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
