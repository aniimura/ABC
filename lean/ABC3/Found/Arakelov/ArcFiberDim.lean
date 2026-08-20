import ABC3.Found.Arakelov.ArcUnitCoord
import ABC3.Found.Arakelov.ArcMetricExists

/-!
# Arakelov (C3) 第 288 ブロック —— **★★★★★★どの複素点でもファイバーは 1 次元**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★第 287 を局所自明な層へ運ぶ

第 287 で「単位加群のファイバーでは 0 でない元が生成元」を得た。
★同型で運べばそのまま任意の自明化された層に移る。

## ★★機構

| 段 | 使うもの |
|---|---|
| 1 | `exists_smul_eq_of_iso` —— 同型で性質を運ぶ(一般形) |
| 2 | 第 253 `trivFiberIso` —— 自明化からファイバーの同型 |
| 3 | 第 271 `arcFiberAt` + 第 269 `bridgeIso` —— `X` の点を `V` の点へ |
| 4 | 第 284 `preimage_eq_top_of_mem` —— 篩から自明化開集合を取る |

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_smul_eq_of_iso` | ★同型で運ぶ(一般形) |
| `exists_smul_eq_triv` | ★★自明化された層のファイバー |
| `exists_smul_eq_arcFiber` | ★★★★★★**局所自明な層の任意の複素点で** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite

/-- ★**「0 でない元は生成元」は同型で運べる**。 -/
theorem exists_smul_eq_of_iso {R : CommRingCat.{0}} {M N : ModuleCat.{0} R} (φ : M ≅ N)
    (h : ∀ x y : (N : Type), x ≠ 0 → ∃ c : (R : Type), c • x = y)
    (x y : (M : Type)) (hx : x ≠ 0) : ∃ c : (R : Type), c • x = y := by
  have hb : Function.Bijective φ.hom.hom := ConcreteCategory.bijective_of_isIso _
  have hx' : φ.hom.hom x ≠ 0 := by
    intro hz
    exact hx (hb.1 (by rw [hz, map_zero]))
  obtain ⟨c, hc⟩ := h (φ.hom.hom x) (φ.hom.hom y) hx'
  refine ⟨c, hb.1 ?_⟩
  rw [map_smul, hc]

variable {V : Scheme.{0}} {L : V.Modules}

/-- ★★**自明化された層のファイバーでも 0 でない元は生成元**。 -/
theorem exists_smul_eq_triv (t : L ≅ unitModules V) (p : Spec (CommRingCat.of ℂ) ⟶ V)
    (x y : ↥(arcFiber p L)) (hx : x ≠ 0) :
    ∃ c : (CommRingCat.of ℂ : Type), c • x = y :=
  exists_smul_eq_of_iso (trivFiberIso t p) (fun a b ha => exists_smul_eq_unit a b ha) x y hx

variable {X : Scheme.{0}}

/-- ★★★★★★**局所自明な層なら、任意の複素点でファイバーは 1 次元**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで「計量の比がベクトルの選び方に依らない」が全点で言える。 -/
theorem exists_smul_eq_arcFiber (F : X.Modules) (hF : IsLocallyTrivial X F.val)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (x y : ↥(arcFiber p F)) (hx : x ≠ 0) :
    ∃ c : (CommRingCat.of ℂ : Type), c • x = y := by
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
  refine exists_smul_eq_of_iso
    (arcFiberAt W F p htop ≪≫ trivFiberIso (bridgeIso W F e) (liftToOpenOfTop W p htop))
    (fun a b ha => exists_smul_eq_unit a b ha) x y hx

/-! ## ★出典の紐付け(`.src`) -/

def exists_smul_eq_arcFiber.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C3——任意の複素点でファイバーが 1 次元であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
