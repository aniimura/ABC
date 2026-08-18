import ABC3.Found.Arakelov.PicCoverUnit

/-!
# Arakelov (B1) 第 139 ブロック —— **`f` 倍は `D(f)` の切断上で全単射**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★層で貼る——`existsUnique_gluing'` の初使用

全射性の筋:

| 段 | 内容 |
|---|---|
| 1 | 各 `D(f·gᵢ)` で `f` 倍は全単射(第 138)なので `zᵢ` が取れる |
| 2 | 交わり `D(f·gᵢ) ⊓ D(f·gⱼ) = D(f·(f·gᵢ·gⱼ))` でも `f` 倍は**単射** |
| 3 | よって `zᵢ` は両立し、層で**貼れる** |
| 4 | 貼った `s` は `f • s = y` を局所的に満たすので大域でも満たす |

★★第 2 段が要点である——交わりも `D(f · h)` の形なので同じ補題が使える。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `bijective_smul_of_eq` | ★開集合を**等式で**指定できる形(`subst` で運ぶ) |
| `inf_specD_eq` | ★交わりも `D(f·h)` の形 |
| `specD_mul_le` | ★`D(f·a·b) ≤ D(a)` |
| `surjective_smul_specD` | ★★★★**`f` 倍は全射** |
-/

universe u v

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable (R : CommRingCat.{u}) (F : (Spec R).Modules)


/-- ★開集合を等式で指定できる形。 -/
theorem bijective_smul_of_eq (V : (Spec R).Opens) (g f : (R : Type u))
    (hV : V = specD R (f * g))
    (e : (restrictPresheafFunctor (Spec R) (specD R g)).obj F.val
      ≅ 𝟙_ (PresheafModulesOn (Spec R) (specD R g))) :
    Function.Bijective (fun s : (((modulesSpecToSheaf.obj F).obj.obj (op V)) : Type u) =>
      f • s) := by
  subst hV
  exact (Module.End.isUnit_iff _).1 (isUnit_smul_of_trivial R F g f e)

/-- ★交わりも `D(f · h)` の形である。 -/
theorem inf_specD_eq (f a b : (R : Type u)) :
    specD R (f * a) ⊓ specD R (f * b) = specD R (f * (f * a * b)) := by
  rw [show f * (f * a * b) = (f * a) * (f * b) from by ring]
  show _ = PrimeSpectrum.basicOpen ((f * a) * (f * b))
  rw [PrimeSpectrum.basicOpen_mul]
  rfl

/-- ★`D(f·a·b) ≤ D(a)`。 -/
theorem specD_mul_le (f a b : (R : Type u)) : specD R (f * a * b) ≤ specD R a := by
  show PrimeSpectrum.basicOpen (f * a * b) ≤ PrimeSpectrum.basicOpen a
  rw [PrimeSpectrum.basicOpen_mul, PrimeSpectrum.basicOpen_mul]
  exact le_trans inf_le_left inf_le_right

/-- ★★★★`f` 倍は `Γ(F, D f)` 上で全射である。 -/
theorem surjective_smul_specD {ι : Type v} (g : ι → (R : Type u))
    (hspan : Ideal.span (Set.range g) = ⊤)
    (htriv : ∀ i, Nonempty ((restrictPresheafFunctor (Spec R) (specD R (g i))).obj F.val
      ≅ 𝟙_ (PresheafModulesOn (Spec R) (specD R (g i))))) (f : (R : Type u)) :
    Function.Surjective (fun s : (((modulesSpecToSheaf.obj F).obj.obj (op (specD R f))) : Type u) =>
      f • s) := by
  intro y
  have hbij : ∀ i, Function.Bijective
      (fun s : (((modulesSpecToSheaf.obj F).obj.obj (op (specD R (f * g i)))) : Type u) =>
        f • s) := fun i => bijective_smul_of_eq R F _ (g i) f rfl (htriv i).some
  choose z hz using fun i => (hbij i).2
    (((modulesSpecToSheaf.obj F).obj.map (homOfLE (specDle_left R f (g i))).op).hom y)
  have hz' : ∀ i, f • z i =
      ((modulesSpecToSheaf.obj F).obj.map (homOfLE (specDle_left R f (g i))).op).hom y := hz
  have hcomp : TopCat.Presheaf.IsCompatible (modulesSpecToSheaf.obj F).obj
      (fun i => specD R (f * g i)) z := by
    intro i j
    have hb := bijective_smul_of_eq R F (specD R (f * g i) ⊓ specD R (f * g j))
      (f * g i * g j) f (inf_specD_eq R f (g i) (g j))
      (trivialOfLe F.val (specD_mul_le R f (g i) (g j)) (htriv i).some)
    refine hb.1 ?_
    show f • ((modulesSpecToSheaf.obj F).obj.map _).hom (z i)
      = f • ((modulesSpecToSheaf.obj F).obj.map _).hom (z j)
    rw [← map_smul, ← map_smul, hz' i, hz' j, ← ConcreteCategory.comp_apply,
      ← ConcreteCategory.comp_apply, ← Functor.map_comp, ← Functor.map_comp]
    rfl
  obtain ⟨s, hs, -⟩ := TopCat.Sheaf.existsUnique_gluing' (modulesSpecToSheaf.obj F)
    (fun i => specD R (f * g i)) (specD R f) (fun i => homOfLE (specDle_left R f (g i)))
    (specD_le_iSup R g hspan f) z hcomp
  refine ⟨s, ?_⟩
  refine TopCat.Sheaf.eq_of_locally_eq' (modulesSpecToSheaf.obj F)
    (fun i => specD R (f * g i)) (specD R f) (fun i => homOfLE (specDle_left R f (g i)))
    (specD_le_iSup R g hspan f) _ y (fun i => ?_)
  show ((modulesSpecToSheaf.obj F).obj.map _).hom (f • s) = _
  rw [map_smul, hs i]
  exact hz' i


/-! ## ★出典の紐付け(`.src`) -/

def surjective_smul_specD.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——f 倍は D(f) の切断上で全射)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
