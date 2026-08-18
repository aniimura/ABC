import ABC3.Found.Arakelov.PicSmulSurj

/-!
# Arakelov (B1) 第 140 ブロック —— **`IsLocalizing` の第 1・3 欄**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★`IsLocalizedModule` の 3 欄のうち 2 つ

| 欄 | 本ブロック |
|---|---|
| `map_units` | ★`isUnit_pow_smul_specD`(第 138・139 から) |
| `surj'` | (第 141 ブロック) |
| `exists_of_eq` | ★★★`exists_pow_smul_eq` |

## ★★`exists_of_eq` の筋——**最大値**

各 `D(gᵢ)` で第 137 が `∃ mᵢ, f^{mᵢ} x|ᵢ = f^{mᵢ} y|ᵢ` を与える。
★`N := Finset.univ.sup m` を取れば全ての `i` で `f^N x|ᵢ = f^N y|ᵢ`、
★★層は分離的だから `f^N x = f^N y`。

★★★ここで**有限**被覆(第 135)が効く——`Finset.sup` が取れる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `top_le_iSup_specD` | ★`⊤ ≤ ⨆ D(gᵢ)` |
| `exists_pow_smul_eq` | ★★★**`exists_of_eq` の中身** |
| `isUnit_pow_smul_specD` | ★★**`map_units` の中身** |
-/

universe u v

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable (R : CommRingCat.{u}) (F : (Spec R).Modules)


/-- ★`⊤ ≤ ⨆ D(gᵢ)`。 -/
theorem top_le_iSup_specD {ι : Type v} (g : ι → (R : Type u))
    (hspan : Ideal.span (Set.range g) = ⊤) :
    (⊤ : (Spec R).Opens) ≤ ⨆ i, specD R (g i) := by
  intro x _
  have htop : (⨆ i, PrimeSpectrum.basicOpen (g i)) = ⊤ :=
    PrimeSpectrum.iSup_basicOpen_eq_top_iff.2 hspan
  have hx2 : x ∈ (⨆ i, PrimeSpectrum.basicOpen (g i)) := by rw [htop]; trivial
  obtain ⟨i, hi⟩ := Opens.mem_iSup.1 hx2
  exact Opens.mem_iSup.2 ⟨i, hi⟩

/-- ★★★`res x = res y` なら `∃ N, fᴺ x = fᴺ y`。 -/
theorem exists_pow_smul_eq {n : ℕ} (g : Fin n → (R : Type u))
    (hspan : Ideal.span (Set.range g) = ⊤)
    (htriv : ∀ i, Nonempty ((restrictPresheafFunctor (Spec R) (specD R (g i))).obj F.val
      ≅ 𝟙_ (PresheafModulesOn (Spec R) (specD R (g i))))) (f : (R : Type u))
    (x y : (((modulesSpecToSheaf.obj F).obj.obj (op (⊤ : (Spec R).Opens))) : Type u))
    (h : ((modulesSpecToSheaf.obj F).obj.map (homOfLE (le_top (a := specD R f))).op).hom x
       = ((modulesSpecToSheaf.obj F).obj.map (homOfLE (le_top (a := specD R f))).op).hom y) :
    ∃ N : ℕ, f ^ N • x = f ^ N • y := by
  have key : ∀ i : Fin n, ∃ m : ℕ,
      f ^ m • ((modulesSpecToSheaf.obj F).obj.map (homOfLE (le_top (a := specD R (g i)))).op).hom x
        = f ^ m •
          ((modulesSpecToSheaf.obj F).obj.map (homOfLE (le_top (a := specD R (g i)))).op).hom y := by
    intro i
    haveI := isLocalizedModule_secRes R F (g i) f (htriv i).some
    have hres : ((modulesSpecToSheaf.obj F).obj.map (homOfLE (specDle R (g i) f)).op).hom
          (((modulesSpecToSheaf.obj F).obj.map (homOfLE (le_top (a := specD R (g i)))).op).hom x)
        = ((modulesSpecToSheaf.obj F).obj.map (homOfLE (specDle R (g i) f)).op).hom
          (((modulesSpecToSheaf.obj F).obj.map (homOfLE (le_top (a := specD R (g i)))).op).hom y) := by
      rw [← ConcreteCategory.comp_apply, ← ConcreteCategory.comp_apply, ← Functor.map_comp]
      have e1 : ((modulesSpecToSheaf.obj F).obj.map
          ((homOfLE (le_top (a := specD R (g i)))).op ≫ (homOfLE (specDle R (g i) f)).op))
          = ((modulesSpecToSheaf.obj F).obj.map
            ((homOfLE (le_top (a := specD R f))).op ≫ (homOfLE (specDle_left R f (g i))).op)) := by
        congr 1
      rw [e1, Functor.map_comp, ConcreteCategory.comp_apply, ConcreteCategory.comp_apply, h]
    obtain ⟨c, hc⟩ := IsLocalizedModule.exists_of_eq
      (S := Submonoid.powers f)
      (f := ((modulesSpecToSheaf.obj F).obj.map (homOfLE (specDle R (g i) f)).op).hom) hres
    obtain ⟨m, hm⟩ := c.2
    refine ⟨m, ?_⟩
    have := hc
    rw [Submonoid.smul_def, Submonoid.smul_def, ← hm] at this
    exact this
  choose m hm using key
  refine ⟨Finset.univ.sup m, ?_⟩
  refine TopCat.Sheaf.eq_of_locally_eq' (modulesSpecToSheaf.obj F)
    (fun i => specD R (g i)) ⊤ (fun i => homOfLE le_top)
    (top_le_iSup_specD R g hspan) _ _ (fun i => ?_)
  show ((modulesSpecToSheaf.obj F).obj.map _).hom (f ^ Finset.univ.sup m • x) = _
  rw [map_smul, map_smul]
  have hle : m i ≤ Finset.univ.sup m := Finset.le_sup (Finset.mem_univ i)
  rw [show Finset.univ.sup m = (Finset.univ.sup m - m i) + m i from
    (Nat.sub_add_cancel hle).symm, pow_add, mul_smul, mul_smul, hm i]

/-- ★★`powers f` の元は `Γ(F, D f)` に可逆に作用する。 -/
theorem isUnit_pow_smul_specD {ι : Type v} (g : ι → (R : Type u))
    (hspan : Ideal.span (Set.range g) = ⊤)
    (htriv : ∀ i, Nonempty ((restrictPresheafFunctor (Spec R) (specD R (g i))).obj F.val
      ≅ 𝟙_ (PresheafModulesOn (Spec R) (specD R (g i))))) (f : (R : Type u))
    (s : Submonoid.powers f) :
    IsUnit ((algebraMap (R : Type u) (Module.End (R : Type u)
      (((modulesSpecToSheaf.obj F).obj.obj (op (specD R f))) : Type u))) (s : (R : Type u))) := by
  obtain ⟨k, hk⟩ := s.2
  rw [← hk, map_pow]
  exact IsUnit.pow k ((Module.End.isUnit_iff _).2
    ⟨injective_smul_specD R F g hspan htriv f, surjective_smul_specD R F g hspan htriv f⟩)


/-! ## ★出典の紐付け(`.src`) -/

def exists_pow_smul_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——IsLocalizing の exists_of_eq)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
