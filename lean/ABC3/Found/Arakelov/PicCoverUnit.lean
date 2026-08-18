import ABC3.Found.Arakelov.PicSecRes

/-!
# Arakelov (B1) 第 138 ブロック —— **`f` は `D(f)` の切断に可逆に作用する**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★`IsLocalizedModule` の第 1 欄(`map_units`)へ向かう

    Γ(F,⊤) → Γ(F, D f) が powers f の局所化

の 3 欄のうち、第 1 欄は「`f` が `Γ(F,D f)` に**可逆に**作用する」である。

★★`D(f)` は `D(f·gᵢ)` たちで覆われ、各々では第 137 が
「`f` は可逆に作用する」を与える。★あとは層の分離性で貼るだけ。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `specD_le_iSup` | ★`D(f) ≤ ⨆ D(f·gᵢ)`(被覆) |
| `isUnit_smul_of_trivial` | ★★`f` は `Γ(F,D(f·g))` に可逆に作用 |
| `injective_smul_specD` | ★★★**`f` 倍は `Γ(F,D f)` 上で単射** |
-/

universe u v

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable (R : CommRingCat.{u}) (F : (Spec R).Modules)

/-- ★**`D(f)` は `D(f·gᵢ)` たちで覆われる**(`gᵢ` が単位イデアルを生成するとき)。 -/
theorem specD_le_iSup {ι : Type v} (g : ι → (R : Type u))
    (hspan : Ideal.span (Set.range g) = ⊤) (f : (R : Type u)) :
    specD R f ≤ ⨆ i, specD R (f * g i) := by
  intro x hx
  have htop : (⨆ i, PrimeSpectrum.basicOpen (g i)) = ⊤ :=
    PrimeSpectrum.iSup_basicOpen_eq_top_iff.2 hspan
  have hx2 : x ∈ (⨆ i, PrimeSpectrum.basicOpen (g i)) := by rw [htop]; trivial
  obtain ⟨i, hi⟩ := Opens.mem_iSup.1 hx2
  refine Opens.mem_iSup.2 ⟨i, ?_⟩
  show x ∈ PrimeSpectrum.basicOpen (f * g i)
  rw [PrimeSpectrum.basicOpen_mul]
  exact ⟨hx, hi⟩

/-- ★**`D(f·g) ≤ D(f)`**。 -/
theorem specDle_left (f g : (R : Type u)) : specD R (f * g) ≤ specD R f := by
  rw [specD, PrimeSpectrum.basicOpen_mul]; exact inf_le_left

/-- ★★**`f` は `Γ(F, D(f·g))` に可逆に作用する**(`F` が `D(g)` 上で自明なとき)。 -/
theorem isUnit_smul_of_trivial (g f : (R : Type u))
    (e : (restrictPresheafFunctor (Spec R) (specD R g)).obj F.val
      ≅ 𝟙_ (PresheafModulesOn (Spec R) (specD R g))) :
    IsUnit ((algebraMap (R : Type u)
      (Module.End (R : Type u)
        (((modulesSpecToSheaf.obj F).obj.obj (op (specD R (f * g)))) : Type u))) f) := by
  haveI := isLocalizedModule_secRes R F g f e
  exact IsLocalizedModule.map_units
    ((modulesSpecToSheaf.obj F).obj.map (homOfLE (specDle R g f)).op).hom
    ⟨f, Submonoid.mem_powers f⟩

/-- ★★★**`f` 倍は `Γ(F, D f)` 上で単射である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★各 `D(f·gᵢ)` で可逆、層は分離的だから大域でも単射。 -/
theorem injective_smul_specD {ι : Type v} (g : ι → (R : Type u))
    (hspan : Ideal.span (Set.range g) = ⊤)
    (htriv : ∀ i, Nonempty ((restrictPresheafFunctor (Spec R) (specD R (g i))).obj F.val
      ≅ 𝟙_ (PresheafModulesOn (Spec R) (specD R (g i))))) (f : (R : Type u)) :
    Function.Injective (fun s : (((modulesSpecToSheaf.obj F).obj.obj (op (specD R f))) : Type u) =>
      f • s) := by
  intro s t hst
  refine TopCat.Sheaf.eq_of_locally_eq' (modulesSpecToSheaf.obj F)
    (fun i => specD R (f * g i)) (specD R f)
    (fun i => homOfLE (specDle_left R f (g i)))
    (specD_le_iSup R g hspan f) s t (fun i => ?_)
  have hbij := (Module.End.isUnit_iff _).1 (isUnit_smul_of_trivial R F (g i) f (htriv i).some)
  refine hbij.1 ?_
  show f • ((modulesSpecToSheaf.obj F).obj.map _ s) = f • ((modulesSpecToSheaf.obj F).obj.map _ t)
  rw [← map_smul, ← map_smul]
  exact congrArg _ hst

/-! ## ★出典の紐付け(`.src`) -/

def injective_smul_specD.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——f 倍は D(f) の切断上で単射)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
