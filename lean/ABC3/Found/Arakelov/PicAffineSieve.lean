import ABC3.Found.Arakelov.PicDivisorTop
import ABC3.Found.Arakelov.PicLocalBij

/-!
# Arakelov (B2) 第 189 ブロック —— **アフィン開で全単射なら局所全単射**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★再利用できる形にした

    (∀ アフィン開 A, f.app A が全単射)  ⟹  f は局所全単射

★★★これは `ofDivisor_mul` だけでなく、**イデアル層を扱うすべての段**で効く
——イデアル層の情報はアフィン開に載っているからである。

## ★★鍵は「像の篩は下向き閉」

篩の元がアフィンである必要は**ない**。アフィン開で覆う篩

    affineSieve U := { V ⟶ U | ∃ A アフィン開, V ≤ A ≤ U }

は下向き閉であり、その元 `V ≤ A` に対しては
**`A` での全単射性を制限で運べばよい**(自然性)。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `affineSieve` / `affineSieve_mem` | ★★アフィン開で覆う篩は被覆篩 |
| `presheafRes_trans` | ★制限の推移律(アーベル群の前層) |
| `locallyBijective_of_bijective_on_affine` | ★★★★**アフィンで全単射なら局所全単射** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{u}}

/-- ★アフィン開で覆う篩。 -/
def affineSieve (U : X.Opens) : Sieve U where
  arrows V _ := ∃ A : X.affineOpens, V ≤ A.1 ∧ A.1 ≤ U
  downward_closed := by
    rintro V W i ⟨A, hVA, hAU⟩ j
    exact ⟨A, le_trans (leOfHom j) hVA, hAU⟩

theorem affineSieve_mem (U : X.Opens) :
    affineSieve U ∈ (Opens.grothendieckTopology X) U := by
  intro x hx
  obtain ⟨A, hA, hxA, hAU⟩ :=
    Opens.isBasis_iff_nbhd.1 X.isBasis_affineOpens (U := U) (x := x) hx
  exact ⟨A, homOfLE hAU, ⟨⟨A, hA⟩, le_rfl, hAU⟩, hxA⟩



variable {X : Scheme.{u}}

theorem presheafRes_trans {F : (X.Opens)ᵒᵖ ⥤ AddCommGrpCat.{u}} {A B C : X.Opens}
    (h1 : C ≤ B) (h2 : B ≤ A) (z : (F.obj (op A) : Type u)) :
    F.map (homOfLE h1).op (F.map (homOfLE h2).op z) = F.map (homOfLE (h1.trans h2)).op z := by
  rw [← ConcreteCategory.comp_apply, ← F.map_comp]
  rfl

theorem locallyBijective_of_bijective_on_affine
    {F G : (X.Opens)ᵒᵖ ⥤ AddCommGrpCat.{u}} (f : F ⟶ G)
    (h : ∀ A : X.affineOpens, Function.Bijective (f.app (op A.1))) :
    Presheaf.IsLocallyInjective (Opens.grothendieckTopology X) f ∧
      Presheaf.IsLocallySurjective (Opens.grothendieckTopology X) f := by
  constructor
  · refine isLocallyInjective_of_coverSieve _ f (fun U x y hxy => ?_)
    refine ⟨affineSieve U, affineSieve_mem U, ?_⟩
    rintro V i ⟨A, hVA, hAU⟩
    have hA : F.map (homOfLE hAU).op x = F.map (homOfLE hAU).op y := by
      refine (h A).1 ?_
      rw [NatTrans.naturality_apply f (homOfLE hAU).op x,
        NatTrans.naturality_apply f (homOfLE hAU).op y, hxy]
    have h2 := congrArg (fun z => F.map (homOfLE hVA).op z) hA
    simp only [presheafRes_trans] at h2
    exact h2
  · refine isLocallySurjective_of_cover _ f (fun U s => ?_)
    refine ⟨affineSieve U, affineSieve_mem U, ?_⟩
    rintro V i ⟨A, hVA, hAU⟩
    obtain ⟨tA, htA⟩ := (h A).2 (G.map (homOfLE hAU).op s)
    refine ⟨F.map (homOfLE hVA).op tA, ?_⟩
    rw [NatTrans.naturality_apply f (homOfLE hVA).op tA, htA, presheafRes_trans]
    exact rfl

/-! ## ★出典の紐付け(`.src`) -/

def locallyBijective_of_bijective_on_affine.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——アフィン開で全単射なら局所全単射)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
