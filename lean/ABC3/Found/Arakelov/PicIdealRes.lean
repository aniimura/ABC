import ABC3.Found.Arakelov.PicIdealBij

/-!
# Arakelov (B2) 第 159 ブロック —— **制限写像と `powers g` の可逆性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★可換図式(第 120 の版)へ向けた道具

第 120(`tildeAwayEquiv_res`)の `idealSections` 版を作るには、
`IsLocalizedModule.ext` を使うので

    ∀ s ∈ powers g, s は Γ(X, D(t·g)) の切断に可逆に作用する

が要る。★本ブロックはそれを用意する。

| 定理 | 内容 |
|---|---|
| `boMul_le` | ★`X.basicOpen (t·g) ≤ X.basicOpen g` |
| `idealResLin` | ★制限の `Γ(X,A)` 線型版 |
| `isUnit_res_g` | ★★`g` は `D(t·g)` 上で可逆(mathlib の `isUnit_res_basicOpen`) |
| `isUnit_smul_pow` | ★★★`powers g` は切断に可逆に作用する |

## ★★★逃げ道——**作用は `rfl` で一致する**

`r • z`(`Γ(X,A)` 作用)と `(res r) • z`(`Γ(X,D f)` 作用)は
`resAlg` が `RingHom.toAlgebra` なので **`rfl` で一致する**(実測)。
★これで `algebraMap_smul` を `rw` する必要が消える
——`rw` は `Module.End` を `A` に取ってしまい当たらなかった。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{u}} (D : X.IdealSheafData) (A : X.affineOpens)

/-- ★`X.basicOpen (t·g) ≤ X.basicOpen g`。 -/
theorem boMul_le (g t : (Γ(X, A.1) : Type u)) :
    X.basicOpen (t * g) ≤ X.basicOpen g := by
  rw [Scheme.basicOpen_mul]; exact inf_le_right

/-- ★制限の `Γ(X,A)` 線型版。 -/
noncomputable def idealResLin (g t : (Γ(X, A.1) : Type u)) :
    letI := resAlg A g
    letI := resAlg A (t * g)
    (idealSections D (X.basicOpen g)) →ₗ[(Γ(X, A.1) : Type u)]
      (idealSections D (X.basicOpen (t * g))) :=
  letI := resAlg A g
  letI := resAlg A (t * g)
  { toFun := fun s => ⟨(X.presheaf.map (homOfLE (boMul_le A g t)).op).hom s.1,
      idealSections_res D (boMul_le A g t) s.1 s.2⟩
    map_add' := fun a b => Subtype.ext (map_add _ _ _)
    map_smul' := fun r a => Subtype.ext (by
      show (X.presheaf.map (homOfLE (boMul_le A g t)).op).hom
          ((X.presheaf.map (homOfLE (X.basicOpen_le g)).op).hom r * a.1) = _
      rw [map_mul, ← CommRingCat.comp_apply, ← Functor.map_comp, ← op_comp]
      rfl) }

/-- ★`g` は `D(t·g)` 上で可逆。 -/
theorem isUnit_res_g (g t : (Γ(X, A.1) : Type u)) :
    IsUnit ((X.presheaf.map (homOfLE (le_trans (boMul_le A g t) (X.basicOpen_le g))).op).hom g) := by
  have h1 : IsUnit ((X.presheaf.map (homOfLE (X.basicOpen_le g)).op).hom g) :=
    X.toRingedSpace.isUnit_res_basicOpen g
  have := h1.map (X.presheaf.map (homOfLE (boMul_le A g t)).op).hom
  rwa [← CommRingCat.comp_apply, ← Functor.map_comp, ← op_comp] at this

/-- ★★`powers g` は切断に可逆に作用する。 -/
theorem isUnit_smul_pow (g t : (Γ(X, A.1) : Type u)) (s : Submonoid.powers g) :
    letI := resAlg A (t * g)
    IsUnit ((algebraMap (Γ(X, A.1) : Type u)
      (Module.End (Γ(X, A.1) : Type u)
        (idealSections D (X.basicOpen (t * g))))) (s : (Γ(X, A.1) : Type u))) := by
  letI := resAlg A (t * g)
  obtain ⟨n, hn⟩ := s.2
  rw [← hn, map_pow]
  refine IsUnit.pow n ?_
  refine (Module.End.isUnit_iff _).2 ?_
  obtain ⟨u, hu⟩ := isUnit_res_g A g t
  have hact : ∀ z : (idealSections D (X.basicOpen (t * g))),
      (algebraMap (Γ(X, A.1) : Type u)
        (Module.End (Γ(X, A.1) : Type u)
          (idealSections D (X.basicOpen (t * g))))) g z
      = ((u : (Γ(X, X.basicOpen (t * g)) : Type u)ˣ) : (Γ(X, X.basicOpen (t * g)) : Type u)) • z := by
    intro z
    show g • z = _
    rw [hu]
    rfl
  constructor
  · intro a b hab
    rw [hact, hact] at hab
    have h2 := congrArg (fun z => ((u⁻¹ : (Γ(X, X.basicOpen (t * g)) : Type u)ˣ)
      : (Γ(X, X.basicOpen (t * g)) : Type u)) • z) hab
    simpa [smul_smul, Units.inv_mul] using h2
  · intro y
    refine ⟨((u⁻¹ : (Γ(X, X.basicOpen (t * g)) : Type u)ˣ)
      : (Γ(X, X.basicOpen (t * g)) : Type u)) • y, ?_⟩
    rw [hact, smul_smul, Units.mul_inv, one_smul]


/-! ## ★出典の紐付け(`.src`) -/

def isUnit_smul_pow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——制限写像と powers g の可逆性)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
