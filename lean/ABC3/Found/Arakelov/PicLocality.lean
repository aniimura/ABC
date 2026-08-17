import ABC3.Found.Arakelov.PicLocalInj

/-!
# Arakelov (B1) 第 8 ブロック —— **局所単射性は局所的である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★これが「貼り合わせ」の核である

第 7 ブロックで局所単射性の**基底**(`M = 𝟙_`)と**輸送**(`M ≅ N`)を取った。
残るのは「**被覆の各員で成り立てば全体で成り立つ**」である。

★★★mathlib に無かったので自作する(2026-08-17 実測)。
★機構は `GrothendieckTopology.transitive` と、
「等化篩の引き戻しは制限の等化篩」という補題だけである。

## ★なぜ局所性が要るか

可逆層 `M` は**局所的にしか**構造層と同型でない。
★したがって `f ▷ M` の局所単射性も、まず局所的に示してから貼り合わせるほかない。
-/

universe w v u

namespace ABC3.Found.Arakelov

open CategoryTheory Opposite

variable {C : Type u} [Category.{v} C] (J : GrothendieckTopology C)
  {D : Type*} [Category* D] {FD : D → D → Type*} {CD : D → Type w}
  [∀ X Y, FunLike (FD X Y) (CD X) (CD Y)] [ConcreteCategory.{w} D FD]

/-- ★**等化篩の引き戻しは、制限した切断の等化篩である**。 -/
theorem pullback_equalizerSieve {F : Cᵒᵖ ⥤ D} {X : C} (x y : ToType (F.obj (op X)))
    {V : C} (i : V ⟶ X) :
    (Presheaf.equalizerSieve x y).pullback i
      = Presheaf.equalizerSieve (F.map i.op x) (F.map i.op y) := by
  ext W j
  show F.map (j ≫ i).op x = F.map (j ≫ i).op y ↔ _
  show F.map (j ≫ i).op x = F.map (j ≫ i).op y ↔
    F.map j.op (F.map i.op x) = F.map j.op (F.map i.op y)
  rw [op_comp, F.map_comp]
  simp only [ConcreteCategory.comp_apply]

/-- ★★★**局所単射性は局所的である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★各切断の対について、**被覆篩を 1 つ見つけて、その各員の上で等化篩が被覆**
であればよい。★機構は `GrothendieckTopology.transitive` である。

★★これで「可逆層は局所的に構造層」を使う道が開く。 -/
theorem isLocallyInjective_of_cover {F₁ F₂ : Cᵒᵖ ⥤ D} (φ : F₁ ⟶ F₂)
    (h : ∀ (X : C) (x y : ToType (F₁.obj (op X))), φ.app (op X) x = φ.app (op X) y →
      ∃ S : Sieve X, S ∈ J X ∧ ∀ ⦃V : C⦄ (i : V ⟶ X), S i →
        Presheaf.equalizerSieve (F₁.map i.op x) (F₁.map i.op y) ∈ J V) :
    Presheaf.IsLocallyInjective J φ where
  equalizerSieve_mem {X} x y hxy := by
    obtain ⟨S, hS, hcov⟩ := h X.unop x y hxy
    refine J.transitive hS _ ?_
    intro V i hi
    rw [pullback_equalizerSieve]
    exact hcov i hi

/-! ## ★出典の紐付け(`.src`) -/

def pullback_equalizerSieve.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——等化篩の引き戻し)",
    sectionId := "genell-def-1-1-i" }

def isLocallyInjective_of_cover.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——局所単射性が局所的であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
