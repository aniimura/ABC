/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.HomColim
import Mathlib.Data.Setoid.Basic

/-!
# filtered な帰納極限は商と交換する

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.105。

原文 (FrdI p.105):
> erality that C is of isotropic type. Then it follows immediately from the definition of

## ★★何のためか

`Proposition 5.5, (ii)` の**骨**である。原文は

> it suffices to obtain natural bijections between the respective sets of morphisms
> between the images of two given objects of C

と言う。★`(𝒞^pf)^un-tr` の射の集合は「**帰納極限をとってから商をとった**もの」、
`(𝒞^un-tr)^pf` の射の集合は「**商をとってから帰納極限をとった**もの」である。
★★**その 2 つが一致する** —— これが「follows immediately from the definitions」の中身。

## ★層を分ける理由

`Proposition 5.5, (ii)` は `un-tr` と `birat` の 2 つを主張する。
★どちらも「Hom の商」であり、添字圏(`𝒞^{bi-Fr}` / `𝒞^{coa-pre}`)だけが違う。
★★そこで**添字圏に依らない形**でここに置く。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `QuotData.functor` | 商をとった帰納系 |
| `toQuot` | `colim F → colim (F/∼)` |
| `toQuot_surjective` | 全射 |
| `toQuot_eq_iff` | ★**核の完全な記述** —— 「ある段で同値になる」 |
| `quotEquiv` | ★★**`(colim F)/∼ ≃ colim (F/∼)`** |
-/

namespace ABC3.Found.FrdI

open CategoryTheory CategoryTheory.Limits

namespace HomColim

universe uJ vJ w

variable {J : Type uJ} [Category.{vJ} J]

/-! ## ★1. 遷移写像と両立する同値関係の族 -/

/-- ★★**遷移写像と両立する同値関係の族**。

★`compat` が「遷移写像が商へ降りる」ことそのもの。 -/
structure QuotData (F : J ⥤ Type w) where
  /-- 各段の同値関係。 -/
  setoid : ∀ j : J, Setoid (F.obj j)
  /-- 遷移写像との両立。 -/
  compat : ∀ {i j : J} (f : i ⟶ j) {x y : F.obj i},
    (setoid i).r x y → (setoid j).r (F.map f x) (F.map f y)

variable {F : J ⥤ Type w}

/-- ★**商をとった帰納系** —— 各段を同値関係で割る。 -/
def QuotData.functor (q : QuotData F) : J ⥤ Type w where
  obj j := Quotient (q.setoid j)
  map {i j} f := TypeCat.ofHom (Quotient.map (F.map f) (fun _ _ h => q.compat f h))
  map_id j := by
    ext x
    refine Quotient.inductionOn x (fun a => ?_)
    show Quotient.mk (q.setoid j) (F.map (𝟙 j) a) = Quotient.mk (q.setoid j) a
    rw [F.map_id]
    rfl
  map_comp {i j k} f g := by
    ext x
    refine Quotient.inductionOn x (fun a => ?_)
    show Quotient.mk (q.setoid k) (F.map (f ≫ g) a)
      = Quotient.mk (q.setoid k) (F.map g (F.map f a))
    rw [F.map_comp]
    rfl

@[simp] theorem QuotData.functor_map_mk (q : QuotData F) {i j : J} (f : i ⟶ j) (x : F.obj i) :
    q.functor.map f (Quotient.mk (q.setoid i) x) = Quotient.mk (q.setoid j) (F.map f x) := rfl

/-! ## ★2. 自然な写像 `colim F → colim (F/∼)` -/

variable [HasColimit F]

/-- ★商の帰納系への余錐。 -/
noncomputable def quotCocone (q : QuotData F) [HasColimit q.functor] : Cocone F :=
  Cocone.mk (HomColim q.functor)
    { app := fun j => TypeCat.ofHom fun x => mk q.functor j (Quotient.mk (q.setoid j) x)
      naturality := by
        intro i j f
        ext x
        show mk q.functor j (Quotient.mk (q.setoid j) (F.map f x))
          = mk q.functor i (Quotient.mk (q.setoid i) x)
        rw [← q.functor_map_mk f x]
        exact mk_map q.functor f _ }

/-- ★★**自然な写像 `colim F → colim (F/∼)`**。 -/
noncomputable def toQuot (q : QuotData F) [HasColimit q.functor] :
    HomColim F → HomColim q.functor :=
  fun z => colimit.desc F (quotCocone q) z

@[simp] theorem toQuot_mk (q : QuotData F) [HasColimit q.functor] (j : J) (x : F.obj j) :
    toQuot q (mk F j x) = mk q.functor j (Quotient.mk (q.setoid j) x) := by
  show colimit.desc F (quotCocone q) (colimit.ι F j x) = _
  rw [← types_comp_apply (colimit.ι F j) (colimit.desc F (quotCocone q)), colimit.ι_desc]
  rfl

/-- ★**全射** —— 代表元を 2 段取るだけ。 -/
theorem toQuot_surjective (q : QuotData F) [HasColimit q.functor] :
    Function.Surjective (toQuot q) := by
  refine induction (P := fun z => ∃ w, toQuot q w = z) q.functor (fun j y => ?_)
  refine Quotient.inductionOn y (fun x => ?_)
  exact ⟨mk F j x, toQuot_mk q j x⟩

/-! ## ★3. 核の完全な記述 -/

/-- ★★★**核の完全な記述** —— 帰納極限で同じ像を持つのは、
**ある段で同値になる**とき、ちょうどそのとき。

★★これが「filtered な帰納極限は商と交換する」の中身である。 -/
theorem toQuot_eq_iff [IsFiltered J] (q : QuotData F) [HasColimit q.functor]
    {i j : J} {x : F.obj i} {y : F.obj j} :
    toQuot q (mk F i x) = toQuot q (mk F j y)
      ↔ ∃ (k : J) (f : i ⟶ k) (g : j ⟶ k), (q.setoid k).r (F.map f x) (F.map g y) := by
  rw [toQuot_mk, toQuot_mk, eq_iff]
  constructor
  · rintro ⟨k, f, g, h⟩
    refine ⟨k, f, g, ?_⟩
    exact Quotient.exact
      (a := F.map f x) (b := F.map g y) (s := q.setoid k) h
  · rintro ⟨k, f, g, h⟩
    exact ⟨k, f, g, Quotient.sound h⟩

/-- ★同じ段での形。 -/
theorem toQuot_eq_iff_same [IsFiltered J] (q : QuotData F) [HasColimit q.functor]
    {i : J} {x y : F.obj i} :
    toQuot q (mk F i x) = toQuot q (mk F i y)
      ↔ ∃ (k : J) (f : i ⟶ k), (q.setoid k).r (F.map f x) (F.map f y) := by
  rw [toQuot_eq_iff]
  constructor
  · rintro ⟨k, f, g, h⟩
    refine ⟨IsFiltered.coeq f g, f ≫ IsFiltered.coeqHom f g, ?_⟩
    have hfg : f ≫ IsFiltered.coeqHom f g = g ≫ IsFiltered.coeqHom f g :=
      IsFiltered.coeq_condition f g
    have h1 : F.map (f ≫ IsFiltered.coeqHom f g) x
        = F.map (IsFiltered.coeqHom f g) (F.map f x) := by
      rw [F.map_comp]; rfl
    have h2 : F.map (f ≫ IsFiltered.coeqHom f g) y
        = F.map (IsFiltered.coeqHom f g) (F.map g y) := by
      rw [hfg, F.map_comp]; rfl
    rw [h1, h2]
    exact q.compat _ h
  · rintro ⟨k, f, h⟩
    exact ⟨k, f, f, h⟩

/-- ★同値関係になっていることを明示する —— `toQuot` の核。 -/
noncomputable abbrev quotKer (q : QuotData F) [HasColimit q.functor] : Setoid (HomColim F) :=
  Setoid.ker (toQuot q)

/-- ★★★★**`(colim F)/∼ ≃ colim (F/∼)`** —— 帰納極限は商と交換する。 -/
noncomputable def quotEquiv (q : QuotData F) [HasColimit q.functor] :
    Quotient (quotKer q) ≃ HomColim q.functor :=
  Setoid.quotientKerEquivOfSurjective _ (toQuot_surjective q)

@[simp] theorem quotEquiv_apply (q : QuotData F) [HasColimit q.functor] (z : HomColim F) :
    quotEquiv q (Quotient.mk (quotKer q) z) = toQuot q z := rfl

end HomColim

end ABC3.Found.FrdI
