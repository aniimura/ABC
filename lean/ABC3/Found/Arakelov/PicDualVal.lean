import ABC3.Found.Arakelov.PicDualGlue

/-!
# Arakelov (B2) 第 177 ブロック —— **貼り合わせた値は加法的・線型・自然**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★双対が層であることの第 2 歩

第 176 で「局所値は貼り合わせ可能」まで出た。本ブロックは

    dualGlueVal : F(V) → 𝒪(V)   (V ≤ ⨆ Uᵢ)

を作り、★**加法的・`𝒪(V)` 線型・制限と可換**であることを示す。
これで第 178 で「射 `s : F|_{⨆U} ⟶ 𝟙_`」が組める。

## ★★証明はすべて「一意性で殴る」形

`X.sheaf.eq_of_locally_eq'` により、**`V ⊓ Uᵢ` へ制限して一致すれば等しい**。
制限した先では `dualGlueVal_res` で局所値に落ちるので、
★**局所値についての命題**(`dualLocVal_add` / `dualLocVal_smul` / `dualLocVal_res`)
に帰着する。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `resTrans` | ★`𝒪` の制限の推移律 |
| `modResTrans` | ★`F` の制限の推移律 |
| `dualLocVal_apply` | ★局所値の展開(`rfl`) |
| `dualLocVal_add` / `dualLocVal_smul` | ★局所値は加法的・線型 |
| `dualGlueVal` | ★★**貼り合わせた値** |
| `dualGlueVal_res` | ★★貼り合わせた値の特徴づけ |
| `dualGlue_ext` | ★局所的に一致すれば等しい |
| `dualGlueVal_add` | ★★★**加法的** |
| `dualGlueVal_smul` | ★★★**`𝒪(V)` 線型** |
| `dualGlueVal_map` | ★★★★**制限と可換**(自然性) |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} {F : X.PresheafOfModules}

/-- ★**`𝒪` の制限の推移律**。 -/
theorem resTrans {A B C : X.Opens} (h1 : C ≤ B) (h2 : B ≤ A) (a : (Γ(X, A) : Type u)) :
    (X.presheaf.map (homOfLE h1).op).hom ((X.presheaf.map (homOfLE h2).op).hom a)
      = (X.presheaf.map (homOfLE (h1.trans h2)).op).hom a := by
  rw [← ConcreteCategory.comp_apply, ← Functor.map_comp]
  rfl

/-- ★**`F` の制限の推移律**。 -/
theorem modResTrans {A B C : X.Opens} (h1 : C ≤ B) (h2 : B ≤ A) (a : (F.obj (op A) : Type u)) :
    (F.map (homOfLE h1).op).hom ((F.map (homOfLE h2).op).hom a)
      = (F.map (homOfLE (h1.trans h2)).op).hom a :=
  (congrArg (fun (m : _ ⟶ _) => (ModuleCat.Hom.hom m) a)
    (F.map_comp (homOfLE h2).op (homOfLE h1).op)).symm

variable {ι : Type u} {U : ι → X.Opens} (φ : ∀ i, ((dualPresheaf F).obj (op (U i)) : Type u))

/-- ★**局所値の展開**(`rfl`)。 -/
theorem dualLocVal_apply (V : X.Opens) (x : (F.obj (op V) : Type u)) (i : ι) :
    dualLocVal φ V x i
      = ((φ i).app (op (infOver U V i))).hom
          ((F.map (homOfLE (inf_le_left : V ⊓ U i ≤ V)).op).hom x) := rfl

/-- ★**局所値は加法的**。 -/
theorem dualLocVal_add (V : X.Opens) (x y : (F.obj (op V) : Type u)) (i : ι) :
    dualLocVal φ V (x + y) i = dualLocVal φ V x i + dualLocVal φ V y i :=
  (congrArg ((φ i).app (op (infOver U V i))).hom
    (map_add (F.map (homOfLE (inf_le_left : V ⊓ U i ≤ V)).op).hom x y)).trans
    (map_add _ _ _)

/-- ★**局所値は線型**——係数は `V ⊓ Uᵢ` へ制限される。 -/
theorem dualLocVal_smul (V : X.Opens) (c : (Γ(X, V) : Type u)) (x : (F.obj (op V) : Type u))
    (i : ι) :
    dualLocVal φ V (c • x) i
      = (X.presheaf.map (homOfLE (inf_le_left : V ⊓ U i ≤ V)).op).hom c * dualLocVal φ V x i :=
  (congrArg ((φ i).app (op (infOver U V i))).hom
    (F.map_smul (homOfLE (inf_le_left : V ⊓ U i ≤ V)).op c x)).trans
    (((φ i).app (op (infOver U V i))).hom.map_smul _ _)

variable (hφ : TopCat.Presheaf.IsCompatible (dualPresheaf F).presheaf U φ)

/-- ★★**貼り合わせた値**——`V ⊓ Uᵢ` の上の値を `X.sheaf` で貼る。 -/
noncomputable def dualGlueVal {V : X.Opens} (hV : V ≤ ⨆ i, U i) (x : (F.obj (op V) : Type u)) :
    (Γ(X, V) : Type u) :=
  (X.sheaf.existsUnique_gluing' (U := fun i => V ⊓ U i) V (fun _ => homOfLE inf_le_left)
    (le_iSup_inf U V hV) (dualLocVal φ V x) (dualLocVal_compat φ hφ V x)).choose

/-- ★★**貼り合わせた値の特徴づけ**。 -/
theorem dualGlueVal_res {V : X.Opens} (hV : V ≤ ⨆ i, U i) (x : (F.obj (op V) : Type u)) (i : ι) :
    (X.presheaf.map (homOfLE (inf_le_left : V ⊓ U i ≤ V)).op).hom (dualGlueVal φ hφ hV x)
      = dualLocVal φ V x i :=
  (X.sheaf.existsUnique_gluing' (U := fun i => V ⊓ U i) V (fun _ => homOfLE inf_le_left)
    (le_iSup_inf U V hV) (dualLocVal φ V x) (dualLocVal_compat φ hφ V x)).choose_spec.1 i

/-- ★**局所的に一致すれば等しい**。 -/
theorem dualGlue_ext {V : X.Opens} (hV : V ≤ ⨆ i, U i) (a b : (Γ(X, V) : Type u))
    (h : ∀ i, (X.presheaf.map (homOfLE (inf_le_left : V ⊓ U i ≤ V)).op).hom a
             = (X.presheaf.map (homOfLE (inf_le_left : V ⊓ U i ≤ V)).op).hom b) : a = b :=
  X.sheaf.eq_of_locally_eq' (fun i => V ⊓ U i) V (fun _ => homOfLE inf_le_left)
    (le_iSup_inf U V hV) a b h

/-- ★★★**貼り合わせた値は加法的である**。 -/
theorem dualGlueVal_add {V : X.Opens} (hV : V ≤ ⨆ i, U i) (x y : (F.obj (op V) : Type u)) :
    dualGlueVal φ hφ hV (x + y) = dualGlueVal φ hφ hV x + dualGlueVal φ hφ hV y := by
  refine dualGlue_ext hV _ _ (fun i => ?_)
  rw [dualGlueVal_res, map_add, dualGlueVal_res, dualGlueVal_res, dualLocVal_add]

/-- ★★★**貼り合わせた値は `𝒪(V)` 線型である**。 -/
theorem dualGlueVal_smul {V : X.Opens} (hV : V ≤ ⨆ i, U i) (c : (Γ(X, V) : Type u))
    (x : (F.obj (op V) : Type u)) :
    dualGlueVal φ hφ hV (c • x) = c * dualGlueVal φ hφ hV x := by
  refine dualGlue_ext hV _ _ (fun i => ?_)
  rw [dualGlueVal_res, map_mul, dualGlueVal_res, dualLocVal_smul]

/-- ★★★★**貼り合わせた値は制限と可換である**(自然性)。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが「射になる」ことの最後の条件である。 -/
theorem dualGlueVal_map {V W : X.Opens} (hV : V ≤ ⨆ i, U i) (hWV : W ≤ V)
    (x : (F.obj (op V) : Type u)) :
    (X.presheaf.map (homOfLE hWV).op).hom (dualGlueVal φ hφ hV x)
      = dualGlueVal φ hφ (hWV.trans hV) ((F.map (homOfLE hWV).op).hom x) := by
  refine dualGlue_ext (hWV.trans hV) _ _ (fun i => ?_)
  have hWi : W ⊓ U i ≤ V ⊓ U i := inf_le_inf_right _ hWV
  refine Eq.trans (resTrans (inf_le_left : W ⊓ U i ≤ W) hWV _) ?_
  refine Eq.trans (resTrans hWi (inf_le_left : V ⊓ U i ≤ V) _).symm ?_
  rw [dualGlueVal_res, dualGlueVal_res, dualLocVal_res φ x i hWi, dualLocVal_apply]
  exact (congrArg _ (modResTrans (inf_le_left : W ⊓ U i ≤ W) hWV x)).symm

/-! ## ★出典の紐付け(`.src`) -/

def dualGlueVal_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——貼り合わせた値が加法的・線型・自然であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
