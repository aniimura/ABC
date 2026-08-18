import ABC3.Found.Arakelov.PicEvalHom

/-!
# Arakelov (B2) 第 176 ブロック —— **双対の局所値は貼り合わせ可能**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★双対が層であることの第 1 歩

§9-200 で **A の道**(双対が層であることを示す)を採ると決めた。その筋は

    1. `V ≤ ⨆ Uᵢ` なら `V` は `V ⊓ Uᵢ` で覆われる
    2. 各 `Uᵢ` の切断 `φᵢ` は `V ⊓ Uᵢ` の上で値を与える
    3. ★その値たちは**重なりで一致する**
    4. `X.sheaf` で貼り合わせて `𝒪(V)` の元を得る

★本ブロックは **1–3** を取る。

## ★★鍵は「切断の自然性」だけ

双対の切断 `ψ : F|_A ⟶ 𝟙_` は `Over A` の上の射だから、
**制限と可換**である(`dual_app_res`)。これは `ψ.naturality` そのもので、

    congrArg (fun m => m.hom x) h.symm

の 1 行で出る。★★重なりでの一致も、仮定 `hφ i j`(`Uᵢ ⊓ Uⱼ` での一致)を
`W = (V ⊓ Uᵢ) ⊓ (V ⊓ Uⱼ)` で評価するだけで出る——
`Over` の対象は `left` が同じなら**射が一意**(`Opens` の `Hom` が Prop)なので、
`Over.map` で運んだ対象と直接作った対象は**定義から等しい**からである。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `le_iSup_inf` | ★`V ≤ ⨆ Uᵢ` なら `V ≤ ⨆ (V ⊓ Uᵢ)` |
| `infOver` | ★`V ⊓ Uᵢ` を `Over (Uᵢ)` の対象として見る |
| `dual_app_res` | ★★**双対の切断は制限と可換** |
| `dualLocVal` | ★各 `i` の上での値 |
| `dualLocVal_res` | ★★局所値の制限は「もっと小さい所での値」 |
| `dualLocVal_compat` | ★★★**局所値は貼り合わせ可能** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

/-- ★**`V ≤ ⨆ Uᵢ` なら `V` は `V ⊓ Uᵢ` たちで覆われる**。 -/
theorem le_iSup_inf {X : Scheme.{u}} {ι : Type u} (U : ι → X.Opens) (V : X.Opens)
    (hV : V ≤ ⨆ i, U i) : V ≤ ⨆ i, (V ⊓ U i) := by
  intro x hx
  obtain ⟨i, hi⟩ := Opens.mem_iSup.1 (hV hx)
  exact Opens.mem_iSup.2 ⟨i, ⟨hx, hi⟩⟩

/-- ★**`V ⊓ Uᵢ` を `Over (Uᵢ)` の対象として見る**。 -/
abbrev infOver {X : Scheme.{u}} {ι : Type u} (U : ι → X.Opens) (V : X.Opens) (i : ι) :
    Over (U i) :=
  Over.mk (homOfLE (inf_le_right : V ⊓ U i ≤ U i))

/-- ★★**双対の切断は制限と可換である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★これは `ψ.naturality` そのものである——双対の切断は `Over A` の上の**射**だから。 -/
theorem dual_app_res {X : Scheme.{u}} {F : X.PresheafOfModules} {A : X.Opens}
    (ψ : ((dualPresheaf F).obj (op A) : Type u))
    {Z W : Over A} (f : W ⟶ Z) (x : (F.obj (op Z.left) : Type u)) :
    (X.presheaf.map (homOfLE (leOfHom f.left)).op).hom (((ψ.app (op Z))).hom x)
      = ((ψ.app (op W))).hom ((F.map (homOfLE (leOfHom f.left)).op).hom x) :=
  congrArg (fun (m : _ ⟶ _) => (ModuleCat.Hom.hom m) x) (ψ.naturality f.op).symm

variable {X : Scheme.{u}} {F : X.PresheafOfModules} {ι : Type u} {U : ι → X.Opens}
  (φ : ∀ i, ((dualPresheaf F).obj (op (U i)) : Type u))

/-- ★**各 `i` の上での値**——`x` を `V ⊓ Uᵢ` へ制限して `φᵢ` を当てる。 -/
noncomputable def dualLocVal (V : X.Opens) (x : (F.obj (op V) : Type u)) (i : ι) :
    (Γ(X, V ⊓ U i) : Type u) :=
  ((φ i).app (op (infOver U V i))).hom
    (F.map (homOfLE (inf_le_left : V ⊓ U i ≤ V)).op x)

/-- ★★**局所値の制限は「もっと小さい所での値」である**。 -/
theorem dualLocVal_res {V : X.Opens} (x : (F.obj (op V) : Type u)) (i : ι)
    {W : X.Opens} (hW : W ≤ V ⊓ U i) :
    (X.presheaf.map (homOfLE hW).op).hom (dualLocVal φ V x i)
      = ((φ i).app (op (Over.mk (homOfLE (hW.trans inf_le_right) : W ⟶ U i)))).hom
          ((F.map (homOfLE (hW.trans inf_le_left)).op).hom x) := by
  refine (dual_app_res (φ i)
    (Over.homMk (homOfLE hW) : Over.mk (homOfLE (hW.trans inf_le_right) : W ⟶ U i) ⟶ infOver U V i)
    ((F.map (homOfLE (inf_le_left : V ⊓ U i ≤ V)).op).hom x)).trans ?_
  congr 1
  exact congrArg (fun (m : _ ⟶ _) => (ModuleCat.Hom.hom m) x) (F.map_comp _ _).symm

/-- ★★★**局所値は貼り合わせ可能である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★仮定 `hφ i j`(`Uᵢ ⊓ Uⱼ` での一致)を `W = (V ⊓ Uᵢ) ⊓ (V ⊓ Uⱼ)` で評価するだけ。
`Over` の対象は `left` が同じなら射が一意だから、`Over.map` で運んだ対象と
直接作った対象は**定義から等しい**。 -/
theorem dualLocVal_compat (hφ : TopCat.Presheaf.IsCompatible (dualPresheaf F).presheaf U φ)
    (V : X.Opens) (x : (F.obj (op V) : Type u)) :
    TopCat.Presheaf.IsCompatible X.presheaf (fun i => V ⊓ U i) (dualLocVal φ V x) := by
  intro i j
  have hWi : (V ⊓ U i) ⊓ (V ⊓ U j) ≤ V ⊓ U i := inf_le_left
  have hWj : (V ⊓ U i) ⊓ (V ⊓ U j) ≤ V ⊓ U j := inf_le_right
  have hWij : (V ⊓ U i) ⊓ (V ⊓ U j) ≤ U i ⊓ U j :=
    le_inf (hWi.trans inf_le_right) (hWj.trans inf_le_right)
  refine (dualLocVal_res φ x i hWi).trans (Eq.trans ?_ (dualLocVal_res φ x j hWj).symm)
  exact congrArg
    (fun (ψ : ((dualPresheaf F).obj (op (U i ⊓ U j)) : Type u)) =>
      ((ψ.app (op (Over.mk (homOfLE hWij)))).hom
        ((F.map (homOfLE (hWi.trans inf_le_left)).op).hom x))) (hφ i j)

/-! ## ★出典の紐付け(`.src`) -/

def dualLocVal_compat.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——双対の局所値が貼り合わせ可能であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
