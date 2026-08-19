import ABC3.Found.Arakelov.ArcUnitComp

/-!
# Arakelov (C3) 第 279 ブロック —— **`Γ`・押し出し・`pushforwardComp` は元のレベルで素通し**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★「型は違うが元は同じ」——第 255 ブロックの型を再利用する

第 278 ブロックの `unit_conj_compat` を `⊤` の切断に落とすには、3 つの橋が要る。
★★どれも **`ModuleCat` の係数環が違うので射としては `rfl` にならない**が、
**`.hom` を噛ませて元のレベルに落とすと 3 つとも `rfl`** である。

| 橋 | 射として | 元として |
|---|---|---|
| `moduleSpecΓFunctor.map ψ` と `ψ.val.app ⊤` | ✗(係数環が違う) | ★`rfl` |
| `(pushforward p).map φ` の `⊤` と `φ` の `⊤` | ✗(底空間が違う) | ★`rfl` |
| `pushforwardComp` の `⊤` | ✗ | ★★`rfl`(恒等) |

★★★これは §9-309 で得た型——**担い手が `rfl` で一致するなら元のレベルで書け**——
そのものである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `gammaMap_apply` | ★`Γ` は `⊤` での評価(元のレベル) |
| `pushforwardMap_apply` | ★押し出しは `⊤` で素通し |
| `pushforwardComp_apply` | ★★`pushforwardComp` は `⊤` で恒等 |
| `gamma_arcEval_naturality` | ★★★★**`Γ` レベルでの `arcEval` の自然性** |

★残るのは `restrictFunctorIsoPullback` と `j` の単位の関係のみ(§9-323)。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite

variable {U X : Scheme.{0}} (p : Spec (CommRingCat.of ℂ) ⟶ X)

/-- ★**`Γ` は `⊤` での評価である**(元のレベルでは `rfl`)。 -/
theorem gammaMap_apply {M N : (Spec (CommRingCat.of ℂ)).Modules} (ψ : M ⟶ N)
    (v : (M.val.obj (op ⊤) : Type)) :
    ((moduleSpecΓFunctor (R := CommRingCat.of ℂ)).map ψ).hom v = (ψ.val.app (op ⊤)).hom v := rfl

/-- ★**押し出しは `⊤` の切断では素通しである**(`p ⁻¹ᵁ ⊤ = ⊤` だから)。 -/
theorem pushforwardMap_apply {M N : (Spec (CommRingCat.of ℂ)).Modules} (φ : M ⟶ N)
    (v : (((Scheme.Modules.pushforward p).obj M).val.obj (op ⊤) : Type)) :
    (((Scheme.Modules.pushforward p).map φ).val.app (op ⊤)).hom v = (φ.val.app (op ⊤)).hom v := rfl

/-- ★★**`pushforwardComp` は `⊤` の切断では恒等である**。 -/
theorem pushforwardComp_apply (q : Spec (CommRingCat.of ℂ) ⟶ U) (j : U ⟶ X)
    (M : (Spec (CommRingCat.of ℂ)).Modules)
    (v : (((Scheme.Modules.pushforward q ⋙ Scheme.Modules.pushforward j).obj M).val.obj
      (op ⊤) : Type)) :
    (((Scheme.Modules.pushforwardComp q j).hom.app M).val.app (op ⊤)).hom v = v := rfl

/-- ★★★★**`Γ` レベルでの `arcEval` の自然性**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★第 250 ブロックの `arcEval_naturality` を、ファイバー(`Γ`)の言葉に翻訳したもの。
★これで「切断を先に動かす」と「ファイバーで動かす」が入れ替えられる。 -/
theorem gamma_arcEval_naturality {A B : X.Modules} (θ : A ⟶ B)
    (t : (A.val.obj (op ⊤) : Type)) :
    ((moduleSpecΓFunctor (R := CommRingCat.of ℂ)).map
        ((Scheme.Modules.pullback p).map θ)).hom (arcEval p A t)
      = arcEval p B ((θ.val.app (op ⊤)).hom t) := by
  rw [gammaMap_apply, ← pushforwardMap_apply p ((Scheme.Modules.pullback p).map θ)]
  exact (arcEval_naturality p θ t).symm

/-! ## ★出典の紐付け(`.src`) -/

def gamma_arcEval_naturality.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C3——Γ レベルでの評価の自然性)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
