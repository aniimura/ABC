/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55BiratPf

/-!
# [FrdI] Proposition 5.5, (ii) の birat の側を**関手**として組む

原文 (FrdI p.104):
> (ii) There is a natural equivalence of categories [compatible with the functors to the

★★`Prop55BiratPf.lean` は **射の全単射**

  `Hom_{(𝒞^birat)^pf}(A,B) ≃ Hom_{(𝒞^pf)^birat}(A,B)`

(`biratPfHomEquiv`)と、**恒等射の保存**(`biratPfHom_id`)・
**合成との両立**(`biratPfHom_comp`)まで作ってある。
★本ファイルはそれを**関手**に束ね、**充満忠実**であることを言う。

## ★★対象について(記録)

`(𝒞^birat)^pf` の**根 1 の部分**(`PfCat (biratPre P G) F'`)から
`(𝒞^pf)^birat`(`BiratCat (pfRootPre P F) Gpf`、対象は対 `(A,n)`)への関手になる。
★像は根 1 の対象だけなので**本質的全射ではない** ——
原文の「equivalence of categories」に届かせるには
根 `n` を一般にした版(`scaleRootEquiv` で根を揃える段)が要る。
★un-tr の側(`Prop55UntrFun.lean`)では両側の対象が一致していたので
そこは無料だったが、birat の側はそうではない。

## ★★★★測って分かったこと(2026-08-25)—— 源を広げるしかない

★**この関手(源が根 1)は原理的に本質的全射にならない**:
`⟨A,n⟩ ≅ ⟨B,1⟩` は Frobenius 次数が不変量なので起きない。
`pfKappa A n : ⟨A,n⟩ ⟶ ⟨A,1⟩` は次数 `n` であり、
★**birat 化は次数 `n` の射を可逆にしない**(可逆にするのは底同型・単元の側)。
したがって**源を `PfRootObj (biratPre P G) F'` 全体へ広げる**しかない。

★★段取りは `Prop55BiratPf.lean` の ★(1425 行あたり)にある 4 段そのもの:

```
Hom_{(𝒞^birat)^pf}((A,n),(B,m)) = Hom^pf_{𝒞^birat}(A^{(m)}, B^{(n)})   -- 定義
    ≃ Hom_{(𝒞^pf)^birat}(⟨A^{(m)},1⟩, ⟨B^{(n)},1⟩)   -- biratPfHomEquiv(★済)
    ≃ Hom_{(𝒞^pf)^birat}(⟨A^{(m)},k⟩, ⟨B^{(n)},k⟩)   -- ★★Σ_k を birat へ降ろす(残)
    ≃ Hom_{(𝒞^pf)^birat}(⟨A,n⟩, ⟨B,m⟩)               -- pfRoot_exists_iso_root(★済)
```

* 2 段目 …… `biratPfHomEquiv`(本ファイルの `biratPfFunctor`)—— **済**
* 4 段目 …… `⟨A,n⟩ ≅ ⟨A^{(m)}, m·n⟩`(在庫 `pfRoot_exists_iso_root`)を
  `toBiratCat` で送るだけ —— **済**(関手は同型を同型へ送る)
* ★★3 段目だけが残る —— `scaleRootEquiv`(`𝒞^pf ≌ 𝒞^pf`、**在庫**)を
  **birat 化へ降ろす**こと。`Σ_k` は `obj` を変えず `root` を `k` 倍し、
  射は `rtRootIso` の共役なので**次数と `Div` を保つ** ——
  したがって Frobenioid の自己同値であり、birat 化は関手的なはずである。
  ★要るのは「Frobenioid の自己同値が birat 化へ降りる」1 本
  (在庫 `psiBiratCor`(`Cor411Birat.lean`)が同じ形をしている)。

★あわせて必要: `F'` の根の選び方を `F` のものと揃えること
(`rtObj (biratPre P G) F' A d = rtObj P F A d`)。
★★これは `Corollary 5.4` の継ぎ目で基準切断を揃えたのと**同じ型の手当て**である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section BiratPfFun

variable {D : Type u} [Category.{v} D] [IsConnected D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}
  {G : Frobenioid P}

/-- ★★★★★★**[FrdI] Proposition 5.5, (ii)(birat の側)の関手** ——

  `(𝒞^birat)^pf` の根 1 の部分 ⥤ `(𝒞^pf)^birat`

★対象は `A ↦ ⟨A, 1⟩`、射は `biratPfHom`。
★恒等・合成は `biratPfHom_id` / `biratPfHom_comp` そのもの。 -/
noncomputable def biratPfFunctor (hfi : IsOfFrobeniusIsotropicType P)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) :
    PfCat (biratPre P G) F' ⥤ BiratCat (pfRootPre P F) Gpf where
  obj A := (⟨biratDown P G A, 1⟩ : PfRootObj P F)
  map {A B} f := biratPfHom hfi Gpf F' (biratDown P G A) (biratDown P G B) f
  map_id A := biratPfHom_id hfi Gpf F' (biratDown P G A)
  map_comp {A B E} f g := biratPfHom_comp hfi Gpf F'
    (biratDown P G A) (biratDown P G B) (biratDown P G E) f g

/-- ★★★★★★**[FrdI] Proposition 5.5, (ii)(birat の側)** —— 関手は**忠実**。 -/
theorem biratPfFunctor_faithful (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) :
    (biratPfFunctor hfi Gpf F').Faithful where
  map_injective {A B} {f g} h :=
    (biratPfHom_bijective hfi hiso Gpf F'
      (biratDown P G A) (biratDown P G B)).1 h

/-- ★★★★★★**[FrdI] Proposition 5.5, (ii)(birat の側)** —— 関手は**充満**。 -/
theorem biratPfFunctor_full (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) :
    (biratPfFunctor hfi Gpf F').Full where
  map_surjective {A B} h :=
    (biratPfHom_bijective hfi hiso Gpf F'
      (biratDown P G A) (biratDown P G B)).2 h

/-- ★★★★★**[FrdI] Proposition 5.5, (ii)(birat の側)** ——
関手が根 1 の対象のあいだで与える全単射は `biratPfHomEquiv` そのもの。 -/
theorem biratPfFunctor_map_eq (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G))
    {A B : PfCat (biratPre P G) F'} (f : A ⟶ B) :
    (biratPfFunctor hfi Gpf F').map f
      = biratPfHomEquiv hfi hiso Gpf F' (biratDown P G A) (biratDown P G B) f :=
  rfl

end BiratPfFun

/-! ### ★出典の紐付け -/

/-- ★★★★★★locator —— `Proposition 5.5, (ii)` の birat の側を関手として組んだもの
(充満忠実まで。本質的全射は根を一般にする段が残る)。 -/
def biratPfFunctor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Proposition 5.5, (ii) — (𝒞^birat)^pf ⥤ (𝒞^pf)^birat(充満忠実)",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
