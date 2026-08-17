import ABC3.Found.Arakelov.PicSheafGroup
import Mathlib.AlgebraicGeometry.Modules.Sheaf
import Mathlib.Algebra.Category.ModuleCat.Sheaf.PullbackFree

/-!
# Arakelov (B1) 第 21 ブロック —— **層の段の引き戻し**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★測ったら在庫があった

第 18–20 ブロックでは**前層**の段で `pullback` の oplax monoidal 性を作った。
★★しかし `PicardData` が要るのは**層**の段である。測ってみると:

| 在庫 | 場所 |
|---|---|
| `Scheme.Modules.pullback f : Y.Modules ⥤ X.Modules` | `AlgebraicGeometry/Modules/Sheaf.lean` |
| `Scheme.Modules.pullbackId : pullback (𝟙 X) ≅ 𝟭 _` | 同上 |
| `Scheme.Modules.pullbackComp : pullback g ⋙ pullback f ≅ pullback (f ≫ g)` | 同上 |
| `SheafOfModules.pullbackObjUnitToUnit` + `IsIso`(`F.Final` のとき) | `Sheaf/PullbackFree.lean` |
| `SheafOfModules.pullbackObjFreeIso : pullback (free I) ≅ free I` | 同上 |

★★★**`pullback_id` と `pullback_comp` は既に在る。**
`PicSheaf` は同型による商なので、**自然同型がそのまま等式になる**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `opensMapFinal` | ★★`Opens.map f.base` は Final(6 行) |
| `pullbackUnitIso` | ★★★★**構造層の引き戻しは構造層** |
| `pullbackId_app` / `pullbackComp_app` | ★対象ごとの同型 |

★★★`Final` が要るのは「切断が押し出しで保たれる」ためである。
`X.Opens` は前順序であり `⊤` を持つので、
**すべての構造射が `⊤` へ向かう**——それだけで連結になる。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X Y Z : Scheme.{u}} (f : X ⟶ Y) (g : Y ⟶ Z)

/-! ## ★★底の関手は Final -/

/-- ★★★**開集合の逆像関手は Final である**。

★★`StructuredArrow U (Opens.map f.base)` は「`U ⊆ f⁻¹V` なる `V`」の圏である。
★★★`V = ⊤` は常にこの条件を満たし、**すべての対象がそこへ射を持つ**。
前順序なので三角形は自動で可換であり、連結性が出る。

★これが `SheafOfModules.pullbackObjUnitToUnit` の同型性の前提である。 -/
instance opensMapFinal : (Opens.map f.base).Final where
  out U := by
    letI t : StructuredArrow U (Opens.map f.base) :=
      StructuredArrow.mk (Y := (⊤ : Y.Opens)) (homOfLE le_top)
    haveI : Nonempty (StructuredArrow U (Opens.map f.base)) := ⟨t⟩
    refine zigzag_isConnected fun j₁ j₂ => ?_
    have h₁ : j₁ ⟶ t := StructuredArrow.homMk (homOfLE le_top) (Subsingleton.elim _ _)
    have h₂ : j₂ ⟶ t := StructuredArrow.homMk (homOfLE le_top) (Subsingleton.elim _ _)
    exact (Zigzag.of_hom h₁).trans (Zigzag.of_inv h₂)

/-! ## ★★★★構造層の引き戻し -/

/-- ★★★★★**構造層の引き戻しは構造層である** `f^* 𝒪_Y ≅ 𝒪_X`。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `Pic` の**単位元が引き戻しで保たれる**ことを与える。
★機構は mathlib の `pullbackObjUnitToUnit` であり、
その同型性の前提が上の `opensMapFinal` である。 -/
noncomputable def pullbackUnitIsoAux :
    (SheafOfModules.pullback f.toRingCatSheafHom).obj (SheafOfModules.unit Y.ringCatSheaf)
      ≅ SheafOfModules.unit X.ringCatSheaf :=
  haveI : IsIso (SheafOfModules.pullbackObjUnitToUnit f.toRingCatSheafHom) := inferInstance
  asIso (SheafOfModules.pullbackObjUnitToUnit f.toRingCatSheafHom)

/-- ★★★★★**構造層の引き戻しは構造層である**——`Scheme.Modules` の語彙で。

★★★**実装の罠(2026-08-18 実測)**: これを 1 段で書くと
`IsIso` のインスタンス探索が**落ちる**。
`Scheme.Modules.pullback` は `def` なので、期待型

    (Scheme.Modules.pullback f).obj (unitModules Y) ≅ unitModules X

は mathlib 側の

    (SheafOfModules.pullback φ).obj (unit S) ≅ unit R

と**構文的に一致せず**、`φ` の型にメタ変数が残ったまま
インスタンス探索に入るためである。

★★★**塞ぎ方: まず自然な型で書き(`pullbackUnitIsoAux`)、
次にそれを別の `def` で移す。**
2 段目では 1 段目が完全に確定しているので、defeq 判定だけで済む。 -/
noncomputable def pullbackUnitIso :
    (Scheme.Modules.pullback f).obj (unitModules Y) ≅ unitModules X :=
  pullbackUnitIsoAux f

/-! ## ★★関手性(対象ごと) -/

/-- ★**恒等射の引き戻しは何もしない**。 -/
noncomputable def pullbackIdApp (M : X.Modules) :
    (Scheme.Modules.pullback (𝟙 X)).obj M ≅ M :=
  (Scheme.Modules.pullbackId X).app M

/-- ★★**合成の引き戻しは引き戻しの合成**。 -/
noncomputable def pullbackCompApp (M : Z.Modules) :
    (Scheme.Modules.pullback (f ≫ g)).obj M
      ≅ (Scheme.Modules.pullback f).obj ((Scheme.Modules.pullback g).obj M) :=
  ((Scheme.Modules.pullbackComp f g).app M).symm

/-! ## ★出典の紐付け(`.src`) -/

def opensMapFinal.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——開集合の逆像関手が Final であること)",
    sectionId := "genell-def-1-1-i" }

def pullbackUnitIso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——構造層の引き戻しが構造層になること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
