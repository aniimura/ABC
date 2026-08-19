import ABC3.Found.Arakelov.ArcGenSection

/-!
# Arakelov (C3) 第 265 ブロック —— **`Over V` 側への型付き橋**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★§9-297 の壁を**型付き恒等関数**で削る

`restrict F V.ι`(mathlib)と `Over V` 側(我々)は台が同じ(第 255 で `rfl` を実測)
だが、**加群構造の綴りが違う**ので元素レベルの証明が書けなかった(§9-297)。

★★第 164 ブロックが既に答えを書いていた:

> **instance を足すより、型付き恒等関数で橋を架ける方が安全**
> ——instance は既存の経路と競合しうるが、恒等関数は競合しない。

★★★本ブロックはその方針を `Over V` 側に適用したものである。

| 橋 | 何を移すか |
|---|---|
| `coefOverV` | ★係数 `Γ(X, ι''ᵁW)` → `Over V` の係数環 |
| `secOverV` | ★切断 `Γ(restrict F ι, W)` → `Over V` の前層の値 |
| `unitValAt` | ★単位対象の値 → `Γ(X, ι''ᵁW)` |
| `unitValAt_smul` | ★★★スカラー倍が積に移る(`rfl`) |

★★★★これで `e.hom.app` の線形性(`map_smul`)が**書けるようになった**
——`coefOverV` と `secOverV` を噛ませればよい。

## ★残り

`bridgeApp_smul`(自然性・スカラー両立)の**算術の段**にもう 1–2 個の橋が要る。
★見通しは立った——**壁ではなく、橋を数える作業**になった。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{0}} (V : X.Opens) (F : X.Modules)

/-- ★`Over V` の係数環への橋。 -/
def coefOverV (W : V.toScheme.Opens) (c : ((X.presheaf.obj (op (V.ι ''ᵁ W))) : Type)) :
    (((((Over.forget V).op ⋙ X.presheaf) ⋙ forget₂ CommRingCat.{0} RingCat.{0}).obj
      (op (overObj V W))) : Type) := c

/-- ★切断の橋。 -/
def secOverV (W : V.toScheme.Opens)
    (s : ((Scheme.Modules.restrict F V.ι).val.obj (op W) : Type)) :
    (((((restrictPresheafFunctor X V).obj F.val).obj (op (overObj V W)))) : Type) := s

/-- ★単位対象の値の橋。 -/
def unitValAt (W : V.toScheme.Opens)
    (x : (((𝟙_ (PresheafModulesOn X V)).obj (op (overObj V W))) : Type)) :
    ((X.presheaf.obj (op (V.ι ''ᵁ W))) : Type) := x

/-- ★★橋はスカラー倍を積に移す。 -/
theorem unitValAt_smul (W : V.toScheme.Opens)
    (c : ((X.presheaf.obj (op (V.ι ''ᵁ W))) : Type))
    (x : (((𝟙_ (PresheafModulesOn X V)).obj (op (overObj V W))) : Type)) :
    unitValAt V W (coefOverV V W c • x) = c * unitValAt V W x := rfl


/-! ## ★出典の紐付け(`.src`) -/

def coefOverV.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——Over V 側への型付き橋)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
