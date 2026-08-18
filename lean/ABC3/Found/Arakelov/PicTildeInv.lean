import ABC3.Found.Arakelov.PicTildeTensorIso

/-!
# Arakelov (B1) 第 92 ブロック —— **`tilde` は可逆加群を可逆層へ送る(逆の側)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★第 91 ブロックの直接の配当

    tilde (Mᵛ) ⊗ tilde M  ≅  tilde (Mᵛ ⊗ M)  ≅  tilde R  ≅  𝒪_{Spec R}

★★mathlib の `Module.Invertible` は「`Mᵛ ⊗ M → R` が全単射」と定義されている
——★★★したがって**双対がそのまま逆になる**。

## ★★機構 —— 3 つの同型の合成

| 段 | 出典 |
|---|---|
| `tilde Mᵛ ⊗ tilde M ≅ tilde (Mᵛ ⊗ M)` | ★★★第 91 ブロック |
| `Mᵛ ⊗ M ≅ R` | ★mathlib `Module.Invertible.linearEquiv` |
| `tilde R ≅ 𝒪` | ★mathlib `tildeSelf`(**`.refl`**) |

★**一発で通った。**

## ★★本ブロックで取れるもの

| 定義 | 内容 |
|---|---|
| `tildeDualTensorIso` | ★★★★**`tilde Mᵛ ⊗ tilde M ≅ 𝒪`** |

## ★★★残り —— 局所自明性

`InvSheaf` の残るフィールドは `trivial` / `invTrivial`(局所自明性)である。
★第 76 ブロック(可逆加群は基本開集合の上で自由)から

    D(r) の上で M_r ≅ R_r

は出ている。★★あとは `(tilde M)|_{D r} ≅ 𝟙_` に運ぶ段が要る。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TensorProduct

variable (R : CommRingCat.{u}) (M : ModuleCat.{u} (R : Type u))
  [Module.Invertible (R : Type u) M]

/-- ★★★★**`tilde Mᵛ ⊗ tilde M ≅ 𝒪_{Spec R}`**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★可逆加群の双対が、そのまま可逆層の逆を与える。 -/
noncomputable def tildeDualTensorIso :
    tensorModules (tilde (ModuleCat.of (R : Type u) (Module.Dual (R : Type u) M))) (tilde M)
      ≅ unitModules (Spec R) :=
  tildeTensorIso R _ M
    ≪≫ (tilde.functor R).mapIso
      (LinearEquiv.toModuleIso (Module.Invertible.linearEquiv (R : Type u) M))
    ≪≫ tildeSelf

/-! ## ★出典の紐付け(`.src`) -/

def tildeDualTensorIso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——tilde は可逆加群の双対を逆へ送る)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
