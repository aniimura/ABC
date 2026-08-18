import ABC3.Found.Arakelov.PicTrivialSheaf

/-!
# Arakelov (B1) 第 112 ブロック —— **局所化は全単射を保つ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★生成元が生成元のまま残ることの中身

`M_g ≅ R_g`(生成元 `s` による乗法)を `D(t·g)` へ制限すると
`M_{t·g} ≅ R_{t·g}`(`s` の制限による乗法)になる。

★★これは「**局所化は全単射を保つ**」だけである
——mathlib に `IsLocalizedModule.map_injective` / `map_surjective` がある(実測)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `bijective_localizedMap` | ★★★**局所化は全単射を保つ** |

## ★★★これで局所自明性が閉じる

    M_g ≅ R_g(第 76)→ 局所化(本ブロック)→ M_h ≅ R_h(生成元の像で)
      → 第 103(生成元の乗法は全単射)→ 第 111(層加群版の同型)
-/

universe u

namespace ABC3.Found.Arakelov

open TensorProduct

variable {R : Type u} [CommRing R] (S : Submonoid R)
  {M N M' N' : Type u} [AddCommGroup M] [Module R M] [AddCommGroup N] [Module R N]
  [AddCommGroup M'] [Module R M'] [AddCommGroup N'] [Module R N']
  (f : M →ₗ[R] M') (g : N →ₗ[R] N')
  [IsLocalizedModule S f] [IsLocalizedModule S g]

/-- ★★★**局所化は全単射を保つ**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★これで「生成元は制限しても生成元」が出る。 -/
theorem bijective_localizedMap (h : M →ₗ[R] N) (hb : Function.Bijective h) :
    Function.Bijective (IsLocalizedModule.map S f g h) :=
  ⟨IsLocalizedModule.map_injective S f g h hb.1,
    IsLocalizedModule.map_surjective S f g h hb.2⟩

/-! ## ★出典の紐付け(`.src`) -/

def bijective_localizedMap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——局所化は全単射を保つ)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
