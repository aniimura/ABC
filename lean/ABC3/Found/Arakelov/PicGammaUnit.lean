import ABC3.Found.Arakelov.PicGammaSpecIso

/-!
# Arakelov (B1) 第 68 ブロック —— **比較射の後段(層化の unit)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★比較射の 2 段のうち後段

第 66 ブロックで比較射の形が

    F.val(⊤) ⊗_R G.val(⊤)  →  (F.val ⊗ G.val)(⊤)  →  (層化 (F.val ⊗ G.val))(⊤)
    ────────── 前段(第 67 の係数環の読み替え) ──   ── 後段(本ブロック) ──

と分かった。★**後段は層化の unit を `⊤` で評価したもの**である。

## ★★機構

第 11 ブロックで作った `sheafifyUnit`(層化の単位)を `op ⊤` で評価するだけである。
★★★局所単射・局所全射であることも第 11 で示してある。

## ★★本ブロックで取れるもの

| 定義 | 内容 |
|---|---|
| `gammaSheafifyUnit` | ★★**層化の unit の `⊤` 成分** |

## ★★★残り(B1)

| # | 主張 |
|---|---|
| 1 | 比較射の前段(`M ⊗_R N → M ⊗_{𝒪(⊤)} N`、★第 67 で係数環は同型と確認) |
| 2 | `tilde` がテンソルを保つ(★★第 64・65 + 第 29) |
| 3 | 可逆層は `tilde` の本質像に入る |
| 4 | `equivPicRing` / 5 `PicardData` の組み立て |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable (R : CommRingCat.{u}) (F G : (Spec R).Modules)

/-- ★★**層化の unit の `⊤` 成分**——比較射の後段である。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★第 11 ブロックの `sheafifyUnit` を `op ⊤` で評価しただけである。 -/
noncomputable abbrev gammaSheafifyUnit :
    (F.val ⊗ G.val).obj (op ⊤)
      ⟶ ((PresheafOfModules.sheafification (R := (Spec R).ringCatSheaf)
          (𝟙 (Spec R).ringCatSheaf.obj)).obj (F.val ⊗ G.val)).val.obj (op ⊤) :=
  (sheafifyUnit (Spec R) (F.val ⊗ G.val)).app (op ⊤)

/-! ## ★出典の紐付け(`.src`) -/

def gammaSheafifyUnit.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——比較射の後段(層化の unit))",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
