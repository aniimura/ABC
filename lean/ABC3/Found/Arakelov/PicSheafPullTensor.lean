import ABC3.Found.Arakelov.PicDeltaIso
import ABC3.Found.Arakelov.PicSchemeDelta

/-!
# Arakelov (B1) 第 41 ブロック —— **層の段で引き戻しはテンソル積を保つ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★前層の結果を層へ運ぶ

第 40 ブロックで**前層の段**の `δ` が同型になった。★層の段へ運ぶ道:

    f^*_sh (tensorModules L M)
      ≅ 層化 (f^*_pre (L.val ⊗ M.val))              ★第 22 + `sheafifyValIso`
      ≅ 層化 (f^*_pre L.val ⊗ f^*_pre M.val)        ★★第 40(δ が同型)
      ≅ 層化 ((層化 (f^*_pre L.val)).val ⊗ (層化 (f^*_pre M.val)).val)
                                                    ★第 11(層化はテンソルを変えない)
      ≅ tensorModules (f^*_sh L) (f^*_sh M)         ★`pullbackValIso` を逆に

★★★第 11 の橋にだけ**局所階数 1** の仮定が要る。可逆層では満たされる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `pullbackValIso` | ★層の引き戻しは「前層で引き戻して層化」に等しい |
| `pullbackSheafTensorIso` | ★★★★★★**`f^*(L ⊗ M) ≅ f^*L ⊗ f^*M`(層の段)** |

★★★★これで `PicardData.pullback_mul` が書ける。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-! ## ★層の引き戻しの前層表示 -/

/-- ★**層の引き戻しは「前層で引き戻して層化」に等しい**。

★★第 22 ブロックの `sheafifyPullbackApp` に `sheafifyValIso`(層は層化で変わらない)
を継いだものである。 -/
noncomputable def pullbackValIso (L : Y.Modules) :
    (Scheme.Modules.pullback f).obj L
      ≅ (sheafifyFunctor X).obj ((pullbackPre f).obj L.val) :=
  (Scheme.Modules.pullback f).mapIso (sheafifyValIso L).symm
    ≪≫ sheafifyPullbackApp f L.val

/-! ## ★★★★★★層の段でテンソル積を保つ -/

/-- ★★★★★★**層の段で引き戻しはテンソル積を保つ**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★仮定は「引き戻して層化したものが局所階数 1」だけである
——可逆層ではこれが成り立つ。

★機構は 4 本の橋(第 22・第 40・第 11・第 22 逆)。 -/
noncomputable def pullbackSheafTensorIso (L M : Y.Modules)
    (hM : IsLocallyRankOne X ((pullbackPre f).obj M.val))
    (hLs : IsLocallyRankOne X
      ((sheafifyFunctor X).obj ((pullbackPre f).obj L.val)).val) :
    (Scheme.Modules.pullback f).obj (tensorModules L M)
      ≅ tensorModules ((Scheme.Modules.pullback f).obj L)
          ((Scheme.Modules.pullback f).obj M) :=
  sheafifyPullbackApp f (L.val ⊗ M.val)
    ≪≫ (sheafifyFunctor X).mapIso (pullbackTensorIso f L.val M.val)
    ≪≫ sheafifyTensorRight X ((pullbackPre f).obj L.val) ((pullbackPre f).obj M.val) hM
    ≪≫ sheafifyTensorLeft X
        ((sheafifyFunctor X).obj ((pullbackPre f).obj L.val)).val
        ((pullbackPre f).obj M.val) hLs
    ≪≫ (sheafifyFunctor X).mapIso (MonoidalCategory.tensorIso
        ((SheafOfModules.forget X.ringCatSheaf).mapIso (pullbackValIso f L)).symm
        ((SheafOfModules.forget X.ringCatSheaf).mapIso (pullbackValIso f M)).symm)

/-! ## ★出典の紐付け(`.src`) -/

def pullbackSheafTensorIso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——層の段で引き戻しがテンソル積を保つこと)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
