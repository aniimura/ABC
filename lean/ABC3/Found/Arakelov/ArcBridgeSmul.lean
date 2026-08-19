import ABC3.Found.Arakelov.ArcOverCoef

/-!
# Arakelov (C3) 第 266 ブロック —— ★★★★★★**§9-297 の壁を越えた**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★「元素レベルでは書けない」は**誤りだった**

§9-297 で「台が同じでも加群構造が 2 つあれば元素レベルの証明は書けない、
関手レベルで 5–15 ブロック」と測った。★★**越えられた。**

★★★**型付き恒等関数の橋を 5 本架けるだけ**だった:

| 橋 | 何を移すか |
|---|---|
| `coefOverV` | 係数 → `Over V` の係数環 |
| `secOverV` | 切断 → `Over V` の前層の値 |
| `unitValAt` | 単位対象の値 → `Γ(X, ι''ᵁW)` |
| `unitValV` | `V.toScheme` 側の単位対象の値 |
| `coefV` | `V.toScheme` 側の係数 |

★★★★橋はすべて `:= x`(恒等関数)であり、法則はすべて `rfl` である。

## ★★★教訓の訂正 —— 「構造が 2 つある」は**行き止まりではない**

§9-297 では「綴りが 2 つなら寄せられるが、構造が 2 つなら書けない」と結論した。
★★正しくは:**構造が 2 つでも、台が同じなら型付き恒等関数で橋を架けられる**。
★★★見分け方は「台が `rfl` で一致するか」であり、
第 255 でそれを実測していたのに、活かせていなかった。

★★★★**測定してあった事実を、壁の判定に使えていなかった**——これが本当の失敗である。

| 定理 | 内容 |
|---|---|
| `unitValV` / `coefV` / `unitValV_smul` | ★`V.toScheme` 側の橋 |
| `bridgeApp_smul` | ★★★★★★**スカラー両立** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{0}} (V : X.Opens) (F : X.Modules)

/-- ★`V.toScheme` 側の単位対象の値の橋。 -/
def unitValV (W : (V.toScheme.Opens)ᵒᵖ)
    (y : (((unitModules V.toScheme).val.obj W) : Type)) :
    ((V.toScheme.presheaf.obj W) : Type) := y

/-- ★係数の橋。 -/
def coefV (W : (V.toScheme.Opens)ᵒᵖ) (c : ((V.toScheme.ringCatSheaf.obj.obj W) : Type)) :
    ((V.toScheme.presheaf.obj W) : Type) := c

/-- ★★スカラー倍は積になる。 -/
theorem unitValV_smul (W : (V.toScheme.Opens)ᵒᵖ)
    (c : ((V.toScheme.ringCatSheaf.obj.obj W) : Type))
    (y : (((unitModules V.toScheme).val.obj W) : Type)) :
    unitValV V W (c • y) = coefV V W c * unitValV V W y := rfl


variable (e : (restrictPresheafFunctor X V).obj F.val ≅ 𝟙_ (PresheafModulesOn X V))

/-- ★★★★スカラー両立。 -/
theorem bridgeApp_smul (W : (V.toScheme.Opens)ᵒᵖ)
    (c : ((V.toScheme.ringCatSheaf.obj.obj W) : Type))
    (s : ((Scheme.Modules.restrict F V.ι).val.obj W : Type)) :
    bridgeApp F V e W (c • s) = c • bridgeApp F V e W s := by
  have hms := map_smul (e.hom.app (op (overObj V W.unop))).hom
    (coefOverV V W.unop ((V.ι.appIso W.unop).inv.hom c)) (secOverV V F W.unop s)
  have h1 : bridgeApp F V e W (c • s)
      = (V.ι.appIso W.unop).hom.hom
          (unitValAt V W.unop (coefOverV V W.unop ((V.ι.appIso W.unop).inv.hom c) •
            (e.hom.app (op (overObj V W.unop))).hom (secOverV V F W.unop s))) :=
    congrArg (fun z => (V.ι.appIso W.unop).hom.hom (unitValAt V W.unop z)) hms
  have h2 : (V.ι.appIso W.unop).hom.hom
        (unitValAt V W.unop (coefOverV V W.unop ((V.ι.appIso W.unop).inv.hom c) •
          (e.hom.app (op (overObj V W.unop))).hom (secOverV V F W.unop s)))
      = coefV V W c * unitValV V W (bridgeApp F V e W s) := by
    rw [unitValAt_smul, map_mul]
    exact congrArg (fun z => z * (V.ι.appIso W.unop).hom.hom
      (unitValAt V W.unop ((e.hom.app (op (overObj V W.unop))).hom (secOverV V F W.unop s))))
      (congrArg (fun (m : _ ⟶ _) => (CommRingCat.Hom.hom m) c)
        (V.ι.appIso W.unop).inv_hom_id)
  exact (h1.trans h2).trans (unitValV_smul V W c (bridgeApp F V e W s)).symm


/-! ## ★出典の紐付け(`.src`) -/

def bridgeApp_smul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——橋のスカラー両立)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
