import ABC3.Found.Arakelov.PicCartierWitness
import Mathlib.AlgebraicGeometry.Modules.Tilde
import Mathlib.Analysis.Complex.Basic

/-!
# Arakelov (C3) 第 244 ブロック —— ★★★★★**複素点でのファイバーと切断の評価**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★§9-282 の結論を訂正する

§9-282 で「(C3) の欠けている欄 `|s|_L` を書くには `X^arc` の点で切断を**評価**でき
なければならず、それは可逆層の**解析化**そのものである」と書いた。

★★★**これは誤りである。**複素点 `p : Spec ℂ ⟶ X` における切断の評価は
**引き戻しの単位射**そのものであり、**完全に代数的**である:

    Γ(X, L) →[η]→ Γ(X, p_* p^* L) = Γ(Spec ℂ, p^* L) = arcFiber p L

★★本ブロックはそれを Lean で確かめたものである——**摩擦ゼロで通った**。

## ★★残る本当の障害は 2 つに絞られた

| 何 | 種類 |
|---|---|
| `p ↦ ‖s‖(p)` の連続性 | ★**条件**であって構成ではない(欄として書けばよい) |
| 連続計量の**存在**(`metric_nonempty`) | ★★★局所自明性 + 1 の分割 ⟹ `X^arc` のパラコンパクト性 |

★★★つまり (C3) の障害は「複素解析空間」ではなく
**点集合位相(1 の分割)**である。★見積もりの桁が変わる。

| 定義・定理 | 内容 |
|---|---|
| `arcFiber` | ★★複素点における層のファイバー(`ℂ` 上の加群) |
| `arcEval` | ★★★★**大域切断を複素点で評価する**(代数的) |
| `arcEval_add` | ★加法的 |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite

variable {X : Scheme.{0}}

/-- ★★**複素点における層のファイバー**——`ℂ` 上の加群。

★`p^* L` の `Spec ℂ` 上の大域切断である。 -/
noncomputable def arcFiber (p : Spec (CommRingCat.of ℂ) ⟶ X) (L : X.Modules) :
    ModuleCat (CommRingCat.of ℂ) :=
  moduleSpecΓFunctor (R := CommRingCat.of ℂ) |>.obj ((Scheme.Modules.pullback p).obj L)

/-- ★★★★**大域切断を複素点で評価する**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★引き戻し ⊣ 押し出しの**単位射**を `⊤` で見るだけである——解析化は要らない。
★★これが `|s|_L` を書くための材料であり、§9-282 の結論の訂正である。 -/
noncomputable def arcEval (p : Spec (CommRingCat.of ℂ) ⟶ X) (L : X.Modules) :
    (L.val.obj (op ⊤) : Type) → (arcFiber p L : Type) :=
  fun s => (((Scheme.Modules.pullbackPushforwardAdjunction p).unit.app L).val.app (op ⊤)).hom s

/-- ★評価は加法的である。 -/
theorem arcEval_add (p : Spec (CommRingCat.of ℂ) ⟶ X) (L : X.Modules)
    (s t : (L.val.obj (op ⊤) : Type)) :
    arcEval p L (s + t) = arcEval p L s + arcEval p L t :=
  map_add _ _ _

/-! ## ★出典の紐付け(`.src`) -/

def arcEval.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——複素点での切断の評価)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
