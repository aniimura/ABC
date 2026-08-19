import ABC3.Found.Arakelov.ArcBridgeIso

/-!
# Arakelov (C3) 第 270 ブロック —— ★★★★★★**橋から直接、連続な局所ノルム**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★迂回路が**不要になった**

第 257–261 で作った `genNorm`(基準ノルムの**比**で選択を消す)は、
**橋が無いとき**の迂回路だった。★第 269 で橋(`restrict F V.ι ≅ 𝒪_V`)が出たので、
**第 253 の `trivNorm` がそのまま使える**:

| 法則 | 出どころ |
|---|---|
| 非負・`0 ↔ 0`・`‖c·v‖ = \|c\|‖v‖` | ★第 253 `trivNorm_*` |
| ★**連続性** | ★★第 253 `continuous_trivNorm` |

★★★**生成切断の非消滅は要らなくなった**——第 261 で仮定として受けたものが、
迂回路ごと消えた。

## ★★迂回路は無駄ではなかった

| 第 257–261 で得たもの | その後 |
|---|---|
| 比で選択を消す発想 | ★使わない(橋があるので) |
| `arcEvalOnTop`(開集合上の評価) | ★★**使う**(§9-298 の迂回) |
| 半線形性・`evalOn` の連続性 | ★★使う |
| 貼り合わせ(第 262–263) | ★★★**`localNorm` 版に書き直しが要る** |

★★★★教訓: **壁があるときの迂回路は、壁が崩れたら作り直す**——
ただし迂回の途中で作った**部品**(評価・連続性)は残る。

| 定義・定理 | 内容 |
|---|---|
| `localNorm` | ★★★★★橋から得る局所ノルム |
| `localNorm_nonneg` / `_eq_zero_iff` / `_smul` | ★3 法則 |
| `continuous_localNorm` | ★★★★★★**連続性**(仮定なし) |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{0}} (V : X.Opens) (F : X.Modules)
  (e : (restrictPresheafFunctor X V).obj F.val ≅ 𝟙_ (PresheafModulesOn X V))

/-- ★★★★★橋から得る局所ノルム(`V.toScheme` の点で)。 -/
noncomputable def localNorm (q : Spec (CommRingCat.of ℂ) ⟶ V.toScheme)
    (w : ↥(arcFiber q (Scheme.Modules.restrict F V.ι))) : ℝ :=
  trivNorm (bridgeIso V F e) q w

theorem localNorm_nonneg (q : Spec (CommRingCat.of ℂ) ⟶ V.toScheme)
    (w : ↥(arcFiber q (Scheme.Modules.restrict F V.ι))) : 0 ≤ localNorm V F e q w :=
  trivNorm_nonneg _ _ _

theorem localNorm_eq_zero_iff (q : Spec (CommRingCat.of ℂ) ⟶ V.toScheme)
    (w : ↥(arcFiber q (Scheme.Modules.restrict F V.ι))) :
    localNorm V F e q w = 0 ↔ w = 0 :=
  trivNorm_eq_zero_iff _ _ _

theorem localNorm_smul (q : Spec (CommRingCat.of ℂ) ⟶ V.toScheme)
    (c : (CommRingCat.of ℂ : Type))
    (w : ↥(arcFiber q (Scheme.Modules.restrict F V.ι))) :
    localNorm V F e q (c • w) = ‖c‖ * localNorm V F e q w :=
  trivNorm_smul _ _ _ _

/-- ★★★★★★連続性——第 253 がそのまま効く。 -/
theorem continuous_localNorm
    (s : (((Scheme.Modules.restrict F V.ι).val.obj (op ⊤)) : Type)) :
    @Continuous _ ℝ (arcTopology V.toScheme) _
      (fun q => localNorm V F e q (arcEval q (Scheme.Modules.restrict F V.ι) s)) :=
  continuous_trivNorm (bridgeIso V F e) s


/-! ## ★出典の紐付け(`.src`) -/

def continuous_localNorm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——橋から直接、連続な局所ノルム)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
