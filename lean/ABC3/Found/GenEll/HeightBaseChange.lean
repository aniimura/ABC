import ABC3.Found.GenEll.ArchADivBase
import ABC3.Found.GenEll.HeightConstruction

/-!
# [GenEll] Definition 1.2, (i) —— **高さは定義体の取り方に依らない**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★4 段を組み上げる

`Definition 1.2, (i)` は高さを **`ht_M̄ : X(ℚ̄) → ℝ`** として定める。
★`X(ℚ̄)` の点は定義体 `F` を選んで `x_F : Spec 𝓞_F → X` で表すので、
**`F` の取り方に依らない**ことを示さないと well-defined でない。

★★本日その 4 段を取った:

| 段 | 場所 |
|---|---|
| 1. 有限素点側を `ADiv` に繋ぐ(`ord_w(J·𝓞_K) = e·ord_v(J)`) | `IdealADivBase.lean` |
| 2. アルキメデス側を `ADiv` に繋ぐ | `ArchADivBase.lean` |
| 3. 正規化次数は底変換で不変 | `BaseChange.lean`(既存) |
| 4. **組み上げ** | ★★★**本ファイル** |

## ★★仮定は「どの射で底を上げるか」の 2 本だけ

有限素点側とアルキメデス側で、それぞれ前段の結論を仮定として受ける。
★どちらも前 2 ファイルが供給する——**未証明の事実ではない**。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

variable (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K]
  [Algebra F K]

/-! ## ★★引き戻した算術因子の底変換 -/

open scoped Classical in
/-- ★★★**引き戻した算術因子は底変換で移る**。

    `x_K^* D̄ = baseChange F K (x_F^* D̄)`

★★2 つの仮定は前 2 ファイル(`IdealADivBase` / `ArchADivBase`)が供給する。 -/
theorem pullbackADiv_baseChange {X : Scheme.{0}} (D : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X)
    (φ : CommRingCat.of (𝓞 F) ⟶ CommRingCat.of (𝓞 K))
    (hJ : pullbackIdeal F D.divisor xF ≠ 0)
    (hfin : pullbackIdeal K D.divisor (Spec.map φ ≫ xF)
      = (pullbackIdeal F D.divisor xF).map (algebraMap (𝓞 F) (𝓞 K)))
    (harch : archADiv K D.green (Spec.map φ ≫ xF)
      = baseChangeArc F K (pullbackADiv F D xF))
    (hlies : ∀ W : FinitePlace K, W.asIdeal.LiesOver (hosComap F K W).asIdeal) :
    pullbackADiv K D (Spec.map φ ≫ xF)
      = baseChange F K (pullbackADiv F D xF) := by
  refine Prod.ext ?_ ?_
  · -- 有限素点側
    show (pullbackADiv K D (Spec.map φ ≫ xF)).fin
      = (baseChange F K (pullbackADiv F D xF)).fin
    rw [pullbackADiv_fin, baseChange_fin, hfin]
    ext W
    haveI := hlies W
    rw [idealADiv_fin_map F K (pullbackIdeal F D.divisor xF) hJ W]
    show _ = baseChangeFin F K (pullbackADiv F D xF) W
    rw [baseChangeFin]
    rfl
  · -- アルキメデス側
    show (pullbackADiv K D (Spec.map φ ≫ xF)).arc
      = (baseChange F K (pullbackADiv F D xF)).arc
    rw [pullbackADiv_arc, baseChange_arc]
    exact harch

/-! ## ★★★高さの底変換不変性 -/

/-- ★★★**[GenEll] Definition 1.2, (i)** —— 高さは定義体の取り方に依らない。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★★★**これで `ht` が `X(ℚ̄)` の上で well-defined になる。**

★機構は `pullbackADiv_baseChange` + `degNormalized_baseChange` の 2 行。
★★本日取った 4 段がここで 1 つになる。 -/
theorem htArith_baseChange {X : Scheme.{0}} (D : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X)
    (φ : CommRingCat.of (𝓞 F) ⟶ CommRingCat.of (𝓞 K))
    (hJ : pullbackIdeal F D.divisor xF ≠ 0)
    (hfin : pullbackIdeal K D.divisor (Spec.map φ ≫ xF)
      = (pullbackIdeal F D.divisor xF).map (algebraMap (𝓞 F) (𝓞 K)))
    (harch : archADiv K D.green (Spec.map φ ≫ xF)
      = baseChangeArc F K (pullbackADiv F D xF))
    (hlies : ∀ W : FinitePlace K, W.asIdeal.LiesOver (hosComap F K W).asIdeal) :
    htArith K D (Spec.map φ ≫ xF) = htArith F D xF := by
  rw [htArith, htArith,
    pullbackADiv_baseChange F K D xF φ hJ hfin harch hlies,
    degNormalized_baseChange]

/-! ## ★出典の紐付け(`.src`)

★条つきである。`Definition 1.2` 全体には `X(ℚ̄)` の型そのものの構成
(数体についての colimit)が残っている——★それは設計の仕事である。 -/

def htArith_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2, (i)(高さの底変換不変性——型は AlgPointClass.lean)",
    sectionId := "genell-def-1-2-i" }

end ABC3.Found.GenEll
