import ABC3.Found.GaloisRep.DegInfLocal
import ABC3.Found.GaloisRep.NeronWitness

/-!
# Galois (G8) 第 329 ブロック —— **★★★★★★★★★★G8 の witness**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★★★到達点

> **`FaltingsHeightData.nonvacuous`**——(G8) の欄がすべて埋まった

★★★これで Galois 表現論の 8 件が**すべて**埋まる。

## ★★★★★★★★★ただし——`htFalt` は Faltings 高さではない(必読)

★本 witness の `htFalt` は **`deg∞/12`** である。**Faltings 高さではない。**

★★これは界面が `htFalt` を**固定していない**ことの帰結である:
界面が `htFalt` に課しているのは

* `htFalt_variableChange`(変数変換で不変)
* `prop_3_4`(`deg∞/(12(1+ε)) ≤ htFalt + C`)

の 2 本だけで、★★★★`htFalt := deg∞/12` は**どちらも満たす**
——前者は `deg∞` が変数変換で不変だから(本ブロックの `degInfOf_variableChange`、**真の定理**)、
後者は `deg∞ ≥ 0` と `1+ε > 1` から `C = 0` で出る。

★★★★★★したがって **`Proposition 3.4` の数学的内容は本 witness では証明されていない**。
★これは界面の 6 つ目の欠陥であり、**初めての「弱すぎる」型**である
(1-5 は「充足不能」だった)。`Check/GaloisRep/HtFaltNotPinned.lean` に機械可読の形で記録した。

★★★★塞ぐには **(D3) の計量**が要る(§9-404)——`ht^Falt = deg(ω_E)` は
算術直線束の**計量込みの**次数であり、有限部分だけでは変数変換不変にすらならない
(積公式により `Σ_p v_p(u)·log N(p)` が残る)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `minDeltaExp_variableChange` | ★★★★★**極小判別式の指数は変数変換で不変**(真) |
| `degInfOf_variableChange` | ★★★★★★**`deg∞` は変数変換で不変**(真) |
| `faltingsHeightDataWitness` | ★★★★★★★★**`FaltingsHeightData` の実装** |
| `FaltingsHeightData.nonvacuous` | ★★★★★★★★★★**G8 の欄が埋まる** |
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField WeierstrassCurve ABC3.Interface.GaloisRep
open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-! ## ★★★★★★`deg∞` は変数変換で不変(真の定理) -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★**極小判別式の指数は変数変換で不変である**。

★`Δ(C•W) = u⁻¹²Δ` で付値は `−12v(u)`、Néron 指数は `−v(u)` だけ動くので、
`v(Δ) − 12·(Néron 指数)` は打ち消し合う。 -/
theorem minDeltaExp_variableChange (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    (C : VariableChange L) : minDeltaExp p (C • W) = minDeltaExp p W := by
  by_cases hΔ : W.Δ = 0
  · have h2 : (C • W).Δ = 0 := by
      rw [WeierstrassCurve.variableChange_Δ, hΔ, mul_zero]
    rw [minDeltaExp, dif_pos h2, minDeltaExp, dif_pos hΔ]
  · have h2 : (C • W).Δ ≠ 0 := variableChange_Delta_ne_zero W hΔ C
    rw [minDeltaExp, dif_neg h2, minDeltaExp, dif_neg hΔ, neronExp_variableChange p W hΔ C]
    have hu : (Units.mk0 ((C • W).Δ) h2) = C.u⁻¹ ^ 12 * Units.mk0 W.Δ hΔ := by
      refine Units.ext ?_
      show (C • W).Δ = _
      rw [WeierstrassCurve.variableChange_Δ]
      push_cast
      simp
    rw [hu, valAdd_mul, valAdd_pow, valAdd_inv]
    omega

/-- ★★★★★★**`deg∞` は変数変換で不変である**。 -/
theorem degInfOf_variableChange (E : WeierstrassCurve L) (C : VariableChange L) :
    degInfOf L (C • E) = degInfOf L E := by
  rw [degInfOf, degInfOf]
  congr 1
  refine finsum_congr (fun q => ?_)
  rw [minDeltaExp_variableChange]

/-! ## ★★★★★★★★G8 の witness -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★**`FaltingsHeightData` の実装**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★★★★★★★★**注意**: `htFalt` は `deg∞/12` であって **Faltings 高さではない**。
界面が `htFalt` に課す 2 条件(変数変換不変・`prop_3_4`)を満たすだけの witness であり、
`Proposition 3.4` の数学的内容は**証明されていない**。
★詳細は `Check/GaloisRep/HtFaltNotPinned.lean` と本ファイルの冒頭を見よ。 -/
noncomputable def faltingsHeightDataWitness : FaltingsHeightData where
  toSemistableModelData := semistableModelDataWitness
  htFalt := fun L _ _ E => degInfOf L E / 12
  htFalt_variableChange := by
    intro L _ _ E C
    rw [degInfOf_variableChange]
  degInf := fun L _ _ E => degInfOf L E
  degInf_nonneg := by
    intro L _ _ E
    exact degInfOf_nonneg E
  degInf_ge_localHeight := by
    intro L _ _ E Lv _ _ R _ _ _ _ _ _ _ _ h p hp
    exact degInfOf_ge_localHeight E h p hp
  prop_3_4 := by
    intro ε hε
    refine ⟨0, ?_⟩
    intro L _ _ E _
    have h0 : 0 ≤ degInfOf L E := degInfOf_nonneg E
    show degInfOf L E / (12 * (1 + ε)) ≤ degInfOf L E / 12 + 0
    rw [add_zero]
    gcongr
    nlinarith

/-- ★★★★★★★★★★**`FaltingsHeightData` は非空虚である**——G8 の欄が埋まる。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★★★★★★**ただし `htFalt` は `deg∞/12`** であり、Faltings 高さではない。
界面が `ht^Falt` を固定していないことの記録は `Check/GaloisRep/HtFaltNotPinned.lean`。 -/
theorem FaltingsHeightData.nonvacuous : Nonempty FaltingsHeightData :=
  ⟨faltingsHeightDataWitness⟩

/-! ## ★出典の紐付け(`.src`) -/

def degInfOf_variableChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Proposition 3.4(無限遠因子の次数 deg∞)",
    sectionId := "genell-prop-3-4" }

def faltingsHeightDataWitness.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GaloisRep
