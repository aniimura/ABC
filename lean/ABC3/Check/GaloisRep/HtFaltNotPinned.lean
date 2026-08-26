import ABC3.Found.GaloisRep.FaltingsWitness

/-!
# 退化封じの検査 —— **G8 の `ht^Falt` は界面に固定されていない**(`Check`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★2026-08-26 に見つかった 6 つ目の穴——**初めての「弱すぎる」型**

`FaltingsHeightData` が `htFalt` に課しているのは 2 本だけである:

* `htFalt_variableChange`(変数変換で不変)
* `prop_3_4`(`deg∞/(12(1+ε)) ≤ htFalt + C`)

★★★**`htFalt := deg∞/12` はどちらも満たす**
——前者は `deg∞` が変数変換で不変だから(第 329 の `degInfOf_variableChange`、**真の定理**)、
後者は `deg∞ ≥ 0` と `1+ε > 1` から `C = 0` で出る。

★★★★★★したがって **`Proposition 3.4` の数学的内容は witness では証明されていない**
——不等式が**恒等的に成り立つ形**で埋まってしまう。

## ★★これで界面の欠陥は 6 つ目、ただし種類が違う

| # | 場所 | 欠陥の型 | 塞いだ |
|---|---|---|---|
| 1 | G6 `localHeight` | **充足不能**(付値が任意) | 第 302 |
| 2 | G6 全体 | **充足不能**(`Δ = 0`) | 第 304 |
| 3 | G7 `omega` | **弱すぎる**(曲線と結ばれていない) | 第 311 |
| 4 | G7 `omegaFrac` | **充足不能**(`Δ = 0`) | 第 317 |
| 5 | G8 `degInf_ge_localHeight` | **充足不能**(分岐で偽) | 第 318 |
| 6 | **G8 `htFalt`** | **弱すぎる**(`deg∞/12` で満たせる) | ★**未** |

★★#3 は塞げた(`omegaFrac_variableChange` を足した)。
★★★★#6 が塞げないのは、**塞ぐ条件が (D3) の計量に依存する**からである(§9-404):
`ht^Falt = deg(ω_E)` は算術直線束の**計量込みの**次数であり、
有限部分だけでは変数変換不変にすらならない
(積公式により `Σ_p v_p(u)·log N(p)` が残り、アルキメデス側が打ち消す)。

## ★★★★★本ファイルが示すこと(限界の明示)

★示すのは「**現在の界面は `ht^Falt` を固定していない**」ことだけである。
★★「`FaltingsHeightData` が意味を持たない」ことでも、
「原文の `Proposition 3.4` が偽」でもない。
★★★`Check/GaloisRep/OmegaNondegenerate.lean` と同じ立場である。
-/

namespace ABC3.Check.GaloisRep

open ABC3.Interface.GaloisRep ABC3.Found.GaloisRep

/-! ## ★★★★★★`htFalt` は `deg∞/12` で満たせる -/

/-- ★★★★★★**界面は `ht^Falt` を固定していない**。

`FaltingsHeightData` の witness であって、`htFalt` が `deg∞/12` に一致するものが存在する。
★したがって `htFalt` の欄は「無限遠因子の次数の定数倍」で埋まってしまい、
`prop_3_4` は**恒等的に成り立つ形**になる。 -/
theorem htFalt_not_pinned :
    ∃ D : FaltingsHeightData,
      ∀ (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L),
        D.htFalt L E = D.degInf L E / 12 :=
  ⟨faltingsHeightDataWitness, fun _ _ _ _ => rfl⟩

/-- ★★★★★**`prop_3_4` はその witness では `C = 0` で成り立つ**——内容が残らない。 -/
theorem prop_3_4_trivial_for_witness (ε : ℝ) (hε : 0 < ε)
    (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L) :
    faltingsHeightDataWitness.degInf L E / (12 * (1 + ε))
      ≤ faltingsHeightDataWitness.htFalt L E + 0 := by
  have h0 : 0 ≤ degInfOf L E := degInfOf_nonneg E
  show degInfOf L E / (12 * (1 + ε)) ≤ degInfOf L E / 12 + 0
  rw [add_zero]
  gcongr
  nlinarith

/-! ## ★出典の紐付け(`.src`) -/

def htFalt_not_pinned.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Check.GaloisRep
