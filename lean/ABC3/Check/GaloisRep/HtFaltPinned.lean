import ABC3.Found.GaloisRep.FaltingsWitness

/-!
# 退化封じの検査 —— **G8 の `ht^Falt` は界面に固定された**(`Check`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★2026-08-26 —— **欠陥 #6 が塞がった**

以前(第 329)は `FaltingsHeightData` が `htFalt` に課しているのが 2 本だけだった:

* `htFalt_variableChange`(変数変換で不変)
* `prop_3_4`(`deg∞/(12(1+ε)) ≤ htFalt + C`)

★★★**`htFalt := deg∞/12` はどちらも満たしてしまい**、
`Proposition 3.4` が**恒等的に成り立つ形**で埋まった——界面の欠陥 #6 である。

★★★★★★塞ぐのに要ったのは**アルキメデス素点での計量**であり、
そこに至る鎖(一意化 → 共体積 → `(g₂,g₃)` が束を決める → `archNorm` → 有界性)は
第 332-356 の 25 ブロックで閉じた。

## ★★★★界面に足した 4 つ

| 欄・条件 | 役割 |
|---|---|
| `archNorm` | アルキメデス素点でのノルム `‖Δ‖_arch(E^σ)` |
| `archNorm_pos` | 楕円曲線では正 |
| **`archNorm_eq`** | **モジュラー判別式で同定**——`j` の全射性により一意に決まる |
| **`htFalt_eq`** | **`12·ht^Falt = deg∞ − (1/d)Σ_σ log((2π)¹²‖Δ‖_arch)`** |

## ★★★★★★本ファイルが示すこと

★`htFalt_determined`: **`degInf` が一致すれば `htFalt` も一致する**
——`htFalt` の欄はもはや自由でない。
★★`htFalt_ne_degInf_div_twelve`: アルキメデス和が `0` でない曲線が 1 つでもあれば、
**`htFalt = deg∞/12` は排除される**。

★★★示さないこと: 「アルキメデス和が `0` でない曲線が実在すること」は
具体的な数値評価になるので、ここでは仮定の形で残す。
★`Check/GaloisRep/OmegaNondegenerate.lean` と同じ立場である。

## ★★界面の欠陥の一覧(2026-08-26 現在)

| # | 場所 | 欠陥の型 | 塞いだ |
|---|---|---|---|
| 1 | G6 `localHeight` | **充足不能**(付値が任意) | 第 302 |
| 2 | G6 全体 | **充足不能**(`Δ = 0`) | 第 304 |
| 3 | G7 `omega` | **弱すぎる**(曲線と結ばれていない) | 第 311 |
| 4 | G7 `omegaFrac` | **充足不能**(`Δ = 0`) | 第 317 |
| 5 | G8 `degInf_ge_localHeight` | **充足不能**(分岐で偽) | 第 318 |
| 6 | G8 `htFalt` | **弱すぎる**(`deg∞/12` で満たせる) | ★**第 357** |
-/

namespace ABC3.Check.GaloisRep

open NumberField WeierstrassCurve ABC3.Interface.GaloisRep ABC3.Found.GenEll

/-! ## ★★★★★★★`ht^Falt` は界面に固定されている -/

/-- ★★★★★★★**界面は `ht^Falt` を固定する**——`degInf` が一致すれば `htFalt` も一致する。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`archNorm_eq` と `j` の全射性(第 348)から `archNorm` が一致し、
`htFalt_eq` から `htFalt` が一致する。 -/
theorem htFalt_determined (D D' : FaltingsHeightData) (L : Type) [Field L] [NumberField L]
    (E : WeierstrassCurve L) [E.IsElliptic] (hdeg : D.degInf L E = D'.degInf L E) :
    D.htFalt L E = D'.htFalt L E := by
  have harch : ∀ σ : (L →+* ℂ), D.archNorm L E σ = D'.archNorm L E σ := by
    intro σ
    obtain ⟨τ, hτ⟩ := jFun_surjective ((E.map σ).j)
    have h : ModularForm.E₄ τ ^ 3 / ModularForm.discriminant τ = (E.map σ).j := hτ
    rw [D.archNorm_eq L E σ τ h, D'.archNorm_eq L E σ τ h]
  have h1 := D.htFalt_eq L E
  have h2 := D'.htFalt_eq L E
  have hsum : (∑ σ : (L →+* ℂ), Real.log ((2 * Real.pi) ^ 12 * D.archNorm L E σ))
      = ∑ σ : (L →+* ℂ), Real.log ((2 * Real.pi) ^ 12 * D'.archNorm L E σ) :=
    Finset.sum_congr rfl (fun σ _ => by rw [harch σ])
  rw [hdeg, hsum] at h1
  linarith

/-- ★★★★★★**`htFalt := deg∞/12` は排除される**——アルキメデス和が `0` でない曲線があれば。

★これが欠陥 #6 の塞ぎの内容である。 -/
theorem htFalt_ne_degInf_div_twelve (D : FaltingsHeightData) (L : Type) [Field L] [NumberField L]
    (E : WeierstrassCurve L) [E.IsElliptic]
    (hne : (∑ σ : (L →+* ℂ), Real.log ((2 * Real.pi) ^ 12 * D.archNorm L E σ)) ≠ 0) :
    D.htFalt L E ≠ D.degInf L E / 12 := by
  intro heq
  have hd : ((Module.finrank ℚ L : ℝ)) ≠ 0 := by
    have hpos : (0:ℝ) < (Module.finrank ℚ L : ℝ) := by exact_mod_cast Module.finrank_pos
    exact ne_of_gt hpos
  have h1 := D.htFalt_eq L E
  rw [heq] at h1
  refine hne ?_
  have hz : (∑ σ : (L →+* ℂ), Real.log ((2 * Real.pi) ^ 12 * D.archNorm L E σ))
      / (Module.finrank ℚ L : ℝ) = 0 := by linarith
  exact (div_eq_zero_iff.1 hz).resolve_right hd

/-! ## ★出典の紐付け(`.src`) -/

def htFalt_determined.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def htFalt_ne_degInf_div_twelve.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Check.GaloisRep
