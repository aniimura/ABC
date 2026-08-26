import ABC3.Found.GenEll.PeterssonBound

/-!
# GenEll 第 356 ブロック —— **★★★★★`archNorm` を全曲線へ延ばす**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★なぜ楕円性を外すのか——界面の `htFalt_variableChange` と `prop_3_4`

`Interface/GaloisRep/Reduction.lean` の `htFalt_variableChange` と `prop_3_4` は
**楕円性を仮定していない**(全 `WeierstrassCurve` について量化している)。
★したがって witness を組むには、`archNorm` の変数変換不変性と上界も
**楕円性なしの形**で要る。

★★`Δ = 0` の曲線では `curveArchInv = 0`(定義の `dif_neg` 側)であり、
変数変換は `Δ` の零性を保つ(`Δ(C•W) = u⁻¹²Δ`)ので、両方とも `0` になって一致する。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `isElliptic_variableChange_iff` | ★変数変換は楕円性を保つ |
| `curveArchInv_eq_zero`・`curveArchInv_nonneg` | ★非楕円では `0`、つねに非負 |
| `curveArchInv_variableChange'` | ★★★★★**変数変換不変(仮定なし)** |
| `archNorm_nonneg`・`archNorm_variableChange'` | ★★★★★同(埋め込み版) |
| `exists_bound_curveArchInv'` | ★★★★★★**一様な上界 `1 ≤ M`** |
-/

namespace ABC3.Found.GenEll

open Complex WeierstrassCurve

/-! ## ★変数変換は楕円性を保つ -/

theorem isElliptic_variableChange_iff (W : WeierstrassCurve ℂ) (C : VariableChange ℂ) :
    (C • W).IsElliptic ↔ W.IsElliptic := by
  rw [WeierstrassCurve.isElliptic_iff, WeierstrassCurve.isElliptic_iff,
    WeierstrassCurve.variableChange_Δ]
  simp [isUnit_iff_ne_zero]

theorem curveArchInv_eq_zero (W : WeierstrassCurve ℂ) (h : ¬ W.IsElliptic) :
    curveArchInv W = 0 := by
  classical
  rw [curveArchInv, dif_neg h]

theorem curveArchInv_nonneg (W : WeierstrassCurve ℂ) : 0 ≤ curveArchInv W := by
  by_cases h : W.IsElliptic
  · exact (curveArchInv_pos h).le
  · rw [curveArchInv_eq_zero W h]

/-! ## ★★★★★楕円性を仮定しない不変性と上界 -/

/-- ★★★★★**変数変換で変わらない**(楕円性を仮定しない形)。 -/
theorem curveArchInv_variableChange' (W : WeierstrassCurve ℂ) (C : VariableChange ℂ) :
    curveArchInv (C • W) = curveArchInv W := by
  by_cases h : W.IsElliptic
  · exact curveArchInv_variableChange W h C
  · rw [curveArchInv_eq_zero _ (fun hc => h ((isElliptic_variableChange_iff W C).1 hc)),
      curveArchInv_eq_zero W h]

theorem archNorm_nonneg {L : Type*} [Field L] (E : WeierstrassCurve L) (σ : L →+* ℂ) :
    0 ≤ archNorm E σ := curveArchInv_nonneg _

/-- ★★★★★**変数変換で変わらない**(埋め込み版、楕円性を仮定しない形)。 -/
theorem archNorm_variableChange' {L : Type*} [Field L] (E : WeierstrassCurve L)
    (C : WeierstrassCurve.VariableChange L) (σ : L →+* ℂ) :
    archNorm (C • E) σ = archNorm E σ := by
  rw [archNorm, archNorm, ← WeierstrassCurve.map_variableChange E C σ]
  exact curveArchInv_variableChange' _ _

/-- ★★★★★★**一様な上界**(楕円性を仮定しない形、`1 ≤ M`)。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`1 ≤ M` にしておくと `log((2π)¹²M) ≥ 0` が言えて、非楕円の `log 0 = 0` も抑えられる。 -/
theorem exists_bound_curveArchInv' :
    ∃ M : ℝ, 1 ≤ M ∧ ∀ W : WeierstrassCurve ℂ, curveArchInv W ≤ M := by
  obtain ⟨M, hMpos, hM⟩ := exists_bound_curveArchInv
  refine ⟨max M 1, le_max_right _ _, fun W => ?_⟩
  by_cases h : W.IsElliptic
  · exact le_trans (hM W h) (le_max_left _ _)
  · rw [curveArchInv_eq_zero W h]
    exact le_trans zero_le_one (le_max_right _ _)

/-! ## ★出典の紐付け(`.src`) -/

def curveArchInv_variableChange'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def exists_bound_curveArchInv'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
