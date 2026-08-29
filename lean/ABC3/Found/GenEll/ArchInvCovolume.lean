/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.CurveArchInv
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★`curveArchInv` の共体積表示（`Found`、無条件）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★★★★★これは何か

    **`log(curveArchInv W) = log‖Δ_W‖ − 12·log‖u‖ + 6·log(covol P)`**

（`C • W = latticeCurve P`、`u = C.u`）。

★`curveArchInv` は `‖Δ_lat‖·covol⁶` として定義されている（`§9-353`）。
本ファイルは**それを曲線 `W` 自身の判別式で書き直す**——
変数変換の因子 `u` が明示的に現れる形にする。

## ★★なぜ要るのか —— `htFaltOf` の共体積表示へ

`§9-1020`（第 576）の測定:

    `12·htFaltOf = −12·log(2π) − (12/d)·Σ_p neronExp_p·log N(p) − (6/d)·Σ_σ log covol(Λ_{E^σ})`

★この式を Lean で証明するには 2 つの道具が要る:

| 道具 | 状態 |
|---|---|
| 対数版の積公式 `Σ_σ log‖σx‖ = Σᶠ_p v_p(x)·log N(p)` | ★**済**（`§9-1021`、第 577） |
| `log(archNorm)` を `log‖σΔ_E‖`・`log‖u_σ‖`・`log covol` に分解 | ★★**本ファイル** |

★★★この 2 つを合わせると `Σ_σ log‖σΔ_E‖` が `deg∞` の側と**打ち消し合い**、
`neronExp` と `covol` だけが残る。

## ★★★これが同種写像評価に効く道

☆残っている `ht^Falt(E/H) ≤ ht^Falt(E) + 2·log(l)` は、上の表示のもとで

    `12(htFaltOf(E′) − htFaltOf(E))`
      `= −(12/d)Σ_p [neronExp_p(E′) − neronExp_p(E)]·log N(p) + (12/d)Σ_σ log‖u′_σ/u_σ‖`
        `− (6/d)Σ_σ log(covol P′_σ / covol P_σ)`

になる。★アルキメデスの共体積の比は指数 `l` で `1/l`（`§9-1017`・`§9-1019`、**済**）。
☆残るのは `neronExp` と `u` の項——それが `[DelSB616]` の段 1・段 2 である。
-/

namespace ABC3.Found.GenEll

open Complex WeierstrassCurve

/-! ## ★★★★★★★★★★★★共体積表示 -/

/-- ★★★★★★★★★★★★**`log(curveArchInv W) = log‖Δ_W‖ − 12·log‖u‖ + 6·log(covol P)`**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`curveArchInv W = ‖Δ_lat(P)‖·covol(P)⁶`（`§9-354` の `curveArchInv_eq`）と
`Δ_lat(P) = (u⁻¹)¹²·Δ_W`（`WeierstrassCurve.variableChange_Δ`）から出る。
★★変数変換の因子 `u` が**明示される**のが要点である
——`htFaltOf` の共体積表示ではこの項が `neronExp` と対になる。 -/
theorem log_curveArchInv_eq {W : WeierstrassCurve ℂ} (h : W.IsElliptic)
    {P : PeriodPair} {C : VariableChange ℂ} (hC : C • W = latticeCurve P) :
    Real.log (curveArchInv W)
      = Real.log ‖W.Δ‖ - 12 * Real.log ‖(C.u : ℂ)‖ + 6 * Real.log (covol P) := by
  haveI := h
  have hΔ : W.Δ ≠ 0 := h.isUnit.ne_zero
  have hu : ((C.u : ℂ)) ≠ 0 := C.u.ne_zero
  have hcov : (0:ℝ) < covol P := covol_pos P
  have hdisc : latticeDisc P = ((C.u⁻¹ : ℂˣ) : ℂ) ^ 12 * W.Δ := by
    rw [← latticeCurve_Δ, ← hC, WeierstrassCurve.variableChange_Δ]
  have hnorm : ‖latticeDisc P‖ = ‖(C.u : ℂ)‖⁻¹ ^ 12 * ‖W.Δ‖ := by
    rw [hdisc, norm_mul, norm_pow]
    congr 2
    simp
  rw [curveArchInv_eq h hC, archInv, hnorm]
  have h1 : (0:ℝ) < ‖(C.u : ℂ)‖ := norm_pos_iff.2 hu
  have h2 : (0:ℝ) < ‖W.Δ‖ := norm_pos_iff.2 hΔ
  rw [Real.log_mul (by positivity) (by positivity), Real.log_mul (by positivity) (by positivity),
    Real.log_pow, Real.log_inv, Real.log_pow]
  push_cast
  ring

/-- ★★★★★★★★**埋め込みごとの形**——`archSum` が使う `archNorm` で。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`archNorm E σ = curveArchInv (E.map σ)` と `(E.map σ).Δ = σ(E.Δ)`
（`WeierstrassCurve.map_Δ`）を通しただけである。
★★これを `σ` について足すと `archSum` の分解になる。 -/
theorem log_archNorm_eq {L : Type*} [Field L] (E : WeierstrassCurve L) [E.IsElliptic]
    (σ : L →+* ℂ) {P : PeriodPair} {C : VariableChange ℂ}
    (hC : C • (E.map σ) = latticeCurve P) :
    Real.log (archNorm E σ)
      = Real.log ‖σ E.Δ‖ - 12 * Real.log ‖(C.u : ℂ)‖ + 6 * Real.log (covol P) := by
  rw [archNorm, log_curveArchInv_eq (inferInstance : (E.map σ).IsElliptic) hC,
    WeierstrassCurve.map_Δ]

/-! ## ★出典の紐付け(`.src`) -/

def log_curveArchInv_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(curveArchInv の共体積表示——log‖Δ_W‖ − 12log‖u‖ + 6log covol)",
    sectionId := "genell-prop-3-4" }

def log_archNorm_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(archNorm の共体積表示——埋め込みごと。★無条件)",
    sectionId := "genell-prop-3-4" }

def log_archNorm_eq.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "curveArchInv_eq(curveArchInv = archInv P、§9-354)"
      (.inProject "ABC3" "ABC3.Found.GenEll.curveArchInv_eq") 3,
    .citation "[mathlib]" "WeierstrassCurve.variableChange_Δ(Δ は u⁻¹² 倍)"
      (.inMathlib "WeierstrassCurve.variableChange_Δ") 1,
    .implicitStep
      ("★★★★★★★★到達点(2026-08-29、第 578): §9-1020(第 576)の測定" ++
       "『12·htFaltOf = −12log(2π) − (12/d)Σ_p neronExp_p·log N(p) " ++
       "− (6/d)Σ_σ log covol(Λ_{E^σ})』を Lean で証明するのに要る**2 つ目の道具**。" ++
       "★1 つ目(対数版の積公式)は §9-1021(第 577)で済んでいる。" ++
       "★★両者を合わせると Σ_σ log‖σΔ_E‖ が deg∞ の側と打ち消し合い、" ++
       "neronExp と covol だけが残る") 8,
    .implicitStep
      ("☆同種写像評価への道: 上の表示のもとで " ++
       "12(htFaltOf(E′) − htFaltOf(E)) は neronExp の差・u の比・covol の比の 3 項になる。" ++
       "★covol の比は指数 l で 1/l(§9-1017・§9-1019、**済**)。" ++
       "☆残るのは neronExp と u の項——[DelSB616] の段 1・段 2 である") 9 ]

end ABC3.Found.GenEll
