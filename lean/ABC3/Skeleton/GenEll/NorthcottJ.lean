/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.EllModuliObjects
import ABC3.Found.GaloisRep.HtJBound
import ABC3.Meta.Claim

/-!
# `j` の高さの Northcott 性（`Skeleton`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★★★これは何か

`Interface/GenEll/EllModuli.lean` の `EllModuliData` の

    northcott : ∀ (C : ℝ) (d : ℕ), 0 < d → {x ∈ degLe d | faltingsHeight x ≤ C}.Finite

欄を埋めるのに要る**唯一の外部入力**を、本プロジェクトの語彙で固定する。

原文は `Proposition 1.4, (iv)` から取る。★実質は**古典的な Northcott の定理**

> 次数と高さを抑えた代数的数は有限個しかない

である。

## ★★★★なぜ `Skeleton` に置くか——mathlib の測定（2026-08-31）

mathlib には `Mathlib/Order/Northcott.lean`（`Northcott` 型クラス）と
`Mathlib/NumberTheory/Height/Northcott.lean` があるが、**具体的な体に対する
`Northcott` の instance は 1 つも無い**（`Northcott (mulHeight₁ (K := K))` は
どこでも証明されていない）。`Mathlib/NumberTheory/Height/EllipticCurve.lean` も
「TODO: Define the naïve height」の段階である。

★したがってこれは**mathlib の欠落**であり、依存グラフの節点 1 つとして立てる。

## ★★★何を証明すれば終わりか

1. `htJ L E` は `E.j` の**絶対 Weil 高さ**である（`htFinJ + htArchJ`、どちらも `[L:ℚ]` で割ってある）
2. 古典的 Northcott: `{x : ℚ̄ | [ℚ(x):ℚ] ≤ d ∧ h(x) ≤ B}` は有限
   （最小多項式の係数が Mahler 測度で抑えられる ⟹ 多項式が有限個）
3. `[ℚ(E.j) : ℚ] ≤ [L : ℚ] = E.deg`

☆mathlib には `Mathlib/NumberTheory/MahlerMeasure.lean` があるので、2 の素材はある。

## ★これを消費する側

`Found/GenEll/EllModuliObjects.lean` の `northcottJ_of_northcott_htJ`（`§9-1170`、第 743）が
本命題を仮説として受け、`EllModuliData.northcott` 欄の形にする。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Found.GenEll ABC3.Found.GaloisRep

/-- **[GenEll] `Proposition 3.4` が使う `Proposition 1.4, (iv)`**
——`j` の高さの Northcott 性。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★★★★これが `EllModuliData.northcott` 欄に要る唯一の外部入力である。 -/
theorem northcott_htJ (B : ℝ) (d : ℕ) :
    {x : ℂ | ∃ E : SSCurve, E.j = x ∧ E.deg ≤ d ∧ htJ E.fld E.W ≤ B}.Finite := by
  sorry

def northcott_htJ.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Proposition 1.4, (iv) の中身——j の高さの Northcott 性)",
    sectionId := "genell-prop-3-4" }

def northcott_htJ.needs : List ProofObligation :=
  [ .citation "[GenEll]" "Proposition 1.4, (iv)(高さと次数で抑えた類は有限)"
      (.absent "mathlib: Mathlib/Order/Northcott.lean と Mathlib/NumberTheory/Height/ を 'Northcott (mulHeight' で grep して instance 0 件(2026-08-31)") 3,
    .implicitStep
      ("★★★古典的な Northcott の定理: 次数 ≤ d かつ絶対 Weil 高さ ≤ B の代数的数は有限個。" ++
       "最小多項式の係数が Mahler 測度で抑えられることから出る" ++
       "(mathlib の Mathlib/NumberTheory/MahlerMeasure.lean が素材)") 12,
    .implicitStep
      ("★htJ L E = htFinJ + htArchJ は E.j の絶対 Weil 高さであること" ++
       "(どちらも [L:ℚ] で割ってあるので基底変換で不変)") 5,
    .implicitStep "★[ℚ(E.j) : ℚ] ≤ [L : ℚ] = E.deg" 2 ]

end ABC3.Skeleton.GenEll
