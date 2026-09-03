/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.GalActVc
import ABC3.Found.GenEll.GalPointCoord
import ABC3.Meta.Claim

/-!
# 第 1379 ブロック —— **`galPoint` は変数変換と可換**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——`C • E` の局所データを `E` に戻す

第 1372（`exists_h2_h1_of_bad_prime_ram`）は
`[IsMinimal R (E.baseChange Lv)]` を要求するが、
`SSCurve` の `E.W` 自体は `p`-極小とは限らない。
☆`SemistableAt` が与えるのは **`C • E.W`** の極小性である。

★★★したがって局所データ（`P₀`・`hcard`・`Pfix`・`Qmov`）は `C • E.W` の側で出る。
それを `E.W` に戻すのが本ブロックである。

☆第 1288（`galAct_vcPoint`）が「`galAct` は変数変換と可換」を持っているので、
**`galPoint` と `galAct` が一致する**ことだけ示せばよい
——どちらも座標に `σ` を当てるだけだからである（第 1298・第 1285）。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Found.GaloisRep ABC3.Meta
open ABC3.Interface.GaloisRep

open scoped Classical

section GalPoint

variable {K L : Type} [Field K] [Field L] [Algebra K L]

/-- ★★★★**`galPoint` は原点でない点を原点でない点に送る**（第 1379）。 -/
theorem galPoint_ne_zero (W : WeierstrassCurve K) (σ : L ≃ₐ[K] L)
    {P : (W.baseChange L).toAffine.Point} (hP : P ≠ 0) : galPoint W σ P ≠ 0 := by
  cases P with
  | zero => exact absurd rfl hP
  | some x y h =>
      intro hcon
      exact absurd hcon (by simp [galPoint, WeierstrassCurve.Affine.Point.map])

/-- ★★★★★★**底変換した曲線は `σ` で固定される**（第 1379）。 -/
theorem baseChange_map_algEquiv (W : WeierstrassCurve K) (σ : L ≃ₐ[K] L) :
    (W.baseChange L).map (σ : L →+* L) = W.baseChange L := by
  show (W.map (algebraMap K L)).map (σ : L →+* L) = W.map (algebraMap K L)
  rw [WeierstrassCurve.map_map]
  congr 1
  exact RingHom.ext (fun x => σ.commutes x)

/-- ★★★★★★★★**`galPoint` は `galAct` と一致する**——★**無条件**（第 1379）。

☆どちらも座標に `σ` を当てるだけである。 -/
theorem galPoint_eq_galAct (W : WeierstrassCurve K) (σ : L ≃ₐ[K] L)
    (P : (W.baseChange L).toAffine.Point) :
    galPoint W σ P
      = galAct (σ : L →+* L) (W.baseChange L) (baseChange_map_algEquiv W σ) P := by
  by_cases hP : P = 0
  · subst hP
    rw [map_zero, map_zero]
  · refine pointCoords_injective (galPoint_ne_zero W σ hP)
      (galAct_ne_zero (σ : L →+* L) (W.baseChange L) _ hP) ?_
    rw [pointCoords_galPoint, pointCoords_galAct]
    rfl

end GalPoint

section GalPointVc

variable {K L : Type} [Field K] [Field L] [Algebra K L]

/-- ★★★★★★**変数変換の底変換**（第 1379）。 -/
theorem baseChange_variableChange (W : WeierstrassCurve K) (C : VariableChange K) :
    (C • W).baseChange L = (C.map (algebraMap K L)) • (W.baseChange L) :=
  (WeierstrassCurve.map_variableChange W C (algebraMap K L)).symm

/-- ★★★★**底変換した変数変換は `σ` で固定される**（第 1379）。 -/
theorem map_variableChange_algEquiv (C : VariableChange K) (σ : L ≃ₐ[K] L) :
    (C.map (algebraMap K L)).map (σ : L →+* L) = C.map (algebraMap K L) := by
  rw [VariableChange.map_map]
  congr 1
  exact RingHom.ext (fun x => σ.commutes x)

/-- ★★★★★★★★★★★★★★★★
**`galPoint` は変数変換と可換**——★**無条件**（第 1379）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これで `C • E` の局所データを `E` に戻せる。 -/
theorem galAct_vcPoint_baseChange (W : WeierstrassCurve K) (C : VariableChange K)
    (σ : L ≃ₐ[K] L) (P : (W.baseChange L).toAffine.Point) :
    galAct (σ : L →+* L) ((C.map (algebraMap K L)) • (W.baseChange L))
        (by rw [← WeierstrassCurve.map_variableChange,
          map_variableChange_algEquiv C σ, baseChange_map_algEquiv W σ])
        (vcPoint (C.map (algebraMap K L)) (W.baseChange L) P)
      = vcPoint (C.map (algebraMap K L)) (W.baseChange L) (galPoint W σ P) := by
  rw [galPoint_eq_galAct W σ P]
  exact galAct_vcPoint (σ : L →+* L) (W.baseChange L) (C.map (algebraMap K L))
    (map_variableChange_algEquiv C σ) (baseChange_map_algEquiv W σ) P

end GalPointVc

/-! ## ★出典の紐付け(`.src`) -/

def galPoint_eq_galAct.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(galPoint は galAct と一致する。★無条件)",
    sectionId := "genell-thm-3-8" }

def galAct_vcPoint_baseChange.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(galPoint は変数変換と可換。★無条件)",
    sectionId := "genell-thm-3-8" }

def baseChange_variableChange.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(変数変換の底変換。★無条件)",
    sectionId := "genell-thm-3-8" }

def galAct_vcPoint_baseChange.needs : List ProofObligation :=
  [ .citation "[ABC3]" "galAct_vcPoint(第 1288、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.galAct_vcPoint") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1379）**——`C • E` の局所データを `E` に戻す段である。" ++
       "☆第 1372 は `[IsMinimal R (E.baseChange Lv)]` を要求するが、" ++
       "`SemistableAt` が与えるのは `C • E.W` の極小性である。") 19 ]

end ABC3.Found.GenEll
