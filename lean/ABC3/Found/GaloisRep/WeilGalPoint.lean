import ABC3.Found.GaloisRep.WeilGalois
import ABC3.Found.GaloisRep.GalRep

/-!
# Galois (G5) 第 194 ブロック —— **★★★★★★★★★`galPoint` の言葉での Galois 同変性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★半線型輸送を `galPoint` に繋ぐ

第 193 の `weilPairingVal_semi` は **`σ : L ≃+* L` が係数を固定する**という形で述べてある。
★スケルトンと (G3) の表現は `σ : L ≃ₐ[K] L` と `galPoint W σ = Point.map σ.toAlgHom` を使う。

### ★★2 つの橋

| 補題 | 内容 |
|---|---|
| `fixesCoeffs_baseChange` | `σ : L ≃ₐ[K] L` は `W.baseChange L` の係数を固定する(`σ.commutes`) |
| `semiPoint_eq_galPoint` | `semiPoint = galPoint`——どちらも `(x,y) ↦ (σx, σy)` |

★★`semiPoint` は加法公式を直接輸送して作った(第 193)、`galPoint` は mathlib の
`Point.map`。★★★**同じ写像である**ことを `Point.some` の成分比較で言う。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `fixesCoeffs_baseChange` | 基底変換した曲線の係数は `σ` で固定される |
| `semiPoint_eq_galPoint` | ★★★**`semiPoint = galPoint`** |
| `weilPairingVal_galPoint` | ★★★★★★★★★**`σ(e_n(P,Q)) = e_n(σP, σQ)`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine

variable {K L : Type} [Field K] [DecidableEq K] [Field L] [DecidableEq L] [Algebra K L]

/-- `σ : L ≃ₐ[K] L` は基底変換した曲線の係数を固定する。 -/
theorem fixesCoeffs_baseChange (W : WeierstrassCurve K) (σ : L ≃ₐ[K] L) :
    FixesCoeffs (W.baseChange L).toAffine (σ.toRingEquiv) := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  all_goals
    show σ.toRingEquiv (algebraMap K L _) = algebraMap K L _
    exact σ.commutes _

/-- ★★★**`semiPoint = galPoint`**。 -/
theorem semiPoint_eq_galPoint (W : WeierstrassCurve K) [((W.baseChange L).toAffine).IsElliptic]
    (σ : L ≃ₐ[K] L) (P : (W.baseChange L).toAffine.Point) :
    semiPoint (W.baseChange L).toAffine (fixesCoeffs_baseChange W σ) P
      = ABC3.Interface.GaloisRep.galPoint W σ P := by
  match P with
  | 0 => rw [map_zero, map_zero]
  | Point.some x y hns =>
      rw [semiPoint_some]
      show _ = WeierstrassCurve.Affine.Point.map (σ.toAlgHom) (Point.some x y hns)
      rw [WeierstrassCurve.Affine.Point.map_some]
      exact point_some_congr rfl rfl _ _

/-- ★★★★★★★★★**Galois 同変性**——`σ(e_n(P,Q)) = e_n(σP, σQ)`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 193 の半線型同変性を `galPoint` の言葉に翻訳しただけ。 -/
theorem weilPairingVal_galPoint [IsAlgClosed L] [CharZero L] (W : WeierstrassCurve K)
    [((W.baseChange L).toAffine).IsElliptic]
    (n : ℕ) (hn : 1 ≤ n) (σ : L ≃ₐ[K] L) (P Q : (W.baseChange L).toAffine.Point) :
    σ (weilPairingVal (W.baseChange L).toAffine n P Q)
      = weilPairingVal (W.baseChange L).toAffine n
          (ABC3.Interface.GaloisRep.galPoint W σ P)
          (ABC3.Interface.GaloisRep.galPoint W σ Q) := by
  rw [← semiPoint_eq_galPoint W σ P, ← semiPoint_eq_galPoint W σ Q]
  exact weilPairingVal_semi (W.baseChange L).toAffine (Ne.isUnit (two_ne_zero (α := L)))
    (fixesCoeffs_baseChange W σ) hn P Q

/-! ## ★出典の紐付け(`.src`) -/

def weilPairingVal_galPoint.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の Galois 同変性)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
