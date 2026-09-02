/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.GalTateUnipFull
import ABC3.Found.GaloisRep.TateUnipotent
import ABC3.Meta.Claim

/-!
# 第 1303 ブロック —— **`h2` と `h1` の一括**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★★★★★これは何か——`alpha_mem_map_of_galTate` の入力そのもの

第 1237（`alpha_mem_map_of_galTate`）が要求する `h2`・`h1` を、
**1 つの `σ` について同時に**出す形にまとめる。

☆受け取るのは 3 つだけ:

| # | 仮説 | Tate 曲線での出どころ |
|---|---|---|
| 1 | `σ` が原始 `l` 乗根 `ζ` を固定 | `ζ_l ∈ L`（大域で添加してよい） |
| 2 | `σ` が位数 `l` の点を固定 | `μ_l ⊂ E(K₀)`（第 1297）＋ 第 1302 で大域へ |
| 3 | `σ` が動かす `l`-捩れ点がある | `E[l] ⊄ E(K₀)`（第 1279・1298）＋ 第 1271 で大域へ |

★★★これで `α` が像に入る段（II 側）は**この 3 つの入力だけ**に帰着した。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Interface.GaloisRep ABC3.Found.GaloisRep
open ABC3.Meta

open scoped Classical

variable {K L : Type} [Field K] [DecidableEq K] [CharZero K]
  [Field L] [DecidableEq L] [Algebra K L] [IsAlgClosed L]

set_option maxHeartbeats 800000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**`h2` と `h1` を 1 つの `σ` について同時に出す**——★**無条件**（第 1303）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆第 1296（`h2`）と第 1270（`h1`）を並べただけである。 -/
theorem galTate_h2_h1_of_fixed_moved
    (W : WeierstrassCurve K) [WeierstrassCurve.IsElliptic (W.baseChange L).toAffine]
    (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) (σ : L ≃ₐ[K] L)
    (ζ : L) (hζ : IsPrimitiveRoot ζ l) (hfixζ : σ ζ = ζ)
    (Pfix : (W.baseChange L).toAffine.Point) (hPfix : addOrderOf Pfix = l)
    (hfixed : galPoint W σ Pfix = Pfix)
    (Qmov : (W.baseChange L).toAffine.Point) (hQ : l • Qmov = 0)
    (hmoved : galPoint W σ Qmov ≠ Qmov) :
    (∀ x : tateModule (W.baseChange L) l, ∃ u : tateModule (W.baseChange L) l,
        galTate W l σ (galTate W l σ x) + x
          = galTate W l σ x + galTate W l σ x + l • u) ∧
      (∃ x : tateModule (W.baseChange L) l, ∀ u : tateModule (W.baseChange L) l,
        galTate W l σ x ≠ x + l • u) := by
  haveI : CharZero L := charZero_of_injective_algebraMap (algebraMap K L).injective
  refine ⟨galTate_unipotent_of_fixed_root_point W l e σ ζ hζ hfixζ Pfix hPfix hfixed, ?_⟩
  exact exists_galTate_ne_of_galPoint W l σ Qmov hQ hmoved

/-! ## ★出典の紐付け(`.src`) -/

def galTate_h2_h1_of_fixed_moved.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(h2 と h1 を 1 つの σ について同時に出す。★無条件)",
    sectionId := "genell-thm-3-8" }

def galTate_h2_h1_of_fixed_moved.needs : List ProofObligation :=
  [ .citation "[ABC3]" "galTate_unipotent_of_fixed_root_point(第 1296、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.galTate_unipotent_of_fixed_root_point") 1,
    .citation "[ABC3]" "exists_galTate_ne_of_galPoint(第 1270、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_galTate_ne_of_galPoint") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1303）**——`alpha_mem_map_of_galTate`（第 1237）の" ++
       "入力そのものである。☆受け取るのは 3 つだけ:" ++
       "(1) `σ` が `ζ_l` を固定、(2) `σ` が位数 `l` の点を固定、(3) `σ` が動かす `l`-捩れ点がある。" ++
       "★★★これで II 側はこの 3 つの入力だけに帰着した。") 3 ]

end ABC3.Found.GenEll
