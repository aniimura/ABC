/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.FixedVecFromPoint
import ABC3.Found.GenEll.DetModOne
import ABC3.Meta.Claim

/-!
# 第 1296 ブロック —— **`h2` は「`ζ_l` と固定点」だけから出る**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★★★★★これは何か——`h2` の**最短の形**

第 1292-1295 を繋ぐ:

| 入力 | 出どころ |
|---|---|
| `σ` が原始 `l` 乗根 `ζ` を固定 | `ζ_l ∈ L`（大域の体に添加してよい） |
| `σ` が位数 `l` の点 `P` を固定 | `μ_l ⊂ E(K₀)`（Tate 一意化、基礎局所体の上だけでよい） |

☆この 2 つから `det = 1`（第 1295）と固定ベクトル（第 1294）が出て、
第 1292（`2×2` の Cayley–Hamilton）で行列が幂単、第 1293 で `h2` になる。

★★★**Tate 一意化の同変性も局所体の拡大も要らない**——
`alpha_mem_map_of_galTate`（第 1237）の `h2` がこれで完全に出る。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve ABC3.Interface.GaloisRep ABC3.Found.GaloisRep ABC3.Meta

open scoped Matrix

variable {K L : Type} [Field K] [DecidableEq K] [CharZero K]
  [Field L] [DecidableEq L] [Algebra K L] [IsAlgClosed L]

set_option maxHeartbeats 800000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**`σ` が `ζ_l` と位数 `l` の点を固定すれば `T_l E` で `mod l` 幂単**——★**無条件**（第 1296）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆第 1292（`det = 1` ＋ 固定ベクトル ⇒ 幂単）に第 1294・1295 を渡し、第 1293 で降ろす。 -/
theorem galTate_unipotent_of_fixed_root_point
    (W : WeierstrassCurve K) [WeierstrassCurve.IsElliptic (W.baseChange L).toAffine]
    (l : ℕ) [Fact l.Prime] (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l]))
    (σ : L ≃ₐ[K] L)
    (ζ : L) (hζ : IsPrimitiveRoot ζ l) (hfixζ : σ ζ = ζ)
    (P : (W.baseChange L).toAffine.Point) (hP : addOrderOf P = l)
    (hfixP : galPoint W σ P = P) :
    ∀ x : tateModule (W.baseChange L) l, ∃ u : tateModule (W.baseChange L) l,
      galTate W l σ (galTate W l σ x) + x
        = galTate W l σ x + galTate W l σ x + l • u := by
  haveI : CharZero L := charZero_of_injective_algebraMap (algebraMap K L).injective
  obtain ⟨v, hv, hAv⟩ := exists_fixed_vec_of_galPoint_eq W l e σ P hP hfixP
  have hdet := det_glRed_eq_one_of_fixed_root W l e σ ζ hζ hfixζ
  have hA := sq_sub_one_eq_zero_of_det_one_of_fixed
    ((glRedPadic l (galRep W l e σ) : GL (Fin 2) (ZMod l)) :
      Matrix (Fin 2) (Fin 2) (ZMod l)) hdet v hv hAv
  exact galTate_unipotent_of_matrix W l e σ hA

/-! ## ★出典の紐付け(`.src`) -/

def galTate_unipotent_of_fixed_root_point.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(σ が ζ_l と位数 l の点を固定すれば T_l E で mod l 幂単。★無条件)",
    sectionId := "genell-thm-3-8" }

def galTate_unipotent_of_fixed_root_point.needs : List ProofObligation :=
  [ .citation "[ABC3]" "sq_sub_one_eq_zero_of_det_one_of_fixed(第 1292、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.sq_sub_one_eq_zero_of_det_one_of_fixed") 1,
    .citation "[ABC3]" "exists_fixed_vec_of_galPoint_eq(第 1294、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_fixed_vec_of_galPoint_eq") 1,
    .citation "[ABC3]" "det_glRed_eq_one_of_fixed_root(第 1295、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.det_glRed_eq_one_of_fixed_root") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1296）**——`h2` の**最短の形**である。" ++
       "☆入力は「`ζ_l` を固定」と「位数 `l` の点を固定」の 2 つだけ。" ++
       "★★★**Tate 一意化の同変性も局所体の拡大も要らない**。") 3 ]

end ABC3.Found.GenEll
