/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.LocalDataVc
import ABC3.Found.GaloisRep.GalRepWitness
import Mathlib.RingTheory.RootsOfUnity.AlgebraicallyClosed
import ABC3.Meta.Claim

/-!
# 第 1381 ブロック —— **`M`・`ι`・`ζ`・`z` を消す**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★これは何か——最後の配管

第 1380（`exists_h2_h1_of_bad_prime_ram_vc`）は
`M`（代数閉体）・`ι`（`L̄ ↪ M`）・`ζ`（原始 `l` 乗根）・`z`（その局所体での持ち上げ）・
`eq`（`T_l E ≃ ℤ_l²`）を**入力**として受ける。
★本ブロックはそれらを**すべて構成して消す**。

☆道:

| 入力 | 構成 |
|---|---|
| `M` | `AlgebraicClosure L_v`（`IsGalois L_v M` は標数 0 なので自動） |
| `ι` | `IsAlgClosed.lift` |
| `ζ` | `HasEnoughRootsOfUnity`（代数閉体は原始根を持つ、在庫） |
| `z` | ★`ι ζ` も `φ ζ_v` も `M` の原始 `l` 乗根なので、`ι ζ = (φ ζ_v)^i`。`z ≔ ζ_v^i` |
| `eq` | `nonempty_tate_basis`（在庫） |

★★★これで残るのは**分裂性**（`HasSplitMultiplicativeReduction`）だけになる。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine IsDedekindDomain NumberField
open ABC3.Interface.GaloisRep ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**分裂乗法還元から `h2`・`h1` が出る**——★**無条件**（第 1381）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★第 1380 の `M`・`ι`・`ζ`・`z`・`eq` をすべて構成して消した形である。
☆残る仮説は**分裂乗法還元**と `hpe`・`l ∤ e`・`ζ_l ∈ L_v` だけ。 -/
theorem exists_h2_h1_of_split_local
    (L : Type) [Field L] [NumberField L] {Lv : Type} [Field Lv] [CharZero Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {e : ℕ} (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (W : WeierstrassCurve L) [W.IsElliptic] (C : WeierstrassCurve.VariableChange L)
    [(C • W).IsElliptic]
    [((C • W).baseChange Lv).IsElliptic] [WeierstrassCurve.IsMinimal R ((C • W).baseChange Lv)]
    (h : WeierstrassCurve.HasSplitMultiplicativeReduction R ((C • W).baseChange Lv))
    (hj : jExp p W < 0) {l : ℕ} [Fact l.Prime] (hle : ¬ (l ∣ e))
    (hcop : ¬ ((l : ℤ) ∣ jExp p W))
    {ζv : Lvˣ} (hζv : IsPrimitiveRoot ((ζv : Lv)) l) :
    ∃ σbar : AlgebraicClosure L ≃ₐ[L] AlgebraicClosure L,
      (∀ x : tateModule (W.baseChange (AlgebraicClosure L)) l,
          ∃ u : tateModule (W.baseChange (AlgebraicClosure L)) l,
          galTate W l σbar (galTate W l σbar x) + x
            = galTate W l σbar x + galTate W l σbar x + l • u) ∧
        (∃ x : tateModule (W.baseChange (AlgebraicClosure L)) l,
          ∀ u : tateModule (W.baseChange (AlgebraicClosure L)) l,
          galTate W l σbar x ≠ x + l • u) := by
  have hl : l.Prime := Fact.out
  haveI : NeZero l := ⟨hl.ne_zero⟩
  haveI hlL : NeZero ((l : ℕ) : L) := ⟨Nat.cast_ne_zero.2 hl.ne_zero⟩
  -- ★`M` を作る
  haveI : IsGalois Lv (AlgebraicClosure Lv) := ⟨⟩
  haveI : WeierstrassCurve.IsElliptic
      (((W.baseChange Lv).baseChange (AlgebraicClosure Lv)).toAffine) :=
    isElliptic_baseChange_affine (W.baseChange Lv) inferInstance
  -- ★`ι` を作る
  let ι : AlgebraicClosure L →ₐ[L] AlgebraicClosure Lv :=
    IsAlgClosed.lift (R := L) (S := AlgebraicClosure L) (M := AlgebraicClosure Lv)
  -- ★`ζ` を作る
  obtain ⟨ζ, hζ⟩ := HasEnoughRootsOfUnity.exists_primitiveRoot (AlgebraicClosure L) l
  -- ★`z` を作る
  have hprim : IsPrimitiveRoot (algebraMap Lv (AlgebraicClosure Lv) ((ζv : Lv))) l :=
    hζv.map_of_injective (algebraMap Lv (AlgebraicClosure Lv)).injective
  have hpow : (ι ζ) ^ l = 1 := by
    rw [← map_pow, hζ.pow_eq_one, map_one]
  obtain ⟨i, -, hzi⟩ := hprim.eq_pow_of_pow_eq_one hpow
  -- ★`eq` を作る
  obtain ⟨eq⟩ := nonempty_tate_basis (L := AlgebraicClosure L) W inferInstance l
  exact exists_h2_h1_of_bad_prime_ram_vc L (M := AlgebraicClosure Lv) p hpe W C h hj hle hcop
    hζv eq ι ζ hζ (((ζv : Lv)) ^ i) (by rw [map_pow]; exact hzi)

/-! ## ★出典の紐付け(`.src`) -/

def exists_h2_h1_of_split_local.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(分裂乗法還元から h2・h1 が出る。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_h2_h1_of_split_local.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_h2_h1_of_bad_prime_ram_vc(第 1380、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_h2_h1_of_bad_prime_ram_vc") 1,
    .citation "[ABC3]" "nonempty_tate_basis(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.nonempty_tate_basis") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1381）**——`M`・`ι`・`ζ`・`z`・`eq` をすべて構成して消した。" ++
       "☆`z` は「`ι ζ` も `φ ζ_v` も `M` の原始 `l` 乗根なので `ι ζ = (φ ζ_v)^i`」から出る。" ++
       "★★★これで残るのは**分裂性**だけである。") 19 ]

end ABC3.Found.GenEll
