/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.AlphaGlobalTransport
import ABC3.Meta.Claim

/-!
# 第 1312 ブロック —— **`ζ_l` は基礎局所体にあれば足りる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★これは何か——**大域の底変換が要らなくなる**

第 1308 は `ζ₀ ∈ L`（大域の体に原始 `l` 乗根がある）を仮定していた。
★しかし必要なのは **`σ` が `ζ` を固定する**ことだけであり、
それは **`ζ` が基礎局所体 `L_v` にある**ことから出る（`σ` は `L_v` を固定するから）。

☆したがって `L_v` を「`L(ζ_l)` の素点での完備化」に取れば、
**`SSCurve` を底変換する必要はまったく無い**——
`restrictLocalHom` は `L` の上で働くからである。

★★★これで第 1311 の入力 (1)（`ζ₀ ∈ L`）は
「`ζ` の像が `L_v` に入る」に弱まる。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Interface.GaloisRep ABC3.Found.GaloisRep
open ABC3.Meta

open scoped Classical

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★
**`ζ` の像が基礎局所体にあれば `σ` はそれを固定する**——★**無条件**（第 1312）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`σM` は `L_v` を点ごとに固定するので、`hcomm` と `ι` の単射性で戻せる。 -/
theorem galEquiv_fixes_of_mem_local
    (L Lv M : Type) [Field L] [Field Lv] [Field M]
    [Algebra L Lv] [Algebra Lv M] [Algebra L M] [IsScalarTower L Lv M]
    (ι : AlgebraicClosure L →ₐ[L] M) (σM : M ≃ₐ[Lv] M)
    (σbar : AlgebraicClosure L ≃ₐ[L] AlgebraicClosure L)
    (hcomm : ∀ x : AlgebraicClosure L, ι (σbar x) = (σM.restrictScalars L) (ι x))
    (ζ : AlgebraicClosure L) (z : Lv) (hz : algebraMap Lv M z = ι ζ) :
    σbar ζ = ζ := by
  have h1 : ι (σbar ζ) = ι ζ := by
    rw [hcomm]
    show σM (ι ζ) = ι ζ
    rw [← hz]
    exact σM.commutes z
  exact ι.injective h1

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**局所の `σ` から大域の `h2`・`h1`（`ζ` は局所体にあれば足りる）**——★**無条件**（第 1312）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆第 1308 の一般化——`ι`・`σbar`・`hcomm` を仮説として受け取る形にした。 -/
theorem exists_h2_h1_global_of_local'
    (L Lv M : Type) [Field L] [CharZero L] [Field Lv] [Field M]
    [Algebra L Lv] [Algebra Lv M] [Algebra L M] [IsScalarTower L Lv M] [IsAlgClosed M]
    (W : WeierstrassCurve L) [W.IsElliptic] (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange (AlgebraicClosure L)) l ≃+ (Fin 2 → ℤ_[l]))
    (ι : AlgebraicClosure L →ₐ[L] M) (σM : M ≃ₐ[Lv] M)
    (σbar : AlgebraicClosure L ≃ₐ[L] AlgebraicClosure L)
    (hcomm : ∀ x : AlgebraicClosure L, ι (σbar x) = (σM.restrictScalars L) (ι x))
    (ζ : AlgebraicClosure L) (hζ : IsPrimitiveRoot ζ l) (hfixζ : σbar ζ = ζ)
    (Pfix : (W.baseChange M).toAffine.Point) (hPfixOrd : addOrderOf Pfix = l)
    (hfixed : galPoint W (σM.restrictScalars L) Pfix = Pfix)
    (Qmov : (W.baseChange M).toAffine.Point) (hQ : l • Qmov = 0)
    (hmoved : galPoint W (σM.restrictScalars L) Qmov ≠ Qmov) :
    (∀ x : tateModule (W.baseChange (AlgebraicClosure L)) l,
        ∃ u : tateModule (W.baseChange (AlgebraicClosure L)) l,
        galTate W l σbar (galTate W l σbar x) + x
          = galTate W l σbar x + galTate W l σbar x + l • u) ∧
      (∃ x : tateModule (W.baseChange (AlgebraicClosure L)) l,
        ∀ u : tateModule (W.baseChange (AlgebraicClosure L)) l,
        galTate W l σbar x ≠ x + l • u) := by
  have hl : l.Prime := Fact.out
  haveI : CharZero M := charZero_of_injective_algebraMap (algebraMap L M).injective
  haveI : CharZero (AlgebraicClosure L) :=
    charZero_of_injective_algebraMap (algebraMap L (AlgebraicClosure L)).injective
  haveI hEF : (W.baseChange (AlgebraicClosure L)).IsElliptic := by
    unfold WeierstrassCurve.baseChange
    infer_instance
  haveI hEM : (W.baseChange M).IsElliptic := by
    unfold WeierstrassCurve.baseChange
    infer_instance
  have hΔF : (W.baseChange (AlgebraicClosure L)).Δ ≠ 0 :=
    (W.baseChange (AlgebraicClosure L)).isUnit_Δ.ne_zero
  have hΔK : (W.baseChange M).Δ ≠ 0 := (W.baseChange M).isUnit_Δ.ne_zero
  have hcF : ∀ k : ℕ, 1 ≤ k → k ≤ l → ((k : AlgebraicClosure L) ≠ 0) :=
    fun k hk _ => Nat.cast_ne_zero.2 (by omega)
  have hcK : ∀ k : ℕ, 1 ≤ k → k ≤ l → ((k : M) ≠ 0) :=
    fun k hk _ => Nat.cast_ne_zero.2 (by omega)
  obtain ⟨P, hPord, hPfix⟩ := exists_galPoint_fixed_of_map W hΔF hΔK ι σbar
    (σM.restrictScalars L) hcomm l hl.pos hcF hcK Pfix hPfixOrd hfixed
  obtain ⟨Q, hQ0, hQne⟩ := exists_galPoint_ne_of_map W hΔF hΔK ι σbar
    (σM.restrictScalars L) hcomm l hl.pos hcF hcK ⟨Qmov, hQ, hmoved⟩
  exact galTate_h2_h1_of_fixed_moved W l e σbar ζ hζ hfixζ P hPord hPfix Q hQ0 hQne

/-! ## ★出典の紐付け(`.src`) -/

def galEquiv_fixes_of_mem_local.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(ζ の像が基礎局所体にあれば σ はそれを固定する。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_h2_h1_global_of_local'.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(局所の σ から大域の h2・h1——ζ は局所体にあれば足りる。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_h2_h1_global_of_local'.needs : List ProofObligation :=
  [ .citation "[ABC3]" "galTate_h2_h1_of_fixed_moved(第 1303、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.galTate_h2_h1_of_fixed_moved") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1312）**——**大域の底変換が要らなくなる**。" ++
       "☆必要なのは `σ` が `ζ` を固定することだけで、それは `ζ` が基礎局所体 `L_v` に" ++
       "あることから出る。★`L_v` を「`L(ζ_l)` の素点での完備化」に取ればよく、" ++
       "`restrictLocalHom` は `L` の上で働くので `SSCurve` を底変換する必要はない。") 3 ]

end ABC3.Found.GenEll
