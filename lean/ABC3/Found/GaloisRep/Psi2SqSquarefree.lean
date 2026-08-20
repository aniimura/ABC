import ABC3.Found.GaloisRep.TwoTorsionDiscr
import Mathlib.FieldTheory.IsAlgClosed.Basic

/-!
# Galois (G5) 第 131 ブロック —— **★★★★★★`Ψ₂Sq` は squarefree**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★欠けていた 1 件を埋めた

§9-441 で在庫を測った結果、`discr ≠ 0` から `Squarefree Ψ₂Sq` を出すのに
**欠けているのは `roots.Nodup ⟹ Squarefree` の向きだけ**と分かった。★それを埋める。

## ★★★★機構

分解可能で根が相異なるなら、`p = C(lc)·∏_{r ∈ roots}(X − r)` と書ける。
★`separable_prod_X_sub_C_iff'`(mathlib)は「相異なる根の積は separable」を与え、
`Separable.squarefree` で squarefree になる。
★★`C(lc)` は単元なので `squarefree_mul_iff` で落ちる。

## ★★★これで `Ψ₂Sq` が squarefree になる

| 段 | 出所 |
|---|---|
| `discr = 16Δ ≠ 0` | mathlib(`twoTorsionPolynomial_discr_ne_zero`) |
| `discr ≠ 0 ⟹ roots.Nodup` | mathlib(`Cubic.discr_ne_zero_iff_roots_nodup`) |
| **`roots.Nodup ⟹ Squarefree`** | **本ブロック** |

★★★★これで `IsIntegrallyClosed` への経路 2 は
「`F[X][√d]`(`d` squarefree)が整閉」の 1 段を残すのみになった。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `squarefree_of_roots_nodup` | ★★★★**根が相異なる分解可能多項式は squarefree** |
| `squarefree_Psi2Sq` | ★★★★★★**`Ψ₂Sq` は squarefree** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial

variable {F : Type} [Field F] [DecidableEq F]

/-- ★★★★**根が相異なる分解可能多項式は squarefree** である。

★mathlib には `Separable ⟹ roots.Nodup`(`Polynomial.nodup_roots`)はあるが、
**逆向きは無い**(2026-08-20 実測)。 -/
theorem squarefree_of_roots_nodup {p : Polynomial F} (hp : p ≠ 0)
    (hsplit : p.Splits) (hnd : p.roots.Nodup) : Squarefree p := by
  have hval : p.roots.toFinset.val = p.roots := by
    rw [Multiset.toFinset_val, Multiset.Nodup.dedup hnd]
  have hsep : (∏ r ∈ p.roots.toFinset, (Polynomial.X - Polynomial.C r)).Separable :=
    (Polynomial.separable_prod_X_sub_C_iff' (f := fun r : F => r)
      (s := p.roots.toFinset)).2 (fun x _ y _ h => h)
  have hprodeq : (∏ r ∈ p.roots.toFinset, (Polynomial.X - Polynomial.C r))
      = (Multiset.map (fun x => Polynomial.X - Polynomial.C x) p.roots).prod := by
    rw [Finset.prod_eq_multiset_prod, hval]
  have hsqprod : Squarefree
      ((Multiset.map (fun x => Polynomial.X - Polynomial.C x) p.roots).prod) := by
    rw [← hprodeq]
    exact hsep.squarefree
  have hlc : IsUnit (Polynomial.C p.leadingCoeff) :=
    Polynomial.isUnit_C.2 (Polynomial.leadingCoeff_ne_zero.2 hp).isUnit
  rw [Polynomial.Splits.eq_prod_roots hsplit, squarefree_mul_iff]
  exact ⟨IsUnit.isRelPrime_left hlc, hlc.squarefree, hsqprod⟩

/-- ★★★★★★**`Ψ₂Sq` は squarefree である**(楕円曲線・代数閉・標数 ≠ 2)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`discr = 16Δ ≠ 0`(mathlib)から根が相異なり、本ブロックの補題で squarefree になる。
★★これが `IsIntegrallyClosed W.CoordinateRing` への経路 2 の中核である。 -/
theorem squarefree_Psi2Sq (W : WeierstrassCurve F) [W.IsElliptic] [IsAlgClosed F]
    (h2 : IsUnit (2 : F)) : Squarefree W.Ψ₂Sq := by
  have hPsi : W.Ψ₂Sq = W.twoTorsionPolynomial.toPoly := Psi2Sq_eq_twoTorsionPolynomial W
  have h4 : (4 : F) ≠ 0 := by
    have h44 : (4 : F) = 2 * 2 := by norm_num
    rw [h44]
    exact mul_ne_zero h2.ne_zero h2.ne_zero
  have ha : W.twoTorsionPolynomial.a ≠ 0 := h4
  have hsplit : (Polynomial.map (RingHom.id F) W.twoTorsionPolynomial.toPoly).Splits := by
    rw [Polynomial.map_id]
    exact IsAlgClosed.splits _
  have hdiscr : W.twoTorsionPolynomial.discr ≠ 0 :=
    twoTorsionPolynomial_discr_ne_zero_of_isElliptic W h2
  have hnd : (Cubic.map (RingHom.id F) W.twoTorsionPolynomial).roots.Nodup :=
    (Cubic.discr_ne_zero_iff_roots_nodup ha hsplit).1 hdiscr
  have hmapeq : Cubic.map (RingHom.id F) W.twoTorsionPolynomial = W.twoTorsionPolynomial := by
    ext <;> rfl
  rw [hmapeq] at hnd
  have hne0 : W.Ψ₂Sq ≠ 0 := by
    intro h
    refine h4 ?_
    have hc := W.coeff_Ψ₂Sq
    rw [h, Polynomial.coeff_zero] at hc
    exact hc.symm
  rw [hPsi]
  refine squarefree_of_roots_nodup ?_ (IsAlgClosed.splits _) hnd
  rw [← hPsi]
  exact hne0

/-! ## ★出典の紐付け(`.src`) -/

def squarefree_Psi2Sq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——Ψ₂Sq が squarefree であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
