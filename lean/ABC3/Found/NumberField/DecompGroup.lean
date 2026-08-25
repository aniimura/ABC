/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NumberField.SplitDensity
import Mathlib.NumberTheory.RamificationInertia.Galois

/-!
# 完全分解 ⟺ 分解群が自明(鎖 `cheb` の `cheb-spl-det` の前段)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.116。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

## ★★★測り直し —— 分解群は mathlib にある

`decompositionGroup` という名前は無いが、**`MulAction.stabilizer` がそれである**。
`Mathlib/NumberTheory/RamificationInertia/Galois.lean` に

| 在庫 | 中身 |
|---|---|
| `MulAction G (primesOver p B)` | Galois 群の素点への作用 |
| `Ideal.isPretransitive_of_isGaloisGroup` | ★**推移的**(素点は 1 つの軌道) |
| `Ideal.card_stabilizer_eq` | ★★**`#(stabilizer G P) = e·f`** |
| `Ideal.ramificationIdx_eq_of_isGaloisGroup` | Galois なら `e` は素点によらない |

がある。★「`decompositionGroup` が 0 件」から「分解群が無い」と読むのは**誤り**だった。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `finite_int_quotient` | `ℤ/(p)` は有限(`ZMod p` と同型) |
| `isSeparable_residue` | 剰余体の拡大は分離的(有限体は完全体) |
| `card_stabilizer_eq_int` | `#(stabilizer Gal(L/ℚ) P) = e·f`(底が `ℤ`) |
| `splitsCompletely_iff_stabilizer_trivial` | ★★★★**完全分解 ⟺ 分解群が自明** |

★★これが `Spl(LM) = Spl(L) ∩ Spl(M)` の道具である ——
`σ ∈ stabilizer_N(𝔓)` を `L`・`M` へ制限して両方 `1` なら、`N = L ⊔ M` より `σ = 1`。
-/

namespace ABC3.Found.NF

open _root_.NumberField Ideal Pointwise
open scoped _root_.NumberField

variable {L : Type*} [Field L] [NumberField L]

/-! ## ★1. 剰余体は有限体 -/

theorem finite_int_quotient {p : ℕ} (hp : p.Prime) : Finite (ℤ ⧸ Ideal.span {(p : ℤ)}) := by
  haveI : NeZero p := ⟨hp.ne_zero⟩
  exact Finite.of_equiv (ZMod p) (by
    have h := (Int.quotientSpanEquivZMod ((p : ℤ))).symm
    rw [Int.natAbs_natCast] at h
    exact h.toEquiv)

/-- ★剰余体の拡大は分離的 —— 有限体は完全体だから。 -/
theorem isSeparable_residue {p : ℕ} (hp : p.Prime) (P : Ideal (𝓞 L)) (hP : P.IsPrime)
    (hLO : P.LiesOver (Ideal.span {(p : ℤ)})) :
    Algebra.IsSeparable (ℤ ⧸ Ideal.span {(p : ℤ)}) (𝓞 L ⧸ P) := by
  haveI : (Ideal.span {(p : ℤ)}).IsMaximal := span_intCast_isMaximal hp
  haveI := hP
  haveI := hLO
  haveI : P.IsMaximal := hP.isMaximal (ne_bot_of_liesOver_span hp P hLO)
  haveI : Finite (ℤ ⧸ Ideal.span {(p : ℤ)}) := finite_int_quotient hp
  letI : Field (ℤ ⧸ Ideal.span {(p : ℤ)}) := Ideal.Quotient.field _
  letI : Field (𝓞 L ⧸ P) := Ideal.Quotient.field _
  haveI : Algebra.IsAlgebraic (ℤ ⧸ Ideal.span {(p : ℤ)}) (𝓞 L ⧸ P) :=
    Algebra.IsAlgebraic.of_finite _ _
  exact Algebra.IsAlgebraic.isSeparable_of_perfectField
    (K := ℤ ⧸ Ideal.span {(p : ℤ)}) (L := 𝓞 L ⧸ P)

theorem mem_primesOverFinset_iff' {p : ℕ} (hp : p.Prime) (P : Ideal (𝓞 L)) :
    P ∈ IsDedekindDomain.primesOverFinset (Ideal.span {(p : ℤ)}) (𝓞 L)
      ↔ P.IsPrime ∧ P.LiesOver (Ideal.span {(p : ℤ)}) := by
  haveI : (Ideal.span {(p : ℤ)}).IsMaximal := span_intCast_isMaximal hp
  rw [← Finset.mem_coe, IsDedekindDomain.coe_primesOverFinset (span_intCast_ne_bot hp) (𝓞 L)]
  exact Iff.rfl

/-! ## ★2. 分解群の位数は `e·f` -/

variable [IsGalois ℚ L]

/-- ★★**`#(stabilizer Gal(L/ℚ) P) = e·f`**(底が `ℤ` の形)。 -/
theorem card_stabilizer_eq_int {p : ℕ} (hp : p.Prime) (P : Ideal (𝓞 L)) (hP : P.IsPrime)
    (hLO : P.LiesOver (Ideal.span {(p : ℤ)})) :
    Nat.card (MulAction.stabilizer (L ≃ₐ[ℚ] L) P)
      = ramificationIdx (Ideal.span {(p : ℤ)}) P * inertiaDeg (Ideal.span {(p : ℤ)}) P := by
  haveI : (Ideal.span {(p : ℤ)}).IsMaximal := span_intCast_isMaximal hp
  haveI := hP
  haveI := hLO
  haveI : P.IsMaximal := hP.isMaximal (ne_bot_of_liesOver_span hp P hLO)
  haveI := isSeparable_residue hp P hP hLO
  have h := Ideal.card_stabilizer_eq (G := (L ≃ₐ[ℚ] L)) (Ideal.span {(p : ℤ)})
    (span_intCast_ne_bot hp) P
  rwa [Ideal.ramificationIdxIn_eq_ramificationIdx _ P (L ≃ₐ[ℚ] L),
    Ideal.inertiaDegIn_eq_inertiaDeg _ P (L ≃ₐ[ℚ] L)] at h

/-- ★★★★**完全分解 ⟺ 分解群(素点の固定部分群)が自明**。

★→ は `#(stabilizer) = e·f = 1`。
★← は Galois なので `e`・`f` が素点によらないこと(`ramificationIdx_eq_of_isGaloisGroup`)。 -/
theorem splitsCompletely_iff_stabilizer_trivial {p : ℕ} (hp : p.Prime)
    (P : Ideal (𝓞 L)) (hP : P.IsPrime) (hLO : P.LiesOver (Ideal.span {(p : ℤ)})) :
    (primesOver (Ideal.span {(p : ℤ)}) (𝓞 L)).ncard = Module.finrank ℚ L
      ↔ ∀ σ : L ≃ₐ[ℚ] L, σ • P = P → σ = 1 := by
  haveI : (Ideal.span {(p : ℤ)}).IsMaximal := span_intCast_isMaximal hp
  haveI := hP
  haveI := hLO
  have hcard := card_stabilizer_eq_int hp P hP hLO
  constructor
  · intro hsplit σ hσ
    have hef := (ncard_primesOver_eq_finrank_iff_int hp).mp hsplit P
      ((mem_primesOverFinset_iff' hp P).mpr ⟨hP, hLO⟩)
    rw [hef] at hcard
    have hbot : MulAction.stabilizer (L ≃ₐ[ℚ] L) P = ⊥ :=
      Subgroup.eq_bot_of_card_eq _ hcard
    have hmem : σ ∈ MulAction.stabilizer (L ≃ₐ[ℚ] L) P := hσ
    rw [hbot, Subgroup.mem_bot] at hmem
    exact hmem
  · intro htriv
    refine (ncard_primesOver_eq_finrank_iff_int hp).mpr (fun Q hQ => ?_)
    obtain ⟨hQp, hQlo⟩ := (mem_primesOverFinset_iff' hp Q).mp hQ
    haveI := hQp
    haveI := hQlo
    have h1 : ramificationIdx (Ideal.span {(p : ℤ)}) Q
        = ramificationIdx (Ideal.span {(p : ℤ)}) P :=
      Ideal.ramificationIdx_eq_of_isGaloisGroup _ Q P (L ≃ₐ[ℚ] L)
    have h2 : inertiaDeg (Ideal.span {(p : ℤ)}) Q = inertiaDeg (Ideal.span {(p : ℤ)}) P :=
      Ideal.inertiaDeg_eq_of_isGaloisGroup _ Q P (L ≃ₐ[ℚ] L)
    rw [h1, h2, ← hcard]
    have hbot : MulAction.stabilizer (L ≃ₐ[ℚ] L) P = ⊥ :=
      Subgroup.eq_bot_iff_forall _ |>.mpr (fun σ hσ => htriv σ hσ)
    rw [hbot]
    simp

/-- ★`p ∈ SplQ L` の分解群による言い換え。 -/
theorem mem_SplQ_iff_stabilizer_trivial {p : Nat.Primes}
    (P : Ideal (𝓞 L)) (hP : P.IsPrime) (hLO : P.LiesOver (Ideal.span {((p : ℕ) : ℤ)})) :
    p ∈ SplQ L ↔ ∀ σ : L ≃ₐ[ℚ] L, σ • P = P → σ = 1 :=
  splitsCompletely_iff_stabilizer_trivial p.2 P hP hLO

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `cheb-spl-det` の道具(完全分解 ⟺ 分解群が自明)。 -/
def splitsCompletely_iff_stabilizer_trivial.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — 完全分解 ⟺ 分解群が自明",
    sectionId := "frdi-thm-6-4" }

end ABC3.Found.NF
