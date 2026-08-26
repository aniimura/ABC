/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NumberField.SplitDensity
import ABC3.Found.NumberField.GaloisFaithful

/-!
# 完全分解 ⟺ 分解群が自明(鎖 `cheb` の `cheb-spl-det` の前段)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.116。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

## ★★★見立て違いの訂正(2026-08-25、**8 回目**)

「`decompositionGroup` が mathlib に 0 件」から「分解群が無い」と読むのは誤りだった ——
`MulAction.stabilizer` がそのまま分解群であり、
`Mathlib/NumberTheory/RamificationInertia/Galois.lean` に
`Ideal.card_stabilizer_eq`(`#(stabilizer G P) = e·f`)・
`Ideal.isPretransitive_of_isGaloisGroup`(推移的)がある。

★★さらに **`Found/NumberField/GaloisFaithful.lean` に順方向は既にあった**:

| 在庫 | 中身 |
|---|---|
| `ramificationIdxIn_mul_inertiaDegIn_eq_one` | 完全分解 ⟹ `e·f = 1` |
| `isSeparable_residue` | 剰余体拡大は分離的 |
| `subsingleton_stabilizer_of_splitsCompletely` | 完全分解 ⟹ 分解群は自明 |
| `eq_one_of_fixes_prime` | ★素点を 1 つ固定する自己同型は恒等 |

★教訓は毎回同じ: **「残っている」と書いた条は、着手前に必ず Lean の側を読む。**

## ★本ファイルで閉じること —— **逆向き**

| 定理 | 中身 |
|---|---|
| `card_stabilizer_eq_ef` | `#(stabilizer Gal(L/ℚ) P) = e·f`(等式そのもの) |
| `splitsCompletely_of_stabilizer_trivial` | ★★★**分解群が自明 ⟹ 完全分解** |
| `splitsCompletely_iff_stabilizer_trivial` | ★★★★両向き |

★★逆向きが要るのは `Spl(LM) = Spl(L) ∩ Spl(M)` の「⊇」である ——
`σ ∈ stabilizer_N(𝔓)` を `L`・`M` へ制限して両方 `1` なら `N = L ⊔ M` より `σ = 1`、
ゆえに `N` でも完全分解する。
-/

namespace ABC3.Found.NF

open _root_.NumberField Ideal
open scoped _root_.NumberField Pointwise

attribute [local instance] Ideal.Quotient.field

variable {L : Type*} [Field L] [NumberField L]

theorem mem_primesOverFinset_iff' {p : ℕ} (hp : p.Prime) (P : Ideal (𝓞 L)) :
    P ∈ IsDedekindDomain.primesOverFinset (Ideal.span {(p : ℤ)}) (𝓞 L)
      ↔ P.IsPrime ∧ P.LiesOver (Ideal.span {(p : ℤ)}) := by
  haveI : (Ideal.span {(p : ℤ)}).IsMaximal := span_intCast_isMaximal hp
  rw [← Finset.mem_coe, IsDedekindDomain.coe_primesOverFinset (span_intCast_ne_bot hp) (𝓞 L)]
  exact Iff.rfl

variable [IsGalois ℚ L]

/-- ★★**`#(stabilizer Gal(L/ℚ) P) = e·f`**(等式そのもの)。

★`GaloisFaithful.lean` は「完全分解 ⟹ `= 1`」だけを使うが、
逆向きにはこの等式が要る。 -/
theorem card_stabilizer_eq_ef {p : ℕ} (hp : p.Prime) (P : Ideal (𝓞 L)) (hP : P.IsPrime)
    (hLO : P.LiesOver (Ideal.span {(p : ℤ)})) :
    Nat.card (MulAction.stabilizer (L ≃ₐ[ℚ] L) P)
      = ramificationIdx (Ideal.span {(p : ℤ)}) P * inertiaDeg (Ideal.span {(p : ℤ)}) P := by
  haveI : Fact p.Prime := ⟨hp⟩
  haveI : (Ideal.span {(p : ℤ)}).IsMaximal := span_intCast_isMaximal hp
  haveI := hP
  haveI := hLO
  have hPne : P ≠ ⊥ := ne_bot_of_liesOver_span hp P hLO
  haveI : P.IsMaximal := Ideal.IsPrime.isMaximal hP hPne
  haveI := isSeparable_residue (p := p) P hPne
  have h := Ideal.card_stabilizer_eq (G := (L ≃ₐ[ℚ] L)) (Ideal.span {(p : ℤ)})
    (span_intCast_ne_bot hp) P
  rwa [Ideal.ramificationIdxIn_eq_ramificationIdx _ P (L ≃ₐ[ℚ] L),
    Ideal.inertiaDegIn_eq_inertiaDeg _ P (L ≃ₐ[ℚ] L)] at h

/-- ★★★**分解群が自明なら完全分解する**(`Spl(LM) ⊇ Spl(L) ∩ Spl(M)` に要る)。

★Galois なので `e`・`f` は `p` の上のどの素でも等しい
(`Ideal.ramificationIdx_eq_of_isGaloisGroup`)。 -/
theorem splitsCompletely_of_stabilizer_trivial {p : ℕ} (hp : p.Prime)
    (P : Ideal (𝓞 L)) (hP : P.IsPrime) (hLO : P.LiesOver (Ideal.span {(p : ℤ)}))
    (htriv : ∀ σ : L ≃ₐ[ℚ] L, σ • P = P → σ = 1) :
    (primesOver (Ideal.span {(p : ℤ)}) (𝓞 L)).ncard = Module.finrank ℚ L := by
  haveI : (Ideal.span {(p : ℤ)}).IsMaximal := span_intCast_isMaximal hp
  haveI := hP
  haveI := hLO
  have hcard := card_stabilizer_eq_ef hp P hP hLO
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

/-- ★★★★**完全分解 ⟺ 分解群(素点の固定部分群)が自明**。 -/
theorem splitsCompletely_iff_stabilizer_trivial {p : ℕ} (hp : p.Prime)
    (P : Ideal (𝓞 L)) (hP : P.IsPrime) (hLO : P.LiesOver (Ideal.span {(p : ℤ)})) :
    (primesOver (Ideal.span {(p : ℤ)}) (𝓞 L)).ncard = Module.finrank ℚ L
      ↔ ∀ σ : L ≃ₐ[ℚ] L, σ • P = P → σ = 1 := by
  haveI : Fact p.Prime := ⟨hp⟩
  haveI := hP
  haveI := hLO
  exact ⟨fun hsplit σ hσ => eq_one_of_fixes_prime hsplit P σ hσ,
    splitsCompletely_of_stabilizer_trivial hp P hP hLO⟩

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
