/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NumberField.MaxDeg
import Mathlib.NumberTheory.RamificationInertia.Galois

/-!
# 素点への作用は忠実(鎖 `sec6items` の `thm64-i-nondil` の Div-slim の側)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.114。

原文 (FrdI p.114):
> (Arithmetic Frobenioids) For i = 1, 2, let Fi be a number field;

## ★★原文が「clearly」で畳んだ所

`Theorem 6.4, (i)` は「数体の自己同型で全素点を固定するものは恒等」だから
`𝒟` は **Div-slim** である、と述べる。

★★**中身は完全分解する素数 1 つで足りる** ——
`p` が完全分解するなら分解群 `D_P` は自明（`e·f = 1`）なので、
`P` を固定する自己同型は恒等である。

| 段 | 根拠 |
|---|---|
| `g·(e·f) = |Gal|` | `Ideal.ncard_primesOver_mul_ramificationIdxIn_mul_inertiaDegIn` |
| 完全分解 ⟹ `e·f = 1` | 上の等式と `g = [K:ℚ]` |
| `#(stabilizer) = e·f` | `Ideal.card_stabilizer_eq` |

★完全分解する素数の存在は `SplitSeparable.lean` の
`infinite_splitsCompletely'`（**Chebotarev を使わない**）が与える。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `ramificationIdxIn_mul_inertiaDegIn_eq_one` | 完全分解なら `e·f = 1` |
| `subsingleton_stabilizer_of_splitsCompletely` | 分解群は自明 |
| `eq_one_of_fixes_prime` | ★★**素点を 1 つ固定する自己同型は恒等** |
★★**条つき**: 剰余体拡大の分離性 `Algebra.IsSeparable (ℤ/p) (𝓞 K/P)` を
instance 引数として受けている。★これは**真**である（剰余体は有限、有限体は完全）が、
`Field (ℤ ⧸ span {p})` の instance 探索が届かず、`PerfectField.ofFinite` を
そのまま繋げなかった（2026-08-21 実測）。★中身ではなく配線の問題である。
-/

namespace ABC3.Found.NF

open _root_.NumberField Ideal
open scoped _root_.NumberField Pointwise

variable {K : Type*} [Field K] [NumberField K] [IsGalois ℚ K]

/-! ## ★1. 完全分解なら `e·f = 1` -/

/-- ★★**完全分解する素数では `e·f = 1`**。

★Galois の基本等式 `g·(e·f) = |Gal(K/ℚ)| = [K:ℚ]` と `g = [K:ℚ]` から。 -/
theorem ramificationIdxIn_mul_inertiaDegIn_eq_one {p : ℕ} [Fact p.Prime]
    (hsplit : (primesOver (span {(p : ℤ)}) (𝓞 K)).ncard = Module.finrank ℚ K) :
    ramificationIdxIn (span {(p : ℤ)}) (𝓞 K)
      * inertiaDegIn (span {(p : ℤ)}) (𝓞 K) = 1 := by
  have hG := Ideal.ncard_primesOver_mul_ramificationIdxIn_mul_inertiaDegIn
    (p := span {(p : ℤ)}) (by simp [NeZero.ne p]) (𝓞 K) (K ≃ₐ[ℚ] K)
  rw [IsGalois.card_aut_eq_finrank ℚ K, hsplit] at hG
  have hpos : 0 < Module.finrank ℚ K := Module.finrank_pos
  refine Nat.eq_of_mul_eq_mul_left hpos ?_
  rw [mul_one]
  exact hG

/-! ## ★2. 分解群は自明 -/

/-- ★★★**完全分解する素数の分解群は自明**。

★`Ideal.card_stabilizer_eq` が `#(stabilizer) = e·f` を与えるので、
`e·f = 1` から出る。 -/
theorem subsingleton_stabilizer_of_splitsCompletely {p : ℕ} [Fact p.Prime]
    (hsplit : (primesOver (span {(p : ℤ)}) (𝓞 K)).ncard = Module.finrank ℚ K)
    (P : Ideal (𝓞 K)) [P.IsPrime] [P.LiesOver (span {(p : ℤ)})]
    [hsep : Algebra.IsSeparable (ℤ ⧸ span {(p : ℤ)}) (𝓞 K ⧸ P)] :
    Nat.card (MulAction.stabilizer (K ≃ₐ[ℚ] K) P) = 1 := by
  haveI : (span {(p : ℤ)}).IsMaximal := Int.ideal_span_isMaximal_of_prime p
  have hspan : (span {(p : ℤ)} : Ideal ℤ) ≠ ⊥ := by simp [NeZero.ne p]
  have hPne : P ≠ ⊥ := ne_bot_of_liesOver_of_ne_bot hspan P
  haveI : P.IsMaximal := Ideal.IsPrime.isMaximal inferInstance hPne
  rw [Ideal.card_stabilizer_eq (span {(p : ℤ)}) hspan P]
  exact ramificationIdxIn_mul_inertiaDegIn_eq_one hsplit

/-! ## ★3. 忠実性 -/

/-- ★★★★**完全分解する素数の上の素イデアルを 1 つ固定する自己同型は恒等**。 -/
theorem eq_one_of_fixes_prime {p : ℕ} [Fact p.Prime]
    (hsplit : (primesOver (span {(p : ℤ)}) (𝓞 K)).ncard = Module.finrank ℚ K)
    (P : Ideal (𝓞 K)) [P.IsPrime] [P.LiesOver (span {(p : ℤ)})]
    [Algebra.IsSeparable (ℤ ⧸ span {(p : ℤ)}) (𝓞 K ⧸ P)]
    (σ : K ≃ₐ[ℚ] K) (hσ : σ • P = P) : σ = 1 := by
  have hmem : σ ∈ MulAction.stabilizer (K ≃ₐ[ℚ] K) P := hσ
  have hcard := subsingleton_stabilizer_of_splitsCompletely hsplit P
  haveI : Subsingleton (MulAction.stabilizer (K ≃ₐ[ℚ] K) P) :=
    Nat.card_eq_one_iff_unique.mp hcard |>.1
  exact congrArg Subtype.val (Subsingleton.elim (⟨σ, hmem⟩ : MulAction.stabilizer _ P) 1)

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Theorem 6.4, (i)` の「全素点を固定する自己同型は恒等」。 -/
def eq_one_of_fixes_prime.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 全素点を固定する自己同型は恒等(𝒟 は Div-slim)",
    sectionId := "frdi-thm-6-4" }

end ABC3.Found.NF
