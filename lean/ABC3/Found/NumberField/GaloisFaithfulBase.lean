/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NumberField.GaloisFaithful

/-!
# 素点への作用は忠実(★**底が一般の数体**の版)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.114。

原文 (FrdI p.115):
> since any automorphism of a number field that fixes all of the valuations of the
> number field is clearly equal to the identity automorphism

## ★★なぜ底を一般化するのか

`GaloisFaithful.lean` は同じ議論を**底 `ℚ`** で書いており、`[IsGalois ℚ K]` を要求する。
★しかし `Theorem 6.4, (i)` の `𝒟 = B(G)⁰` の対象 `L` は
**`F` の有限拡大**であって、**ℚ 上 Galois とは限らない**。

★★そこで「全素点を固定する自己同型 `τ` は恒等」を示すには

    M := L^{⟨τ⟩}

と置いて **`L/M` を巡回 Galois 拡大として見る**のが道筋である
(節点 `thm64-i-divslim` / `nf-spl-base`)。★本ファイルはその
**底が一般の数体の場合**の 3 本を用意する。

## ★★★中身は `ℚ` の場合と同じ

| 段 | 根拠 |
|---|---|
| `g·(e·f) = #Gal` | `Ideal.ncard_primesOver_mul_ramificationIdxIn_mul_inertiaDegIn` |
| 完全分解 ⟹ `e·f = 1` | 上の等式と `g = [L:M]` |
| `#(stabilizer) = e·f` | `Ideal.card_stabilizer_eq`(mathlib、`IsGaloisGroup G R S` だけ要求) |

★`IsGaloisGroup (L ≃ₐ[M] L) (𝓞 M) (𝓞 L)` は **mathlib の instance がそのまま効く**
(底が `ℚ` のときのように `𝓞 ℚ ≅ ℤ` で移す必要が無い)。

## ☆残るもの —— 節点 `nf-spl-base`

本ファイルが要求する `hsplit`(`q` が `L` で完全分解する)の**存在**は、
底が `ℚ` なら在庫の `infinite_splitsCompletely`(★Chebotarev を使わない)が与える。
★一般の底 `M` では

1. `L̃` を `L/ℚ` の Galois 閉包とし、`L̃` で完全分解する有理素数 `p` を取る、
2. `p` は `L` でも完全分解する(`e`・`f` の乗法性)、
3. `P | p`、`q := P ∩ 𝓞 M` とすると `e(P/q) = f(P/q) = 1`

という降下が要る。★そこが節点 `nf-spl-base` である。

## ★実装上の注意(在庫と同じ)

* `Ideal.Quotient.field` は mathlib で**大域 instance ではない**ので
  `attribute [local instance]` が要る。
* `Ideal` への群作用は **`Pointwise` スコープ**にある(`open scoped Pointwise`)。
-/

namespace ABC3.Found.NF

open _root_.NumberField Ideal
open scoped _root_.NumberField Pointwise

-- ★`Ideal.Quotient.field` は mathlib で**大域 instance ではない**(local instance)。
attribute [local instance] Ideal.Quotient.field

variable {M L : Type} [Field M] [Field L] [NumberField M] [NumberField L]
  [Algebra M L] [IsGalois M L]

/-! ## ★1. 完全分解なら `e·f = 1` -/

/-- ★★**完全分解なら `e·f = 1`**(底が一般の数体の版)。

★Galois の基本等式 `g·e·f = #Gal(L/M) = [L:M]` を当てるだけ。 -/
theorem ramificationIdxIn_mul_inertiaDegIn_eq_one_base
    (q : Ideal (𝓞 M)) [q.IsMaximal] (hq : q ≠ ⊥)
    (hsplit : (q.primesOver (𝓞 L)).ncard = Module.finrank M L) :
    q.ramificationIdxIn (𝓞 L) * q.inertiaDegIn (𝓞 L) = 1 := by
  have hG := Ideal.ncard_primesOver_mul_ramificationIdxIn_mul_inertiaDegIn
    (p := q) hq (𝓞 L) (L ≃ₐ[M] L)
  rw [IsGalois.card_aut_eq_finrank M L, hsplit] at hG
  have hpos : 0 < Module.finrank M L := Module.finrank_pos
  refine Nat.eq_of_mul_eq_mul_left hpos ?_
  rw [mul_one]
  exact hG

/-! ## ★2. 剰余体拡大は分離的 -/

omit [IsGalois M L] in
/-- ★★**剰余体拡大は分離的**(底が一般の数体の版)—— 剰余体は有限、有限体は完全だから。 -/
theorem isSeparable_residue_base
    (q : Ideal (𝓞 M)) [q.IsMaximal] (hqne : q ≠ ⊥)
    (P : Ideal (𝓞 L)) [P.IsPrime] [P.LiesOver q] (hPne : P ≠ ⊥) :
    Algebra.IsSeparable (𝓞 M ⧸ q) (𝓞 L ⧸ P) := by
  haveI : P.IsMaximal := Ideal.IsPrime.isMaximal inferInstance hPne
  haveI : Finite (𝓞 M ⧸ q) := Ideal.finiteQuotientOfFreeOfNeBot q hqne
  haveI : Finite (𝓞 L ⧸ P) := Ideal.finiteQuotientOfFreeOfNeBot P hPne
  haveI : PerfectField (𝓞 M ⧸ q) := PerfectField.ofFinite
  haveI : Algebra.IsAlgebraic (𝓞 M ⧸ q) (𝓞 L ⧸ P) := Algebra.IsAlgebraic.of_finite _ _
  infer_instance

/-! ## ★3. 分解群は自明、したがって素点を固定する自己同型は恒等 -/

/-- ★★★**完全分解する素点の分解群は自明**(底が一般の数体の版)。 -/
theorem subsingleton_stabilizer_of_splitsCompletely_base
    (q : Ideal (𝓞 M)) [q.IsMaximal] (hqne : q ≠ ⊥)
    (hsplit : (q.primesOver (𝓞 L)).ncard = Module.finrank M L)
    (P : Ideal (𝓞 L)) [P.IsPrime] [P.LiesOver q] (hPne : P ≠ ⊥) :
    Nat.card (MulAction.stabilizer (L ≃ₐ[M] L) P) = 1 := by
  haveI : P.IsMaximal := Ideal.IsPrime.isMaximal inferInstance hPne
  haveI := isSeparable_residue_base q hqne P hPne
  rw [Ideal.card_stabilizer_eq q hqne P]
  exact ramificationIdxIn_mul_inertiaDegIn_eq_one_base q hqne hsplit

/-- ★★★★**完全分解する素点を 1 つ固定する自己同型は恒等**(底が一般の数体の版)。

★★これが `Theorem 6.4, (i)` の「数体の自己同型で全素点を固定するものは恒等」の
**最後の段**である —— `τ ∈ Aut(L/F)` が全素点を固定するとき
`M := L^{⟨τ⟩}` と置いて本定理を当てればよい。
★残るのは `hsplit`(`q` が `L` で完全分解する)の**存在**だけで、
それが節点 `nf-spl-base` である。 -/
theorem eq_one_of_fixes_prime_base
    (q : Ideal (𝓞 M)) [q.IsMaximal] (hqne : q ≠ ⊥)
    (hsplit : (q.primesOver (𝓞 L)).ncard = Module.finrank M L)
    (P : Ideal (𝓞 L)) [P.IsPrime] [P.LiesOver q] (hPne : P ≠ ⊥)
    (σ : L ≃ₐ[M] L) (hσ : σ • P = P) : σ = 1 := by
  have hmem : σ ∈ MulAction.stabilizer (L ≃ₐ[M] L) P := hσ
  have hcard := subsingleton_stabilizer_of_splitsCompletely_base q hqne hsplit P hPne
  haveI : Subsingleton (MulAction.stabilizer (L ≃ₐ[M] L) P) :=
    Nat.card_eq_one_iff_unique.mp hcard |>.1
  exact congrArg Subtype.val
    (Subsingleton.elim (⟨σ, hmem⟩ : MulAction.stabilizer (L ≃ₐ[M] L) P) 1)

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Theorem 6.4, (i)` の「全素点を固定する自己同型は恒等」
(★**底が一般の数体**の版。`ℚ` 上 Galois とは限らない `L` に当てるため)。 -/
def eq_one_of_fixes_prime_base.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 完全分解する素点を固定する自己同型は恒等(一般の底)",
    sectionId := "frdi-thm-6-4" }

def eq_one_of_fixes_prime_base.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Ideal.card_stabilizer_eq(#stabilizer = e·f)"
      (.inMathlib "Ideal.card_stabilizer_eq") 114,
    .citation "[mathlib]" "IsGaloisGroup (L ≃ₐ[M] L) (𝓞 M) (𝓞 L)(instance)"
      (.inMathlib "IsGaloisGroup") 114,
    .citation "[ABC3]" "ℚ を底とする版(同じ議論)"
      (.inProject "ABC3" "ABC3.Found.NF.eq_one_of_fixes_prime") 114,
    .implicitStep
      ("★仮定 hsplit(q が L で完全分解する)の存在が残る。底が ℚ なら在庫の " ++
       "infinite_splitsCompletely(Chebotarev を使わない)が与えるが、" ++
       "一般の底 M では Galois 閉包経由の降下が要る ＝ 節点 nf-spl-base") 116 ]

end ABC3.Found.NF
