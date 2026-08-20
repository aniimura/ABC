/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import Mathlib.NumberTheory.RamificationInertia.Basic
import Mathlib.NumberTheory.RamificationInertia.Inertia
import Mathlib.NumberTheory.RamificationInertia.Ramification
import Mathlib.NumberTheory.NumberField.Basic
import Mathlib.RingTheory.DedekindDomain.Ideal.Lemmas
import Mathlib.NumberTheory.RamificationInertia.Galois
import Mathlib.FieldTheory.Galois.Basic

/-!
# 素点の個数と拡大次数 —— `Theorem 6.4, (iv)` の (b)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.116。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

## ★★何を閉じたか

依存グラフ(`ResearchPaper/frdi-decomposition.json`)の鎖 `cheb` の **`cheb-split-count`** ——
「`p` が `L` で完全分解する ⟺ `deg(L,v) = [L:ℚ]`」の**数え上げの側**である。

★★中身は **`Σ_{P|p} e_P·f_P = [L:K]`**(基本等式、mathlib の
`Ideal.sum_ramification_inertia`)と「`e_P·f_P ≥ 1`」だけである:

* `#{P | P ∣ p} ≤ [L:K]` がつねに成り立つ(`ncard_primesOver_le_finrank`)
* 等号 ⟺ すべての `P` で `e_P·f_P = 1` ⟺ 完全分解(`eq_finrank_iff`)

★**Tchebotarev は要らない** —— 原文が Tchebotarev を引くのは (a) と (c) であって、
(b) は基本等式だけで出る。
-/

namespace ABC3.Found.NF

open NumberField IsDedekindDomain
open scoped NumberField Classical

variable (K L : Type*) [Field K] [Field L] [NumberField K] [NumberField L] [Algebra K L]

/-- ★`e_P · f_P ≥ 1`。 -/
theorem one_le_ramificationIdx_mul_inertiaDeg (p : HeightOneSpectrum (𝓞 K))
    (P : Ideal (𝓞 L)) (hP : P ∈ IsDedekindDomain.primesOverFinset p.asIdeal (𝓞 L)) :
    1 ≤ Ideal.ramificationIdx p.asIdeal P * Ideal.inertiaDeg p.asIdeal P := by
  rw [← Finset.mem_coe, IsDedekindDomain.coe_primesOverFinset p.ne_bot (𝓞 L)] at hP
  haveI : P.IsPrime := hP.1
  haveI : P.LiesOver p.asIdeal := hP.2
  have h2 : 0 < Ideal.inertiaDeg p.asIdeal P := Ideal.inertiaDeg_pos' p.asIdeal P
  have h1 : Ideal.ramificationIdx p.asIdeal P ≠ 0 :=
    Ideal.IsDedekindDomain.ramificationIdx_ne_zero_of_liesOver P p.ne_bot
  exact Nat.one_le_iff_ne_zero.mpr (Nat.mul_ne_zero h1 h2.ne')

/-- ★★**`p` の上にある素イデアルの個数は `[L:K]` 以下**。 -/
theorem ncard_primesOver_le_finrank (p : HeightOneSpectrum (𝓞 K)) :
    (Ideal.primesOver p.asIdeal (𝓞 L)).ncard ≤ Module.finrank K L := by
  have hsum := Ideal.sum_ramification_inertia (𝓞 L) K L p.ne_bot
  rw [← IsDedekindDomain.coe_primesOverFinset p.ne_bot (𝓞 L), Set.ncard_coe_finset, ← hsum]
  have h := Finset.card_nsmul_le_sum (IsDedekindDomain.primesOverFinset p.asIdeal (𝓞 L))
    (fun P => Ideal.ramificationIdx p.asIdeal P * Ideal.inertiaDeg p.asIdeal P) 1
    (fun P hP => one_le_ramificationIdx_mul_inertiaDeg K L p P hP)
  simpa using h

/-- ★★★**完全分解 ⟺ すべての `P` で `e_P·f_P = 1`**。

★★これが原文 p.116 の (b)「`p` が `L` で完全分解 ⟺ `deg(L,v) = [L:ℚ]`」の中身である
—— **Tchebotarev は要らない**。 -/
theorem ncard_primesOver_eq_finrank_iff (p : HeightOneSpectrum (𝓞 K)) :
    (Ideal.primesOver p.asIdeal (𝓞 L)).ncard = Module.finrank K L ↔
      ∀ P ∈ IsDedekindDomain.primesOverFinset p.asIdeal (𝓞 L),
        Ideal.ramificationIdx p.asIdeal P * Ideal.inertiaDeg p.asIdeal P = 1 := by
  have hsum := Ideal.sum_ramification_inertia (𝓞 L) K L p.ne_bot
  rw [← IsDedekindDomain.coe_primesOverFinset p.ne_bot (𝓞 L), Set.ncard_coe_finset]
  constructor
  · intro hcard P hP
    by_contra hne
    have hge : 2 ≤ Ideal.ramificationIdx p.asIdeal P * Ideal.inertiaDeg p.asIdeal P := by
      have := one_le_ramificationIdx_mul_inertiaDeg K L p P hP
      omega
    have hlt : (IsDedekindDomain.primesOverFinset p.asIdeal (𝓞 L)).card
        < ∑ Q ∈ IsDedekindDomain.primesOverFinset p.asIdeal (𝓞 L),
            Ideal.ramificationIdx p.asIdeal Q * Ideal.inertiaDeg p.asIdeal Q := by
      have hle : ∀ Q ∈ IsDedekindDomain.primesOverFinset p.asIdeal (𝓞 L),
          1 ≤ Ideal.ramificationIdx p.asIdeal Q * Ideal.inertiaDeg p.asIdeal Q :=
        fun Q hQ => one_le_ramificationIdx_mul_inertiaDeg K L p Q hQ
      calc (IsDedekindDomain.primesOverFinset p.asIdeal (𝓞 L)).card
          = ∑ _Q ∈ IsDedekindDomain.primesOverFinset p.asIdeal (𝓞 L), 1 := by simp
        _ < _ := Finset.sum_lt_sum hle ⟨P, hP, by omega⟩
    omega
  · intro hall
    rw [← hsum, Finset.sum_congr rfl hall]
    simp

/-! ## ★★Galois の場合 —— 1 つの素で決まる

★★Galois 拡大では `e` も `f` も `p` の上のすべての素で等しいので、
**1 つの素 `P` で `e_P = f_P = 1` なら完全分解する**。
★mathlib の `Ideal.ncard_primesOver_mul_ramificationIdxIn_mul_inertiaDegIn`
（Galois の場合の基本等式 `g·e·f = |G|`）を当てるだけである。

★★★**これが迂回路 `cheb-spl-infinite` の後半**である ——
前半は `PrimeDivisorsOfValues.lean` の Schur の補題、
残るのは「`f` が mod `p` で根を持つ ⇒ `f_P = 1` なる `P` がある」
（Dedekind–Kummer、mathlib の `NumberField/Ideal/KummerDedekind.lean`）だけである。 -/

/-- ★★★★★**Galois 拡大では 1 つの不分岐で剰余次数 1 の素があれば完全分解する**。

★Galois の場合の基本等式 `g·e·f = |Gal(L/K)| = [L:K]` を使うだけ。 -/
theorem ncard_primesOver_eq_finrank_of_unramified_degOne [IsGalois K L]
    (p : HeightOneSpectrum (𝓞 K)) (P : Ideal (𝓞 L))
    [P.IsPrime] [P.LiesOver p.asIdeal]
    (he : Ideal.ramificationIdx p.asIdeal P = 1)
    (hf : Ideal.inertiaDeg p.asIdeal P = 1) :
    (Ideal.primesOver p.asIdeal (𝓞 L)).ncard = Module.finrank K L := by
  haveI : p.asIdeal.IsPrime := p.isPrime
  have hG := Ideal.ncard_primesOver_mul_ramificationIdxIn_mul_inertiaDegIn
    (p := p.asIdeal) p.ne_bot (𝓞 L) (L ≃ₐ[K] L)
  rw [Ideal.ramificationIdxIn_eq_ramificationIdx p.asIdeal P (L ≃ₐ[K] L),
      Ideal.inertiaDegIn_eq_inertiaDeg p.asIdeal P (L ≃ₐ[K] L), he, hf] at hG
  rw [IsGalois.card_aut_eq_finrank K L] at hG
  simpa using hG

/-! ### ★出典の紐付け -/

/-- ★locator —— `Theorem 6.4, (iv)` の (b)(完全分解の判定)。 -/
def ncard_primesOver_eq_finrank_iff.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — 完全分解の判定(基本等式だけで出る)",
    sectionId := "frdi-thm-6-4" }

/-- ★★locator —— Galois の場合の完全分解の判定。 -/
def ncard_primesOver_eq_finrank_of_unramified_degOne.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — Galois なら 1 つの素で完全分解が決まる",
    sectionId := "frdi-thm-6-4" }

end ABC3.Found.NF
