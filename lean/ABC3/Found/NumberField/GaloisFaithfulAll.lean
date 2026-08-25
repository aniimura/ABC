/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NumberField.GaloisFaithfulBase
import ABC3.Found.NumberField.SplitExponent

/-!
# 全素イデアルを固定する数体の自己同型は恒等（★仮定なし）

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.115。

原文 (FrdI p.115):
> since any automorphism of a number field that fixes all of the valuations of the
> number field is clearly equal to the identity automorphism

## ★★原文の `clearly` の中身

原文はここを `clearly` で畳んでいるが、実際には

1. `L` を含む **ℚ 上 Galois** な数体 `N` を取り（Galois 閉包）、
2. `N` で**完全分解する有理素数 `p`** を取り、
3. `p` の上の素イデアル `R` の分解群が自明であることから、
   `R ∩ 𝓞 L` を固定する自己同型は恒等である

という 3 段を要する。★★2 は普通は **Chebotarev の密度定理**で出すが、
本プロジェクトは `infinite_splitsCompletely`（Schur ＋ Dedekind–Kummer ＋
Galois の基本等式）で**迂回している**。

## ★★★本ファイルで最後の仮定が落ちた

`GaloisFaithfulBase.lean` の `eq_one_of_fixes_prime_of_galoisClosure` は
`hsplit`（完全分解する素数）と `R` を**仮定として受けて**いた。
★`SplitExponent.lean` の `exists_splitsCompletely_prime`（仮定なし）と
本ファイルの `exists_prime_over_ne_bot` を与えることで、

    eq_one_of_fixes_all_primes : 全素イデアルを固定する自己同型は恒等

が**仮定なし**で立つ。これが `Theorem 6.4, (i)` の Div-slim の数論の側である。
-/

namespace ABC3.Found.NF

open _root_.NumberField Ideal
open scoped _root_.NumberField Pointwise

/-- ★**有理素数の上には `0` でない素イデアルがある**。

★在庫の `exists_prime_over`(`SplCompositum.lean`)は節変数 `Ω` 用で `≠ ⊥` を言わないので、
★ここでは体を明示し `≠ ⊥` まで込めた版を別名で置く。 -/
theorem exists_prime_over_ne_bot (N : Type) [Field N] [NumberField N] (p : ℕ) [hp : Fact p.Prime] :
    ∃ R : Ideal (𝓞 N), R.IsPrime ∧ R.LiesOver (span {(p : ℤ)}) ∧ R ≠ ⊥ := by
  haveI : (span {(p : ℤ)}).IsPrime := by
    rw [Ideal.span_singleton_prime (by exact_mod_cast hp.out.ne_zero)]
    exact Nat.prime_iff_prime_int.mp hp.out
  obtain ⟨⟨R, hR1, hR2⟩⟩ := (inferInstance : Nonempty (primesOver (span {(p : ℤ)}) (𝓞 N)))
  refine ⟨R, hR1, hR2, ?_⟩
  intro hbot
  have hover := hR2.over
  rw [hbot] at hover
  rw [Ideal.under, Ideal.comap_bot_of_injective] at hover
  · exact hp.out.ne_zero (by exact_mod_cast (Ideal.span_singleton_eq_bot.mp hover))
  · exact RingHom.injective_int _

/-- ★★★★★★**全素イデアルを固定する数体の自己同型は恒等**（★仮定なし）。

原文 (FrdI p.115):
> since any automorphism of a number field that fixes all of the valuations of the
> number field is clearly equal to the identity automorphism

★★★**Chebotarev の密度定理は 1 度も使わない**。 -/
theorem eq_one_of_fixes_all_primes (L : Type) [Field L] [NumberField L]
    (τ : L ≃ₐ[ℚ] L)
    (hτ : ∀ 𝔭 : Ideal (𝓞 L), 𝔭.IsPrime → 𝔭 ≠ ⊥ → τ • 𝔭 = 𝔭) : τ = 1 := by
  obtain ⟨N, _, _, _, _⟩ := exists_galoisClosure L
  obtain ⟨p, hp, hsplit⟩ := exists_splitsCompletely_prime N
  haveI : Fact p.Prime := ⟨hp⟩
  obtain ⟨R, hR1, hR2, hRne⟩ := exists_prime_over_ne_bot N p
  haveI := hR1
  haveI := hR2
  refine eq_one_of_fixes_prime_of_galoisClosure hsplit R hRne τ (hτ _ inferInstance ?_)
  exact Ideal.IsIntegral.comap_ne_bot (𝓞 L) hRne

/-! ### ★出典の紐付け -/

def exists_prime_over_ne_bot.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 115,
    item := "Theorem 6.4, (i) — 有理素数の上には 0 でない素イデアルがある",
    sectionId := "frdi-thm-6-4" }

def eq_one_of_fixes_all_primes.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 115,
    item := "Theorem 6.4, (i) — 全素イデアルを固定する数体の自己同型は恒等(仮定なし)",
    sectionId := "frdi-thm-6-4" }

def eq_one_of_fixes_all_primes.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "eq_one_of_fixes_prime_of_galoisClosure(分解群が自明)"
      (.inProject "ABC3" "ABC3.Found.NF.eq_one_of_fixes_prime_of_galoisClosure") 115,
    .citation "[ABC3]" "exists_splitsCompletely_prime(完全分解する素数。Chebotarev 不使用)"
      (.inProject "ABC3" "ABC3.Found.NF.exists_splitsCompletely_prime") 116,
    .derivation "原文「clearly equal to the identity automorphism」の中身(3 段)" 115 ]

end ABC3.Found.NF
