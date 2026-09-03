/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.SigmaBound
import ABC3.Found.GenEll.CartierPullback

/-!
# LogCondSigma —— `[GenEll] Definition 1.5` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open NumberField IsDedekindDomain

/-! ## ★★★根基はモノイドの意味でも radical -/

/-- ★★★**`rad(I)` はモノイドの意味でも `IsRadical`**。

★イデアル論の `IsRadical`(`r^n ∈ I → r ∈ I`)ではなく、
モノイドの `IsRadical`(`x ∣ y^n → x ∣ y`)である。`Squarefree` へ渡すのはこちら。 -/
theorem isRadical_ideal_radical {F : Type*} [Field F] [NumberField F] (I : Ideal (𝓞 F)) :
    IsRadical (I.radical) := by
  intro n K hdvd
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · rw [pow_zero] at hdvd
    exact dvd_trans hdvd (one_dvd K)
  rw [Ideal.dvd_iff_le] at hdvd ⊢
  calc K ≤ K.radical := Ideal.le_radical
    _ = (K ^ n).radical := (Ideal.radical_pow K (by omega)).symm
    _ ≤ (I.radical).radical := Ideal.radical_mono hdvd
    _ = I.radical := Ideal.radical_idem I

theorem radical_ne_zero {F : Type*} [Field F] [NumberField F] {I : Ideal (𝓞 F)} (hI : I ≠ 0) :
    I.radical ≠ 0 := by
  intro h
  apply hI
  have hle : I ≤ I.radical := Ideal.le_radical
  rw [h] at hle
  exact le_bot_iff.1 hle

/-- ★★★★**根基イデアルの分解の重複度は `≤ 1`**。 -/
theorem count_radical_le_one {F : Type*} [Field F] [NumberField F] (I : Ideal (𝓞 F)) (hI : I ≠ 0)
    (v : FinitePlace F) :
    (Associates.mk v.asIdeal).count (Associates.mk I.radical).factors ≤ 1 := by
  classical
  have hrad0 : I.radical ≠ 0 := radical_ne_zero hI
  have hsq : Squarefree (I.radical) := (isRadical_ideal_radical I).squarefree hrad0
  by_contra hcon
  have h2 : 2 ≤ (Associates.mk v.asIdeal).count (Associates.mk I.radical).factors := by omega
  have hirr : Irreducible (Associates.mk v.asIdeal) := by
    rw [Associates.irreducible_mk]
    exact (Ideal.prime_of_isPrime v.ne_bot v.isPrime).irreducible
  have hne : (Associates.mk I.radical) ≠ 0 := by simpa using hrad0
  have hdvd : (Associates.mk v.asIdeal) ^ 2 ≤ Associates.mk I.radical :=
    (Associates.prime_pow_dvd_iff_le hne hirr).2 h2
  have hd2 : v.asIdeal * v.asIdeal ∣ I.radical := by
    rw [← pow_two, ← Associates.mk_le_mk_iff_dvd, Associates.mk_pow]
    exact hdvd
  exact (Ideal.IsPrime.ne_top v.isPrime) (Ideal.isUnit_iff.1 (hsq _ hd2))

/-! ## ★出典の紐付け(`.src`) -/

def count_radical_le_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (ii)((−)_red の係数が ≤ 1 であること)",
    sectionId := "genell-def-1-5" }

end ABC3.Found.GenEll
