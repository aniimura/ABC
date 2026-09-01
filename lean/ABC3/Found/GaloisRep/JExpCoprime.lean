/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.DegInfBaseChange
import ABC3.Meta.Claim

/-!
# 第 1208 ブロック —— **`l` と局所高さの素性は分岐指数が `l` と素なら底変換で保たれる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——第 1199 の最後の仮説

第 1199（`lemma_3_5_height_ineq_over_extension`）の仮説のうち

    ∀ P, jExp P (E ⊗ L'') < 0 → ¬ (l ∣ jExp P (E ⊗ L''))

だけが、`L` の側の `PrimeToLocalHeights` から**自動では出ない**。

☆`jExp_baseChange`（在庫）が `v_P(j) = e(P|p) · v_p(j)` を与えるので、
`l` が素数なら `l ∤ e(P|p)` と `l ∤ v_p(j)` から出る。

## ★★★★測ったこと——**新しい義務が 1 つ見えた**

したがって残るのは **`l ∤ e(P|p)`** である。原文の道では
`L'' = L(H)` の次数が `l−1` を割るのでこれは自動だが、
本形式化ではまだその段（`Gal(L(H)/L) ↪ 𝔽_l^×`）を取っていない。

☆`M` は `{k • Q}` の座標で生成するので `Gal`-安定（第 1205 の直線の安定性）であり、
標数 0 だから `M/L` は Galois、よって `e ∣ [M:L]` である。
★そこから `[M:L] ∣ l−1` を出せば `l ∤ e` が従う——**これが次の義務**である。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField ABC3.Meta

variable (L L' : Type) [Field L] [NumberField L] [Field L'] [NumberField L']
  [Algebra L L'] [IsScalarTower ℚ L L']

/-- ★★★★★★★★★★★★
**分岐指数が `l` と素なら `l ∤ v_P(j)` は底変換で保たれる**——★**無条件**（第 1208）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`jExp_baseChange`（在庫）が `v_P(j) = e(P|p) · v_p(j)` を与えるので、
`l` が素数なら積の一方を割ることになり、両方の仮定に反する。

★★★これが第 1199 の最後の仮説
（`∀ P, jExp P < 0 → ¬ (l ∣ jExp P)`）を `L` の側へ降ろす段である。 -/
theorem not_dvd_jExp_baseChange (p : HeightOneSpectrum (𝓞 L))
    (P : HeightOneSpectrum (𝓞 L')) [P.asIdeal.LiesOver p.asIdeal]
    (E : WeierstrassCurve L) [E.IsElliptic]
    {l : ℕ} (hl : l.Prime)
    (he : ¬ ((l : ℤ) ∣ (p.asIdeal.ramificationIdx P.asIdeal : ℤ)))
    (hj : ¬ ((l : ℤ) ∣ jExp p E)) :
    ¬ ((l : ℤ) ∣ jExp P (E.baseChange L')) := by
  rw [jExp_baseChange L L' p P E]
  intro hdvd
  have hp : Prime (l : ℤ) := Nat.prime_iff_prime_int.mp hl
  rcases hp.dvd_mul.1 hdvd with h | h
  · exact he h
  · exact hj h

/-! ## ★★★★★★★★★★★★次数が `l` より小さければ `l ∤ e` -/

section RamBound

variable [Algebra (𝓞 L) L'] [IsScalarTower (𝓞 L) L L']
  [IsScalarTower (𝓞 L) (𝓞 L') L'] [Module.Finite (𝓞 L) (𝓞 L')]

/-- ★★★★★★★★★★★★
**`[L':L] < l` なら `l ∤ e(P|p)`**——★**無条件**（第 1209）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`e ≥ 1`（`ramificationIdx_ne_zero_of_liesOver`）かつ
`e ≤ [L':L]`（`ramificationIdx_le_finrank`）なので、`e < l` であり、
`l ∣ e` なら `l ≤ e` となって矛盾する。

★★★これで第 1208 の仮説 `l ∤ e(P|p)` が**次数の不等式ひとつ**に落ちた
——原文の道では `[L(H):L] ∣ l−1 < l` である。 -/
theorem not_dvd_ramificationIdx_of_finrank_lt (p : HeightOneSpectrum (𝓞 L))
    (P : HeightOneSpectrum (𝓞 L')) [P.asIdeal.LiesOver p.asIdeal]
    {l : ℕ} (hdeg : Module.finrank L L' < l) :
    ¬ ((l : ℤ) ∣ (p.asIdeal.ramificationIdx P.asIdeal : ℤ)) := by
  have hne : p.asIdeal.ramificationIdx P.asIdeal ≠ 0 :=
    Ideal.IsDedekindDomain.ramificationIdx_ne_zero_of_liesOver P.asIdeal p.ne_bot
  haveI : p.asIdeal.IsMaximal := p.isMaximal
  have hinj : Function.Injective (algebraMap (𝓞 L) (𝓞 L')) := by
    have h3 : Function.Injective (algebraMap (𝓞 L) L') := by
      rw [IsScalarTower.algebraMap_eq (𝓞 L) L L']
      exact (algebraMap L L').injective.comp (IsFractionRing.injective (𝓞 L) L)
    intro a b hab
    apply h3
    rw [IsScalarTower.algebraMap_eq (𝓞 L) (𝓞 L') L']
    simp only [RingHom.comp_apply, hab]
  haveI : NoZeroSMulDivisors (𝓞 L) (𝓞 L') := by
    refine ⟨fun {c x} hcx => ?_⟩
    rcases eq_or_ne c 0 with rfl | hc
    · exact Or.inl rfl
    · right
      have hmul : algebraMap (𝓞 L) (𝓞 L') c * x = 0 := by
        rwa [Algebra.smul_def] at hcx
      rcases mul_eq_zero.1 hmul with h | h
      · exact absurd (hinj (by rw [h, map_zero])) hc
      · exact h
  have hle : p.asIdeal.ramificationIdx P.asIdeal ≤ Module.finrank L L' :=
    Ideal.ramificationIdx_le_finrank (𝓞 L') L L' P.asIdeal
  intro hdvd
  have hnat : l ∣ p.asIdeal.ramificationIdx P.asIdeal := by exact_mod_cast hdvd
  have hlel := Nat.le_of_dvd (Nat.pos_of_ne_zero hne) hnat
  omega

/-- ★★★★★★★★★★★★★★★★
**`[L':L] < l` なら `l ∤ v_P(j)` は底変換で保たれる**——★**無条件**（第 1209）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★これが第 1199 の最後の仮説を `L` の側へ完全に降ろした形である
——**残るのは `[L'':L] < l` を示すことだけ**になった。 -/
theorem not_dvd_jExp_baseChange_of_finrank_lt (p : HeightOneSpectrum (𝓞 L))
    (P : HeightOneSpectrum (𝓞 L')) [P.asIdeal.LiesOver p.asIdeal]
    (E : WeierstrassCurve L) [E.IsElliptic]
    {l : ℕ} (hl : l.Prime) (hdeg : Module.finrank L L' < l)
    (hj : ¬ ((l : ℤ) ∣ jExp p E)) :
    ¬ ((l : ℤ) ∣ jExp P (E.baseChange L')) :=
  not_dvd_jExp_baseChange L L' p P E hl
    (not_dvd_ramificationIdx_of_finrank_lt L L' p P hdeg) hj

end RamBound

/-! ## ★出典の紐付け(`.src`) -/

def not_dvd_ramificationIdx_of_finrank_lt.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5([L':L] < l なら l ∤ e(P|p)。★無条件)",
    sectionId := "genell-lemma-3-5" }

def not_dvd_jExp_baseChange_of_finrank_lt.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5([L':L] < l なら l ∤ v_P(j) は底変換で保たれる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def not_dvd_jExp_baseChange.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(分岐指数が l と素なら l ∤ v_P(j) は底変換で保たれる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def not_dvd_jExp_baseChange.needs : List ProofObligation :=
  [ .citation "[ABC3]" "jExp_baseChange(証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.jExp_baseChange") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1208）**——第 1199 の最後の仮説を `L` の側へ降ろす段である。" ++
       "☆残るのは **`l ∤ e(P|p)`**。原文の道では `L'' = L(H)` の次数が `l−1` を" ++
       "割るので自動だが、本形式化ではまだ `Gal(L(H)/L) ↪ 𝔽_l^×` を取っていない。" ++
       "★`M` は `{k • Q}` の座標で生成するので `Gal`-安定（第 1205）であり、" ++
       "標数 0 だから `M/L` は Galois、よって `e ∣ [M:L]` である——" ++
       "そこから `[M:L] ∣ l−1` を出せばよい。") 3 ]

end ABC3.Found.GaloisRep
