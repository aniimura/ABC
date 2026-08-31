/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.Prop34
import ABC3.Found.GaloisRep.DegInfLocal
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★`Lemma 3.7` の (a) が要る 3 つ（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.18。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

## ★★★★★★★★★★★★★★これは何か

`Skeleton/GenEll/Section3.lean` の `lemma_3_7` の**(a) の証明**は、
`EllModuliData` の欄を 3 つ使う:

| 欄 | 内容 | ★本ファイル |
|---|---|---|
| `faltingsHeight_bddBelow` | `ht^Falt` は下に有界 | ★**証明する**（無条件） |
| `prop_3_4` の第 1・第 2 の `≲` | `deg∞ ≤ 14·ht^Falt + A` | ★**より強い `12` で証明する**（無条件） |
| `primeToLocalHeights_of_lt` | `l·log 2 > d·deg∞` なら `l` は局所高さと素 | ★**証明する** |

★★どれも**受けた仮定ではなく、具体的な量から出る定理**である。

## ★機構

`12·ht^Falt = deg∞ − archSum/d`（`§9-670`）と
`archSum ≤ d·log((2π)¹²M)`（`§9-355` の一様上界、`archSum_le`）から

    `deg∞ = 12·ht^Falt + archSum/d ≤ 12·ht^Falt + log((2π)¹²M)`

★`deg∞ ≥ 0` と合わせて `ht^Falt ≥ −log((2π)¹²M)/12`——★★下に有界である。

★★★局所高さの側は `d·deg∞ = Σ_p v_p(Δ_min)·log N(p) ≥ v_p(Δ_min)·log 2`
（`N(p) ≥ 2`）。`l·log 2 > d·deg∞` なら `v_p(Δ_min) < l` となり、
`l` が素数で `v_p(Δ_min) > 0` なら `l ∤ v_p(Δ_min)`。

## ☆`Lemma 3.7` に残るもの

☆(b)・(c) が要求する**例外集合**（`lcyclicExc`・`noMultRedExc`）。
★これらの Galois-finite 性は `Lemma 3.5`（さらに `Lemma 3.2, (ii)` の Tate 曲線と
Faltings 高さの同種写像評価）を経由する。★★**受けて済ませてはならない**
——存在そのものが `Lemma 3.7` の内容だからである（`check.mjs` B6）。
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField WeierstrassCurve ABC3.Found.GenEll
open scoped Classical

/-! ## ★★★★★★★★★★`deg∞ ≤ 12·ht^Falt + A` -/

/-- ★★★★★★★★★★**`deg∞ ≤ 12·ht^Falt + A`**——★**無条件**、`A` は普遍定数。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

★`Skeleton` は `Proposition 3.4` の 2 本の `≲` を合成して `14` を得ていたが、
★★`htFaltOf` の定義に直接入れば **`12`** で済む（`archSum` の一様上界だけ）。 -/
theorem exists_degInfOf_le_htFalt :
    ∃ A : ℝ, 0 ≤ A ∧ ∀ (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L),
      degInfOf L E ≤ 12 * htFaltOf L E + A := by
  obtain ⟨M, hM1, hMb⟩ := exists_bound_curveArchInv'
  have hpi : (1:ℝ) ≤ (2 * Real.pi) ^ 12 := by
    have h2 : (1:ℝ) ≤ 2 * Real.pi := by nlinarith [Real.two_le_pi]
    exact one_le_pow₀ h2
  have hbig : (1:ℝ) ≤ (2 * Real.pi) ^ 12 * M := by nlinarith
  refine ⟨Real.log ((2 * Real.pi) ^ 12 * M), Real.log_nonneg hbig, fun L _ _ E => ?_⟩
  have hd : (0:ℝ) < (Module.finrank ℚ L : ℝ) := by exact_mod_cast Module.finrank_pos
  have hS := archSum_le L E M hM1 hMb
  have h12 : 12 * htFaltOf L E = degInfOf L E - archSum L E / (Module.finrank ℚ L : ℝ) := by
    rw [htFaltOf]; field_simp
  rw [h12]
  have hdiv : archSum L E / (Module.finrank ℚ L : ℝ) ≤ Real.log ((2 * Real.pi) ^ 12 * M) := by
    rw [div_le_iff₀ hd]
    linarith [hS]
  linarith

/-! ## ★★★★★★★★`ht^Falt` は下に有界 -/

/-- ★★★★★★★★**`ht^Falt` は下に有界**——★**無条件**。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

★`Interface/GenEll/EllModuli.lean` の `faltingsHeight_bddBelow` 欄の中身である。
★★`deg∞ ≥ 0` と上の評価から出る。 -/
theorem exists_htFalt_bddBelow :
    ∃ B : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L),
      B ≤ htFaltOf L E := by
  obtain ⟨A, hA0, hA⟩ := exists_degInfOf_le_htFalt
  refine ⟨-A / 12, fun L _ _ E => ?_⟩
  have h0 : 0 ≤ degInfOf L E := degInfOf_nonneg E
  have h := hA L E
  linarith

/-! ## ★★★★★★★局所高さと素であること -/

variable {L : Type} [Field L] [NumberField L]

/-- ★★★★★**`v_p(Δ_min)·log 2 ≤ d·deg∞`**——原文の観察そのもの。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

★`N(p) ≥ 2` と、和の 1 項が全体以下であることによる。
★★`Found/GaloisRep/DegInfLocal.lean` の `degInfOf_ge_localHeight` の中核を
**完備化を経由しない形**で取り出したものである。 -/
theorem minDeltaExp_log_two_le (E : WeierstrassCurve L) (hΔ : E.Δ ≠ 0)
    (p : HeightOneSpectrum (𝓞 L)) :
    (minDeltaExp p E : ℝ) * Real.log 2 ≤ (Module.finrank ℚ L : ℝ) * degInfOf L E := by
  rw [finrank_mul_degInfOf]
  calc (minDeltaExp p E : ℝ) * Real.log 2
      ≤ (minDeltaExp p E : ℝ) * Real.log (Ideal.absNorm p.asIdeal) := by
        refine mul_le_mul_of_nonneg_left (log_two_le_log_absNorm p) ?_
        exact_mod_cast minDeltaExp_nonneg p E
    _ ≤ ∑ᶠ q : HeightOneSpectrum (𝓞 L),
          (minDeltaExp q E : ℝ) * Real.log (Ideal.absNorm q.asIdeal) := by
        refine single_le_finsum p (hasFiniteSupport_degInf E hΔ) (fun q => ?_)
        exact mul_nonneg (by exact_mod_cast minDeltaExp_nonneg q E) (log_absNorm_nonneg q)

/-- ★★★★★★**`l·log 2 > d·deg∞` なら `l` はすべての局所高さより大きい**。 -/
theorem minDeltaExp_lt_of_lt (E : WeierstrassCurve L) (hΔ : E.Δ ≠ 0) (l : ℕ)
    (h : (Module.finrank ℚ L : ℝ) * degInfOf L E < (l : ℝ) * Real.log 2)
    (p : HeightOneSpectrum (𝓞 L)) : minDeltaExp p E < (l : ℤ) := by
  have hlog : (0:ℝ) < Real.log 2 := Real.log_pos (by norm_num)
  have h1 := minDeltaExp_log_two_le E hΔ p
  have h2 : (minDeltaExp p E : ℝ) * Real.log 2 < (l : ℝ) * Real.log 2 := by linarith
  have h3 : (minDeltaExp p E : ℝ) < (l : ℝ) := lt_of_mul_lt_mul_right (by linarith) hlog.le
  exact_mod_cast h3

/-- ★★★★★★★**`l` は局所高さと素である**。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

★`Interface/GenEll/EllModuli.lean` の `primeToLocalHeights_of_lt` 欄の中身である。
★★局所高さは**正**の `v_p(Δ_min)`（悪い還元の素点）を指すので `hne` を置く
——`v_p = 0` では `Nat.Coprime l 0 = (l = 1)` となり素数では偽になる。 -/
theorem coprime_minDeltaExp (E : WeierstrassCurve L) (hΔ : E.Δ ≠ 0) (l : ℕ) (hl : Nat.Prime l)
    (h : (Module.finrank ℚ L : ℝ) * degInfOf L E < (l : ℝ) * Real.log 2)
    (p : HeightOneSpectrum (𝓞 L)) (hne : minDeltaExp p E ≠ 0) :
    Nat.Coprime l (minDeltaExp p E).toNat := by
  have hlt := minDeltaExp_lt_of_lt E hΔ l h p
  have hnn := minDeltaExp_nonneg p E
  have hpos : 0 < (minDeltaExp p E).toNat := by omega
  rw [Nat.Prime.coprime_iff_not_dvd hl]
  intro hdvd
  have := Nat.le_of_dvd hpos hdvd
  omega

/-! ## ★★★★★★★★悪い還元の素点が 1 つでもあれば `deg∞ > 0` -/

/-- ★★★★★★★★**`v_p(Δ_min) ≠ 0` なる `p` があれば `deg∞ > 0`**——★**無条件**。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

★★`Check/GenEll/EllModuliDegInfPos.lean` の測定（界面は `deg∞ > 0` を強制する）に
応えるために要る。★`Curve` 欄を「乗法還元を持つ半安定曲線」に取る根拠である。 -/
theorem degInfOf_pos_of_minDeltaExp_ne_zero (E : WeierstrassCurve L) (hΔ : E.Δ ≠ 0)
    (p : HeightOneSpectrum (𝓞 L)) (hp : minDeltaExp p E ≠ 0) : 0 < degInfOf L E := by
  have hle : (minDeltaExp p E : ℝ) * Real.log (Ideal.absNorm p.asIdeal)
      ≤ ∑ᶠ q : HeightOneSpectrum (𝓞 L),
          (minDeltaExp q E : ℝ) * Real.log (Ideal.absNorm q.asIdeal) := by
    refine single_le_finsum p (hasFiniteSupport_degInf E hΔ) (fun q => ?_)
    exact mul_nonneg (by exact_mod_cast minDeltaExp_nonneg q E) (log_absNorm_nonneg q)
  have hn : (0:ℝ) < (minDeltaExp p E : ℝ) := by
    have := lt_of_le_of_ne (minDeltaExp_nonneg p E) (Ne.symm hp)
    exact_mod_cast this
  have hl : (0:ℝ) < Real.log (Ideal.absNorm p.asIdeal) :=
    lt_of_lt_of_le (Real.log_pos (by norm_num)) (log_two_le_log_absNorm p)
  have hterm : (0:ℝ) < (minDeltaExp p E : ℝ) * Real.log (Ideal.absNorm p.asIdeal) :=
    mul_pos hn hl
  rw [← finrank_mul_degInfOf] at hle
  have hd : (0:ℝ) < (Module.finrank ℚ L : ℝ) := by exact_mod_cast Module.finrank_pos
  by_contra hc
  push_neg at hc
  nlinarith

/-! ## ★出典の紐付け(`.src`) -/

def degInfOf_pos_of_minDeltaExp_ne_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(悪い還元の素点が 1 つでもあれば deg∞ > 0。★無条件)",
    sectionId := "genell-lemma-3-7" }

def exists_degInfOf_le_htFalt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(deg∞ ≤ 12·ht^Falt + A——★無条件。Skeleton の 14 より強い)",
    sectionId := "genell-lemma-3-7" }

def exists_htFalt_bddBelow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(ht^Falt は下に有界——★無条件。界面の faltingsHeight_bddBelow の中身)",
    sectionId := "genell-lemma-3-7" }

def coprime_minDeltaExp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(l·log 2 > d·deg∞ なら l は局所高さと素——界面の primeToLocalHeights_of_lt の中身)",
    sectionId := "genell-lemma-3-7" }

def coprime_minDeltaExp.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "archSum_le(アルキメデス和の一様上界、§9-355)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.archSum_le") 2,
    .citation "[ABC3]" "finrank_mul_degInfOf(d·deg∞ の和の表示)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.finrank_mul_degInfOf") 2,
    .implicitStep
      ("★★★★★★到達点(2026-08-29、第 561): Skeleton の lemma_3_7 の**(a) の証明が" ++
       "使う 3 つの欄**——faltingsHeight_bddBelow・deg∞ ≤ 14ht^Falt + A・" ++
       "primeToLocalHeights_of_lt——を、どれも**受けた仮定ではなく定理として**取った。" ++
       "★deg∞ の側は Skeleton の 14 ではなく **12** で済む" ++
       "(htFaltOf の定義に直接入れば archSum の一様上界だけでよい)") 8,
    .folklore
      ("☆Lemma 3.7 の (b)・(c) が要求する**例外集合**(lcyclicExc・noMultRedExc)の" ++
       "Galois-finite 性。★Lemma 3.5(さらに Lemma 3.2, (ii) の Tate 曲線と " ++
       "Faltings 高さの同種写像評価)を経由する。" ++
       "★★**受けて済ませてはならない**——存在そのものが Lemma 3.7 の内容である") 9 ]

end ABC3.Found.GaloisRep
