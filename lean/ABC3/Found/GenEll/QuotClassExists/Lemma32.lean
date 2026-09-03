/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.VeluSetCard
import ABC3.Found.GenEll.EllModuliObjects
import ABC3.Found.GaloisRep.Lemma35Concrete
import ABC3.Found.GaloisRep.UniformFamily
import ABC3.Found.GaloisRep.BadPrimeData
import ABC3.Meta.Claim
import ABC3.Found.GenEll.QuotClassExists.Lemma35

/-!
# QuotClassExists —— `[GenEll] Lemma 3.2` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep IsDedekindDomain NumberField ABC3.Meta
open scoped Classical

/-- ★★★★★★★★★★★★
**`Δ_min` の `l` 倍関係は全素点で成り立つ**——★**無条件**（第 1247）。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

☆悪い素点（`v_p(j) < 0`）では `minDeltaExp_eq_mul_of_jExp_mul`（在庫）、
良い素点（`0 ≤ v_p(j)`）では両辺とも `0` である。

★★★これが第 1244（`quotLCyclicJ` の `deg∞`）が要る `hloc` の形であり、
残るのは `Lemma 3.2, (ii)` が与える `v_p(j(E′)) = l·v_p(j(E))` だけである。 -/
theorem minDeltaExp_eq_mul_of_jExp_all {L : Type} [Field L] [NumberField L]
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    (hss : ∀ p : HeightOneSpectrum (𝓞 L), SemistableAt p E)
    (hss' : ∀ p : HeightOneSpectrum (𝓞 L), SemistableAt p E') {l : ℕ}
    (hbad : ∀ p : HeightOneSpectrum (𝓞 L), jExp p E < 0 → jExp p E' = (l : ℤ) * jExp p E)
    (hgood : ∀ p : HeightOneSpectrum (𝓞 L), 0 ≤ jExp p E → 0 ≤ jExp p E') :
    ∀ p : HeightOneSpectrum (𝓞 L), minDeltaExp p E' = l * minDeltaExp p E := by
  intro p
  rcases lt_or_ge (jExp p E) 0 with hneg | hnn
  · exact minDeltaExp_eq_mul_of_jExp_mul p E E' (hss p) (hss' p) hneg (hbad p hneg)
  · have h1 : minDeltaExp p E = 0 := by
      rw [minDeltaExp_eq_maxJ_of_semistable p E (hss p)]
      exact max_eq_left (by omega)
    have h2 : minDeltaExp p E' = 0 := by
      rw [minDeltaExp_eq_maxJ_of_semistable p E' (hss' p)]
      have hg := hgood p hnn
      exact max_eq_left (by omega)
    rw [h1, h2, mul_zero]

/-- ★★★★★★★★★★★★
**悪い素点の関係を全素点へ束ねる**——★**無条件**（第 1258）。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

☆`hbad` は `minDeltaExp_eq_mul_at_bad_prime`
（`Skeleton/GenEll/TateIsogeny.lean`、証明済み）が与える形である。
★良い素点（`0 ≤ v_p(j)`）では両辺とも `0`。

★★★これが第 1244（`quotLCyclicJ` の `deg∞`）が要る `hloc` を
**局所の結果から直接**作る段である。 -/
theorem minDeltaExp_eq_mul_all {L : Type} [Field L] [NumberField L]
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    (hss : ∀ p : HeightOneSpectrum (𝓞 L), SemistableAt p E)
    (hss' : ∀ p : HeightOneSpectrum (𝓞 L), SemistableAt p E') {l : ℕ}
    (hbad : ∀ p : HeightOneSpectrum (𝓞 L), jExp p E < 0 →
      minDeltaExp p E' = l * minDeltaExp p E)
    (hgood : ∀ p : HeightOneSpectrum (𝓞 L), 0 ≤ jExp p E → 0 ≤ jExp p E') :
    ∀ p : HeightOneSpectrum (𝓞 L), minDeltaExp p E' = l * minDeltaExp p E := by
  intro p
  rcases lt_or_ge (jExp p E) 0 with hneg | hnn
  · exact hbad p hneg
  · have h1 : minDeltaExp p E = 0 := by
      rw [minDeltaExp_eq_maxJ_of_semistable p E (hss p)]
      exact max_eq_left (by omega)
    have h2 : minDeltaExp p E' = 0 := by
      rw [minDeltaExp_eq_maxJ_of_semistable p E' (hss' p)]
      have hg := hgood p hnn
      exact max_eq_left (by omega)
    rw [h1, h2, mul_zero]

def minDeltaExp_eq_mul_all.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(悪い素点の関係を全素点へ束ねる。★無条件)",
    sectionId := "genell-lemma-3-2" }

def minDeltaExp_eq_mul_of_jExp_all.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(Δ_min の l 倍関係は全素点で成り立つ。★無条件)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GenEll
