/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.SemistableFin
import ABC3.Found.GaloisRep.HtFaltBounds
import ABC3.Found.GaloisRep.VeluNormalized
import ABC3.Found.GaloisRep.DegInfBaseChange
import ABC3.Meta.Claim

/-!
# Lemma35Concrete —— `[GenEll] Lemma 3.2` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField WeierstrassCurve ABC3.Found.GenEll
open scoped Classical
variable {L : Type} [Field L] [NumberField L]

/-! ## ★★★★★★★`deg∞(E′) = l·deg∞(E)` の大域化 -/

/-- ★★★★★★★**局所の関係を大域の `deg∞` へ足し上げる**。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★`Lemma 3.2, (ii)`（`Found/GenEll/Lemma32.lean`、§9-1011）が各素点で与える
`v_p(Δ_min(E′)) = l·v_p(Δ_min(E))` を、`deg∞` の定義（有限台の和）に流すだけである。 -/
theorem degInfOf_eq_of_local (E E' : WeierstrassCurve L) (l : ℕ)
    (hloc : ∀ p : HeightOneSpectrum (𝓞 L), minDeltaExp p E' = l * minDeltaExp p E) :
    degInfOf L E' = (l : ℝ) * degInfOf L E := by
  rw [degInfOf, degInfOf, ← mul_div_assoc]
  congr 1
  rw [mul_finsum]
  refine finsum_congr (fun p => ?_)
  rw [hloc p]
  push_cast
  ring

/-! ## ★出典の紐付け(`.src`)——★★**すべて条つき。項目全体の `.src` は置かない** -/

def degInfOf_eq_of_local.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(deg∞(E′) = l·deg∞(E) の大域化——局所の関係を足し上げる)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GaloisRep
