/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.BadPrimeData
import Mathlib.FieldTheory.KummerExtension

/-!
# 第 1007 ブロック —— **★★★★★★★★不分岐 2 次拡大の葉**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★これは何か

第 1005／1006 に残った 2 本の仮説のうち、`hsplit`（完備化で**分裂**乗法還元）は
**`p ∣ 2` で非分裂の場合**だけが未処理である。

☆古典的な処理は「不分岐 2 次拡大に上げれば分裂になる」であり、
第 992／993 が捻り `d`（完備化の整数環の**単元**）を与えているので、
上げる先は **`Lv(√d)`** でよい。

★★その第一歩が本ファイルである——`d` が `R` で平方でなければ
**`K` でも平方でない**（したがって `X² − d` は既約で、`K(√d)` は 2 次拡大）。

☆証明は付値環の二分法だけである: `a² = d` なら `a` か `a⁻¹` が `R` の元で、
どちらの場合も `d` の平方根が `R` に取れてしまう。
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField WeierstrassCurve

/-! ## ★★★★★★★★`R` で平方でなければ `K` でも平方でない -/

/-- ★★★★★★★★**離散付値環で平方でない元は、その商体でも平方でない**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1007）**——不分岐 2 次拡大の最初の葉である。
☆`a² = d` とすると、付値環の二分法より `a ∈ R` か `a⁻¹ ∈ R`。

* `a = b ∈ R` なら `b² = d` でそのまま矛盾。
* `a⁻¹ = c ∈ R` なら `c²d = 1` なので `(cd)² = (c²d)d = d` でやはり矛盾。

★これで `X² − d` が `K[X]` で既約（根を持たない 2 次式）であることが言え、
`K(√d)` が 2 次拡大であることの土台になる。 -/
theorem not_isSquare_in_fractionField {R : Type} [CommRing R] [IsDomain R]
    [IsDiscreteValuationRing R] {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]
    (d : R) (hns : ∀ b : R, b * b ≠ d) (a : K) :
    a * a ≠ algebraMap R K d := by
  intro hcon
  have hd0 : d ≠ 0 := fun h => hns 0 (by simp [h])
  have hdne : algebraMap R K d ≠ 0 :=
    (map_ne_zero_iff _ (IsFractionRing.injective R K)).2 hd0
  have hane : a ≠ 0 := by
    intro h; rw [h, zero_mul] at hcon; exact hdne hcon.symm
  rcases ValuationRing.isInteger_or_isInteger R a with ⟨b, hb⟩ | ⟨c, hc⟩
  · exact hns b (IsFractionRing.injective R K (by rw [map_mul, hb, hcon]))
  · have hcdR : c * c * d = 1 := by
      refine IsFractionRing.injective R K ?_
      rw [map_mul, map_mul, hc, ← hcon, map_one]
      field_simp
    refine hns (c * d) ?_
    calc c * d * (c * d) = (c * c * d) * d := by ring
      _ = d := by rw [hcdR, one_mul]

/-! ## ★★★★★★★★`X² − d` は既約 -/

open Polynomial in
/-- ★★★★★★★★**`d` が `R` で平方でなければ `X² − d` は `K[X]` で既約**。

★★★★**2026-09-01（第 1008）**——mathlib の
`X_pow_sub_C_irreducible_of_prime`（`p = 2`）に第 1007 を当てるだけである。
☆これで `AdjoinRoot (X² − d)` が**体**になり、`K(√d)` が建つ。 -/
theorem irreducible_X_sq_sub_C_fractionField {R : Type} [CommRing R] [IsDomain R]
    [IsDiscreteValuationRing R] {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]
    (d : R) (hns : ∀ b : R, b * b ≠ d) :
    Irreducible (X ^ 2 - C (algebraMap R K d)) :=
  X_pow_sub_C_irreducible_of_prime Nat.prime_two
    (fun b hb => not_isSquare_in_fractionField d hns b (by rw [← hb]; ring))

/-! ## ★出典の紐付け(`.src`) -/

def not_isSquare_in_fractionField.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(離散付値環で平方でない元は商体でも平方でない。★無条件)",
    sectionId := "genell-lemma-3-5" }

def irreducible_X_sq_sub_C_fractionField.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(d が平方でなければ X² − d は既約。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
