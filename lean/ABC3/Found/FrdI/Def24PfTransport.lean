/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Def24SuppElt
import ABC3.Found.FrdI.MonoidTransport

/-!
# 台の理論を **`M` から `M^pf` へ**移す —— `Proposition 5.5, (iii)` の `𝒞^pf` の側

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.105。

原文 (FrdI p.105):
> if C is of standard (respectively, rationally standard) type, then so are Cun-tr, Crlf.

## ★★なぜ要るか

`Theorem 6.4, (i)` は `C, C^pf, C^rlf, C^un-tr, (C^pf)^un-tr` の 5 圏が
rationally standard だと言う。★`𝒞^pf` の因子の単系は **`Φ^pf`** なので、
`Definition 4.5, (ii)` の `ι`(因子分解写像)も **`Φ^pf` の上に立て直す**必要がある。

## ★★★橋は 2 本

| 橋 | 条件 | 宣言 |
|---|---|---|
| `Prime M ≃ Prime (M^pf)` | torsion-free | `primeEquivPf`(在庫、`MonoidPrime.lean`) |
| `M^pf ≃ (M^pf)^pf` | perfect ＋ divisorial | ★`pfOfEquiv`(本ファイル) |

★★2 本目は「`M^pf` は perfect(`Pf.isPerfectMonoid_pf`、無条件)かつ
divisorial(`Pf.isDivisorial'`)」から、在庫の
`Pf.of_injective_of_divisorial` / `Pf.of_surjective_of_perfect` で出る。

## ★★★★鍵は `pCarrierPf` の対応

`factorMap` は `boundSup (ι p) (pCarrierPf M p)` なので、移送の中身は

    pCarrierPf (M^pf) (primeToPf p) = Pf.of '' (pCarrierPf M p)

ただ 1 本である(`pCarrierPf_primeToPf`)。
★`⊇` は `isPrimaryElt_pf_of_iff` ＋ `primeToPf_mk`、
`⊆` は在庫の **`exists_nsmul_of_mem_primeCarrier_pf`**(`MonoidPrime.lean`)がそのまま効く。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped NNReal

universe w

variable {M : Type w} [AddCommMonoid M]

/-! ## ★1. divisorial なら torsion-free -/

/-- ★**divisorial ⟹ torsion-free**(`n •` が単射だから)。 -/
theorem isTorsionFreeNaive_of_divisorial (hdiv : IsDivisorial M) : IsTorsionFreeNaive M := by
  intro a n hn h
  have h0 : n • a = n • (0 : M) := by rw [h, smul_zero]
  exact nsmul_inj_of_divisorial hdiv hn h0

/-! ## ★2. `M ≃ M^pf`(perfect ＋ divisorial) -/

/-- ★★**`Pf.of : M → M^pf` は perfect ＋ divisorial なら全単射**。

★これを `M := Φ^pf` に当てて `Φ^pf ≃ (Φ^pf)^pf` を得る
(`Φ^pf` は無条件で perfect、`Φ` が divisorial なら `Φ^pf` も divisorial)。 -/
noncomputable def pfOfEquiv (hdiv : IsDivisorial M) (hperf : IsPerfectMonoid M) : M ≃+ Pf M :=
  AddEquiv.ofBijective (Pf.of : M →+ Pf M)
    ⟨Pf.of_injective_of_divisorial hdiv, Pf.of_surjective_of_perfect hperf⟩

@[simp] theorem pfOfEquiv_apply (hdiv : IsDivisorial M) (hperf : IsPerfectMonoid M) (a : M) :
    pfOfEquiv hdiv hperf a = Pf.of a := rfl

/-! ## ★3. `pCarrierPf` の対応 -/

/-- ★★★★★**`p` 成分は `Pf.of` でちょうど移る**。

★これが `factorMap`(したがって `SuppElt`)の移送の中身である。 -/
theorem pCarrierPf_primeToPf (hdiv : IsDivisorial M) (p : Prime M) :
    pCarrierPf (Pf M) (primeToPf (isTorsionFreeNaive_of_divisorial hdiv) p)
      = (Pf.of : Pf M →+ Pf (Pf M)) '' (pCarrierPf M p) := by
  have htf := isTorsionFreeNaive_of_divisorial hdiv
  have hdivPf : IsDivisorial (Pf M) := Pf.isDivisorial' hdiv
  ext Y
  constructor
  · rintro (⟨n, B, hB, hn⟩ | hY)
    · obtain ⟨X, rfl⟩ := Pf.of_surjective_of_perfect (Pf.isPerfectMonoid_pf (M := M)) Y
      have hnX : ((n : ℕ+) : ℕ) • X = B := by
        refine Pf.of_injective_of_divisorial hdivPf ?_
        rw [map_nsmul]
        exact hn
      obtain ⟨k, b, hb, hk⟩ := exists_nsmul_of_mem_primeCarrier_pf htf p hB
      refine ⟨X, Or.inl ⟨k * n, b, hb, ?_⟩, rfl⟩
      have hstep : ((k : ℕ+) : ℕ) • (((n : ℕ+) : ℕ) • X) = Pf.of b := by
        rw [hnX]; exact hk
      rw [← hstep, smul_smul]
      congr 1
    · have hY0 : Y = 0 := hY
      exact ⟨0, Or.inr rfl, by rw [hY0, map_zero]⟩
  · rintro ⟨X, (⟨n, b, hb, hn⟩ | hX), rfl⟩
    · refine Or.inl ⟨n, Pf.of b, ?_, ?_⟩
      · obtain ⟨hbp, hbe⟩ := hb
        refine ⟨(isPrimaryElt_pf_of_iff htf).mpr hbp, ?_⟩
        rw [← primeToPf_mk htf hbp, hbe]
      · rw [← map_nsmul, hn]
    · have hX0 : X = 0 := hX
      refine Or.inr ?_
      rw [hX0, map_zero]
      rfl

/-! ## ★4. `factorMap` と `SuppElt` の移送 -/

/-- ★**`M^pf ≃ (M^pf)^pf`** —— 上の `pfOfEquiv` を `M^pf` に当てたもの。 -/
noncomputable def pfPfEquiv (hdiv : IsDivisorial M) : Pf M ≃+ Pf (Pf M) :=
  pfOfEquiv (Pf.isDivisorial' hdiv) Pf.isPerfectMonoid_pf

@[simp] theorem pfPfEquiv_apply (hdiv : IsDivisorial M) (x : Pf M) :
    pfPfEquiv hdiv x = Pf.of x := rfl

/-- ★★**`Φ^pf` の上の因子分解写像 `ι^pf`** —— 2 本の橋で `ι` を運んだもの。

★`Definition 4.5, (ii)` の `IsStrictlyRational` は `ι` を**自由な引数として受ける**だけで
`IsPerfFactorialWith` を要求しない。★したがって移送に要るのは
`pCarrierPf` の対応(`pCarrierPf_primeToPf`)ただ 1 本である。 -/
noncomputable def iotaPf (hdiv : IsDivisorial M) (ι : Prime M → Pf M → ℝ≥0) :
    Prime (Pf M) → Pf (Pf M) → ℝ≥0 :=
  fun q Y => ι ((primeEquivPf (isTorsionFreeNaive_of_divisorial hdiv)).symm q)
    ((pfPfEquiv hdiv).symm Y)

section BoundImage

variable {N : Type*} [AddCommMonoid N]

/-- ★**`Bound` は加法的同型で像に移る**(`MLe` が移るから)。 -/
theorem bound_image (e : M ≃+ N) (S : Set M) (X : M) :
    Bound N (e '' S) (e X) = e '' (Bound M S X) := by
  ext A
  constructor
  · rintro ⟨⟨a, ha, rfl⟩, hle⟩
    exact ⟨a, ⟨ha, (mLe_addEquiv e).mp hle⟩, rfl⟩
  · rintro ⟨a, ⟨ha, hle⟩, rfl⟩
    exact ⟨⟨a, ha, rfl⟩, (mLe_addEquiv e).mpr hle⟩

end BoundImage

/-- ★★★★★**因子分解写像は `Pf.of` でちょうど移る**。 -/
theorem factorMap_iotaPf (hdiv : IsDivisorial M) (ι : Prime M → Pf M → ℝ≥0)
    (X : Pf M) (p : Prime M) :
    factorMap (iotaPf hdiv ι) (pfPfEquiv hdiv X)
        (primeToPf (isTorsionFreeNaive_of_divisorial hdiv) p)
      = factorMap ι X p := by
  have htf := isTorsionFreeNaive_of_divisorial hdiv
  have hp : (primeEquivPf htf).symm (primeToPf htf p) = p :=
    (primeEquivPf htf).symm_apply_apply p
  show boundSup (iotaPf hdiv ι (primeToPf htf p))
      (pCarrierPf (Pf M) (primeToPf htf p)) (pfPfEquiv hdiv X)
    = boundSup (ι p) (pCarrierPf M p) X
  rw [pCarrierPf_primeToPf hdiv p]
  show sSup (iotaPf hdiv ι (primeToPf htf p) ''
      Bound (Pf (Pf M)) ((Pf.of : Pf M →+ Pf (Pf M)) '' pCarrierPf M p) (pfPfEquiv hdiv X))
    = sSup (ι p '' Bound (Pf M) (pCarrierPf M p) X)
  have him : (Pf.of : Pf M →+ Pf (Pf M)) '' pCarrierPf M p
      = (pfPfEquiv hdiv) '' pCarrierPf M p := rfl
  rw [him, bound_image (pfPfEquiv hdiv) (pCarrierPf M p) X, Set.image_image]
  congr 1
  refine Set.image_congr ?_
  intro a _
  show ι ((primeEquivPf htf).symm (primeToPf htf p))
    ((pfPfEquiv hdiv).symm ((pfPfEquiv hdiv) a)) = ι p a
  rw [hp, AddEquiv.symm_apply_apply]

/-- ★★★★★★**台も `Pf.of` でちょうど移る** —— これが `𝒞^pf` の
`IsStrictlyRational` に要る形である。 -/
theorem mem_suppElt_iotaPf (hdiv : IsDivisorial M) (ι : Prime M → Pf M → ℝ≥0)
    (a : M) (p : Prime M) :
    primeToPf (isTorsionFreeNaive_of_divisorial hdiv) p
        ∈ SuppElt (iotaPf hdiv ι) (Pf.of a : Pf M)
      ↔ p ∈ SuppElt ι a := by
  have htf := isTorsionFreeNaive_of_divisorial hdiv
  have hkey := factorMap_iotaPf hdiv ι (Pf.of a) p
  show factorMap (iotaPf hdiv ι) (Pf.of (Pf.of a : Pf M)) (primeToPf htf p) ≠ 0
    ↔ factorMap ι (Pf.of a) p ≠ 0
  rw [show (Pf.of (Pf.of a : Pf M) : Pf (Pf M)) = pfPfEquiv hdiv (Pf.of a) from rfl, hkey]

/-! ## ★4. `M^gp → (M^pf)^gp` は単射

★`Definition 4.5, (iii), (b)`(`((𝒞^pf)^un-tr)^birat` の Frobenius-compact 対象)を
`𝒞^pf` で使うとき、`Φ^birat` の元が **`M^pf` の側でも無限位数**であることが要る。
★`Pf.of` が単射(`Pf.of_injective_of_divisorial`)で両側が簡約的なので、
`gpMap_injective` がそのまま効く。 -/

/-- ★★`M` が divisorial なら `M^gp → (M^pf)^gp` は単射。 -/
theorem gpMap_pfOf_injective (hdiv : IsDivisorial M) :
    Function.Injective (gpMap M (Pf.of : M →+ Pf M)) := by
  haveI := isCancelAdd_of_isIntegralMonoid M hdiv.1.1
  haveI := isCancelAdd_of_isIntegralMonoid (Pf M) (Pf.isDivisorial' hdiv).1.1
  exact gpMap_injective M (Pf.of_injective_of_divisorial hdiv)

/-- ★★★したがって「無限位数」は `M^pf` の側でも保たれる。 -/
theorem nsmul_gpMap_pfOf_ne_zero (hdiv : IsDivisorial M) {x : Gp M} (n : ℕ)
    (hx : n • x ≠ 0) : n • gpMap M (Pf.of : M →+ Pf M) x ≠ 0 := by
  rw [← map_nsmul]
  intro h
  exact hx (gpMap_pfOf_injective hdiv (h.trans (map_zero _).symm))

/-! ### ★出典の紐付け -/

def gpMap_pfOf_injective.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — M^gp → (M^pf)^gp は単射",
    sectionId := "frdi-prop-5-5" }

def iotaPf.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — Φ^pf の上の因子分解写像 ι^pf",
    sectionId := "frdi-prop-5-5" }

def factorMap_iotaPf.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — 因子分解写像は Pf.of で移る",
    sectionId := "frdi-prop-5-5" }

def mem_suppElt_iotaPf.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — 台は Pf.of で移る",
    sectionId := "frdi-prop-5-5" }

def pfOfEquiv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — M^pf ≃ (M^pf)^pf(perfect ＋ divisorial)",
    sectionId := "frdi-prop-5-5" }

def pCarrierPf_primeToPf.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — pCarrierPf は Pf.of でちょうど移る",
    sectionId := "frdi-prop-5-5" }

def pCarrierPf_primeToPf.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_nsmul_of_mem_primeCarrier_pf(⊆ の側)"
      (.inProject "ABC3" "ABC3.Found.FrdI.exists_nsmul_of_mem_primeCarrier_pf") 105,
    .citation "[ABC3]" "primeToPf / isPrimaryElt_pf_of_iff(⊇ の側)"
      (.inProject "ABC3" "ABC3.Found.FrdI.primeToPf") 105 ]

end ABC3.Found.FrdI
