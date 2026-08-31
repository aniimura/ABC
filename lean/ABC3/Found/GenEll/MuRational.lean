/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.MuDescend
import ABC3.Found.GaloisRep.TateSetupDvr
import Mathlib.Data.Nat.Prime.Int

/-!
# 第 946 ブロック —— **★★★★★★★★★★★★★★★★★★★★有理点なら Galois は不要**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

## ★★★★★★★★★★★★★★★★★★★★これは何か——**層 (D) の鍵**

第 945 は「有理点なら `ζ` は下の体に降りる」を示した。
★★**だが、点がはじめから有理なら、Galois を一度も通さなくてよい**。

☆`Lemma 3.2, (i)`（第 887）が Galois 安定性を使うのは、点が拡大体にしか
ない場合に `x` を握るためである。点が `K`-有理なら `x ∈ Kˣ` そのものであり、

    `x^l = q^k`  ⇒  `l·v(x) = k·v(q)`

★`l` は素で `l ∤ v(q)` なので `l ∣ k`——**付値の計算 1 行**である。

★★★★これで `Lemma 3.5` の組み立てから
`IsGalois`・拡大体 `L`・`hσv`・`q₀`・`v` が**すべて消える**。
☆残るのは `q`・`Δ ≠ 0`・`l` が素・`l ∤ v(q)`・`l • P = 0` だけである。

| 定理 | 内容 |
|---|---|
| `dvd_of_pow_eq_zpow_of_coprime` | ★`x^l = Q^k` と `l ∤ v(Q)` から `l ∣ k` |
| `exists_mu_of_rational` | ★★★★★★★★`[x] = [ζ]`（`ζ^l = 1`） |
| `exists_mu_point_of_rational` | ★★★★★★★★★★★★点の側で述べた形 |
| `exists_mu_point_dvr` | ★★★★★★★★★★★★★★★★★★★★**完備 DVR だけで述べた形** |
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine QuotientGroup ABC3.Found.GaloisRep

/-! ## ★付値の計算 1 行 -/

/-- ★**`x^l = Q^k` と `l ∤ v(Q)` から `l ∣ k`**。

☆`l·v(x) = k·v(Q)` なので `l ∣ k·v(Q)`。`l` は素だから `l ∣ k` か `l ∣ v(Q)`。 -/
theorem dvd_of_pow_eq_zpow_of_coprime {K : Type} [Field K]
    (v : Kˣ →* Multiplicative ℤ) {l : ℕ} (hl : l.Prime)
    (x Q : Kˣ) (k : ℤ) (hxk : x ^ l = Q ^ k)
    (hcop : ¬ ((l : ℤ) ∣ vAdd v Q)) : (l : ℤ) ∣ k := by
  have hv : (l : ℤ) * vAdd v x = k * vAdd v Q := by
    have h1 : vAdd v (x ^ l) = (l : ℤ) * vAdd v x := by
      rw [← zpow_natCast x l, vAdd_zpow]
    have h2 : vAdd v (Q ^ k) = k * vAdd v Q := vAdd_zpow v Q k
    rw [← h1, ← h2, hxk]
  have hdvd : (l : ℤ) ∣ k * vAdd v Q := ⟨vAdd v x, hv.symm⟩
  rcases (Nat.prime_iff_prime_int.1 hl).dvd_mul.1 hdvd with h | h
  · exact h
  · exact absurd h hcop

/-! ## ★★★★★★★★類の側で述べた形 -/

section Rational

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [DecidableEq K] [Algebra R K]

/-- ★★★★★★★★**`K`-有理な位数 `l` の類は `μ_l` の類**。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

★★**Galois 安定性を一切使わない**——`x` が `Kˣ` の元であることがその代わりである。 -/
theorem exists_mu_of_rational (S : TateSetup R I K) {l : ℕ} (hl : l.Prime)
    (hcop : ¬ ((l : ℤ) ∣ vAdd S.v S.Q))
    (x : Kˣ) (k : ℤ) (hxk : x ^ l = S.Q ^ k) :
    ∃ ζ : Kˣ, ζ ^ l = 1 ∧
      (QuotientGroup.mk x : Kˣ ⧸ Subgroup.zpowers S.Q) = QuotientGroup.mk ζ :=
  exists_rootOfUnity_mk_eq S.Q hl.pos x k hxk
    (dvd_of_pow_eq_zpow_of_coprime S.v hl x S.Q k hxk hcop)

/-- ★★★★★★★★★★★★**`K`-有理な位数 `l` の点は `μ_l` の点**。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

☆道は 3 段:

1. `Φ` は全射なので `P = tatePhi([x])` と書ける
2. `l • P = 0` から `x^l = Q^k`（第 916）
3. 付値の計算で `l ∣ k`、したがって `[x] = [ζ]`（第 905） -/
theorem exists_mu_point_of_rational (S : TateSetup R I K) {l : ℕ} (hl : l.Prime)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (hcop : ¬ ((l : ℤ) ∣ vAdd S.v S.Q))
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    (P : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point)
    (hP : l • P = 0) :
    ∃ ζ : Kˣ, ζ ^ l = 1 ∧ P = tatePhi S hΔ (QuotientGroup.mk ζ) := by
  -- ★段 1
  obtain ⟨c, hc⟩ := Φ.surjective P
  obtain ⟨x, hx⟩ := QuotientGroup.mk_surjective (s := Subgroup.zpowers S.Q)
    (Additive.toMul c)
  have hPx : tatePhi S hΔ (QuotientGroup.mk x) = P := by
    rw [← hΦ, hx]
    simpa using hc
  -- ★段 2
  obtain ⟨k, hk⟩ := exists_zpow_of_nsmul_tatePhi_eq_zero S hΔ Φ hΦ x l
    (by rw [hPx]; exact hP)
  -- ★段 3
  obtain ⟨ζ, hζl, hζc⟩ := exists_mu_of_rational S hl hcop x k hk.symm
  exact ⟨ζ, hζl, by rw [← hPx, hζc]⟩

end Rational

/-! ## ★★★★★★★★★★★★★★★★★★★★完備 DVR だけで述べた形 -/

section Dvr

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
  {K : Type} [Field K] [DecidableEq K] [Algebra R K] [IsFractionRing R K]

/-- ★★★★★★★★★★★★★★★★★★★★**完備な DVR だけで、
`K`-有理な位数 `l` の点は `μ_l` の点である**。

原文 (GenEll p.15):
> be a one-dimensional Fl-subspace which is stabilized by GK . Then either vK (qE ) ∈

★★★★**2026-09-01（第 946）**——引数は
`q`・`Δ ≠ 0`・`l` が素・`l ∤ v(q)`・`l • P = 0` の 5 つだけである。
☆Tate 一意化 `Φ` は `dvrTatePhiAddEquiv`（第 899）が無条件に与え、
Galois は一切使わない。 -/
theorem exists_mu_point_dvr (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (hq0 : q ≠ 0)
    (hΔ : ((tateCurveAt (mkTateSetup (K := K) q hq hq0).q
      (mkTateSetup (K := K) q hq hq0).hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    {l : ℕ} (hl : l.Prime)
    (hcop : ¬ ((l : ℤ) ∣ vAdd (mkTateSetup (K := K) q hq hq0).v
      (mkTateSetup (K := K) q hq hq0).Q))
    (P : ((tateCurveAt (mkTateSetup (K := K) q hq hq0).q
      (mkTateSetup (K := K) q hq hq0).hq).map (algebraMap R K)).toAffine.Point)
    (hP : l • P = 0) :
    ∃ ζ : Kˣ, ζ ^ l = 1 ∧
      P = tatePhi (mkTateSetup (K := K) q hq hq0) hΔ (QuotientGroup.mk ζ) :=
  exists_mu_point_of_rational (mkTateSetup q hq hq0) hl hΔ hcop
    (dvrTatePhiAddEquiv q hq hq0 hΔ) (fun _ => rfl) P hP

end Dvr

/-! ## ★出典の紐付け(`.src`) -/

def dvd_of_pow_eq_zpow_of_coprime.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (i)(x^l = Q^k と l ∤ v(Q) から l ∣ k)",
    sectionId := "genell-lemma-3-2" }

def exists_mu_of_rational.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (i)(K-有理な位数 l の類は μ_l の類。★Galois 不要)",
    sectionId := "genell-lemma-3-2" }

def exists_mu_point_of_rational.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (i)(K-有理な位数 l の点は μ_l の点。★Galois 不要)",
    sectionId := "genell-lemma-3-2" }

def exists_mu_point_dvr.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (i)(完備 DVR だけで述べた形。★無条件)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GenEll
