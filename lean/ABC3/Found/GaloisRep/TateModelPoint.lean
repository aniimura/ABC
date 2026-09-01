/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateParamMap
import ABC3.Found.GaloisRep.AdicCompleteIntegers
import ABC3.Found.GaloisRep.TateVeluMu
import ABC3.Found.GaloisRep.TateSetupDvr
import ABC3.Found.GenEll.CyclotomicUnits
import ABC3.Found.GenEll.VeluPointSet
import ABC3.Found.GenEll.PointTransport
import Mathlib.NumberTheory.NumberField.Completion.FinitePlace

/-!
# 第 970 ブロック —— **★★★★★★★★★★★★★★★★★★★★Tate モデルの上の点と `j`**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★これは何か

`minDeltaExp_eq_mul_of_torsion`（第 965）は Tate モデルの上の点 `P` と
`j` の一致 `hW′j` を受ける。★本ブロックは**その 2 つを大域のデータから作る**。

☆道は 5 段:

| 段 | 中身 | 出どころ |
|---|---|---|
| 1 | `(E_q) ⊗ Lv = C₀ • (E ⊗ Lv)` | `tateParamR_spec` ＋ `tateModel_baseChange`（第 944） |
| 2 | Vélu の商の楕円性 | `isElliptic_veluQuotient_vcPoint`（第 969） |
| 3 | `j` の一致（`C₀ • (E ⊗ Lv)` 側） | `j_map_velu_vcPoint`（第 967） |
| 4 | 点の位数は運べる | `addOrderOf_rhPoint` ＋ `addOrderOf_vcPoint`（在庫） |
| 5 | 曲線の等式に沿って点を運ぶ | `exists_point_image_eq`（第 968） |

★★配管の注意（実測）: 最後に `rw [himg, hbase]` を **`.j` の上で**やると
`IsElliptic` のせいで motive が壊れる。☆曲線の等式を先に作って
`j_congr_curve`（第 913）で渡すこと——第 944 と同じ穴である。

★★楕円性を `∃`/`∧` の中で運ぶと、後続の `.j` がそれをインスタンスとして
使えない。☆そこで結論を **`∀ _hell : …, j の等式`** の形にした。
呼ぶ側は `haveI` で入れてから使う。
-/

namespace ABC3.Found.GaloisRep

open ABC3.Found.GenEll IsDedekindDomain NumberField WeierstrassCurve Finset
open scoped Classical

/-- ★★★★★★★★★★★★★★★★★★★★**大域の Vélu の商から、Tate モデルの上の
位数 `l` の点と `j` の一致を作る**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 970）**——これで `minDeltaExp_eq_mul_of_torsion`（第 965）の
`P`・`hP`・`hP0`・`hW′j` がすべて揃う。 -/
theorem exists_point_j_tateModel {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange (p.adicCompletion L)).IsElliptic]
    [(E.baseChange (p.adicCompletion L)).IsMinimal (p.adicCompletionIntegers L)]
    [(E'.baseChange (p.adicCompletion L)).IsElliptic]
    (h : (E.baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    {l : ℕ} (hl : l.Prime) {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (h2K : (2 : (p.adicCompletion L)) ≠ 0)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    ∃ P : ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
        (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).toAffine.Point,
      l • P = 0 ∧ P ≠ 0 ∧
      ∀ _hell : (veluQuotientFull
          ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
            (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
            (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
          (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))).IsElliptic,
        (E'.baseChange (p.adicCompletion L)).j
          = (veluQuotientFull
            ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
              (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
            (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))).j := by
  set Lv := p.adicCompletion L
  set R := p.adicCompletionIntegers L
  set φL : L →+* Lv := algebraMap L Lv with hφL
  set φR : R →+* Lv := algebraMap R Lv with hφR
  -- ★段 1: Tate モデルは `E ⊗ Lv` の変数変換
  obtain ⟨hq, C₀, hne, hCE⟩ := tateParamR_spec (E.baseChange Lv) h
  have hbase : (tateCurveAt (tateParamR (E.baseChange Lv) h)
        (tateParamR_mem (E.baseChange Lv) h)).map φR
      = (C₀.map φR) • (E.map φL) :=
    tateModel_baseChange (E.baseChange Lv) h hCE
  -- ★段 2・3: 楕円性と `j` の一致（`C₀ • (E ⊗ Lv)` 側）
  haveI := isElliptic_veluQuotient_vcPoint φL (C₀.map φR) E E' hQ h2K hE'
  have hj := j_map_velu_vcPoint φL (C₀.map φR) E E' hQ h2K hE'
  -- ★段 4: 点の位数は運べる
  have hQ1 : addOrderOf (rhPoint φL E Q) = l := by rw [addOrderOf_rhPoint φL E Q, hQ]
  have hQ2 : addOrderOf (ABC3.Found.GenEll.vcPoint (C₀.map φR) (E.map φL) (rhPoint φL E Q)) = l := by
    rw [ABC3.Found.GenEll.addOrderOf_vcPoint (C₀.map φR) (E.map φL) (rhPoint φL E Q), hQ1]
  -- ★段 5: 曲線の等式に沿って点を運ぶ
  obtain ⟨P, hP, hP0, himg⟩ := exists_point_image_eq hbase hl _ hQ2
  refine ⟨P, hP, hP0, fun hell => ?_⟩
  haveI := hell
  have hcurve : veluQuotientFull
      ((tateCurveAt (tateParamR (E.baseChange Lv) h)
        (tateParamR_mem (E.baseChange Lv) h)).map φR)
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))
      = veluQuotientFull ((C₀.map φR) • (E.map φL))
        (((range l).erase 0).image
          (fun k : ℕ => pointCoords (k • ABC3.Found.GenEll.vcPoint (C₀.map φR) (E.map φL) (rhPoint φL E Q)))) := by
    rw [himg, hbase]
  rw [ABC3.Found.GenEll.j_congr_curve hcurve]
  exact hj

/-! ## ★★★★★★★★★★★★★★★★第 971 —— Tate モデル側の商の楕円性

★第 965（`minDeltaExp_eq_mul_of_torsion`）は Tate モデル側の Vélu の商の楕円性
`hellQ` を受ける。☆その出どころは `hvw` が持っている `veluCurve` の楕円性である——
`veluQuotientFull_tate_mu`（第 890）が両者を繋ぐ。

★★第 969（`E′` から来る楕円性）は `C • (E ⊗ Lv)` **側**の商についてであった。
こちらは **Tate モデル側**である。☆二つは第 968 で移り合う。

★配管の注意: `(mkTateSetup q hq hq0).q` と `q` は定義上等しいが構文上は違うので、
`rw` の前に `show` で `.q` の形に揃える（第 927 と同じ穴）。 -/

open Finset in
open scoped Classical in
/-- ★★★★★★★★★★★★★★★★**`μ_l` による Vélu の商（Tate モデル側）は楕円**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 971）**——これで第 965 の `hellQ` が
`hvw` から直に出る。 -/
theorem isElliptic_veluQuotient_tate_mu {R : Type} [CommRing R] [IsDomain R] [CharZero R]
    [IsDiscreteValuationRing R] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [CharZero K] [Algebra R K] [IsFractionRing R K]
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (hq0 : q ≠ 0)
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    {l : ℕ} (hl : l.Prime) (hlu : IsUnit ((l : R))) (h2K : (2 : K) ≠ 0)
    (ζ : R) (uζ : Kˣ) (hζ : IsPrimitiveRoot ζ l)
    (hζu : algebraMap R K ζ = (uζ : K)) (hζl : uζ ^ l = 1)
    (hord : ∀ n : ℕ, 0 < n → n < l → uζ ^ n ≠ 1)
    (v w : R)
    (hv : v = ∑ i ∈ (range l).erase 0,
          veluV2 (tateCurveAt q hq)
            (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
    (hw : 2 * w = ∑ i ∈ (range l).erase 0,
          (veluU (tateCurveAt q hq)
              (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            + 2 * (veluV2 (tateCurveAt q hq)
                    (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                    (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                  * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)))
    (hell : ((veluCurve (tateCurveAt q hq) v w).map (algebraMap R K)).IsElliptic) :
    (ABC3.Found.GenEll.veluQuotientFull ((tateCurveAt q hq).map (algebraMap R K))
      (((range l).erase 0).image (fun k : ℕ => ABC3.Found.GenEll.pointCoords
        (k • tatePhi (mkTateSetup (K := K) q hq hq0) hΔ (QuotientGroup.mk uζ))))).IsElliptic := by
  haveI := hell
  have hquot := veluQuotientFull_tate_mu (mkTateSetup q hq hq0) hΔ
    (dvrTatePhiAddEquiv q hq hq0 hΔ) (fun _ => rfl) hl.pos ζ uζ hζu hζl hord
    (ABC3.Found.GenEll.isUnit_one_sub_pow_of_isUnit_natCast hl.pos hζ hlu) v w h2K hv hw
  show (ABC3.Found.GenEll.veluQuotientFull ((tateCurveAt (mkTateSetup (K := K) q hq hq0).q
      (mkTateSetup (K := K) q hq hq0).hq).map (algebraMap R K))
    (((range l).erase 0).image (fun k : ℕ => ABC3.Found.GenEll.pointCoords
      (k • tatePhi (mkTateSetup (K := K) q hq hq0) hΔ (QuotientGroup.mk uζ))))).IsElliptic
  rw [hquot]
  exact inferInstanceAs (((veluCurve (tateCurveAt q hq) v w).map (algebraMap R K)).IsElliptic)

def isElliptic_veluQuotient_tate_mu.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(μ_l による Vélu の商(Tate モデル側)は楕円。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★第 975 —— `j` が同じ相方に付け替える

★★第 974 の測定: mathlib の `IsMinimal R W` は「`W` 自身が極小」であって
変数変換で不変ではない。☆`SemistableAt p E′` が与えるのは `C′ • E′` の極小性なので、
第 972 に渡せるのは `C′ • E′` である。
★一方 Vélu の関係（第 974）が与えるのは `C • E′` の側である。

☆両者は `j` が同じ（変数変換で `j` は不変）。**そこで結論の `E′` を
`j` が同じ相方に付け替えられるようにしておく**。 -/

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★**`j` が同じ相方に付け替えた形**の第 970。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 975）**——`E″` が `E′` と同じ `j` を持てば、
`hW′j` は `E″` についても成り立つ。
☆これで「Vélu の関係は `C • E′` で、極小性は `C′ • E′` で」という
食い違いが吸収できる。 -/
theorem exists_point_j_tateModel' {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E E' E'' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange (p.adicCompletion L)).IsElliptic]
    [(E.baseChange (p.adicCompletion L)).IsMinimal (p.adicCompletionIntegers L)]
    [(E'.baseChange (p.adicCompletion L)).IsElliptic]
    [(E''.baseChange (p.adicCompletion L)).IsElliptic]
    (h : (E.baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    {l : ℕ} (hl : l.Prime) {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (h2K : (2 : (p.adicCompletion L)) ≠ 0)
    (hE' : E' = ABC3.Found.GenEll.veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => ABC3.Found.GenEll.pointCoords (k • Q))))
    (hjj : (E''.baseChange (p.adicCompletion L)).j
      = (E'.baseChange (p.adicCompletion L)).j) :
    ∃ P : ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
        (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).toAffine.Point,
      l • P = 0 ∧ P ≠ 0 ∧
      ∀ _hell : (ABC3.Found.GenEll.veluQuotientFull
          ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
            (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
            (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
          (((range l).erase 0).image
            (fun k : ℕ => ABC3.Found.GenEll.pointCoords (k • P)))).IsElliptic,
        (E''.baseChange (p.adicCompletion L)).j
          = (ABC3.Found.GenEll.veluQuotientFull
            ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
              (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
            (((range l).erase 0).image
              (fun k : ℕ => ABC3.Found.GenEll.pointCoords (k • P)))).j := by
  obtain ⟨P, hP, hP0, hj⟩ := exists_point_j_tateModel p E E' h hl hQ h2K hE'
  exact ⟨P, hP, hP0, fun hell => by rw [hjj]; exact hj hell⟩

def exists_point_j_tateModel'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(j が同じ相方に付け替えた形。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★第 977 —— 第 972 の残りの小さな仮説

★第 972 は `h2 : (2 : R) ≠ 0`・`h2K : (2 : Lv) ≠ 0`・`hql : q^l ∈ 𝔪`・
`hΔ : Δ(E_q ⊗ Lv) ≠ 0` を受ける。☆どれも短いのでまとめて取っておく。

* `2 ≠ 0` は標数 0（第 897 の `charZero_adicCompletion(Integers)`）
* `q^l ∈ 𝔪` は `Ideal.pow_mem_of_mem`
* `Δ ≠ 0` は `tateModel_baseChange`（第 944）——Tate モデルは `E ⊗ Lv` の変数変換だから -/

/-- ★**完備化の体は標数 0 なので `2 ≠ 0`**。 -/
theorem two_ne_zero_adicCompletion (L : Type) [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) : (2 : p.adicCompletion L) ≠ 0 := by
  haveI := charZero_adicCompletion L p
  norm_num

/-- ★**完備化の整数環も標数 0 なので `2 ≠ 0`**。 -/
theorem two_ne_zero_adicCompletionIntegers (L : Type) [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) : (2 : p.adicCompletionIntegers L) ≠ 0 := by
  haveI := charZero_adicCompletionIntegers L p
  norm_num

/-- ★**`q ∈ I` なら `q^l ∈ I`**（`l > 0`）。 -/
theorem pow_mem_of_mem_ideal {R : Type} [CommRing R] {I : Ideal R} {q : R} (hq : q ∈ I)
    {l : ℕ} (hl : 0 < l) : q ^ l ∈ I :=
  Ideal.pow_mem_of_mem I hq l hl

/-- ★★★★★★★★**Tate モデルを `K` に上げた曲線の `Δ` は `0` でない**。

☆`tateModel_baseChange`（第 944）により Tate モデルは `W` の変数変換なので、
`W` が楕円なら `Δ ≠ 0` である。★これが第 972 の `hΔ` である。 -/
theorem tateModel_map_Delta_ne_zero {R : Type} [CommRing R] [IsDomain R]
    [IsDiscreteValuationRing R] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [CharZero K] [Algebra R K] [IsFractionRing R K]
    (W : WeierstrassCurve K) [W.IsElliptic] [WeierstrassCurve.IsMinimal R W]
    (h : WeierstrassCurve.HasSplitMultiplicativeReduction R W) :
    ((tateCurveAt (tateParamR W h) (tateParamR_mem W h)).map
      (algebraMap R K)).toAffine.Δ ≠ 0 := by
  obtain ⟨hq, C₀, hne, hCE⟩ := tateParamR_spec W h
  have hbase := tateModel_baseChange W h hCE
  show ((tateCurveAt (tateParamR W h) (tateParamR_mem W h)).map (algebraMap R K)).Δ ≠ 0
  rw [hbase]
  exact ((C₀.map (algebraMap R K)) • W).isUnit_Δ.ne_zero

def two_ne_zero_adicCompletion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(完備化の体は標数 0 なので 2 ≠ 0。★無条件)",
    sectionId := "genell-lemma-3-5" }

def tateModel_map_Delta_ne_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Tate モデルを K に上げた曲線の Δ は 0 でない。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★出典の紐付け(`.src`) -/

def exists_point_j_tateModel.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(大域の Vélu の商から Tate モデルの上の点と j の一致を作る。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
