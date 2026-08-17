import ABC3.Found.GenEll.LogDiffValue

/-!
# [GenEll] Proposition 1.4, (i) の足場 —— 次数は積を和に送る(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

## ★★★`ht` の加法性の「次数側」

`Proposition 1.4, (i)` は `ht_{L̄⊗M̄} = ht_{L̄} + ht_{M̄}` である。
★構成の側から見ると、この等式は **2 段**に分かれる:

| 段 | 主張 | 状態 |
|---|---|---|
| **引き戻し側** | `x_F^*(D·E) = x_F^*D · x_F^*E` | ★`StalkPullback.lean` で**茎については取れた** |
| **次数側** | `deg(I·J) = deg(I) + deg(J)` | ★★**本ファイル** |

★★2 段が繋がれば `Proposition 1.4, (i)` が構成から出る。
★★★**残る穴は `IdealSheafData.comap` の積の保存**である
(茎では取れているが、大域のイデアルに戻す所が未接続——下記「★残る穴」)。

## ★機構

`deg(idealADiv I) = log(absNorm I)`(`LogDiffValue.lean`)であり、
`Ideal.absNorm` は**モノイド準同型**なので `absNorm(I·J) = absNorm I · absNorm J`。
あとは `Real.log_mul` である。

★より強く**因子の水準で**加法的であること(`idealADiv_mul`)も示す。
これは `Associates.count_mul`(重複度は積で足される)から出る。
★★因子の水準の方が強い——次数を取る前に成り立つからである。

## ★残る穴(B5 を避けるための明示)

本ファイルは `Ideal (𝓞 F)` の積についてのみ述べる。
`X` 上の因子 `D` を `x_F` で引き戻して `Ideal (𝓞 F)` を得る所
(`CartierPullback.lean` の `pullbackIdeal`)は
`Scheme.IdealSheafData.comap` を経由しており、
★**mathlib にも本プロジェクトにも `comap` の積の保存は無い**(2026-08-17 実測)。
★★したがって `Proposition 1.4, (i)` の条なし `.src` は**まだ付けない**。
-/

namespace ABC3.Found.GenEll

open NumberField IsDedekindDomain Ideal

variable (F : Type*) [Field F] [NumberField F]

/-! ## ★★因子の水準での加法性 -/

/-- ★★**イデアルの積は算術因子の和に移る** —— `idealADiv(I·J) = idealADiv I + idealADiv J`。

★これは次数を取る**前**に成り立つ、より強い主張である。
機構は `Associates.count_mul`(素因子の重複度は積で足される)。 -/
theorem idealADiv_fin_apply (I : Ideal (𝓞 F)) (hI : I ≠ 0) (v : FinitePlace F) :
    (idealADiv F I).fin v
      = ((Associates.mk v.asIdeal).count (Associates.mk I).factors : ℤ) := by
  classical
  unfold idealADiv
  rw [dif_neg hI]
  rfl

theorem idealADiv_mul (I J : Ideal (𝓞 F)) (hI : I ≠ 0) (hJ : J ≠ 0) :
    idealADiv F (I * J) = idealADiv F I + idealADiv F J := by
  have hIJ : I * J ≠ 0 := mul_ne_zero hI hJ
  have hmI : (Associates.mk I) ≠ 0 := by
    simpa only [ne_eq, Associates.mk_eq_zero] using hI
  have hmJ : (Associates.mk J) ≠ 0 := by
    simpa only [ne_eq, Associates.mk_eq_zero] using hJ
  refine Prod.ext (Finsupp.ext fun v => ?_) ?_
  · -- 有限素点の成分: 重複度が足される
    have hirr : Irreducible (Associates.mk v.asIdeal) := by
      rw [Associates.irreducible_mk]
      exact (Ideal.prime_of_isPrime v.ne_bot v.isPrime).irreducible
    show (idealADiv F (I * J)).fin v = (idealADiv F I + idealADiv F J).fin v
    rw [ADiv.fin_add, Finsupp.add_apply, idealADiv_fin_apply F _ hIJ,
      idealADiv_fin_apply F I hI, idealADiv_fin_apply F J hJ, ← Nat.cast_add,
      ← Associates.mk_mul_mk, Associates.count_mul hmI hmJ hirr]
  · -- アルキメデス側はどちらも `0`
    show (idealADiv F (I * J)).arc = (idealADiv F I + idealADiv F J).arc
    rw [ADiv.arc_add, idealADiv_arc, idealADiv_arc, idealADiv_arc, add_zero]

/-! ## ★★次数の加法性 -/

/-- ★★★**次数は積を和に送る** —— `deg(idealADiv (I·J)) = deg(...I) + deg(...J)`。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

★`idealADiv_mul` と `deg_add` から出る。
★★これが `Proposition 1.4, (i)` の**次数側**である。 -/
theorem deg_idealADiv_mul (I J : Ideal (𝓞 F)) (hI : I ≠ 0) (hJ : J ≠ 0) :
    deg (idealADiv F (I * J)) = deg (idealADiv F I) + deg (idealADiv F J) := by
  rw [idealADiv_mul F I J hI hJ, deg_add]

/-- ★**正規化した次数も積を和に送る**。

★`ht` は正規化した次数で定義されるので、実際に使うのはこちらである。 -/
theorem degNormalized_idealADiv_mul (I J : Ideal (𝓞 F)) (hI : I ≠ 0) (hJ : J ≠ 0) :
    degNormalized (idealADiv F (I * J))
      = degNormalized (idealADiv F I) + degNormalized (idealADiv F J) := by
  rw [idealADiv_mul F I J hI hJ, degNormalized_add]

/-- ★**ノルムの形での言い換え** —— `log absNorm` が加法的であること。

★`deg_idealADiv` を経由した別証。数値の水準でも確かめておく。 -/
theorem log_absNorm_mul (I J : Ideal (𝓞 F)) (hI : I ≠ 0) (hJ : J ≠ 0) :
    Real.log (Ideal.absNorm (I * J))
      = Real.log (Ideal.absNorm I) + Real.log (Ideal.absNorm J) := by
  rw [← deg_idealADiv F _ (mul_ne_zero hI hJ), ← deg_idealADiv F I hI,
    ← deg_idealADiv F J hJ, deg_idealADiv_mul F I J hI hJ]

/-! ## ★出典の紐付け(`.src`)

★条つきである。`Proposition 1.4, (i)` 全体には**引き戻しの積の保存**が要り、
それは `IdealSheafData.comap` の側に穴が残っている(冒頭の docstring 参照)。 -/

def deg_idealADiv_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (i)(次数側——イデアルの積が次数の和に移ることのみ)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
