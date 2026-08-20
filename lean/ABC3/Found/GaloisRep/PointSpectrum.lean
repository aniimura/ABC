import ABC3.Found.GaloisRep.Dedekind
import Mathlib.RingTheory.DedekindDomain.Ideal.Lemmas

/-!
# Galois (G5) 第 138 ブロック —— **★★★★★★高さ 1 の素イデアル ↔ アフィン点**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★因子計算の座標系を敷く

第 137 で `F[W]` が Dedekind 環になり、mathlib の因子機構
(`IsDedekindDomain.HeightOneSpectrum`)が使えるようになった。
★因子を**幾何の言葉**(曲線上の点)で読むために、両者の対応を付ける。

    HeightOneSpectrum F[W]   ↔   { (c, y₀) : W.Equation c y₀ }

| 向き | 出所 |
|---|---|
| 点 → 素イデアル | mathlib の `quotientXYIdealEquiv`(商が `F` だから極大) |
| 素イデアル → 点 | **第 135**(`P = (x − c, y − y₀)`)+ 方程式の確認 |

★★点が曲線上にあることは、生成点の方程式(第 114 `equation_gen`)を
`P` で割って `F ↪ F[W]/P` の単射性で降ろせば出る。

## ★★★mathlib の `XYIdeal` と第 135 の記述が一致する

    XYIdeal W c (C y₀) = span {XClass W c, YClass W (C y₀)}
                       = span {x − c, y − y₀}

★`XClass W c = mk W (C (X − C c)) = x − c`、`YClass W (C y₀) = mk W (Y − C (C y₀)) = y − y₀`。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `XClass_eq` / `YClass_eq` | ★mathlib の記述と我々の生成元の一致 |
| `xyIdeal_eq_span` | ★★`XYIdeal = (x − c, y − y₀)` |
| `equation_of_prime` | ★★★★**素イデアルが定める点は曲線上にある** |
| `xyIdeal_isMaximal` | ★★★曲線上の点の `XYIdeal` は極大 |
| `xyIdeal_ne_bot` | ★`XYIdeal ≠ 0` |
| `prime_eq_xyIdeal` | ★★★★★**0 でない素イデアルは点の `XYIdeal`** |
| `pointSpectrum` | ★★★★点 → 高さ 1 の素イデアル |
| `exists_point_of_heightOneSpectrum` | ★★★★★★**高さ 1 の素イデアルは点に対応する** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain

variable {F : Type} [Field F]

/-! ## ★mathlib の `XClass` / `YClass` と我々の生成元 -/

theorem XClass_eq (W : WeierstrassCurve.Affine F) (c : F) :
    CoordinateRing.XClass W c = genX W - algebraMap F W.CoordinateRing c := by
  rw [CoordinateRing.XClass, genX, map_sub]; rfl

theorem YClass_eq (W : WeierstrassCurve.Affine F) (y₀ : F) :
    CoordinateRing.YClass W (Polynomial.C y₀) = genY W - algebraMap F W.CoordinateRing y₀ := by
  rw [CoordinateRing.YClass, genY, map_sub]; rfl

/-- ★★`XYIdeal` は第 135 の 2 元生成そのものである。 -/
theorem xyIdeal_eq_span (W : WeierstrassCurve.Affine F) (c y₀ : F) :
    CoordinateRing.XYIdeal W c (Polynomial.C y₀)
      = Ideal.span ({genX W - algebraMap F W.CoordinateRing c,
          genY W - algebraMap F W.CoordinateRing y₀} : Set _) := by
  rw [CoordinateRing.XYIdeal, XClass_eq, YClass_eq]

/-! ## ★★★★素イデアルが定める点は曲線上にある -/

/-- ★★★★**素イデアルが定める点は曲線上にある**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★生成点の方程式(第 114)を `P` で割り、`F ↪ F[W]/P` の単射性で `F` に降ろす。 -/
theorem equation_of_prime (W : WeierstrassCurve.Affine F)
    (P : Ideal W.CoordinateRing) [hPp : P.IsPrime] {c y₀ : F}
    (hx : genX W - algebraMap F W.CoordinateRing c ∈ P)
    (hy : genY W - algebraMap F W.CoordinateRing y₀ ∈ P) : W.Equation c y₀ := by
  haveI : Nontrivial (W.CoordinateRing ⧸ P) := (Ideal.Quotient.isDomain P).toNontrivial
  set π := Ideal.Quotient.mk P with hπ
  set ψ : F →+* (W.CoordinateRing ⧸ P) := π.comp (algebraMap F W.CoordinateRing) with hψ
  have hinj : Function.Injective ψ := ψ.injective
  have hgx : π (genX W) = ψ c := by rw [hψ, RingHom.comp_apply, ← Ideal.Quotient.eq] at *; exact hx
  have hgy : π (genY W) = ψ y₀ := by rw [hψ, RingHom.comp_apply, ← Ideal.Quotient.eq] at *; exact hy
  have heq := (equation_iff _ _).1 (equation_gen W)
  simp only [WeierstrassCurve.map_a₁, WeierstrassCurve.map_a₂, WeierstrassCurve.map_a₃,
    WeierstrassCurve.map_a₄, WeierstrassCurve.map_a₆] at heq
  have hq := congrArg π heq
  simp only [map_add, map_mul, map_pow, hgx, hgy, hψ, RingHom.comp_apply] at hq
  rw [equation_iff]
  refine hinj ?_
  simp only [hψ, RingHom.comp_apply, map_add, map_mul, map_pow]
  exact hq

/-! ## ★★★点 → 素イデアル -/

/-- ★★★**曲線上の点の `XYIdeal` は極大イデアルである**——mathlib の商が `F` になることから。 -/
theorem xyIdeal_isMaximal (W : WeierstrassCurve.Affine F) {c y₀ : F} (h : W.Equation c y₀) :
    (CoordinateRing.XYIdeal W c (Polynomial.C y₀)).IsMaximal := by
  refine Ideal.Quotient.maximal_of_isField _ ?_
  exact (CoordinateRing.quotientXYIdealEquiv (W' := W) h).toMulEquiv.isField
    (Field.toIsField F)

/-- ★`XYIdeal` は 0 でない。 -/
theorem xyIdeal_ne_bot (W : WeierstrassCurve.Affine F) (c y₀ : F) :
    CoordinateRing.XYIdeal W c (Polynomial.C y₀) ≠ ⊥ := by
  intro hbot
  have hmem : CoordinateRing.XClass W c ∈ CoordinateRing.XYIdeal W c (Polynomial.C y₀) :=
    Ideal.subset_span (by simp)
  rw [hbot, Ideal.mem_bot] at hmem
  exact CoordinateRing.XClass_ne_zero c hmem

/-- ★★★★曲線上の点が定める高さ 1 の素イデアル。 -/
noncomputable def pointSpectrum (W : WeierstrassCurve.Affine F) {c y₀ : F}
    (h : W.Equation c y₀) : HeightOneSpectrum W.CoordinateRing where
  asIdeal := CoordinateRing.XYIdeal W c (Polynomial.C y₀)
  isPrime := (xyIdeal_isMaximal W h).isPrime
  ne_bot := xyIdeal_ne_bot W c y₀

/-! ## ★★★★★素イデアル → 点 -/

/-- ★★★★★**0 でない素イデアルは曲線上の点の `XYIdeal` である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 135 の 2 元生成に、点が曲線上にあることを添えたもの。 -/
theorem prime_eq_xyIdeal [IsAlgClosed F] (h2 : IsUnit (2 : F))
    (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    (P : Ideal W.CoordinateRing) (hP : P ≠ ⊥) [hPp : P.IsPrime] :
    ∃ c y₀ : F, W.Equation c y₀ ∧ P = CoordinateRing.XYIdeal W c (Polynomial.C y₀) := by
  obtain ⟨c, s, hs, hz, hspan⟩ := exists_prime_gens h2 W P hP
  refine ⟨c, (s - W.a₁ * c - W.a₃) / 2, ?_, ?_⟩
  · exact equation_of_prime W P (hspan ▸ Ideal.subset_span (by simp))
      (hspan ▸ Ideal.subset_span (by simp))
  · rw [xyIdeal_eq_span]; exact hspan

/-- ★★★★★★**高さ 1 の素イデアルは曲線上の点に対応する**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★これで因子を**幾何の言葉**(曲線上の点)で読めるようになる。 -/
theorem exists_point_of_heightOneSpectrum [IsAlgClosed F] (h2 : IsUnit (2 : F))
    (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    (v : HeightOneSpectrum W.CoordinateRing) :
    ∃ (c y₀ : F) (h : W.Equation c y₀), v = pointSpectrum W h := by
  haveI := v.isPrime
  obtain ⟨c, y₀, heq, hv⟩ := prime_eq_xyIdeal h2 W v.asIdeal v.ne_bot
  exact ⟨c, y₀, heq, HeightOneSpectrum.ext hv⟩

/-! ## ★出典の紐付け(`.src`) -/

def exists_point_of_heightOneSpectrum.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——高さ 1 の素イデアルとアフィン点の対応)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
