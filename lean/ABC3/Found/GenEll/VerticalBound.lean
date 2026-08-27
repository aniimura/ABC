/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.VerticalTwist
import ABC3.Found.GenEll.DegMul
import ABC3.Found.GenEll.LogDiffValue
import ABC3.Found.GenEll.HeightNonneg

/-!
# [GenEll] 垂直な差の**一様上界**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> Now observe that if M is an arithmetic line bundle that arises [by pull-back to X] from an arithmetic line bundle on Spec(Z), then

## ★★★★★★★★★`VerticalTwist.lean` より**弱い仮定**で足りる

`VerticalTwist.lean` は「差が `Spec ℤ` 上の因子の**底変換そのもの**」を要求していた。
★実際にはもっと弱くてよい——**差が有理整数 `n` を含むイデアルであれば足りる**。

★★理由: 垂直因子 `V` は `V ≤ m·X_q` を満たすので、
点 `x : Spec 𝓞_F → X` に沿った引き戻しは `q^m` を**含む**。
★★★したがって「`x^*V` が `q^m` を含む」だけが要る——
**`V` が `Spec ℤ` から来る必要はない**（ファイバーが可約でもよい）。
★★★★**交点数は要らない**。

## ★★★★★★★★鍵は 1 行

> **`n ∈ I`  ⟹  `degNormalized(idealADiv I) ≤ log n`**

★`deg(idealADiv I) = log(absNorm I)`（`LogDiffValue.lean`）と
`absNorm` の単調性（`Ideal.absNorm_dvd_absNorm_of_le`）から出る。

★★`n = q^m` のとき右辺は `m·log q` で、
**点にも定義体にも依らない**——これが原文の `≈` の中身である。

## ★有理整数イデアルの次数

    degNormalized(idealADiv (n·𝓞_F)) = log n     （**`F` に依らない**）

★機構は `absNorm(n·𝓞_F) = n^{[F:ℚ]}`（`Algebra.norm_algebraMap`）で
`[F:ℚ]` が約分されること、それだけである。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

variable (F : Type) [Field F] [NumberField F]

/-! ## ★有理整数が生成するイデアルの正規化次数 -/

/-- ★★★★★**有理整数イデアルの正規化次数は `log n`** —— `F` に依らない。

原文 (GenEll p.6):
> Now observe that if M is an arithmetic line bundle that arises [by pull-back to X] from an arithmetic line bundle on Spec(Z), then

★`absNorm(n·𝓞_F) = |N_{F/ℚ}(n)| = n^{[F:ℚ]}`（`Algebra.norm_algebraMap`）。
★★正規化で `[F:ℚ]` が約分されるので、**定義体に依らない**。 -/
theorem degNormalized_idealADiv_span_natCast (n : ℕ) (hn : n ≠ 0) :
    degNormalized (idealADiv F (Ideal.span {((n : ℕ) : 𝓞 F)})) = Real.log n := by
  have hd : 0 < Module.finrank ℚ F := Module.finrank_pos
  have hne0 : ((n : ℕ) : 𝓞 F) ≠ 0 := Nat.cast_ne_zero.mpr hn
  have hspan : Ideal.span {((n : ℕ) : 𝓞 F)} ≠ 0 := by
    simpa [Ideal.span_singleton_eq_bot] using hne0
  have hnorm : Ideal.absNorm (Ideal.span {((n : ℕ) : 𝓞 F)}) = n ^ Module.finrank ℚ F := by
    rw [Ideal.absNorm_span_singleton]
    have hcast : ((n : ℕ) : 𝓞 F) = algebraMap ℤ (𝓞 F) (n : ℤ) := by push_cast; ring
    rw [hcast, Algebra.norm_algebraMap, RingOfIntegers.rank]
    simp [Int.natAbs_pow]
  rw [degNormalized, deg_idealADiv F _ hspan, hnorm]
  push_cast
  rw [Real.log_pow]
  have : ((Module.finrank ℚ F : ℕ) : ℝ) ≠ 0 := by positivity
  field_simp

/-! ## ★★次数はイデアルについて反単調 -/

/-- ★★**大きいイデアルほど次数は小さい**（`I ≤ J ⟹ deg J ≤ deg I`）。

★イデアルが大きい＝因子が小さい、という向きである。
★★機構は `absNorm` が包含で割り切れること（`Ideal.absNorm_dvd_absNorm_of_le`）。 -/
theorem deg_idealADiv_antitone (I J : Ideal (𝓞 F)) (hI : I ≠ 0) (h : I ≤ J) :
    deg (idealADiv F J) ≤ deg (idealADiv F I) := by
  have hJ : J ≠ 0 := fun hz => hI (le_bot_iff.mp (le_trans h hz.le))
  have hIn : Ideal.absNorm I ≠ 0 := by
    simpa [Ideal.absNorm_eq_zero_iff] using hI
  have hJn : Ideal.absNorm J ≠ 0 := by
    simpa [Ideal.absNorm_eq_zero_iff] using hJ
  have hle : Ideal.absNorm J ≤ Ideal.absNorm I :=
    Nat.le_of_dvd (Nat.pos_of_ne_zero hIn) (Ideal.absNorm_dvd_absNorm_of_le h)
  rw [deg_idealADiv F I hI, deg_idealADiv F J hJ]
  refine Real.log_le_log ?_ ?_
  · exact_mod_cast Nat.pos_of_ne_zero hJn
  · exact_mod_cast hle

/-- ★**正規化した側**の反単調性。 -/
theorem degNormalized_idealADiv_antitone (I J : Ideal (𝓞 F)) (hI : I ≠ 0) (h : I ≤ J) :
    degNormalized (idealADiv F J) ≤ degNormalized (idealADiv F I) := by
  have hd : (0 : ℝ) ≤ (Module.finrank ℚ F : ℝ) := by positivity
  rw [degNormalized, degNormalized]
  exact div_le_div_of_nonneg_right (deg_idealADiv_antitone F I J hI h) hd

/-! ## ★★★★★★★★一様上界 -/

/-- ★★★★★★★★**有理整数 `n` を含むイデアルの正規化次数は `log n` 以下**。

原文 (GenEll p.6):
> Now observe that if M is an arithmetic line bundle that arises [by pull-back to X] from an arithmetic line bundle on Spec(Z), then

★★**上界が `F` にも `I` にも依らず `n` だけで決まる**——これが原文の `≈` の中身である。

★★★垂直因子 `V ≤ m·X_q` の引き戻しは `q^m` を含むので、
本補題が `deg_F(x^*V) ≤ m·log q` を与える。**交点数は要らない**。 -/
theorem degNormalized_idealADiv_le_log_of_natCast_mem
    (I : Ideal (𝓞 F)) (n : ℕ) (hn : n ≠ 0) (hmem : ((n : ℕ) : 𝓞 F) ∈ I) :
    degNormalized (idealADiv F I) ≤ Real.log n := by
  have hne0 : ((n : ℕ) : 𝓞 F) ≠ 0 := Nat.cast_ne_zero.mpr hn
  have hspan0 : Ideal.span {((n : ℕ) : 𝓞 F)} ≠ 0 := by
    simpa [Ideal.span_singleton_eq_bot] using hne0
  have hle : Ideal.span {((n : ℕ) : 𝓞 F)} ≤ I := (Ideal.span_singleton_le_iff_mem _).2 hmem
  calc degNormalized (idealADiv F I)
      ≤ degNormalized (idealADiv F (Ideal.span {((n : ℕ) : 𝓞 F)})) :=
        degNormalized_idealADiv_antitone F _ _ hspan0 hle
    _ = Real.log n := degNormalized_idealADiv_span_natCast F n hn

/-! ## ★★★★★★★★★高さへの適用 -/

/-- ★★★★★★★★★**2 つの `ℤ`-モデルの高さの差は `log n` で抑えられる**。

原文 (GenEll p.6):
> Now observe that if M is an arithmetic line bundle that arises [by pull-back to X] from an arithmetic line bundle on Spec(Z), then

`VerticalTwist.lean` の `htArith_bdeq_of_baseChange'` は
「差が `Spec ℤ` から来る因子の**底変換そのもの**」を要求していた。
★★本定理はそれを **「差が有理整数 `n` を含むイデアルである」** に弱める。

★★★定数 `log n` は**点にも定義体にも依らない**。
★★★★下からは `0`——差が有効因子だからである（`idealADiv_isEffective`）。 -/
theorem htArith_bdeq_of_idealADiv_diff {X X' : Scheme.{0}}
    (D : ArithCartier X) (E : ArithCartier X')
    (ePt : (specRingOfIntegers F ⟶ X) → (specRingOfIntegers F ⟶ X'))
    (n : ℕ) (hn : n ≠ 0)
    (J : (specRingOfIntegers F ⟶ X) → Ideal (𝓞 F))
    (hJn : ∀ xF, ((n : ℕ) : 𝓞 F) ∈ J xF)
    (h : ∀ xF, pullbackADiv F E (ePt xF) - pullbackADiv F D xF = idealADiv F (J xF)) :
    BDeq (fun xF => htArith F D xF) (fun xF => htArith F E (ePt xF)) := by
  refine ⟨Real.log n, fun xF => ?_⟩
  show |htArith F D xF - htArith F E (ePt xF)| ≤ Real.log n
  have hdiff : htArith F E (ePt xF) - htArith F D xF
      = degNormalized (idealADiv F (J xF)) := by
    rw [htArith, htArith, ← degNormalized_sub, h xF]
  have hlo : 0 ≤ degNormalized (idealADiv F (J xF)) :=
    degNormalized_nonneg_of_isEffective F _ (idealADiv_isEffective F _)
  have hhi : degNormalized (idealADiv F (J xF)) ≤ Real.log n :=
    degNormalized_idealADiv_le_log_of_natCast_mem F _ n hn (hJn xF)
  rw [abs_le]
  constructor <;> linarith

/-! ### ★出典の紐付け(`.src`)

★★`Proposition 1.4, (ii)` の証明の中の 1 段の**一般形**である。命題全体ではない。 -/

def degNormalized_idealADiv_span_natCast.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "§1 地の文(正規化次数 deg_F——有理整数イデアルの値は log n で F に依らない)",
    sectionId := "genell-deg" }

def degNormalized_idealADiv_le_log_of_natCast_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (ii)(証明中の段——垂直な差の一様上界)",
    sectionId := "genell-prop-1-4" }

def htArith_bdeq_of_idealADiv_diff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (ii)(証明中の段——ht_{L⊗M} ≈ ht_L の一般形)",
    sectionId := "genell-prop-1-4" }

def htArith_bdeq_of_idealADiv_diff.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "deg_idealADiv(次数はノルムの対数)"
      (.inProject "ABC3" "ABC3.Found.GenEll.deg_idealADiv") 4,
    .citation "[mathlib]" "Ideal.absNorm_dvd_absNorm_of_le(包含でノルムが割り切れる)"
      (.inMathlib "Ideal.absNorm_dvd_absNorm_of_le") 6,
    .citation "[mathlib]" "Algebra.norm_algebraMap(N(n) = n^{[F:ℚ]})"
      (.inMathlib "Algebra.norm_algebraMap") 6,
    .implicitStep
      ("★★★★原文は ht_{L⊗M} ≈ ht_L を 1 行で済ませている。" ++
       "★形式化では (a) 有理整数イデアルの次数が log n であること、" ++
       "(b) 次数がイデアルについて反単調であること、の 2 段に分かれた") 6,
    .implicitStep
      ("★★★VerticalTwist.lean より**弱い仮定**で足りる" ++
       "——差が Spec ℤ から来る必要はなく、有理整数 n を含めばよい。" ++
       "★垂直因子 V ≤ m·X_q の引き戻しは q^m を含むので、" ++
       "ファイバーが可約でも交点数は要らない") 6 ]

end ABC3.Found.GenEll
