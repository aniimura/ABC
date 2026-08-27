import ABC3.Found.GenEll.SigmaBound
import ABC3.Found.GenEll.CartierPullback

/-!
# [GenEll] `log-cond` の `Σ`-限界(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9–p.10。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

## ★★★★★★★★残っていた橋を渡す

`SigmaBound.lean` は「係数が `≤ 1` の算術因子の `Σ` 上の寄与は `Σ_{q∈Σ} log q` 以下」を
取った。★導手 `f_x^D ≝ (D_x)_red` は `CartierPullback.lean` で
**イデアルの根基**として作られており(`idealADiv F I.radical`)、
その**係数が `≤ 1` であること**が繋ぎに要る。

## ★★★「Dedekind 環の根基イデアルは squarefree」は mathlib に無かった

★`Ideal.radical_isRadical` はあるが、それは**イデアル論の `IsRadical`**
(`r^n ∈ I → r ∈ I`)であって、**モノイドの `IsRadical`**
(`x ∣ y^n → x ∣ y`)ではない。★★`Squarefree` へ渡すには後者が要る。

★★★**両者を繋ぐ道は短かった**——`Ideal.dvd_iff_le` でイデアルの割り算を包含に直せば、
`radical_pow` と `radical_idem` だけでモノイド版が出る:

    K^n ≤ rad(I) ⟹ K ≤ rad(K) = rad(K^n) ≤ rad(rad(I)) = rad(I)

★★★★あとは `IsRadical.squarefree`(mathlib)で `Squarefree`、
`Associates.prime_pow_dvd_iff_le` で重複度 `≤ 1` になる。
-/

namespace ABC3.Found.GenEll

open NumberField IsDedekindDomain

/-! ## ★★★根基はモノイドの意味でも radical -/

/-- ★★★**`rad(I)` はモノイドの意味でも `IsRadical`**。

★イデアル論の `IsRadical`(`r^n ∈ I → r ∈ I`)ではなく、
モノイドの `IsRadical`(`x ∣ y^n → x ∣ y`)である。`Squarefree` へ渡すのはこちら。 -/
theorem isRadical_ideal_radical {F : Type*} [Field F] [NumberField F] (I : Ideal (𝓞 F)) :
    IsRadical (I.radical) := by
  intro n K hdvd
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · rw [pow_zero] at hdvd
    exact dvd_trans hdvd (one_dvd K)
  rw [Ideal.dvd_iff_le] at hdvd ⊢
  calc K ≤ K.radical := Ideal.le_radical
    _ = (K ^ n).radical := (Ideal.radical_pow K (by omega)).symm
    _ ≤ (I.radical).radical := Ideal.radical_mono hdvd
    _ = I.radical := Ideal.radical_idem I

theorem radical_ne_zero {F : Type*} [Field F] [NumberField F] {I : Ideal (𝓞 F)} (hI : I ≠ 0) :
    I.radical ≠ 0 := by
  intro h
  apply hI
  have hle : I ≤ I.radical := Ideal.le_radical
  rw [h] at hle
  exact le_bot_iff.1 hle

/-- ★★★★**根基イデアルの分解の重複度は `≤ 1`**。 -/
theorem count_radical_le_one {F : Type*} [Field F] [NumberField F] (I : Ideal (𝓞 F)) (hI : I ≠ 0)
    (v : FinitePlace F) :
    (Associates.mk v.asIdeal).count (Associates.mk I.radical).factors ≤ 1 := by
  classical
  have hrad0 : I.radical ≠ 0 := radical_ne_zero hI
  have hsq : Squarefree (I.radical) := (isRadical_ideal_radical I).squarefree hrad0
  by_contra hcon
  have h2 : 2 ≤ (Associates.mk v.asIdeal).count (Associates.mk I.radical).factors := by omega
  have hirr : Irreducible (Associates.mk v.asIdeal) := by
    rw [Associates.irreducible_mk]
    exact (Ideal.prime_of_isPrime v.ne_bot v.isPrime).irreducible
  have hne : (Associates.mk I.radical) ≠ 0 := by simpa using hrad0
  have hdvd : (Associates.mk v.asIdeal) ^ 2 ≤ Associates.mk I.radical :=
    (Associates.prime_pow_dvd_iff_le hne hirr).2 h2
  have hd2 : v.asIdeal * v.asIdeal ∣ I.radical := by
    rw [← pow_two, ← Associates.mk_le_mk_iff_dvd, Associates.mk_pow]
    exact hdvd
  exact (Ideal.IsPrime.ne_top v.isPrime) (Ideal.isUnit_iff.1 (hsq _ hd2))

/-- ★★★★★**導手の係数は `≤ 1`** —— これが `SigmaBound` へ渡す形。 -/
theorem idealADiv_radical_coeff_le_one {F : Type*} [Field F] [NumberField F]
    (I : Ideal (𝓞 F)) (hI : I ≠ 0) (v : FinitePlace F) :
    (idealADiv F I.radical).fin v ≤ 1 := by
  classical
  have hrad0 : I.radical ≠ 0 := radical_ne_zero hI
  have hval : (idealADiv F I.radical).fin v
      = ((Associates.mk v.asIdeal).count (Associates.mk I.radical).factors : ℤ) := by
    unfold idealADiv
    rw [dif_neg hrad0]
    rfl
  rw [hval]
  exact_mod_cast count_radical_le_one I hI v

/-! ## ★★★★★★★★★`log-cond` の `Σ`-限界 -/

variable {X : AlgebraicGeometry.Scheme.{0}}

/-- ★★★★★★★★★**`Σ` の上に台を持つ `log-cond` は `Σ_{q∈Σ} log q` 以下**。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

★★**右辺は `Σ` だけで決まる定数**であり、点 `x` にも定義体 `F` にも依らない。
★★★これが `Remark 1.5.1`(BD-class が `(X_ℚ, D_ℚ)` だけに依る)と
`Proposition 1.7`(`Σ` 上の寄与は `≈ 0`)が共通して使う形である。 -/
theorem logCond_le_sum_log (F : Type) [Field F] [NumberField F]
    (D : X.IdealSheafData) (xF : specRingOfIntegers F ⟶ X)
    (hI : pullbackIdeal F D xF ≠ 0)
    (Sig : Finset ℕ) (hprime : ∀ q ∈ Sig, q.Prime)
    (ch : FinitePlace F → ℕ)
    (hmem : ∀ v ∈ (conductorADiv F D xF).fin.support, ch v ∈ Sig)
    (hover : ∀ v ∈ (conductorADiv F D xF).fin.support,
      (v.asIdeal).LiesOver (Ideal.span {((ch v : ℕ) : ℤ)})) :
    logCond F D xF ≤ ∑ q ∈ Sig, Real.log q := by
  rw [logCond]
  refine degNormalized_le_sum_log (conductorADiv F D xF) ?_ ?_ Sig hprime ch hmem hover
  · rw [conductorADiv]
    exact idealADiv_arc F _
  · intro v
    rw [conductorADiv]
    exact idealADiv_radical_coeff_le_one _ hI v

/-! ## ★★★★★★★★★★★2 つのモデルの `log-cond` の差 -/

/-- ★★★★★★★★★★**`Σ` の外で導手が一致すれば `log-cond` の差は `Σ_{q∈Σ} log q` 以下**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX(Q) depends only on the pair

★★これが `Remark 1.5.1` の結論の**計算部分**である。
原文の証明は「同型が `ℤ[Σ^{-1}]` へ延びる」で `Σ` の外の一致を与え、
そこから BD-class の一致を結論する。★★★本定理は**後半**を取る
——前半(spreading out)は仮定 `hagree` として受けている。

★★★★**右辺は `Σ` だけで決まる定数**であり、点 `x` にも定義体 `F` にも依らない。
これが BD-class が吸収する「有限個の素数を除けば」の中身である。 -/
theorem logCond_sub_le_sum_log (F : Type) [Field F] [NumberField F]
    {X X' : AlgebraicGeometry.Scheme.{0}}
    (D : X.IdealSheafData) (D' : X'.IdealSheafData)
    (xF : specRingOfIntegers F ⟶ X) (xF' : specRingOfIntegers F ⟶ X')
    (hI : pullbackIdeal F D xF ≠ 0)
    (Sig : Finset ℕ) (hprime : ∀ q ∈ Sig, q.Prime)
    (ch : FinitePlace F → ℕ)
    (hagree : ∀ v, ch v ∉ Sig →
      (conductorADiv F D xF).fin v = (conductorADiv F D' xF').fin v)
    (hover : ∀ v : FinitePlace F,
      (v.asIdeal).LiesOver (Ideal.span {((ch v : ℕ) : ℤ)})) :
    logCond F D xF - logCond F D' xF' ≤ ∑ q ∈ Sig, Real.log q := by
  rw [logCond, logCond]
  refine degNormalized_sub_le_sum_log (conductorADiv F D xF) (conductorADiv F D' xF')
    (idealADiv_arc F _) (idealADiv_arc F _) ?_ ?_ Sig hprime ch hagree hover
  · intro v
    rw [conductorADiv]
    exact idealADiv_radical_coeff_le_one _ hI v
  · intro v
    exact (idealADiv_isEffective F _).1 v

/-- ★★★★★★★★★★★**`log-cond` は `Σ` の外で一致すれば BD-同値**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX(Q) depends only on the pair

★両向きの差を取って絶対値で抑えた形。★★`BDeq` の定数が
**`Σ_{q∈Σ} log q`(点にも定義体にも依らない)**であることが要点である。 -/
theorem abs_logCond_sub_le_sum_log (F : Type) [Field F] [NumberField F]
    {X X' : AlgebraicGeometry.Scheme.{0}}
    (D : X.IdealSheafData) (D' : X'.IdealSheafData)
    (xF : specRingOfIntegers F ⟶ X) (xF' : specRingOfIntegers F ⟶ X')
    (hI : pullbackIdeal F D xF ≠ 0) (hI' : pullbackIdeal F D' xF' ≠ 0)
    (Sig : Finset ℕ) (hprime : ∀ q ∈ Sig, q.Prime)
    (ch : FinitePlace F → ℕ)
    (hagree : ∀ v, ch v ∉ Sig →
      (conductorADiv F D xF).fin v = (conductorADiv F D' xF').fin v)
    (hover : ∀ v : FinitePlace F,
      (v.asIdeal).LiesOver (Ideal.span {((ch v : ℕ) : ℤ)})) :
    |logCond F D xF - logCond F D' xF'| ≤ ∑ q ∈ Sig, Real.log q := by
  have h1 := logCond_sub_le_sum_log F D D' xF xF' hI Sig hprime ch hagree hover
  have h2 := logCond_sub_le_sum_log F D' D xF' xF hI' Sig hprime ch
    (fun v hv => (hagree v hv).symm) hover
  rw [abs_sub_le_iff]
  exact ⟨h1, h2⟩

/-! ## ★出典の紐付け(`.src`) -/

def count_radical_le_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (ii)((−)_red の係数が ≤ 1 であること)",
    sectionId := "genell-def-1-5" }

def logCond_le_sum_log.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i)(Σ 上の log-cond の寄与は Σ_{q∈Σ} log q で抑えられる)",
    sectionId := "genell-prop-1-7" }

def abs_logCond_sub_le_sum_log.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(log-cond の BD-class が (X_ℚ, D_ℚ) だけに依ることの計算部分)",
    sectionId := "genell-rem-1-5-1" }

end ABC3.Found.GenEll
