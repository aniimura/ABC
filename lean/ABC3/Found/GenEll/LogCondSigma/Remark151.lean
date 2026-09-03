/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.SigmaBound
import ABC3.Found.GenEll.CartierPullback
import ABC3.Found.GenEll.LogCondSigma.Proposition17

/-!
# LogCondSigma —— `[GenEll] Remark 1.5.1` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open NumberField IsDedekindDomain

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

def abs_logCond_sub_le_sum_log.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(log-cond の BD-class が (X_ℚ, D_ℚ) だけに依ることの計算部分)",
    sectionId := "genell-rem-1-5-1" }

end ABC3.Found.GenEll
