import ABC3.Found.GaloisRep.LocalHeightDelta
import ABC3.Found.GaloisRep.BadPrimeData
import ABC3.Found.GaloisRep.KernelNormVal
import Mathlib.AlgebraicGeometry.EllipticCurve.ModelsWithJ

/-!
# 潜在的乗法還元 —— `v_p(j) < 0` なら `ofJ j` は `p` で乗法還元

★★★これは Vélu の商の `jExp < 0`(道 C の①)の**下位ノード (i)** である。

## 何を示すか

`Skeleton/GenEll/VeluSemistable.lean` に残る 1 本の `sorry` を対偶・Tate 曲線で攻めるとき、
まず `jExp p X < 0`(＝ `v_p(j(X)) < 0`)から**半安定なモデル**を作らなければならない。
`X` そのものは加法還元でありうる(2 次のねじれが入っている)ので、
在庫の `semistableAt_veluQuot_badPrime_all` はそのままでは当たらない。

★そこで **`j` だけで決まる標準モデル `WeierstrassCurve.ofJ j`** を使う。
`ofJ j` は `j ≠ 0, 1728` のとき

  `ofJNe0Or1728 j = ⟨j − 1728, 0, 0, −36(j−1728)³, −(j−1728)⁵⟩`

であり、`c₄ = j(j−1728)³`・`Δ = j²(j−1728)⁹` である。
`v_p(j) = n < 0` とすると `v_p(j − 1728) = n`(`1728` は整数だから `v_p(1728) ≥ 0`)なので、
**`u = j` の変数変換**を施すと

| 量 | 変換後の値 | 付値 |
|---|---|---|
| `a₁` | `(j−1728)/j` | `0` |
| `a₂`・`a₃` | `0` | —— |
| `a₄` | `−36·(j−1728)³/j⁴` | `≥ −n > 0` |
| `a₆` | `−(j−1728)⁵/j⁶` | `−n > 0` |
| `c₄` | `(j−1728)³/j³` | ★`0`(単元) |
| `Δ` | `(j−1728)⁹/j¹⁰` | `−n > 0` |

となり、**整モデルで `v_p(c₄) = 0`**——`isMinimal_of_c4_valAdd_eq_zero` で極小、
したがって `SemistableAt p (ofJ j)` である。★**体の拡大を一切使わない**のが要点である。

## ★得られるもの

| 定理 | 内容 |
|---|---|
| `valAdd_sub_1728` | ★`v_p(j) < 0` なら `v_p(j − 1728) = v_p(j)` |
| `semistableAt_ofJNe0Or1728_of_valAdd_neg` | ★★★★標準モデルは `p` で乗法還元 |
| `semistableAt_ofJ_of_valAdd_neg` | ★★★★★`SemistableAt p (ofJ j)` |
| `jExp_congr_j` | ★★★`jExp` は `j` だけで決まる |
| `jExp_ofJ` | ★★`jExp p (ofJ j) = v_p(j)` |
| `semistableAt_ofJ_j_of_jExp_neg` | ★★★★★★**潜在的乗法還元の標準モデル** |

## ☆逸脱の記録

無し。原典 (GenEll p.17 Lemma 3.5) は「同種なので自動」と畳んでいるが、
その半安定性の側を埋めるための補助節点であり、前提の追加・読み替えはしていない。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve IsDedekindDomain NumberField
open ABC3.Found.GaloisRep

variable {L : Type} [Field L] [NumberField L]

/-! ## ★`primeSubring` と `valAdd` の小さな橋 -/

/-- ☆整数は `p` の付値環に入る。 -/
theorem mem_primeSubring_intCast (p : HeightOneSpectrum (𝓞 L)) (m : ℤ) :
    ((m : L)) ∈ primeSubring p := by
  rw [mem_primeSubring_iff]
  have h : ((m : L)) = ((m : 𝓞 L) : L) := by norm_cast
  rw [h]
  exact IsDedekindDomain.HeightOneSpectrum.valuation_le_one p _

/-- ☆`valAdd ≥ 0` なら付値環の元。 -/
theorem mem_primeSubring_of_valAdd_nonneg (p : HeightOneSpectrum (𝓞 L)) {x : L} (hx : x ≠ 0)
    (h : 0 ≤ valAdd p (Units.mk0 x hx)) : x ∈ primeSubring p :=
  (mem_primeSubring_iff p x).2 ((valAdd_nonneg_iff p (Units.mk0 x hx)).1 h)

/-- ☆**商の `valAdd`**——`Units.mk0` の中で割り算する形。 -/
theorem valAdd_divL (p : HeightOneSpectrum (𝓞 L)) {a b : L} (ha : a ≠ 0) (hb : b ≠ 0)
    (hab : a / b ≠ 0) :
    valAdd p (Units.mk0 (a / b) hab) = valAdd p (Units.mk0 a ha) - valAdd p (Units.mk0 b hb) := by
  have h1 : a / b * b = a := div_mul_cancel₀ a hb
  have hne : a / b * b ≠ 0 := by rw [h1]; exact ha
  have h2 := valAdd_mulL p hab hb hne
  have h3 : valAdd p (Units.mk0 (a / b * b) hne) = valAdd p (Units.mk0 a ha) :=
    valAdd_eq_of_valuation_eq p _ _ (by simp [h1])
  omega

/-! ## ★`v_p(j) < 0` では `1728` が消える -/

/-- ★`v_p(j) < 0` なら `j ≠ 1728`——`1728` は整数で `v_p(1728) ≥ 0` だから。 -/
theorem sub_1728_ne_zero_of_valAdd_neg (p : HeightOneSpectrum (𝓞 L)) {j : L} (hj0 : j ≠ 0)
    (hneg : valAdd p (Units.mk0 j hj0) < 0) : j - 1728 ≠ 0 := by
  intro h
  have hj : j = ((1728 : ℤ) : L) := by push_cast; exact sub_eq_zero.1 h
  have hmem : j ∈ primeSubring p := by rw [hj]; exact mem_primeSubring_intCast p 1728
  have := (valAdd_nonneg_iff p (Units.mk0 j hj0)).2 ((mem_primeSubring_iff p j).1 hmem)
  omega

/-- ★★**`v_p(j) < 0` なら `v_p(j − 1728) = v_p(j)`**——超距離不等式の等号の場合。 -/
theorem valAdd_sub_1728 (p : HeightOneSpectrum (𝓞 L)) {j : L} (hj0 : j ≠ 0)
    (hneg : valAdd p (Units.mk0 j hj0) < 0) (hd0 : j - 1728 ≠ 0) :
    valAdd p (Units.mk0 (j - 1728) hd0) = valAdd p (Units.mk0 j hj0) := by
  have hb : ValAtLeast p 0 (-1728 : L) := by
    refine valAtLeast_of_mem ?_
    have h : (-1728 : L) = (((-1728 : ℤ)) : L) := by push_cast; ring
    rw [h]
    exact mem_primeSubring_intCast p (-1728)
  have hab : j + (-1728 : L) ≠ 0 := by
    intro h; exact hd0 (by linear_combination h)
  have h := valAdd_add_eq_of_lt (p := p) (n := 0) hj0 hab hb hneg
  rw [← h]
  exact valAdd_eq_of_valuation_eq p _ _ (by simp [sub_eq_add_neg])

/-! ## ★比 `(j−1728)^m / j^k` の付値 -/

/-- ★**`v_p((j−1728)^m / j^k) = (m − k)·v_p(j)`**。 -/
theorem valAdd_ratio (p : HeightOneSpectrum (𝓞 L)) {j : L} (hj0 : j ≠ 0) (hd0 : j - 1728 ≠ 0)
    (hdv : valAdd p (Units.mk0 (j - 1728) hd0) = valAdd p (Units.mk0 j hj0)) (m k : ℕ)
    (h : ((j - 1728) ^ m / j ^ k) ≠ 0) :
    valAdd p (Units.mk0 ((j - 1728) ^ m / j ^ k) h)
      = ((m : ℤ) - (k : ℤ)) * valAdd p (Units.mk0 j hj0) := by
  rw [valAdd_divL p (pow_ne_zero m hd0) (pow_ne_zero k hj0) h,
    valAdd_powL p hd0 m (pow_ne_zero m hd0), valAdd_powL p hj0 k (pow_ne_zero k hj0), hdv]
  ring

/-- ★★**`m ≤ k` なら `(j−1728)^m / j^k` は整**（`v_p(j) < 0` のとき）。 -/
theorem ratio_mem_primeSubring (p : HeightOneSpectrum (𝓞 L)) {j : L} (hj0 : j ≠ 0)
    (hneg : valAdd p (Units.mk0 j hj0) < 0) (hd0 : j - 1728 ≠ 0)
    (hdv : valAdd p (Units.mk0 (j - 1728) hd0) = valAdd p (Units.mk0 j hj0)) {m k : ℕ}
    (hmk : m ≤ k) : ((j - 1728) ^ m / j ^ k) ∈ primeSubring p := by
  refine mem_primeSubring_of_valAdd_nonneg p
    (div_ne_zero (pow_ne_zero m hd0) (pow_ne_zero k hj0)) ?_
  rw [valAdd_ratio p hj0 hd0 hdv m k _]
  have h1 : (0 : ℤ) ≤ (k : ℤ) - (m : ℤ) := by omega
  have h2 : (0 : ℤ) ≤ -valAdd p (Units.mk0 j hj0) := by omega
  have key : ((m : ℤ) - (k : ℤ)) * valAdd p (Units.mk0 j hj0)
      = ((k : ℤ) - (m : ℤ)) * (-valAdd p (Units.mk0 j hj0)) := by ring
  rw [key]
  exact mul_nonneg h1 h2

/-! ## ★★★★標準モデルの半安定性 -/

/-- ★★★**係数がこの形なら極小で `v_p(c₄) = 0`**——変数変換のあとの計算をまとめた形。 -/
theorem semistable_pack (p : HeightOneSpectrum (𝓞 L)) {j : L} (hj0 : j ≠ 0)
    (hneg : valAdd p (Units.mk0 j hj0) < 0) (V : WeierstrassCurve L)
    (h1 : V.a₁ = (j - 1728) ^ 1 / j ^ 1) (h2 : V.a₂ = 0) (h3 : V.a₃ = 0)
    (h4 : V.a₄ = -36 * ((j - 1728) ^ 3 / j ^ 4))
    (h6 : V.a₆ = -((j - 1728) ^ 5 / j ^ 6))
    (hc : V.c₄ = (j - 1728) ^ 3 / j ^ 3)
    (hD : V.Δ = (j - 1728) ^ 9 / j ^ 10) :
    IsMinimal (primeSubring p) V ∧ ∃ h : V.c₄ ≠ 0, valAdd p (Units.mk0 V.c₄ h) = 0 := by
  have hd0 : j - 1728 ≠ 0 := sub_1728_ne_zero_of_valAdd_neg p hj0 hneg
  have hdv : valAdd p (Units.mk0 (j - 1728) hd0) = valAdd p (Units.mk0 j hj0) :=
    valAdd_sub_1728 p hj0 hneg hd0
  have h36 : (-36 : L) ∈ primeSubring p := by
    have e : (-36 : L) = (((-36 : ℤ)) : L) := by push_cast; ring
    rw [e]; exact mem_primeSubring_intCast p (-36)
  haveI hint : WeierstrassCurve.IsIntegral (primeSubring p) V := by
    refine isIntegral_of_mem (primeSubring p) V ?_ ?_ ?_ ?_ ?_
    · rw [h1]; exact ratio_mem_primeSubring p hj0 hneg hd0 hdv (by omega)
    · rw [h2]; exact zero_mem _
    · rw [h3]; exact zero_mem _
    · rw [h4]; exact mul_mem h36 (ratio_mem_primeSubring p hj0 hneg hd0 hdv (by omega))
    · rw [h6]; exact neg_mem (ratio_mem_primeSubring p hj0 hneg hd0 hdv (by omega))
  have hDne : V.Δ ≠ 0 := by
    rw [hD]; exact div_ne_zero (pow_ne_zero 9 hd0) (pow_ne_zero 10 hj0)
  have hcne : V.c₄ ≠ 0 := by
    rw [hc]; exact div_ne_zero (pow_ne_zero 3 hd0) (pow_ne_zero 3 hj0)
  have hcv : valAdd p (Units.mk0 V.c₄ hcne) = 0 := by
    have e : valAdd p (Units.mk0 V.c₄ hcne)
        = valAdd p (Units.mk0 ((j - 1728) ^ 3 / j ^ 3)
            (div_ne_zero (pow_ne_zero 3 hd0) (pow_ne_zero 3 hj0))) :=
      valAdd_eq_of_valuation_eq p _ _ (by simp [hc])
    rw [e, valAdd_ratio p hj0 hd0 hdv 3 3]
    norm_num
  exact ⟨isMinimal_of_c4_valAdd_eq_zero p V hDne hcne hcv, hcne, hcv⟩

/-- ★★★★★★★★★★★★★★★★
**`v_p(j) < 0` なら標準モデル `ofJNe0Or1728 j` は `p` で乗法還元**——★体の拡大は要らない。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`u = j` の変数変換で `v_p(c₄) = 0` になる、というだけの計算である。 -/
theorem semistableAt_ofJNe0Or1728_of_valAdd_neg (p : HeightOneSpectrum (𝓞 L)) {j : L}
    (hj0 : j ≠ 0) (hneg : valAdd p (Units.mk0 j hj0) < 0) :
    SemistableAt p (ofJNe0Or1728 j) := by
  refine Or.inr ⟨(⟨Units.mk0 j hj0, 0, 0, 0⟩ : VariableChange L), ?_⟩
  refine semistable_pack p hj0 hneg _ ?_ ?_ ?_ ?_ ?_ ?_ ?_
  · simp [WeierstrassCurve.variableChange_a₁, ofJNe0Or1728]
    all_goals try field_simp
  · simp [WeierstrassCurve.variableChange_a₂, ofJNe0Or1728]
  · simp [WeierstrassCurve.variableChange_a₃, ofJNe0Or1728]
  · simp [WeierstrassCurve.variableChange_a₄, ofJNe0Or1728]
    all_goals try field_simp
  · simp [WeierstrassCurve.variableChange_a₆, ofJNe0Or1728]
    all_goals try field_simp
  · rw [WeierstrassCurve.variableChange_c₄, ofJNe0Or1728_c₄]
    simp
    all_goals try field_simp
  · rw [WeierstrassCurve.variableChange_Δ, ofJNe0Or1728_Δ]
    simp
    all_goals try field_simp

def semistableAt_ofJNe0Or1728_of_valAdd_neg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(v_p(j) < 0 なら標準モデルは乗法還元。★無条件・拡大不要)",
    sectionId := "genell-lemma-3-5" }

/-- ★★★★★★★★★★★★★★★★★★
**`v_p(j) < 0` なら `ofJ j` は `p` で半安定**。 -/
theorem semistableAt_ofJ_of_valAdd_neg [DecidableEq L] (p : HeightOneSpectrum (𝓞 L)) {j : L}
    (hj0 : j ≠ 0) (hneg : valAdd p (Units.mk0 j hj0) < 0) :
    SemistableAt p (ofJ j) := by
  have hd0 : j - 1728 ≠ 0 := sub_1728_ne_zero_of_valAdd_neg p hj0 hneg
  rw [ofJ_ne_0_ne_1728 j hj0 (sub_ne_zero.mp hd0)]
  exact semistableAt_ofJNe0Or1728_of_valAdd_neg p hj0 hneg

/-- ★★★**`jExp` は `j` だけで決まる**——ねじれで移り合う曲線は同じ `jExp` を持つ。

★これがあるので「`jExp` の計算は標準モデル `ofJ j` の上で行ってよい」と言える。 -/
theorem jExp_congr_j (p : HeightOneSpectrum (𝓞 L)) (E F : WeierstrassCurve L)
    [E.IsElliptic] [F.IsElliptic] (h : E.j = F.j) : jExp p E = jExp p F := by
  by_cases h0 : E.j = 0
  · rw [jExp, dif_pos h0, jExp, dif_pos (h.symm.trans h0)]
  · rw [jExp, dif_neg h0, jExp, dif_neg (fun hx => h0 (h.trans hx))]
    exact valAdd_eq_of_valuation_eq p _ _ (by simp [h])

/-- ★★**`jExp p (ofJ j) = v_p(j)`**——`ofJ_j` から直ちに。 -/
theorem jExp_ofJ [DecidableEq L] (p : HeightOneSpectrum (𝓞 L)) {j : L} (hj0 : j ≠ 0) :
    jExp p (ofJ j) = valAdd p (Units.mk0 j hj0) := by
  have hne : (ofJ j).j ≠ 0 := by rw [ofJ_j]; exact hj0
  rw [jExp, dif_neg hne]
  exact valAdd_eq_of_valuation_eq p _ _ (by simp [ofJ_j])

/-- ★★★★★★★★★★★★★★★★★★★★★★
**潜在的乗法還元の標準モデル**——`jExp p E < 0` なら、`E` と `j` が等しい標準モデル
`ofJ E.j` は `p` で**半安定**(実は乗法還元)であり、`jExp` も一致する。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★これが道 C の①の下位ノード (i) である。
`E` 自身は 2 次のねじれで加法還元でありうるが、**`j` だけで決まるモデルは半安定**であり、
`jExp` は `j` だけで決まるので、`jExp` の計算は `ofJ E.j` の上で行ってよい。 -/
theorem semistableAt_ofJ_j_of_jExp_neg [DecidableEq L] (p : HeightOneSpectrum (𝓞 L))
    (E : WeierstrassCurve L) [E.IsElliptic] (hj : jExp p E < 0) :
    SemistableAt p (ofJ E.j) ∧ jExp p (ofJ E.j) = jExp p E := by
  have hj0 : E.j ≠ 0 := by
    intro h
    rw [jExp, dif_pos h] at hj
    omega
  have hv : valAdd p (Units.mk0 E.j hj0) < 0 := by
    rw [jExp, dif_neg hj0] at hj
    exact hj
  refine ⟨semistableAt_ofJ_of_valAdd_neg p hj0 hv, ?_⟩
  rw [jExp_ofJ p hj0, jExp, dif_neg hj0]

def semistableAt_ofJ_j_of_jExp_neg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(潜在的乗法還元の標準モデルは半安定。★拡大不要)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
