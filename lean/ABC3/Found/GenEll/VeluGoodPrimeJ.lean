/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluGoodPrimeMem
import ABC3.Found.GenEll.PrimeOverL
import ABC3.Found.GaloisRep.GoodPrimeFromJ
import ABC3.Found.GaloisRep.VeluKernelNormVal
import ABC3.Meta.Claim

/-!
# 第 1438 ブロック —— **良い素点：`j` が整なら `p ∣ l` でも半安定**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——**正しい仮定に直す**

第 1436 の `semistableAt_veluQuot_all_of_goodMem` が残していた仮定
「`p ∣ l` かつ良い素点のとき核の座標が `p` で整」は、
★**一般には成り立たない**——局所体の分岐が `e ≥ l−1` なら
形式群の `l`-捩れ `Ê(𝔪)[l]` が非自明になり、核が深い点を含みうる。

☆本ブロックは同じ場合を **`jExp p E′ ≥ 0`（`j(E′)` の整性）**に置き換える。
★★★こちらは**真**である——`E` が良還元なら同種な `E′` の `j` は整
（古典的にはモジュラー多項式 `Φ_l(j, j′) = 0` の単項性）。

☆道:

1. 判別式の恒等式（第 1402）と `3 ∣ v_p(N)`（第 1396）から **`v_p(Δ(E′)) = 12S`**
2. 一様化元の冪 `n = π^S` を取り（本ブロックの `exists_valAdd_eq`）
3. 第 1431 `minDeltaExp_eq_zero_of_jExp_nonneg` で **`minDeltaExp p E′ = 0`**

| 定理 | 内容 |
|---|---|
| `valAdd_zpow` | ☆`valAdd (x^n) = n·valAdd x`（`n : ℤ`） |
| `exists_valAdd_eq` | ★★**任意の `S : ℤ` に対し `v_p(n) = S` な `n ∈ L` がある** |
| `valAdd_2_eq_zero` | ☆`p ∣ l`・`l ≥ 5` なら `v_p(2) = 0` |
| `semistableAt_veluQuot_good_of_jExp` | ★★★★★★★★**良い素点で `p ∣ l` でも半安定** |
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GaloisRep ABC3.Skeleton.GenEll ABC3.Meta

open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-! ## ☆一様化元 -/

/-- ☆`valAdd` は整数冪でも線型。 -/
theorem valAdd_zpow (p : HeightOneSpectrum (𝓞 L)) (x : Lˣ) (n : ℤ) :
    valAdd p (x ^ n) = n * valAdd p x := by
  rcases n with m | m
  · rw [Int.ofNat_eq_natCast, zpow_natCast, valAdd_pow]
  · rw [zpow_negSucc, valAdd_inv, valAdd_pow, Int.negSucc_eq]
    push_cast
    ring

/-- ★★**任意の `S : ℤ` に対し `v_p(n) = S` となる `n ∈ L` がある**——
mathlib の `valuation_exists_uniformizer` の冪。 -/
theorem exists_valAdd_eq (p : HeightOneSpectrum (𝓞 L)) (S : ℤ) :
    ∃ (n : L) (hn : n ≠ 0), valAdd p (Units.mk0 n hn) = S := by
  obtain ⟨π, hπ⟩ := HeightOneSpectrum.valuation_exists_uniformizer L p
  have hπ0 : π ≠ 0 := by
    intro h
    rw [h, map_zero] at hπ
    exact WithZero.exp_ne_zero hπ.symm
  have hv1 : valAdd p (Units.mk0 π hπ0) = 1 := by
    unfold valAdd
    have hne : (HeightOneSpectrum.valuation L p) ((Units.mk0 π hπ0 : Lˣ) : L) ≠ 0 := by
      rw [Units.val_mk0, hπ]
      exact WithZero.exp_ne_zero
    have h : WithZero.unzero hne = Multiplicative.ofAdd (-1 : ℤ) := by
      rw [← WithZero.coe_inj, WithZero.coe_unzero, Units.val_mk0, hπ]
      rfl
    rw [h]
    rfl
  refine ⟨(π : L) ^ S, zpow_ne_zero S hπ0, ?_⟩
  have hU : Units.mk0 ((π : L) ^ S) (zpow_ne_zero S hπ0) = (Units.mk0 π hπ0) ^ S := by
    refine Units.ext ?_
    push_cast
    rfl
  rw [hU, valAdd_zpow, hv1, mul_one]

/-- ☆`p ∣ l`・`l ≥ 5` なら `v_p(2) = 0`。 -/
theorem valAdd_2_eq_zero (p : HeightOneSpectrum (𝓞 L)) {l : ℕ} (hl : l.Prime)
    (hl5 : 5 ≤ l) (hnu : ¬ IsUnit ((l : primeSubring p))) :
    valAdd p (Units.mk0 (2 : L) two_ne_zero) = 0 := by
  have hne : ((2 : ℕ) : L) ≠ 0 := by norm_num
  have h := valAdd_natCast_eq_zero_of_five_le p hl hl5 hnu 1 0 (by norm_num) hne
  have hEq : Units.mk0 (((2 : ℕ) : L)) hne = Units.mk0 (2 : L) two_ne_zero := by
    refine Units.ext ?_
    push_cast
    rfl
  rwa [hEq] at h

/-! ## ★★★★★★★★良い素点で `j` が整なら半安定 -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★**良い素点では `j(E′)` が整なら Vélu の商は半安定**
——★**`p ∣ l` を許す**（第 1438）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆道は 3 段:

1. 恒等式 `Δ(E)^l = Δ(E′)·N⁴`（第 1402）と `3 ∣ v_p(N)`（第 1396）で `v_p(Δ(E′)) = 12S`
2. `v_p(n) = S` な `n ∈ L` を取る（`exists_valAdd_eq`）
3. 第 1431（`jExp ≥ 0` から `v(c₄) ≥ 4S`・`v(c₆) ≥ 6S`）で `minDeltaExp = 0` -/
theorem semistableAt_veluQuot_good_of_jExp [inst : DecidableEq L]
    (p : HeightOneSpectrum (𝓞 L)) (E : WeierstrassCurve L) [E.IsElliptic]
    [WeierstrassCurve.IsIntegral (primeSubring p) E]
    (hΔ0 : valAdd p (Units.mk0 E.Δ E.isUnit_Δ.ne_zero) = 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (h2 : valAdd p (Units.mk0 (2 : L) two_ne_zero) = 0)
    (h48 : valAdd p (Units.mk0 (48 : L) (by norm_num)) = 0)
    (h864 : valAdd p (Units.mk0 (864 : L) (by norm_num)) = 0)
    (h1728 : valAdd p (Units.mk0 (1728 : L) (by norm_num)) = 0)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    [hE' : (veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).IsElliptic]
    (hjE' : 0 ≤ jExp p (veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))) :
    SemistableAt p (veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) := by
  have hinst : inst = fun a b => Classical.propDecidable (a = b) := by
    funext a b
    exact Subsingleton.elim _ _
  subst hinst
  have hΔE : E.Δ ≠ 0 := E.isUnit_Δ.ne_zero
  have hid := disc_pow_eq_veluQuot_mul E hl hodd Q hQ
  obtain ⟨hN, h3⟩ := three_dvd_valAdd_veluKernelNorm p E hΔE hΔ0 h2 hl hodd Q hQ
  have hΔ' : (veluQuotientFull E (((range l).erase 0).image
      (fun k : ℕ => pointCoords (k • Q)))).Δ ≠ 0 :=
    (veluQuotientFull E (((range l).erase 0).image
      (fun k : ℕ => pointCoords (k • Q)))).isUnit_Δ.ne_zero
  have hu : (Units.mk0 E.Δ hΔE) ^ l
      = (Units.mk0 ((veluQuotientFull E (((range l).erase 0).image
            (fun k : ℕ => pointCoords (k • Q)))).Δ) hΔ')
        * (Units.mk0 (veluKernelNorm E (((range l).erase 0).image
            (fun k : ℕ => pointCoords (k • Q)))) hN) ^ 4 := by
    refine Units.ext ?_
    simpa using hid
  have h1 := congrArg (valAdd p) hu
  rw [valAdd_pow, valAdd_mul, valAdd_pow, hΔ0, mul_zero] at h1
  obtain ⟨m, hm⟩ := h3
  rw [hm] at h1
  obtain ⟨n, hn, hvn⟩ := exists_valAdd_eq p (-m)
  have hΔS : valAdd p (Units.mk0 ((veluQuotientFull E (((range l).erase 0).image
        (fun k : ℕ => pointCoords (k • Q)))).Δ) hΔ')
      = 12 * valAdd p (Units.mk0 n hn) := by
    rw [hvn]
    omega
  exact Or.inl (minDeltaExp_eq_zero_of_jExp_nonneg p _ hn hΔS hjE' h48 h864 h1728)

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★
**良い素点では `j(E′)` が整なら Vélu の商は半安定**——★大域の形（第 1438）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1403 と同じく極小モデルへ正規化してから上の補題を当てる。 -/
theorem semistableAt_veluQuot_goodPrime_of_jExp [inst : DecidableEq L]
    (p : HeightOneSpectrum (𝓞 L)) (E : WeierstrassCurve L) [E.IsElliptic]
    (hss : SemistableAt p E) (hj : 0 ≤ jExp p E)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (h2 : valAdd p (Units.mk0 (2 : L) two_ne_zero) = 0)
    (h48 : valAdd p (Units.mk0 (48 : L) (by norm_num)) = 0)
    (h864 : valAdd p (Units.mk0 (864 : L) (by norm_num)) = 0)
    (h1728 : valAdd p (Units.mk0 (1728 : L) (by norm_num)) = 0)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    [hVell : (veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).IsElliptic]
    (hjE' : 0 ≤ jExp p (veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))) :
    SemistableAt p (veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) := by
  have hinst : inst = fun a b => Classical.propDecidable (a = b) := by
    funext a b
    exact Subsingleton.elim _ _
  subst hinst
  have hΔ : E.Δ ≠ 0 := E.isUnit_Δ.ne_zero
  obtain ⟨C, hC⟩ := WeierstrassCurve.exists_isMinimal (primeSubring p) E
  haveI := hC
  have hmd : minDeltaExp p E = 0 := by
    rw [minDeltaExp_eq_maxJ_of_semistable p E hss]
    exact max_eq_left (by omega)
  haveI hCell : (C • E).IsElliptic :=
    ⟨isUnit_iff_ne_zero.2 (variableChange_Delta_ne_zero E hΔ C)⟩
  have hΔC : valAdd p (Units.mk0 ((C • E).Δ) (C • E).isUnit_Δ.ne_zero) = 0 := by
    have h := minDeltaExp_eq p E hΔ C hC
    rw [hmd] at h
    rw [← h]
  have hQ' : addOrderOf (vcPoint C E Q) = l := by rw [addOrderOf_vcPoint, hQ]
  haveI hell' := isElliptic_veluQuotientFull_nsmul_nf' L (C • E) hQ'
  have heq := veluQuotientFull_vcPoint_eq C E _ hQ two_ne_zero rfl
  haveI hCVell : (C • (veluQuotientFull E (((range l).erase 0).image
      (fun k : ℕ => pointCoords (k • Q))))).IsElliptic :=
    ⟨isUnit_iff_ne_zero.2 (variableChange_Delta_ne_zero _
      (veluQuotientFull E (((range l).erase 0).image
        (fun k : ℕ => pointCoords (k • Q)))).isUnit_Δ.ne_zero C)⟩
  have hjC : 0 ≤ jExp p (veluQuotientFull (C • E)
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • vcPoint C E Q)))) := by
    have h1 := jExp_congr_j p _ _ (j_congr_curve heq)
    rw [← h1, jExp_variableChange]
    exact hjE'
  have hgood := semistableAt_veluQuot_good_of_jExp p (C • E) hΔC hl hodd h2 h48 h864 h1728
    (vcPoint C E Q) hQ' hjC
  rw [← heq] at hgood
  have hres := semistableAt_variableChange p _ C⁻¹ hgood
  rwa [inv_smul_smul] at hres

/-! ## ★出典の紐付け(`.src`) -/

def semistableAt_veluQuot_goodPrime_of_jExp.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(良い素点では j(E′) が整なら Vélu の商は半安定——大域の形)",
    sectionId := "genell-lemma-3-5" }


def valAdd_zpow.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(valAdd は整数冪でも線型。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_valAdd_eq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(任意の S に対し v_p(n) = S な n がある。★無条件)",
    sectionId := "genell-lemma-3-5" }

def valAdd_2_eq_zero.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(p ∣ l・l ≥ 5 なら v_p(2) = 0)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuot_good_of_jExp.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(良い素点では j(E′) が整なら Vélu の商は半安定。★p ∣ l を許す)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuot_good_of_jExp.needs : List ProofObligation :=
  [ .citation "[ABC3]" "disc_pow_eq_veluQuot_mul(第 1402、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.disc_pow_eq_veluQuot_mul") 1,
    .citation "[ABC3]" "three_dvd_valAdd_veluKernelNorm(第 1396、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.three_dvd_valAdd_veluKernelNorm") 1,
    .citation "[ABC3]" "minDeltaExp_eq_zero_of_jExp_nonneg(第 1431、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.minDeltaExp_eq_zero_of_jExp_nonneg") 1,
    .citation "[mathlib]" "HeightOneSpectrum.valuation_exists_uniformizer(一様化元)"
      (.inMathlib "IsDedekindDomain.HeightOneSpectrum.valuation_exists_uniformizer") 1,
    .implicitStep
      ("★★★★★**2026-09-02（第 1438）**——第 1436 が残していた仮定" ++
       "「核の座標が `p` で整」は**一般には偽**である" ++
       "（分岐が `e ≥ l−1` なら形式群の `l`-捩れが非自明）。" ++
       "☆本ブロックは同じ場合を **`jExp p E′ ≥ 0`**（`j` の整性）に置き換えた" ++
       "——こちらは真であり、古典的にはモジュラー多項式の単項性から出る。") 17 ]

end ABC3.Found.GenEll
