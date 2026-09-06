/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluDualJ
import ABC3.Found.GenEll.VeluJDescent
import ABC3.Found.GenEll.VeluImage
import ABC3.Found.GenEll.VeluEllipticNF
import ABC3.Found.GenEll.JScale
import ABC3.Found.GenEll.VeluSemistableAll
import ABC3.Skeleton.GenEll.TateLocalModelK
import ABC3.Meta.Claim

/-!
# 第 1445 ブロック —— **道 C(対偶)の組み上げ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か

`Skeleton/GenEll/VeluSemistable.lean` に残る `sorry`
（**`p ∣ l` かつ良い素点のとき `0 ≤ jExp p (E/⟨Q⟩)`**）を
対偶で攻める道（道 C）の**最後の組み上げ**である:

    jExp p E′ < 0 と仮定
      → N1 で有限次拡大の上で E′ は半安定・jExp < 0   …… 第 1444
      → ② で E′ の l-巡回商のひとつが E と同じ j を持つ …… 第 1442
      → ① その商の jExp < 0                            …… ☆残る
      → jExp_congr_j で jExp p E < 0、仮定に矛盾

## ★★★★★★★★☆**拡大は合成体ではなく塔で取る**（測定、2026-09-06）

☆②が返す拡大 `N` と N1 が返す拡大 `M` は別物なので、
素朴には合成体 `N ⊔ M` が要る。★本ファイルは**順序を変えて塔にする**ことで
合成体を避けた:

1. まず②を `L` の上で当てて `N` を取る（`Q′` は `E′ ⊗ N` の上の位数 `l` の点）
2. `N` の中で `p` の上にある素点 `P` を取り、`jExp P (E′ ⊗ N) < 0` を出す
3. `N` を底体として N1 を当てて `M`（`N` の上の有限次拡大）を取る
4. `L ⊆ N ⊆ M` の塔の上で①を当てる

★`veluQuotientFull` は底変換と可換（`veluQuotientFull_baseChange`、第 1197）なので、
②が `N` の上で与えた `j` の等式はそのまま `M` へ上がる。
★`jExp_baseChange`（`v_P(j) = e·v_p(j)`、`e > 0`）は**符号を保存かつ反映**するので、
塔を上げ下げしても不等式の向きは変わらない。

## ★★★★★★★★★★★★①はさらに半分になった

在庫の `minDeltaExp_eq_mul_at_bad_prime_any_K`（第 1140）は
**`l ∤ v_p(j)` なら `minDeltaExp p (E/⟨Q⟩) = l · minDeltaExp p E`** を与える。
☆`semistableAt_veluQuot_badPrime_all`（第 1436）で商の半安定性が出るので、
半安定な曲線では `minDeltaExp = max(0, −jExp)` から `jExp < 0` が出る。

★したがって①のうち **`l ∤ v_P(j)` の場合は閉じた**
（`jExp_neg_veluQuot_of_not_dvd`、本ファイル、sorry 0）。
☆残るのは **`l ∣ v_P(j)`** の場合だけで、これは Tate 母数 `q` の付値が
`l` で割れる場合、すなわち**核が `μ_l` でなく「深い点」で生成されうる**場合である
（`exists_primitiveRoot_of_torsion_point` が `hcop` を要求するのはこのため）。

★★★これが `JExpNegVeluStable`（本ファイルの `def`）である
——`Skeleton/GenEll/VeluSemistable.lean` の `.needs` が
「①は下位 2 つに落ちる: (i) 潜在的乗法還元は有限次拡大で乗法還元になる、
(ii) 深い核の Vélu の商の `jExp < 0`」と書いた**(ii) そのもの**である。
(i) は第 1444 で閉じた。

## ★逸脱の記録

☆**逸脱なし**。原典が「同種なので自動」と括弧で畳んだ段を、
形式化の都合で対偶に組み替えているだけで、前提の追加も読み替えもしていない。
★`JExpNegVeluStable` に `v_p(48) = 0`・`v_p(864) = 0`（すなわち `p ∤ 6`）を
入れてあるのは**仮説を弱める**ためで、呼び出し側（`p ∣ l` かつ `l ≥ 5`）では
`valAdd_48_eq_zero`・`valAdd_864_eq_zero`（第 1435）から無条件に出る。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

/-! ## ★配管——上にある素点と、付値の持ち上げ -/

/-- ☆**有限次中間体には `p` の上にある素点がある**（整拡大での持ち上げ）。

☆`SSCurveExt.exists_heightOneSpectrum_over`（第 1347）の中間体版である。 -/
theorem exists_heightOneSpectrum_over_intermediate
    (K : Type) [Field K] [NumberField K]
    (M : IntermediateField K (AlgebraicClosure K)) [FiniteDimensional K M]
    (p : HeightOneSpectrum (𝓞 K)) :
    letI : NumberField (M : Type) := NumberField.of_module_finite K M
    ∃ P : HeightOneSpectrum (𝓞 (M : Type)), P.asIdeal.LiesOver p.asIdeal := by
  letI : NumberField (M : Type) := NumberField.of_module_finite K M
  obtain ⟨P, -, hPp, hPo⟩ :=
    Ideal.exists_ideal_over_prime_of_isIntegral (R := 𝓞 K) (S := 𝓞 (M : Type))
      p.asIdeal ⊥ (by simp)
  haveI := hPp
  have hPne : P ≠ ⊥ := by
    intro h
    rw [h] at hPo
    exact p.ne_bot (by simpa using hPo.symm)
  exact ⟨⟨P, hPp, hPne⟩, ⟨hPo.symm⟩⟩

/-- ☆**付値 `0` は拡大へ上がる**——`valAdd_algebraMap`（`v_P = e·v_p`）から。 -/
theorem valAdd_eq_zero_of_algebraMap (K K' : Type) [Field K] [NumberField K]
    [Field K'] [NumberField K'] [Algebra K K'] [Algebra (𝓞 K) K']
    [IsScalarTower (𝓞 K) K K'] [IsScalarTower (𝓞 K) (𝓞 K') K']
    [Module.Finite (𝓞 K) (𝓞 K')]
    (p : HeightOneSpectrum (𝓞 K)) (P : HeightOneSpectrum (𝓞 K'))
    [P.asIdeal.LiesOver p.asIdeal]
    {x : K} (hx : x ≠ 0) {y : K'} (hy : y ≠ 0) (hxy : algebraMap K K' x = y)
    (h : valAdd p (Units.mk0 x hx) = 0) :
    valAdd P (Units.mk0 y hy) = 0 := by
  have hu : Units.mk0 y hy = Units.map (algebraMap K K' : K →* K') (Units.mk0 x hx) :=
    Units.ext hxy.symm
  rw [hu, valAdd_algebraMap K K' p P, h, mul_zero]

/-! ## ★★★★★★★★①のうち `l ∤ v_p(j)` の場合 -/

/-- ★★★★★★★★★★★★★★★★
**`l ∤ v_p(j)` なら Vélu の商も乗法還元**（第 1445）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆道は 3 本の在庫を並べるだけである:

| 段 | 在庫 | 第 |
|---|---|---|
| 商は半安定（`p ∣ l` も可） | `semistableAt_veluQuot_badPrime_all` | 1436 |
| `minDeltaExp(E′) = l·minDeltaExp(E)` | `minDeltaExp_eq_mul_at_bad_prime_any_K` | 1140 |
| 半安定なら `minDeltaExp = max(0, −jExp)` | `minDeltaExp_eq_maxJ_of_semistable` | 在庫 |

★`minDeltaExp p E > 0`（`minDeltaExp_pos_of_jExp_neg`）なので
`minDeltaExp p E′ = l·minDeltaExp p E > 0`、したがって `jExp p E′ < 0` である。 -/
theorem jExp_neg_veluQuot_of_not_dvd {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (Y : WeierstrassCurve L) [Y.IsElliptic]
    (hss : SemistableAt p Y) (hj : jExp p Y < 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (h48 : valAdd p (Units.mk0 (48 : L) (by norm_num)) = 0)
    (h864 : valAdd p (Units.mk0 (864 : L) (by norm_num)) = 0)
    (hcop : ¬ ((l : ℤ) ∣ jExp p Y))
    (Q : Y.toAffine.Point) (hQ : addOrderOf Q = l)
    [hVell : (veluQuotientFull Y (((range l).erase 0).image
      (fun k : ℕ => pointCoords (k • Q)))).IsElliptic] :
    jExp p (veluQuotientFull Y
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) < 0 := by
  have hssX := semistableAt_veluQuot_badPrime_all p Y hss hj hl hodd h48 h864 Q hQ
  have hmd := ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_at_bad_prime_any_K p Y
    (veluQuotientFull Y (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))
    hss hssX hj hl hcop hQ rfl
  have hpos : 0 < minDeltaExp p Y := minDeltaExp_pos_of_jExp_neg p Y hss hj
  have hmaxX := minDeltaExp_eq_maxJ_of_semistable p
    (veluQuotientFull Y (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) hssX
  have hlpos : (0 : ℤ) < (l : ℤ) := by exact_mod_cast hl.pos
  rw [hmd] at hmaxX
  have hgt : (0 : ℤ) < max 0 (-jExp p (veluQuotientFull Y
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))) := by
    rw [← hmaxX]
    exact mul_pos hlpos hpos
  have h2 := (lt_max_iff.mp hgt).resolve_left (lt_irrefl (0 : ℤ))
  linarith

/-! ## ★★★★★★★★★★★★①の残り——深い核の場合 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**①の残り**（＝`Skeleton/GenEll/VeluSemistable.lean` の `.needs` の (ii)）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆**乗法還元は `l`-同種で保たれる**という主張の、
`l ∣ v_P(j)`（＝ Tate 母数 `q` の付値が `l` で割れる）場合である。
★この場合だけ核が `μ_l` でなく**深い点**で生成されうるので、
`exists_primitiveRoot_of_torsion_point`（`hcop` を要求する）が使えない。

☆`l ∤ v_P(j)` の側は `jExp_neg_veluQuot_of_not_dvd`（本ファイル）で閉じている。
★`p ∤ 6`（`v_P(48) = v_P(864) = 0`）を仮説に入れてあるのは**仮説を弱める**ためで、
呼び出し側（`p ∣ l` かつ `l ≥ 5`）では第 1435 から無条件に出る。 -/
def JExpNegVeluStable : Prop :=
  ∀ (K : Type) [Field K] [NumberField K] (P : HeightOneSpectrum (𝓞 K))
    (Y : WeierstrassCurve K) [Y.IsElliptic],
    SemistableAt P Y → jExp P Y < 0 →
    ∀ {l : ℕ}, l.Prime → 5 ≤ l →
    valAdd P (Units.mk0 (48 : K) (by norm_num)) = 0 →
    valAdd P (Units.mk0 (864 : K) (by norm_num)) = 0 →
    ((l : ℤ) ∣ jExp P Y) →
    ∀ (Q : Y.toAffine.Point) (_hQ : addOrderOf Q = l)
      [(veluQuotientFull Y (((range l).erase 0).image
        (fun k : ℕ => pointCoords (k • Q)))).IsElliptic],
      jExp P (veluQuotientFull Y
        (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) < 0

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★★★★★★★★★
**①の全体**——`l ∤ v_p(j)` は本ファイルで閉じ、`l ∣ v_p(j)` は `JExpNegVeluStable`。 -/
theorem jExp_neg_veluQuot_of_stable (h1 : JExpNegVeluStable)
    {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (Y : WeierstrassCurve L) [Y.IsElliptic]
    (hss : SemistableAt p Y) (hj : jExp p Y < 0)
    {l : ℕ} (hl : l.Prime) (hl5 : 5 ≤ l)
    (h48 : valAdd p (Units.mk0 (48 : L) (by norm_num)) = 0)
    (h864 : valAdd p (Units.mk0 (864 : L) (by norm_num)) = 0)
    (Q : Y.toAffine.Point) (hQ : addOrderOf Q = l)
    [hVell : (veluQuotientFull Y (((range l).erase 0).image
      (fun k : ℕ => pointCoords (k • Q)))).IsElliptic] :
    jExp p (veluQuotientFull Y
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) < 0 := by
  by_cases hdvd : (l : ℤ) ∣ jExp p Y
  · exact h1 L p Y hss hj hl hl5 h48 h864 hdvd Q hQ
  · exact jExp_neg_veluQuot_of_not_dvd p Y hss hj hl (by omega) h48 h864 hdvd Q hQ

/-! ## ★★★★★★★★★★★★★★★★★★★★組み上げ（対偶） -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**道 C の組み上げ**——①から `0 ≤ jExp p (E/⟨Q⟩)` が出る（第 1445）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`Skeleton/GenEll/VeluSemistable.lean` の `sorry` が要求していたのは
**`p ∣ l` かつ良い素点（`0 ≤ jExp p E`）のとき `0 ≤ jExp p (E/⟨Q⟩)`**、
すなわち `j(E′)` の整性であった。

★★★対偶で `jExp p E′ < 0` と仮定して矛盾を出す:

| 段 | 使うもの | 第 |
|---|---|---|
| ② 双対同種の `j` | `exists_dual_veluQuot_j_numberField` | 1442 |
| `N` の中で `p` の上の素点 | `exists_heightOneSpectrum_over_intermediate` | 本ファイル |
| `jExp` は拡大で符号を保存・反映 | `jExp_baseChange` | 在庫 |
| N1 半安定なモデル | `exists_finite_extension_semistableAt_of_jExp_neg` | 1444 |
| ① `l`-巡回商も `jExp < 0` | `jExp_neg_veluQuot_of_stable` | 本ファイル |
| Vélu の商は底変換と可換 | `veluQuotientFull_baseChange` | 1197 |
| `j` が等しければ `jExp` も等しい | `jExp_congr_j` | 1443 |

★合成体は要らない——②→N1 の順に**塔**（`L ⊆ N ⊆ M`）で取る。 -/
theorem jExp_nonneg_veluQuot_of_stable
    (h1 : JExpNegVeluStable)
    {L : Type} [Field L] [NumberField L] [inst : DecidableEq L]
    (p : HeightOneSpectrum (𝓞 L))
    (E : WeierstrassCurve L) [E.IsElliptic]
    {l : ℕ} (hl : l.Prime) (hl5 : 5 ≤ l)
    (h48 : valAdd p (Units.mk0 (48 : L) (by norm_num)) = 0)
    (h864 : valAdd p (Units.mk0 (864 : L) (by norm_num)) = 0)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    [hVell : (veluQuotientFull E (((range l).erase 0).image
      (fun k : ℕ => pointCoords (k • Q)))).IsElliptic]
    (hj : 0 ≤ jExp p E) :
    0 ≤ jExp p (veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) := by
  have hinst : inst = fun a b => Classical.propDecidable (a = b) := by
    funext a b
    exact Subsingleton.elim _ _
  subst hinst
  set E' := veluQuotientFull E
    (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))) with hE'def
  by_contra hcon
  push_neg at hcon
  -- ★★段 1（②、第 1442）——双対同種を数体の上で取る
  obtain ⟨N, hNfin, Q', hQ'ord, hZj⟩ :=
    exists_dual_veluQuot_j_numberField L E E' hl hQ hE'def
  haveI : FiniteDimensional L (N : Type) := hNfin
  letI : NumberField (N : Type) := NumberField.of_module_finite L N
  letI : IsScalarTower (𝓞 L) L (N : Type) := isScalarTower_ringOfIntegers_base L N
  letI : IsScalarTower (𝓞 L) (𝓞 (N : Type)) (N : Type) := isScalarTower_ringOfIntegers_top L N
  letI : DecidableEq (N : Type) := fun a b => Classical.propDecidable (a = b)
  obtain ⟨P, hPo⟩ := exists_heightOneSpectrum_over_intermediate L N p
  haveI := hPo
  haveI hE'Nell : (E'.baseChange (N : Type)).IsElliptic := by
    show (E'.map (algebraMap L (N : Type))).IsElliptic
    infer_instance
  haveI hENell : (E.baseChange (N : Type)).IsElliptic := by
    show (E.map (algebraMap L (N : Type))).IsElliptic
    infer_instance
  have h48N : valAdd P (Units.mk0 (48 : (N : Type)) (by norm_num)) = 0 :=
    valAdd_eq_zero_of_algebraMap L (N : Type) p P (by norm_num) (by norm_num)
      (map_ofNat _ 48) h48
  have h864N : valAdd P (Units.mk0 (864 : (N : Type)) (by norm_num)) = 0 :=
    valAdd_eq_zero_of_algebraMap L (N : Type) p P (by norm_num) (by norm_num)
      (map_ofNat _ 864) h864
  have hjN : jExp P (E'.baseChange (N : Type)) < 0 := by
    rw [jExp_baseChange L (N : Type) p P E']
    exact mul_neg_of_pos_of_neg
      (by exact_mod_cast Nat.pos_of_ne_zero (ramificationIdx_ne_zero p P)) hcon
  -- ★★段 2（N1、第 1444）——半安定なモデルを `N` の上の有限次拡大で取る
  obtain ⟨M, hMfin, hM⟩ :=
    exists_finite_extension_semistableAt_of_jExp_neg (N : Type) P
      (E'.baseChange (N : Type)) hjN
  haveI : FiniteDimensional (N : Type) (M : Type) := hMfin
  letI : NumberField (M : Type) := NumberField.of_module_finite (N : Type) M
  letI : IsScalarTower (𝓞 (N : Type)) (N : Type) (M : Type) :=
    isScalarTower_ringOfIntegers_base (N : Type) M
  letI : IsScalarTower (𝓞 (N : Type)) (𝓞 (M : Type)) (M : Type) :=
    isScalarTower_ringOfIntegers_top (N : Type) M
  letI : DecidableEq (M : Type) := fun a b => Classical.propDecidable (a = b)
  obtain ⟨P', hP'o⟩ := exists_heightOneSpectrum_over_intermediate (N : Type) M P
  haveI := hP'o
  obtain ⟨hssM, hjM⟩ := hM P' hP'o
  have h48M : valAdd P' (Units.mk0 (48 : (M : Type)) (by norm_num)) = 0 :=
    valAdd_eq_zero_of_algebraMap (N : Type) (M : Type) P P' (by norm_num) (by norm_num)
      (map_ofNat _ 48) h48N
  have h864M : valAdd P' (Units.mk0 (864 : (M : Type)) (by norm_num)) = 0 :=
    valAdd_eq_zero_of_algebraMap (N : Type) (M : Type) P P' (by norm_num) (by norm_num)
      (map_ofNat _ 864) h864N
  haveI hYell : ((E'.baseChange (N : Type)).map
      (algebraMap (N : Type) (M : Type))).IsElliptic := by
    infer_instance
  have hssM' : SemistableAt P' ((E'.baseChange (N : Type)).map
      (algebraMap (N : Type) (M : Type))) := hssM
  have hjM' : jExp P' ((E'.baseChange (N : Type)).map
      (algebraMap (N : Type) (M : Type))) < 0 := hjM
  have hQ2ord : addOrderOf (rhPoint (algebraMap (N : Type) (M : Type))
      (E'.baseChange (N : Type)) Q') = l := by
    rw [addOrderOf_rhPoint]
    exact hQ'ord
  haveI hZell : (veluQuotientFull (E'.baseChange (N : Type))
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q')))).IsElliptic :=
    isElliptic_veluQuotientFull_nsmul_nf' (N : Type) _ hQ'ord
  haveI hZMell : (veluQuotientFull ((E'.baseChange (N : Type)).map
      (algebraMap (N : Type) (M : Type)))
      (((range l).erase 0).image (fun k : ℕ =>
        pointCoords (k • rhPoint (algebraMap (N : Type) (M : Type))
          (E'.baseChange (N : Type)) Q')))).IsElliptic :=
    isElliptic_veluQuotientFull_nsmul_nf' (M : Type) _ hQ2ord
  -- ★★段 3——`M` の上で①を当て、②の `j` の等式へ戻す
  have hZbc := veluQuotientFull_baseChange (algebraMap (N : Type) (M : Type))
    (E'.baseChange (N : Type)) _ hQ'ord rfl
  have h1M := jExp_neg_veluQuot_of_stable h1 P' ((E'.baseChange (N : Type)).map
    (algebraMap (N : Type) (M : Type))) hssM' hjM' hl hl5 h48M h864M _ hQ2ord
  have hjeq : jExp P' (veluQuotientFull ((E'.baseChange (N : Type)).map
      (algebraMap (N : Type) (M : Type)))
      (((range l).erase 0).image (fun k : ℕ =>
        pointCoords (k • rhPoint (algebraMap (N : Type) (M : Type))
          (E'.baseChange (N : Type)) Q'))))
      = jExp P' ((veluQuotientFull (E'.baseChange (N : Type))
        (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q')))).map
          (algebraMap (N : Type) (M : Type))) :=
    jExp_congr_j P' _ _ (j_congr_curve hZbc.symm)
  rw [hjeq] at h1M
  have hZj2 : (veluQuotientFull (E'.baseChange (N : Type))
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q')))).j
      = (E.map (algebraMap L (N : Type))).j := hZj _ rfl
  have hjZ : jExp P (veluQuotientFull (E'.baseChange (N : Type))
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q'))))
      = jExp P (E.baseChange (N : Type)) := jExp_congr_j P _ _ hZj2
  have hjEN : 0 ≤ jExp P (E.baseChange (N : Type)) := by
    rw [jExp_baseChange L (N : Type) p P E]
    exact mul_nonneg (by positivity) hj
  have h1M2 : jExp P' ((veluQuotientFull (E'.baseChange (N : Type))
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q')))).baseChange
        (M : Type)) < 0 := h1M
  rw [jExp_baseChange (N : Type) (M : Type) P P' _, hjZ] at h1M2
  nlinarith [h1M2, hjEN,
    (by positivity : (0 : ℤ) ≤ (P.asIdeal.ramificationIdx P'.asIdeal : ℤ))]

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def jExp_neg_veluQuot_of_not_dvd.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(l ∤ v_p(j) なら Vélu の商も乗法還元。★無条件)",
    sectionId := "genell-lemma-3-5" }

def jExp_neg_veluQuot_of_not_dvd.needs : List ProofObligation :=
  [ .citation "[ABC3]" "semistableAt_veluQuot_badPrime_all(第 1436、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.semistableAt_veluQuot_badPrime_all") 1,
    .citation "[ABC3]" "minDeltaExp_eq_mul_at_bad_prime_any_K(第 1140、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_at_bad_prime_any_K") 1,
    .citation "[ABC3]" "minDeltaExp_eq_maxJ_of_semistable(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.minDeltaExp_eq_maxJ_of_semistable") 1 ]

def JExpNegVeluStable.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(乗法還元は l-同種で保たれる——l ∣ v_P(j)、すなわち深い核の場合)",
    sectionId := "genell-lemma-3-5" }

def jExp_nonneg_veluQuot_of_stable.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(道 C の組み上げ——①から j(E′) の整性が出る)",
    sectionId := "genell-lemma-3-5" }

def jExp_nonneg_veluQuot_of_stable.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_dual_veluQuot_j_numberField(②、第 1442、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_dual_veluQuot_j_numberField") 17,
    .citation "[ABC3]" "exists_finite_extension_semistableAt_of_jExp_neg(N1、第 1444、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_finite_extension_semistableAt_of_jExp_neg") 17,
    .citation "[ABC3]" "veluQuotientFull_baseChange(第 1197、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.veluQuotientFull_baseChange") 1,
    .citation "[ABC3]" "jExp_baseChange(在庫、証明済み——符号を保存かつ反映)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.jExp_baseChange") 1,
    .citation "[ABC3]" "jExp_congr_j(第 1443、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.jExp_congr_j") 1,
    .implicitStep
      ("★★★★★**2026-09-06(第 1445)**——道 C の配管はこれで全部つながった。" ++
       "☆合成体は要らなかった: ②→N1 の順に取れば `L ⊆ N ⊆ M` の**塔**になり、" ++
       "`veluQuotientFull_baseChange` で②の `j` の等式がそのまま `M` へ上がる。" ++
       "★残るのは `JExpNegVeluStable`(＝①の `l ∣ v_P(j)` の場合、深い核)だけで、" ++
       "止まった理由は**数学が足りない**(配管ではない)。" ++
       "☆`l ∤ v_P(j)` の場合は `jExp_neg_veluQuot_of_not_dvd` で閉じている。") 17 ]

end ABC3.Found.GenEll
