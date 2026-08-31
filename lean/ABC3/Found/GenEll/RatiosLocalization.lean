/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.LocalChartRatio
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★局所での比の組を座標から作る（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★★★★★これは何か —— 代数の側を 1 本にまとめる

`§9-943`〜`§9-949` の道具立てを使うには、素点 `Q` ごとに

* 添字 `j`（座標の最小割り切り成分）
* 比の組 `r : Fin (N+1) → (𝓞_F)_Q`（`r_j = 1`、`r_k = x_k/x_j`）
* `𝔞_Q = (x_j)`
* `r_0·x_j = x_0`
* ★**生成点での整合** `σ(r_k)·c_j = c_k`（`c` は生成点の正規化座標）

が要る。★★★本ファイルはそれを **`§9-941` の同次座標から一気に作る**。

## ★★★機構

1. `(𝓞_F)_Q` は離散付値環（`§9-943` と同じ理由）——割り切りが全順序
2. だから `∃ j, ∀ k, x_j ∣ x_k`（`(𝓞_F)_Q` の中で）——`r_k` はその商
3. `r_j = 1` は `x_j ≠ 0` と消去律から（`𝓞_F → (𝓞_F)_Q` は**単射**）
4. ★生成点での整合は `r_k·x_j = x_k` を `F` へ送り、
   `x_k = c_k·x_i`（`§9-941`）を代入して `x_i ≠ 0` で割るだけ

## ★残っている段（明示）

★★残るのは**幾何の側との接続**——`§9-947`（局所化した点の同定）に
本ファイルの `r` を渡し、`§9-948`（`X_{s_j}` を通る）と
`§9-949`（`(s_0/s_j)(y_Q)·x_j = x_0`）を繋いで
`§9-940` の有限素点の整合を丸ごと出すこと——である。
★道具は全部揃っており、残るのは**配線**（`IsScalarTower` で生成点を同定すること）だけである。
-/

namespace ABC3.Found.GenEll

open NumberField

/-! ## ★★★★★★★★★★★★★★★★★★★★局所での比の組 -/

/-- ★★★★★★★★★★★★★★★★★★★★**同次座標から素点ごとの比の組を作る**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★これが `§9-943`〜`§9-949` の道具立てに渡すデータ一式である:
`j`（最小割り切り成分）、`r`（比の組、`r_j = 1`）、`𝔞_Q = (x_j)`、
`r_0·x_j = x_0`、そして★**生成点での整合** `σ(r_k)·c_j = c_k`。

★★機構は「付値環では割り切りが全順序」＋「`𝓞_F → (𝓞_F)_Q` は単射」＋
「`x_k = c_k·x_i` を `F` へ送って `x_i ≠ 0` で割る」だけである。 -/
theorem exists_ratios_localization (F : Type) [Field F] [NumberField F]
    {N : ℕ} (x : Fin (N+1) → 𝓞 F) (i : Fin (N+1)) (hxi : x i ≠ 0)
    (c : Fin (N+1) → F) (hrel : ∀ k, ((x k : F)) = c k * ((x i : F)))
    (Q : Ideal (𝓞 F)) (hQ : Q.IsMaximal) :
    haveI := hQ.isPrime
    ∃ (j : Fin (N+1)) (r : Fin (N+1) → Localization.AtPrime Q), r j = 1 ∧
      c j ≠ 0 ∧
      Ideal.map (algebraMap (𝓞 F) (Localization.AtPrime Q)) (Ideal.span (Set.range x))
        = Ideal.span {algebraMap (𝓞 F) (Localization.AtPrime Q) (x j)} ∧
      r 0 * algebraMap (𝓞 F) (Localization.AtPrime Q) (x j)
        = algebraMap (𝓞 F) (Localization.AtPrime Q) (x 0) ∧
      ∀ k, algebraMap (Localization.AtPrime Q) F (r k) * c j = c k := by
  haveI := hQ.isPrime
  haveI hQ0 : Q ≠ ⊥ := Ring.ne_bot_of_isMaximal_of_not_isField hQ (RingOfIntegers.not_isField F)
  haveI : IsDiscreteValuationRing (Localization.AtPrime Q) :=
    IsLocalization.AtPrime.isDiscreteValuationRing_of_dedekind_domain (𝓞 F) hQ0 _
  have hinj : Function.Injective (algebraMap (𝓞 F) (Localization.AtPrime Q)) :=
    IsLocalization.injective _ Q.primeCompl_le_nonZeroDivisors
  have htower : ∀ a : 𝓞 F,
      algebraMap (Localization.AtPrime Q) F (algebraMap (𝓞 F) (Localization.AtPrime Q) a)
        = ((a : F)) := by
    intro a; rw [← IsScalarTower.algebraMap_apply]
  obtain ⟨j, hj⟩ := exists_dvd_all_of_valuationRing
    (fun k => algebraMap (𝓞 F) (Localization.AtPrime Q) (x k))
  choose r hr using fun k => hj k
  have hxiF : ((x i : F)) ≠ 0 := by
    simpa using (fun h => hxi (by exact_mod_cast h) : ((x i : F)) ≠ 0)
  have hxjne : algebraMap (𝓞 F) (Localization.AtPrime Q) (x j) ≠ 0 := by
    intro h0
    have h := hr i
    rw [h0, zero_mul] at h
    exact hinj.ne hxi (by rw [h, map_zero])
  have hrj : r j = 1 := by
    have h := hr j
    nth_rewrite 1 [← mul_one (algebraMap (𝓞 F) (Localization.AtPrime Q) (x j))] at h
    exact (mul_left_cancel₀ hxjne h).symm
  have hcj : c j ≠ 0 := by
    intro h0
    apply hxjne
    have h := hrel j
    rw [h0, zero_mul] at h
    rw [show (0 : F) = ((0 : 𝓞 F) : F) from by simp] at h
    have hz : x j = 0 := by exact_mod_cast h
    rw [hz, map_zero]
  refine ⟨j, r, hrj, hcj, ?_, ?_, ?_⟩
  · rw [Ideal.map_span, ← Set.range_comp]
    exact span_range_eq_span_singleton_of_dvd _ j hj
  · rw [hr 0]; ring
  · intro k
    have h1 := hr k
    have h2 := congrArg (algebraMap (Localization.AtPrime Q) F) h1
    rw [map_mul, htower, htower] at h2
    rw [hrel k, hrel j] at h2
    have h3 : algebraMap (Localization.AtPrime Q) F (r k) * c j * ((x i : F))
        = c k * ((x i : F)) := by rw [h2]; ring
    exact mul_right_cancel₀ hxiF h3

/-! ## ★出典の紐付け(`.src`) -/

def exists_ratios_localization.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(同次座標から素点ごとの比の組を作る)",
    sectionId := "genell-prop-1-4" }

def exists_ratios_localization.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_dvd_all_of_valuationRing(付値環では割り切りが全順序、§9-943)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_dvd_all_of_valuationRing") 2,
    .citation "[mathlib]" "IsLocalization.injective(𝓞_F → (𝓞_F)_Q は単射)"
      (.inMathlib "IsLocalization.injective") 2,
    .citation "[ABC3]" "exists_homogeneous_coords(同次座標の構成、§9-941)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_homogeneous_coords") 2,
    .implicitStep
      ("★★★★測定(2026-08-29): §9-943〜949 の道具立てに渡すデータ一式" ++
       "(j・r・𝔞_Q = (x_j)・r_0·x_j = x_0・生成点での整合)は" ++
       "**§9-941 の同次座標から一気に作れる**。" ++
       "生成点での整合は r_k·x_j = x_k を F へ送り、x_k = c_k·x_i を代入して " ++
       "x_i ≠ 0 で割るだけである") 5,
    .implicitStep
      ("★残るのは幾何の側との接続——§9-947(局所化した点の同定)に本ファイルの r を渡し、" ++
       "§9-948(X_{s_j} を通る)と §9-949((s_0/s_j)(y_Q)·x_j = x_0)を繋いで " ++
       "§9-940 の有限素点の整合を丸ごと出すこと。" ++
       "★道具は全部揃っており、残るのは配線(IsScalarTower で生成点を同定)だけである") 4 ]

end ABC3.Found.GenEll
