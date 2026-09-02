/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.VeluKernelNorm
import ABC3.Meta.Claim

/-!
# 第 1387 ブロック —— **Vélu の商の判別式の恒等式**（`Skeleton`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★これは何か——残る葉の**根を 1 本の式にした**

`Skeleton/GenEll/VeluSemistable.lean` の `semistableAt_veluQuotientFull` の
良い素点側は、第 1385-1386 によって**次の恒等式ただ 1 本**に落ちた:

    Δ(E)^l = Δ(E/C) · ( ∏_{P ∈ C∖{O}} (2 y_P + a₁ x_P + a₃) )^4

## ★★★★★★★★数値で確かめてある（2026-09-02、第 1384）

☆`tools/velu-disc-check.py` で `l = 3, 5, 7` について
Tate 標準形の族（13 例）すべて**厳密に成立**することを確かめた。

★重みの検算——`Δ` は重み 12、`2y_P + a₁x_P + a₃` は重み 3 なので
右辺第 2 因子は `3(l−1)·4 = 12(l−1)`、合計 `12 + 12(l−1) = 12l` ✓。

★★☆**部分群であることは本質である**——核でない `±`-閉集合
（例: `y² = x³+1` の非捩れ点 `(2,3)`）では成り立たない（比 `1/797121`）。
☆したがって自由な媒介変数の多項式恒等式ではなく、
`genericParam`・`c4DenomFree` の安い手は使えない。

## ★★★★★★★★証明の道

★`ℂ` で証明して埋め込みで降ろす（第 1334 `isElliptic_velu_congr_curve` と同じ型）。
☆`ℂ` では一意化（第 1330-1335）と `latticeCurve_eq_veluQuotientFull`（在庫）があるので、
`Δ(Λ)` の積公式と `℘′(u) = −σ(2u)/σ(u)⁴` から出るはずである。
★★もう 1 つの道は「`ℂ` 上ではモジュラー群の作用で核を `μ_l` 型に直せる」
——そうすれば第 1131（`Δ(veluCurve (tateCurveAt q)) = l¹²·Δ(tateCurveAt (q^l))`、証明済み）
がそのまま当たる。
-/

namespace ABC3.Skeleton.GenEll

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GenEll ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] Vélu の商の判別式の恒等式**（第 1387）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

    Δ(E)^l = Δ(E/C) · ( ∏_{P ∈ C∖{O}} (2 y_P + a₁ x_P + a₃) )^4

☆`l = 3, 5, 7` について数値で確かめてある（`tools/velu-disc-check.py`）。

★★★これが `semistableAt_veluQuotientFull` の良い素点側に残る**ただ 1 本の節点**である。 -/
theorem disc_pow_eq_veluQuot_mul {L : Type} [Field L] [DecidableEq L]
    (E : WeierstrassCurve L) [E.IsElliptic]
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l) :
    E.Δ ^ l
      = (veluQuotientFull E
          (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).Δ
        * (veluKernelNorm E
          (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) ^ 4 := by
  sorry

/-- ★★★★★★★★★★★★★★★★★★★★
**良い素点では Vélu の商は半安定**——★恒等式（第 1387）から（第 1387）。

☆第 1385（恒等式 ⟹ 半安定）と第 1386（`N` の整性）に流すだけである。 -/
theorem semistableAt_veluQuot_good {L : Type} [Field L] [NumberField L]
    [inst : DecidableEq L]
    (p : HeightOneSpectrum (𝓞 L)) (E : WeierstrassCurve L) [E.IsElliptic]
    [WeierstrassCurve.IsIntegral (primeSubring p) E]
    (hΔ0 : valAdd p (Units.mk0 E.Δ E.isUnit_Δ.ne_zero) = 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2) (hlu : IsUnit ((l : primeSubring p)))
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (hΔ' : (veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).Δ ≠ 0) :
    SemistableAt p (veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) := by
  have hinst : inst = fun a b => Classical.propDecidable (a = b) := by
    funext a b
    exact Subsingleton.elim _ _
  subst hinst
  set S := ((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)) with hS
  have hid := disc_pow_eq_veluQuot_mul E hl hodd Q hQ
  haveI hintE' : WeierstrassCurve.IsIntegral (primeSubring p) (veluQuotientFull E S) :=
    isIntegral_veluQuotientFull_of_addOrderOf_prime p E hl hlu Q hQ
  have hNint : veluKernelNorm E S ∈ primeSubring p :=
    veluKernelNorm_mem_primeSubring p E S (veluKernel_coords_mem p E hl hlu Q hQ)
  have hNne : veluKernelNorm E S ≠ 0 := by
    intro h0
    rw [h0] at hid
    simp only [ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, mul_zero] at hid
    exact pow_ne_zero l E.isUnit_Δ.ne_zero hid
  exact semistableAt_of_disc_pow_eq p E (veluQuotientFull E S)
    E.isUnit_Δ.ne_zero hΔ' hintE' hΔ0 hNne hNint hid

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def disc_pow_eq_veluQuot_mul.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商の判別式の恒等式)",
    sectionId := "genell-lemma-3-5" }

def disc_pow_eq_veluQuot_mul.needs : List ProofObligation :=
  [ .citation "[ABC3]" "latticeCurve_eq_veluQuotientFull(ℂ 側の Vélu、在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.latticeCurve_eq_veluQuotientFull") 1,
    .citation "[ABC3]" "vAdd_Delta_veluCurve_tate(Tate 曲線での同型の式、第 1131、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.vAdd_Delta_veluCurve_tate") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1387）**——`l = 3, 5, 7` について数値で確かめてある" ++
       "（`tools/velu-disc-check.py`、13 例すべて厳密に成立）。" ++
       "☆核でない ±閉集合では成り立たないので、部分群構造が本質である" ++
       "——自由な媒介変数の多項式恒等式ではない。" ++
       "★証明の道は `ℂ` で証明して埋め込みで降ろす（第 1334 と同じ型）。") 17 ]

def semistableAt_veluQuot_good.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(良い素点では Vélu の商は半安定——恒等式から)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuot_good.needs : List ProofObligation :=
  [ .citation "[ABC3]" "semistableAt_of_disc_pow_eq(第 1385、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.semistableAt_of_disc_pow_eq") 1,
    .citation "[ABC3]" "veluKernelNorm_mem_primeSubring(第 1386、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.veluKernelNorm_mem_primeSubring") 1,
    .citation "[ABC3]" "isIntegral_veluQuotientFull_of_addOrderOf_prime(第 1074、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isIntegral_veluQuotientFull_of_addOrderOf_prime") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1387）**——良い素点側は**恒等式 1 本**に落ちた。" ++
       "☆`p ∣ l` の場合は `hlu` が無いので別扱いである。") 17 ]

end ABC3.Skeleton.GenEll
