/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.VeluKernelNorm
import ABC3.Found.GaloisRep.VeluDiscDescent
import ABC3.Found.GaloisRep.VeluDiscVarChange
import ABC3.Found.GenEll.JSurjective
import ABC3.Found.GenEll.VeluEllipticNF
import ABC3.Found.GenEll.VeluDiscLattice
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

## ★★★★★★★★★★★★★★★★閉じた（2026-09-02、第 1397-1402）

★★★**σ 函数も q 展開も因子の理論も要らなかった**。道は 6 段:

| 段 | 内容 | 番 |
|---|---|---|
| 1 | `Δ = 16·((e₁−e₂)(e₁−e₃)(e₂−e₃))²`（`e_i` は半周期の `℘` の値） | 第 1397 |
| 2 | `∏_{i≠j}(℘(v_j+w) − e_i) = −D`（**`w` に依らない**） | 第 1398 |
| 3 | 代表系の置換（平行移動と**負号**） | 第 1399 |
| 4 | **同種のノルム** `∏_{w∈T}(℘(z+w) − e_i) = c_i(℘_{Λ′}(z) − e′_i)` | 第 1400 |
| 5 | 帳簿（`D^l = (c₁c₂c₃)²D′`・`N² = 4^{l−1}c₁c₂c₃`） | 第 1401 |
| 6 | 語彙の変換（`⟨Q⟩` の像へ） | 第 1402 |

☆最も重い (4) は **Liouville（第 598、在庫）と `R` の偶関数性（第 1399）だけ**で出た
——`R − c` が原点で 2 位で消えるので極が打ち消し合う。

## ☆参考：当初測った 3 つの道（2026-09-02、第 1396）

☆どれも「新しい理論を 1 つ建てる」規模である。全文は
`ResearchPaper/mathlib-gap.json` の `veluDiscIdentityRoutes20260902`。

| 道 | 中身 | 状態 |
|---|---|---|
| (i) 解析（`ℂ`） | `Δ(Λ)` の η 積公式・`℘′(w) = −σ(2w)/σ(w)⁴` | ★**最短**——`ℂ` 側の足場（`Uniformization.lean` 5700 行）が既にある。mathlib に σ 函数も `Δ` の q 展開も無い |
| (ii) 代数（終結式） | `N = ±Res(ψ_C, 4x³+b₂x²+2b₄x+b₆)`・`Δ′ = Δ^l/Res⁴`（Dewaghe）。`X(x)−X(x′) = ∏_{P∈C}(x−x(x′+P))/ψ(x)²` と `Δ = 16∏_{i<j}(e_i−e_j)²` を使う | ☆分解体と因子の理論が要る |
| (iii) 剰余体の Vélu／NOS／モジュラー多項式 `Φ_l` | いずれも mathlib に無い | ☆丸ごと建てることになる |

★☆**自由な媒介変数の多項式恒等式ではない**ことは実験で確かめてある
（核でない `±` 閉集合では比 `1/797121`）——`genericParam` の安い手は使えない。
-/

namespace ABC3.Skeleton.GenEll

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GenEll ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

/-- ★★★★**曲線が等しければ恒等式は移る**——★**無条件**（第 1392）。 -/
theorem disc_pow_eq_velu_congr_curve {F : Type} [Field F] [DecidableEq F]
    {W₁ W₂ : WeierstrassCurve F} (h : W₁ = W₂) (Q : W₁.toAffine.Point) (l : ℕ) :
    (W₁.Δ ^ l
      = (veluQuotientFull W₁
          (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).Δ
        * (veluKernelNorm W₁
          (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) ^ 4)
      ↔ (W₂.Δ ^ l
      = (veluQuotientFull W₂
          (((range l).erase 0).image
            (fun k : ℕ => pointCoords (k • (h ▸ Q : W₂.toAffine.Point))))).Δ
        * (veluKernelNorm W₂
          (((range l).erase 0).image
            (fun k : ℕ => pointCoords (k • (h ▸ Q : W₂.toAffine.Point))))) ^ 4) := by
  subst h
  exact Iff.rfl

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] Vélu の商の判別式の恒等式——格子曲線の上**（第 1392）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

    Δ(E)^l = Δ(E/C) · ( ∏_{P ∈ C∖{O}} (2 y_P + a₁ x_P + a₃) )^4

★★★★**2026-09-02（第 1402）に閉じた**——`ABC3.Found.GenEll.disc_pow_eq_lattice`。
☆σ 函数も q 展開も要らず、Liouville（第 598）と `R` の偶関数性（第 1399）で出た。 -/
theorem disc_pow_eq_veluQuot_mul_lattice (P : PeriodPair)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (Q : (latticeCurve P).toAffine.Point) (hQ : addOrderOf Q = l) :
    (latticeCurve P).Δ ^ l
      = (veluQuotientFull (latticeCurve P)
          (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).Δ
        * (veluKernelNorm (latticeCurve P)
          (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) ^ 4 := by
  obtain ⟨m, hm⟩ := hl.odd_of_ne_two hodd
  exact ABC3.Found.GenEll.disc_pow_eq_lattice P (show l = 2 * m + 1 by omega) Q hQ

/-- ★★★★★★★★★★★★★★★★★★★★
**`ℂ` の上の恒等式**——格子曲線から降りる（第 1392）。

☆一意化（第 1330-1335）で変数変換して格子曲線に直し、
変数変換不変性（第 1391）で戻す。 -/
theorem disc_pow_eq_veluQuot_mul_complex
    (E : WeierstrassCurve ℂ) [hell : E.IsElliptic]
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l) :
    E.Δ ^ l
      = (veluQuotientFull E
          (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).Δ
        * (veluKernelNorm E
          (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) ^ 4 := by
  obtain ⟨P, C, hCP⟩ := exists_periodPair_of_isElliptic E hell
  haveI hCE : (C • E).IsElliptic := by rw [hCP]; exact isElliptic_latticeCurve' P
  refine ABC3.Found.GaloisRep.disc_pow_eq_veluQuot_mul_of_variableChange C E hl hQ
    (by norm_num) ?_
  rw [disc_pow_eq_velu_congr_curve hCP (ABC3.Found.GenEll.vcPoint C E Q) l]
  refine disc_pow_eq_veluQuot_mul_lattice P hl hodd _ ?_
  rw [addOrderOf_congr_curve hCP, addOrderOf_vcPoint, hQ]

/-- ★★★★★★★★★★★★★★★★★★★★
**数体の上での恒等式**——`ℂ` の側から降りる（第 1390）。

☆第 1334（`isElliptic_veluQuotientFull_nsmul_of_embed`）と同じ型の降下である。 -/
theorem disc_pow_eq_veluQuot_mul {L : Type} [Field L] [NumberField L] [inst : DecidableEq L]
    (E : WeierstrassCurve L) [E.IsElliptic]
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l) :
    E.Δ ^ l
      = (veluQuotientFull E
          (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).Δ
        * (veluKernelNorm E
          (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) ^ 4 := by
  have hinst : inst = fun a b => Classical.propDecidable (a = b) := by
    funext a b
    exact Subsingleton.elim _ _
  subst hinst
  set f := embedComplex L with hf
  haveI hW1 : (E.map f).IsElliptic :=
    ⟨isUnit_iff_ne_zero.2 (by
      rw [WeierstrassCurve.map_Δ]
      exact (map_ne_zero_iff f f.injective).2 E.isUnit_Δ.ne_zero)⟩
  refine ABC3.Found.GaloisRep.disc_pow_eq_of_embed f E _ ?_
  rw [← image_pointCoords_rhPoint_nsmul f E hQ]
  have hQ₁ : addOrderOf (rhPoint f E Q) = l := by rw [addOrderOf_rhPoint, hQ]
  exact disc_pow_eq_veluQuot_mul_complex (E.map f) hl hodd (rhPoint f E Q) hQ₁

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

def disc_pow_eq_velu_congr_curve.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(曲線が等しければ恒等式は移る。★無条件)",
    sectionId := "genell-lemma-3-5" }

def disc_pow_eq_velu_congr_curve.needs : List ProofObligation := []

def disc_pow_eq_veluQuot_mul_lattice.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商の判別式の恒等式——格子曲線の上)",
    sectionId := "genell-lemma-3-5" }

def disc_pow_eq_veluQuot_mul_lattice.needs : List ProofObligation :=
  [ .citation "[ABC3]" "latticeCurve_eq_veluQuotientFull(ℂ 側の Vélu、在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.latticeCurve_eq_veluQuotientFull") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1392）**——★★★**これが残るただ 1 つの節点である**。" ++
       "☆第 1390（埋め込みで降りる）と第 1391（変数変換で不変）で、" ++
       "数体の場合も `ℂ` の場合もここに帰着した。" ++
       "☆`l = 3, 5, 7` について数値で確かめてある（`tools/velu-disc-check.py`、13 例）。" ++
       "★道具は在庫にある——`latticeCurve_eq_veluQuotientFull` が " ++
       "`g₂`・`g₃` の Vélu の式を与えるので、" ++
       "`Δ = g₂³ − 27g₃²` と `∏ ℘′(w)` の関係を取ればよい。") 17 ]

def disc_pow_eq_veluQuot_mul_complex.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商の判別式の恒等式——`ℂ` の上)",
    sectionId := "genell-lemma-3-5" }

def disc_pow_eq_veluQuot_mul_complex.needs : List ProofObligation :=
  [ .citation "[ABC3]" "disc_pow_eq_veluQuot_mul_of_variableChange(第 1391、証明済み。格子曲線への帰着)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.disc_pow_eq_veluQuot_mul_of_variableChange") 1,
    .citation "[ABC3]" "exists_periodPair_of_isElliptic(一意化、在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_periodPair_of_isElliptic") 1,
    .citation "[ABC3]" "latticeCurve_eq_veluQuotientFull(ℂ 側の Vélu、在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.latticeCurve_eq_veluQuotientFull") 1,
    .citation "[ABC3]" "vAdd_Delta_veluCurve_tate(Tate 曲線での同型の式、第 1131、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.vAdd_Delta_veluCurve_tate") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1390）**——第 1390 の降下で**残るのは `ℂ` の側だけ**になった。" ++
       "☆`l = 3, 5, 7` について数値で確かめてある（`tools/velu-disc-check.py`、13 例）。" ++
       "★恒等式は変数変換で不変（`Δ ↦ u⁻¹²Δ`、`N ↦ u^{−3(l−1)}N` で重みが合う）なので、" ++
       "☆一意化（第 1330-1335）で**格子曲線の場合に帰着できる**——変数変換不変性は第 1391 で**証明済み**である。") 17 ]

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
