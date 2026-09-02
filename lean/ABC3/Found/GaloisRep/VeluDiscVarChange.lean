/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.VeluDiscDescent
import ABC3.Found.GenEll.VeluPointSet
import ABC3.Found.GenEll.VeluSetCard
import ABC3.Meta.Claim

/-!
# 第 1391 ブロック —— **判別式の恒等式は変数変換で不変**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——`ℂ` の節点を**格子曲線**まで落とす

第 1390 で節点は `ℂ` の上だけになった。★`ℂ` では一意化（第 1330-1335）が
**変数変換で格子曲線に直す**ので、恒等式が変数変換で不変なら
節点は**格子曲線の場合**まで落ちる。

☆重みの検算がそのまま証明になる:

| 量 | 変数変換 `C = (u,r,s,t)` での変化 |
|---|---|
| `Δ` | `u⁻¹² · Δ` |
| `2y + a₁x + a₃` | `u⁻³ · (2y + a₁x + a₃)` |

★したがって `Δ^l ↦ u^{−12l}Δ^l`、`Δ′·N⁴ ↦ u^{−12}Δ′ · u^{−12(l−1)}N⁴ = u^{−12l}Δ′N⁴`
——**両辺が同じ因子でずれる**ので恒等式は保たれる（`|C∖{O}| = l−1` を使う）。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Finset ABC3.Meta ABC3.Found.GenEll

open scoped Classical

variable {F : Type} [Field F]

/-- ★★★★★★**`2y + a₁x + a₃` は変数変換で `u⁻³` 倍**——★**無条件**（第 1391）。 -/
theorem two_y_add_variableChange (C : WeierstrassCurve.VariableChange F)
    (W : WeierstrassCurve F) (x y : F) :
    2 * (ABC3.Found.GenEll.vcY C x y) + (C • W).a₁ * (ABC3.Found.GenEll.vcX C x) + (C • W).a₃
      = ((C.u⁻¹ : Fˣ) : F) ^ 3 * (2 * y + W.a₁ * x + W.a₃) := by
  rw [ABC3.Found.GenEll.vcX, ABC3.Found.GenEll.vcY, WeierstrassCurve.variableChange_a₁, WeierstrassCurve.variableChange_a₃]
  have hu : ((C.u⁻¹ : Fˣ) : F) ^ 3 = ((C.u⁻¹ : Fˣ) : F) * ((C.u⁻¹ : Fˣ) : F) ^ 2 := by ring
  ring

set_option maxHeartbeats 800000 in
/-- ★★★★★★★★**核のノルムは変数変換で `u^{−3|S|}` 倍**——★**無条件**（第 1391）。 -/
theorem veluKernelNorm_variableChange (C : WeierstrassCurve.VariableChange F)
    (W : WeierstrassCurve F) (S : Finset (F × F)) :
    veluKernelNorm (C • W) (S.image (fun z : F × F => (ABC3.Found.GenEll.vcX C z.1, ABC3.Found.GenEll.vcY C z.1 z.2)))
      = ((C.u⁻¹ : Fˣ) : F) ^ (3 * S.card) * veluKernelNorm W S := by
  have hu0 : ((C.u⁻¹ : Fˣ) : F) ≠ 0 := Units.ne_zero _
  have hinj : ∀ a ∈ S, ∀ b ∈ S,
      (ABC3.Found.GenEll.vcX C a.1, ABC3.Found.GenEll.vcY C a.1 a.2) = (ABC3.Found.GenEll.vcX C b.1, ABC3.Found.GenEll.vcY C b.1 b.2) → a = b := by
    intro a _ b _ hab
    have h1 : ABC3.Found.GenEll.vcX C a.1 = ABC3.Found.GenEll.vcX C b.1 := congrArg Prod.fst hab
    have h2 : ABC3.Found.GenEll.vcY C a.1 a.2 = ABC3.Found.GenEll.vcY C b.1 b.2 := congrArg Prod.snd hab
    have hx : a.1 = b.1 := ABC3.Found.GenEll.vcX_injective C h1
    have hy : a.2 = b.2 := by
      rw [hx] at h2
      exact ABC3.Found.GenEll.vcY_injective C b.1 h2
    exact Prod.ext hx hy
  rw [veluKernelNorm, veluKernelNorm, Finset.prod_image hinj]
  rw [Finset.prod_congr rfl (fun q _ => two_y_add_variableChange C W q.1 q.2),
    Finset.prod_mul_distrib, Finset.prod_const, ← pow_mul]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★
**判別式の恒等式は変数変換で不変**——★**無条件**（第 1391）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★これで `ℂ` の節点は**格子曲線の場合**まで落ちる
——一意化（第 1330-1335）が変数変換で格子曲線に直すからである。 -/
theorem disc_pow_eq_veluQuot_mul_of_variableChange [instF : DecidableEq F]
    (C : WeierstrassCurve.VariableChange F) (E : WeierstrassCurve F)
    [E.IsElliptic] [(C • E).IsElliptic]
    {l : ℕ} (hl : l.Prime) {Q : E.toAffine.Point} (hQ : addOrderOf Q = l) (h2 : (2 : F) ≠ 0)
    (h : (C • E).Δ ^ l
      = (veluQuotientFull (C • E)
          (((range l).erase 0).image (fun k : ℕ => pointCoords (k • ABC3.Found.GenEll.vcPoint C E Q)))).Δ
        * (veluKernelNorm (C • E)
          (((range l).erase 0).image (fun k : ℕ => pointCoords (k • ABC3.Found.GenEll.vcPoint C E Q)))) ^ 4) :
    E.Δ ^ l
      = (veluQuotientFull E
          (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).Δ
        * (veluKernelNorm E
          (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) ^ 4 := by
  have hinstF : instF = fun a b => Classical.propDecidable (a = b) := by
    funext a b
    exact Subsingleton.elim _ _
  subst hinstF
  set K := ((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)) with hK
  have hu0 : ((C.u⁻¹ : Fˣ) : F) ≠ 0 := Units.ne_zero _
  have hK' := image_pointCoords_vcPoint_nsmul C E hQ
  have hquot : veluQuotientFull (C • E)
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • ABC3.Found.GenEll.vcPoint C E Q)))
      = C • (veluQuotientFull E K) :=
    (veluQuotientFull_vcPoint_eq C E _ hQ h2 rfl).symm
  have hcard : K.card + 1 = l := card_image_pointCoords_nsmul hl hQ
  rw [hquot, hK', veluKernelNorm_variableChange C E K,
    WeierstrassCurve.variableChange_Δ, WeierstrassCurve.variableChange_Δ,
    mul_pow, mul_pow, ← pow_mul, ← pow_mul] at h
  -- ★`u^{12l}` を両辺から落とす
  have hexp : 12 * l = 12 + 3 * K.card * 4 := by omega
  have h' : ((C.u⁻¹ : Fˣ) : F) ^ (12 * l) * E.Δ ^ l
      = ((C.u⁻¹ : Fˣ) : F) ^ (12 * l)
        * ((veluQuotientFull E K).Δ * veluKernelNorm E K ^ 4) := by
    rw [h, hexp, pow_add]
    ring
  exact mul_left_cancel₀ (pow_ne_zero _ hu0) h'

/-! ## ★出典の紐付け(`.src`) -/

def two_y_add_variableChange.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(2y + a₁x + a₃ は変数変換で u⁻³ 倍。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluKernelNorm_variableChange.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(核のノルムは変数変換で u^{−3|S|} 倍。★無条件)",
    sectionId := "genell-lemma-3-5" }

def disc_pow_eq_veluQuot_mul_of_variableChange.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(判別式の恒等式は変数変換で不変。★無条件)",
    sectionId := "genell-lemma-3-5" }

def disc_pow_eq_veluQuot_mul_of_variableChange.needs : List ProofObligation :=
  [ .citation "[ABC3]" "veluQuotientFull_vcPoint_eq(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.veluQuotientFull_vcPoint_eq") 1,
    .citation "[ABC3]" "card_image_pointCoords_nsmul(第 1240、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.card_image_pointCoords_nsmul") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1391）**——重みの検算がそのまま証明になる。" ++
       "☆`Δ ↦ u⁻¹²Δ`・`2y+a₁x+a₃ ↦ u⁻³(2y+a₁x+a₃)` なので両辺が同じ因子でずれる" ++
       "（`|C∖{O}| = l−1` を使う）。" ++
       "★★★これで `ℂ` の節点は格子曲線の場合まで落ちる。") 17 ]

end ABC3.Found.GaloisRep
