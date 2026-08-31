/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateDSeries
import ABC3.Found.GaloisRep.TateVelu
import ABC3.Meta.Claim

/-!
# `Skeleton` —— **★★★★★★★★★★Tate 曲線の ODE と `v = ∑_ζ DY`**

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★★★★★★★★★★★★★葉 1 は 2 つの節点に落ちた

第 846 で `v = ∑_ζ DY` が見えた（数値確認済み）。その根拠は ODE 一つである:

    `DY = 3X² − Y + a₄`        …… `tate_ode`（本ファイル）

これを `veluV2 = 3x² + a₄ − y`（`veluV2_tateCurveAt`、証明済み）に入れると

    `veluV2(X_i, Y_i) = DY(ζ^i)`

がそのまま出る（`veluV2_eq_tateDYpair`、★本ファイルで**証明済み**）。
★★★★**`X²` は消え、畳み込み（Besge の恒等式）は要らなくなった**。

## ★★★★★★`tate_ode` はどう出るか

`tate_equation`（証明済み）は `(2Y+X)² = 4X³ + X² + 4a₄X + 4a₆` であり、
`tateDXpair_eq`（第 846、証明済み）は `DX = 2Y + X` である。よって

    `(DX)² = 4X³ + X² + 4a₄X + 4a₆`

万有な環 `TateUniv` の上でこれを `D`（第 845）で微分すると

    `2·DX·D²X = (12X² + 2X + 4a₄)·DX`

★`TateUniv` は整域（UFD `ℤ[A,W]` の局所化）なので `2DX` で割れて

    `D²X = 6X² + X + 2a₄`,  すなわち  `DY = 3X² − Y + a₄`
    （`D²X = D(2Y+X) = 2DY + DX = 2DY + 2Y + X` を使う）

☆級数は `TateUniv` に無く完備化の中にあるので、実際には `q^N` で切り詰めて
`TateUniv` の中で整除性を見る——`tate_equation` を通した道（第 222–229）と同じ形である。

## ★★★★★定数項の指標和

`q = 0` では `DY(u,0) = u²(2+u)/(1−u)⁴` なので、`v` の定数項は

    `240 · ∑_{ζ≠1} ζ²(2+ζ)/(1−ζ)⁴ = l⁴ − 1`     …… `sum_mu_dyterm`（本ファイル）

★`l = 2, 3, 5, 7, 11, 13` で厳密一致（2026-08-31）。
☆第 835 が測った `v(0) = (l⁴−1)/240` と同じ値である。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Found.GaloisRep ABC3.Found.GenEll Finset

/-- **[GenEll] 葉 1 の中核**——**Tate 曲線の ODE**。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★`tate_equation`（証明済み）と `tateDXpair_eq`（第 846、証明済み）から、
万有な環の上で `D` により微分して `2DX` で割ると出る。 -/
theorem tate_ode {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
    (a w q : R) (hq : q ∈ I) (haw : a * w = q) (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) :
    tateDYpair a w q hq
      = 3 * tateXpair a w q hq ^ 2 - tateYpair a w q hq + (tateCurveAt q hq).a₄ := by
  sorry

/-- ★★★★★★★★★★**`veluV2` は `DY` そのものである**——★ODE から**直ちに**出る。

これが第 846 の発見「`v = ∑_ζ DY`」の中身であり、`X²` が消える理由である。 -/
theorem veluV2_eq_tateDYpair {R : Type} [CommRing R] [IsDomain R] {I : Ideal R}
    [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (haw : a * w = q)
    (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) :
    veluV2 (tateCurveAt q hq) (tateXpair a w q hq) (tateYpair a w q hq)
      = tateDYpair a w q hq := by
  rw [veluV2_tateCurveAt, tate_ode a w q hq haw ha hw]
  ring

/-- ★★★★★★★★**`v = ∑_{i≠0} DY(ζ^i)`**——`c4_velu_tate` の左辺の `v` を置き換える形。 -/
theorem sum_veluV2_eq_sum_tateDYpair {R : Type} [CommRing R] [IsDomain R] {I : Ideal R}
    [IsAdicComplete I R] {l : ℕ} (hl : 0 < l) {ζ : R} (hζl : ζ ^ l = 1)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (q : R) (hq : q ∈ I) :
    (∑ i ∈ (range l).erase 0,
        veluV2 (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
          (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
      = ∑ i ∈ (range l).erase 0, tateDYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq := by
  refine Finset.sum_congr rfl (fun i hi => ?_)
  have hpow : (ζ ^ i) * (ζ ^ i) ^ (l - 1) = 1 := by
    rw [← pow_succ']
    rw [Nat.sub_add_cancel hl, ← pow_mul, mul_comm, pow_mul, hζl, one_pow]
  have haw : (ζ ^ i) * (q * (ζ ^ i) ^ (l - 1)) = q := by
    calc (ζ ^ i) * (q * (ζ ^ i) ^ (l - 1)) = q * ((ζ ^ i) * (ζ ^ i) ^ (l - 1)) := by ring
      _ = q := by rw [hpow, mul_one]
  have hwu : IsUnit (1 - q * (ζ ^ i) ^ (l - 1)) :=
    isUnit_one_sub (I := I) (Ideal.mul_mem_right _ _ hq)
  exact veluV2_eq_tateDYpair _ _ q hq haw (hu i hi) hwu

/-- **[GenEll] 葉 1 の定数項**——`240·∑_{ζ≠1} ζ²(2+ζ)/(1−ζ)⁴ = l⁴ − 1`。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★`l = 2, 3, 5, 7, 11, 13` で厳密一致（2026-08-31、第 846）。
☆`MuCharSum.lean` の道具（`P(y) = y^l − (y−1)^l` の Vieta と Newton）で出る。 -/
theorem sum_mu_dyterm {R : Type} [CommRing R] [IsDomain R] {l : ℕ} (hl : l.Prime) {ζ : R}
    (hζ : IsPrimitiveRoot ζ l) (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i)) :
    240 * (∑ i ∈ (range l).erase 0, tateDYterm (ζ ^ i)) = (l : R) ^ 4 - 1 := by
  sorry

/-! ## ★出典の紐付け(`.src`)と証明義務(`.needs`) -/

def tate_ode.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(Tate 曲線の ODE——DY = 3X² − Y + a₄)",
    sectionId := "genell-lemma-3-2" }

def tate_ode.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tate_equation((2Y+X)² = 4X³ + X² + 4a₄X + 4a₆、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tate_equation") 1,
    .citation "[ABC3]" "tateDXpair_eq(DX = 2Y + X、第 846、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateDXpair_eq") 1,
    .citation "[ABC3]" "univD(万有な環の上の微分、第 845、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.univD") 1,
    .implicitStep
      ("★★級数は TateUniv に無く完備化の中にあるので、q^N で切り詰めて " ++
       "TateUniv の中で整除性を見る(tate_equation を通した道・第 222–229 と同じ形)") 8,
    .implicitStep
      ("★★2DX で割る段——TateUniv は UFD ℤ[A,W] の局所化なので整域であり、" ++
       "2DX は q = AW と互いに素(mod q で DX ≡ A(1+A)/(1−A)³ ≠ 0)") 5 ]

def veluV2_eq_tateDYpair.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(veluV2 は DY そのもの。★ODE から直ちに)",
    sectionId := "genell-lemma-3-2" }

def veluV2_eq_tateDYpair.needs : List ProofObligation :=
  [ .citation "[ABC3]" "veluV2_tateCurveAt(Tate 曲線では v_Q = 3x² + a₄ − y、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.veluV2_tateCurveAt") 1,
    .citation "[ABC3]" "tate_ode(本ファイルの ODE)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.tate_ode") 1 ]

def sum_veluV2_eq_sum_tateDYpair.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(v = ∑_{i≠0} DY(ζ^i))",
    sectionId := "genell-lemma-3-2" }

def sum_veluV2_eq_sum_tateDYpair.needs : List ProofObligation :=
  [ .citation "[ABC3]" "veluV2_eq_tateDYpair(各項で veluV2 = DY)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.veluV2_eq_tateDYpair") 1,
    .implicitStep "☆ζ^i · (q·(ζ^i)^{l-1}) = q を ζ^l = 1 から出す段" 1 ]

def sum_mu_dyterm.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(定数項の指標和——240·∑ ζ²(2+ζ)/(1−ζ)⁴ = l⁴ − 1)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_dyterm.needs : List ProofObligation :=
  [ .citation "[ABC3]" "inv_one_sub_isRoot(1/(1−ζ) は y^l = (y−1)^l の根、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.inv_one_sub_isRoot") 1,
    .citation "[ABC3]" "sum_mu_frac(∑ ζ/(1−ζ)² = −(l²−1)/12、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.sum_mu_frac") 1,
    .citation "[ABC3]" "sum_mu_frac_cube(∑ ζ²/(1−ζ)³ = (l²−1)/24、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.sum_mu_frac_cube") 1,
    .implicitStep
      ("★p₄ = ∑ 1/(1−ζ)⁴ を Newton の公式で出す段——" ++
       "P(y) = y^l − (y−1)^l の Vieta は e_k = C(l, k+1)/l") 5 ]

end ABC3.Skeleton.GenEll
