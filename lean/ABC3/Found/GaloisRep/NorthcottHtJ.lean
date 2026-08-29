/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.HtJWeil
import ABC3.Found.GenEll.NorthcottImage
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★`ht^Falt` が有界な曲線の `j` は有限個（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> on Mell(Q). In particular, if C ∈R, then the set of points [E] ∈Mell(Q)≤d such

## ★★★★★★★★★★★★★★★★★★★★これは何か —— `Prop 3.4` の「In particular」

`Proposition 3.4` の後半は

    「`C ∈ ℝ` なら `{[E] ∈ M_ell(ℚ̄)^{≤d} : ht^Falt([E]) ≤ C}` は**有限**」

である。★★`M_ell(ℚ̄)` の点は `j` 不変量で決まるので、これは

    **「`ht^Falt` が有界で定義体の次数が `d` 以下の楕円曲線の `j` の集合は有限」**

と同じことである。★★★★**本ファイルはそれを証明する**——★**無条件**。

## ★段取り —— 2 つの入力がそろった

`Found/GenEll/FiniteFromNorthcott.lean` の `finite_of_le_of_northcott` は
2 つの入力から有限性を出す:

| 入力 | どこ |
|---|---|
| `ht∞ ≤ 12(1+ϵ)·ht^Falt + C` | ★`Found/GaloisRep/HtJBound.lean`（`§9-1003`、無条件） |
| `ht∞` の Northcott 性 | ★`Found/GenEll/NorthcottImage.lean`（`§9-950` の像版） |

★★本ファイルは `ht∞ ≔ htJ`（`j` の Weil 高さ）で**両方をつなぐ**。
`§9-1004` の `htJ_eq_logHeight`（`htJ = Height.logHeight₁(j)/[L:ℚ]`）が
その接着剤である。

## ★機構

`j = a/b`（`a, b ∈ 𝓞_L`、`IsFractionRing.div_surjective`）と書くと

    `htJ = log(mulHeight₁ (a/b))/[L:ℚ] = log(mulHeight ![a,b])/[L:ℚ]`

——★これは `northcott_of_log_mulHeight_image` が受ける形そのものである。
★★出てくる有限性は「座標の組 `![j, 1]` の像が有限」であり、
第 0 成分を見れば `j` の像の有限性になる。

## ☆残るもの

☆`M_ell(ℚ̄)` の点と `j` の 1 対 1 対応（`EllClass` の構成）。
★これは `Interface/GenEll/EllModuli.lean` の witness を作る段であり、
本ファイルはその `northcott` 欄の**中身**を用意したことになる。
-/

namespace ABC3.Found.GaloisRep

open NumberField WeierstrassCurve ABC3.Found.GenEll

/-! ## ★★★★★★★★★★★★★★`htJ` の Northcott 性 -/

/-- ★★★★★★★★★★★★★★**`htJ` が有界な曲線の `j` は有限個**（次数を止めて）。

原文 (GenEll p.17):
> on Mell(Q). In particular, if C ∈R, then the set of points [E] ∈Mell(Q)≤d such

★`j = a/b`（`a, b ∈ 𝓞_L`）と書き、`§9-1004` の `htJ_eq_logHeight` で
`Height.mulHeight ![a,b]` の言葉に直して `northcott_of_log_mulHeight_image` に渡す。 -/
theorem northcott_htJ_image (d : ℕ) {P : Type}
    (fld : P → IntermediateField ℚ ℂ) (hnf : ∀ p, NumberField (fld p))
    (hdeg : ∀ p, Module.finrank ℚ (fld p) ≤ d)
    (E : ∀ p, WeierstrassCurve (fld p)) (hE : ∀ p, (E p).IsElliptic)
    (ht : P → ℝ)
    (hht : ∀ p, haveI := hnf p; haveI := hE p; ht p = htJ (fld p) (E p))
    (C : ℝ) :
    ((fun p : P => haveI := hE p; (((E p).j : fld p) : ℂ)) '' {p | ht p ≤ C}).Finite := by
  classical
  have hchoice : ∀ p, ∃ ab : (𝓞 (fld p)) × (𝓞 (fld p)),
      haveI := hnf p; haveI := hE p;
      ((ab.1 : fld p) / (ab.2 : fld p) = (E p).j) := by
    intro p
    haveI := hnf p; haveI := hE p
    obtain ⟨a, b, _, hab⟩ := IsFractionRing.div_surjective (𝓞 (fld p)) ((E p).j)
    exact ⟨(a, b), by rw [← hab]⟩
  choose ab hab using hchoice
  set x : ∀ p, Fin 2 → 𝓞 (fld p) := fun p => ![(ab p).1, (ab p).2] with hx
  have hht' : ∀ p, haveI := hnf p; ht p
      = Real.log (Height.mulHeight (fun k => ((x p k : fld p))))
        / (Module.finrank ℚ (fld p) : ℝ) := by
    intro p
    haveI := hnf p; haveI := hE p
    rw [hht p, htJ_eq_logHeight (E p), Height.logHeight₁_eq_log_mulHeight₁]
    congr 2
    rw [← hab p, Height.mulHeight₁_div_eq_mulHeight]
    congr 1
    funext k
    fin_cases k <;> simp [hx]
  have hfin := northcott_of_log_mulHeight_image d fld hnf hdeg x ht hht' 1 C
  have himg : (fun p : P => haveI := hE p; (((E p).j : fld p) : ℂ))
      = (fun g : Fin 2 → ℂ => g 0) ∘
        (fun (p : P) (k : Fin 2) => ((((x p k : fld p)) / ((x p 1 : fld p)) : fld p) : ℂ)) := by
    funext p
    haveI := hnf p; haveI := hE p
    simp only [Function.comp_apply, hx]
    rw [← hab p]
    simp
  rw [himg, Set.image_comp]
  exact hfin.image _

/-! ## ★★★★★★★★★★★★★★★★★★★★`Proposition 3.4` の「In particular」 -/

/-- ★★★★★★★★★★★★★★★★★★★★**`ht^Falt` が有界な曲線の `j` は有限個**。

原文 (GenEll p.17):
> on Mell(Q). In particular, if C ∈R, then the set of points [E] ∈Mell(Q)≤d such

★★★**無条件**である——2 つの入力
（`§9-1003` の `htJ ≤ 12(1+ϵ)·ht^Falt + C` と本ファイルの Northcott 性）が
どちらも証明されているため。

★`ϵ = 1`（`12(1+ϵ) = 24`）で使っている——`ϵ` は何でもよい。 -/
theorem finite_j_of_htFalt_le (d : ℕ) {P : Type}
    (fld : P → IntermediateField ℚ ℂ) (hnf : ∀ p, NumberField (fld p))
    (hdeg : ∀ p, Module.finrank ℚ (fld p) ≤ d)
    (E : ∀ p, WeierstrassCurve (fld p)) (hE : ∀ p, (E p).IsElliptic)
    (htF : P → ℝ)
    (hhtF : ∀ p, haveI := hnf p; htF p = htFaltOf (fld p) (E p))
    (C : ℝ) :
    ((fun p : P => haveI := hE p; (((E p).j : fld p) : ℂ)) '' {p | htF p ≤ C}).Finite := by
  obtain ⟨C₀, hC₀⟩ := exists_htJ_le_htFalt' 1 one_pos
  refine Set.Finite.subset
    (northcott_htJ_image d fld hnf hdeg E hE
      (fun p => haveI := hnf p; haveI := hE p; htJ (fld p) (E p)) (fun p => rfl)
      (12 * (1 + 1) * C + C₀)) ?_
  refine Set.image_mono (fun p hp => ?_)
  haveI := hnf p; haveI := hE p
  have h1 : htFaltOf (fld p) (E p) ≤ C := by rw [← hhtF p]; exact hp
  have h2 := hC₀ (fld p) (E p)
  show htJ (fld p) (E p) ≤ 12 * (1 + 1) * C + C₀
  nlinarith [h1, h2]

/-! ## ★出典の紐付け(`.src`) -/

def northcott_htJ_image.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(htJ の Northcott 性——j の像の有限性)",
    sectionId := "genell-prop-3-4" }

def finite_j_of_htFalt_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(In particular——ht^Falt が有界な曲線の j は有限個。★無条件)",
    sectionId := "genell-prop-3-4" }

def finite_j_of_htFalt_le.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]"
      "IsFractionRing.div_surjective(数体の元は整数の商で書ける)"
      (.inMathlib "IsFractionRing.div_surjective") 1,
    .citation "[mathlib]"
      "Height.mulHeight₁_div_eq_mulHeight(商の高さ = 組の高さ)"
      (.inMathlib "Height.mulHeight₁_div_eq_mulHeight") 1,
    .implicitStep
      ("★★★★★★★★★★到達点(2026-08-29、第 556): " ++
       "Prop 3.4 の**「In particular」が無条件で取れた**。" ++
       "★2 つの入力——§9-1003 の htJ ≤ 12(1+ϵ)ht^Falt + C と " ++
       "本ファイルの Northcott 性——がどちらも証明されているため。" ++
       "★★接着剤は §9-1004 の htJ_eq_logHeight(htJ = Height.logHeight₁(j)/[L:ℚ])") 9,
    .implicitStep
      ("☆残るのは M_ell(ℚ̄) の点と j の 1 対 1 対応(EllClass の構成)。" ++
       "★本ファイルは Interface/GenEll/EllModuli.lean の northcott 欄の" ++
       "**中身**を用意したことになる") 7 ]

end ABC3.Found.GaloisRep
