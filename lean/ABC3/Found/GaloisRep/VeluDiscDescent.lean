/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.VeluKernelNorm
import ABC3.Meta.Claim

/-!
# 第 1390 ブロック —— **判別式の恒等式の降下**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——節点を「`ℂ` 側」と「降下」に割る

第 1387 の節点

    Δ(E)^l = Δ(E/C) · ( ∏_{P ∈ C∖{O}} (2 y_P + a₁ x_P + a₃) )^4

は**体の元どうしの等式**なので、単射な環準同型 `f : L →+* A` で
`A` の側から `L` の側へ**そのまま降りる**。

★本ブロックはその降下を作る。☆残るのは `ℂ` の側だけになる
——`ℂ` では一意化（第 1330-1335）と `latticeCurve_eq_veluQuotientFull`（在庫）がある。

☆道具は 2 つとも在庫である:

* `veluQuotientFull_map`——商は底変換と可換
* `WeierstrassCurve.map_Δ`——`Δ` は底変換と可換

★足りないのは `veluKernelNorm` が底変換と可換であること（本ブロックの前半）。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Finset ABC3.Meta ABC3.Found.GenEll

open scoped Classical

/-- ★★★★★★★★**核のノルムは底変換と可換**——★**無条件**（第 1390）。 -/
theorem veluKernelNorm_map {F A : Type} [Field F] [Field A] (f : F →+* A)
    (W : WeierstrassCurve F) (S : Finset (F × F)) :
    f (veluKernelNorm W S)
      = veluKernelNorm (W.map f) (S.image (fun Q => (f Q.1, f Q.2))) := by
  have hinj : Function.Injective (fun Q : F × F => (f Q.1, f Q.2)) := by
    intro a b hab
    have h1 : f a.1 = f b.1 := congrArg Prod.fst hab
    have h2 : f a.2 = f b.2 := congrArg Prod.snd hab
    exact Prod.ext (f.injective h1) (f.injective h2)
  rw [veluKernelNorm, veluKernelNorm, map_prod,
    Finset.prod_image (fun a _ b _ h => hinj h)]
  refine Finset.prod_congr rfl (fun q _ => ?_)
  show f (2 * q.2 + W.a₁ * q.1 + W.a₃)
    = 2 * f q.2 + (W.map f).a₁ * f q.1 + (W.map f).a₃
  rw [WeierstrassCurve.map_a₁, WeierstrassCurve.map_a₃]
  simp only [map_add, map_mul, map_ofNat]

/-- ★★★★★★★★★★★★★★★★
**判別式の恒等式は単射な環準同型で降りる**——★**無条件**（第 1390）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★これで第 1387 の節点は **`ℂ` の側だけ**になる。 -/
theorem disc_pow_eq_of_embed {L A : Type} [Field L] [Field A] (f : L →+* A)
    (E : WeierstrassCurve L) (S : Finset (L × L)) {l : ℕ}
    (h : (E.map f).Δ ^ l
      = (veluQuotientFull (E.map f) (S.image (fun Q => (f Q.1, f Q.2)))).Δ
        * (veluKernelNorm (E.map f) (S.image (fun Q => (f Q.1, f Q.2)))) ^ 4) :
    E.Δ ^ l = (veluQuotientFull E S).Δ * (veluKernelNorm E S) ^ 4 := by
  refine f.injective ?_
  rw [map_pow, map_mul, map_pow, ← WeierstrassCurve.map_Δ, ← WeierstrassCurve.map_Δ,
    veluQuotientFull_map, veluKernelNorm_map]
  exact h

/-! ## ★出典の紐付け(`.src`) -/

def veluKernelNorm_map.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(核のノルムは底変換と可換。★無条件)",
    sectionId := "genell-lemma-3-5" }

def disc_pow_eq_of_embed.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(判別式の恒等式は単射な環準同型で降りる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def disc_pow_eq_of_embed.needs : List ProofObligation :=
  [ .citation "[ABC3]" "veluQuotientFull_map(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.veluQuotientFull_map") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1390）**——第 1387 の節点を「`ℂ` 側」と「降下」に割る段である。" ++
       "☆恒等式は**体の元どうしの等式**なので、単射な環準同型でそのまま降りる" ++
       "（第 1334 `isElliptic_velu_congr_curve` と同じ型）。" ++
       "★★★これで残るのは `ℂ` の側だけ——" ++
       "`ℂ` では一意化（第 1330-1335）と `latticeCurve_eq_veluQuotientFull` がある。") 17 ]

end ABC3.Found.GaloisRep
