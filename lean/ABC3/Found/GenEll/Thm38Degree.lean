/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Thm38KummerExists
import Mathlib.LinearAlgebra.Matrix.GeneralLinearGroup.Card
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★`d₀ = 23040` の正体（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.20。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★これは何か

原文 p.20 の第 1 段:

> there exists a Galois extension L′ of L of degree that divides d₀ = (3²−1)(3²−3)(5²−1)(5²−5)
> = 23040 [i.e., the order of GL₂(F₃)×GL₂(F₅)], so as to render the 3- and 5-torsion points
> of E_L rational over L′

★★★**本ファイルはその「次数が `23040` を割る」段を取る。**

`L′ ≔ L(E[3], E[5])` とすると `Gal(L′/L)` は
`GL₂(𝔽₃) × GL₂(𝔽₅)` に**単射**に入る（3-捩れ・5-捩れへの作用）。
したがって Lagrange から `[L′ : L] ∣ |GL₂(𝔽₃) × GL₂(𝔽₅)| = 23040` である。

★`|GL₂(𝔽_q)| = (q²−1)(q²−q)` は mathlib の `Matrix.card_GL_field` にある
（2026-08-29 に `#check` で確認）。★★`48 × 480 = 23040` は `decide` で出る
——**原文の「目視確認」した数値が機械で確かめられた**。

## ★★残るのは半安定還元の判定

☆「3-・5-捩れが有理的なら半安定還元」——Néron モデルと
Néron–Ogg–Shafarevich の判定（`E[m]` が有理的な `m ≥ 3` があれば半安定）。
★mathlib には `AlgebraicGeometry/EllipticCurve/Reduction.lean` があり
`IsMinimal`・`HasGoodReduction`・`HasMultiplicativeReduction`・三分律まで**ある**
（本日 `#check` で確認）が、**半安定還元の判定そのものは無い**。

## ★これで `Theorem 3.8` に残るのは

| 段 | 状態 |
|---|---|
| 群論（`Lemma 3.1, (iv)`・`§9-992`） | ★済み |
| Galois 表現 `galRep` | ★済み |
| `α` が像に入る | ★済み（`§9-993`・`§9-994`・`§9-995`） |
| `l`-巡回 ⟷ 安定直線 | ★標数 0 では同義（第 543） |
| `torsionExt` の**次数** `∣ 23040` | ★★**本ファイル** |
| ☆`torsionExt` の**半安定性** | ☆残る |
| ☆`Lemma 3.7` | ☆`Prop 3.4` 待ち |

★`.src` は条つき——指標には数えない。
-/

namespace ABC3.Found.GenEll

/-! ## ★★★★★★★★`|GL₂(𝔽₃) × GL₂(𝔽₅)| = 23040` -/

/-- ★★★★★★★★**原文の `d₀ = (3²−1)(3²−3)(5²−1)(5²−5) = 23040`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`|GL₂(𝔽₃)| = 48`、`|GL₂(𝔽₅)| = 480`、積が `23040`。
★★mathlib の `Matrix.card_GL_field`（`|GL_n(𝔽)| = ∏ (q^n − q^i)`）から `decide` で出る
——**原文が「目視確認」した数値が機械で確かめられた**。 -/
theorem gl2_three_five_card :
    Nat.card (GL (Fin 2) (ZMod 3)) * Nat.card (GL (Fin 2) (ZMod 5)) = 23040 := by
  haveI : Fact (Nat.Prime 3) := ⟨Nat.prime_three⟩
  haveI : Fact (Nat.Prime 5) := ⟨Nat.prime_five⟩
  rw [Matrix.card_GL_field, Matrix.card_GL_field]
  decide

/-! ## ★★★★★★★★★★★★★次数は `23040` を割る -/

/-- ★★★★★★★★★★★★★**`GL₂(𝔽₃) × GL₂(𝔽₅)` に単射で入る有限群の位数は `23040` を割る**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`L′ ≔ L(E[3], E[5])` として `Gal(L′/L) ↪ GL₂(𝔽₃) × GL₂(𝔽₅)`
（3-捩れ・5-捩れへの作用）を入れれば、これが原文の
『a Galois extension `L′` of `L` of degree that divides `d₀`』である。

★★機構は Lagrange だけ——`Nat.card G = Nat.card f.range ∣ Nat.card (G₁ × G₂) = 23040`。 -/
theorem card_dvd_23040_of_injective {G : Type} [Group G] [Finite G]
    (f : G →* (GL (Fin 2) (ZMod 3)) × (GL (Fin 2) (ZMod 5)))
    (hf : Function.Injective f) : Nat.card G ∣ 23040 := by
  have h1 : Nat.card G = Nat.card (f.range : Subgroup _) :=
    Nat.card_congr (Equiv.ofInjective f hf)
  have h2 : Nat.card (f.range : Subgroup _)
      ∣ Nat.card ((GL (Fin 2) (ZMod 3)) × (GL (Fin 2) (ZMod 5))) :=
    Subgroup.card_subgroup_dvd_card _
  rw [Nat.card_prod, gl2_three_five_card] at h2
  rw [h1]
  exact h2

/-! ## ★出典の紐付け(`.src`)——★**条つきである。指標には数えない** -/

def gl2_three_five_card.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(d₀ = (3²−1)(3²−3)(5²−1)(5²−5) = 23040)",
    sectionId := "genell-thm-3-8" }

def card_dvd_23040_of_injective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(L′/L の次数は 23040 を割る)",
    sectionId := "genell-thm-3-8" }

def card_dvd_23040_of_injective.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Matrix.card_GL_field(|GL_n(𝔽)| = ∏ (q^n − q^i))"
      (.inMathlib "Matrix.card_GL_field") 2,
    .citation "[mathlib]" "Subgroup.card_subgroup_dvd_card(Lagrange)"
      (.inMathlib "Subgroup.card_subgroup_dvd_card") 1,
    .folklore
      ("Gal(L(E[3],E[5])/L) ↪ GL₂(𝔽₃) × GL₂(𝔽₅)（3-捩れ・5-捩れへの作用）。" ++
       "★本プロジェクトには mod l 表現の構成がある（Found/GaloisRep/）——**接続は未了**") 4,
    .folklore
      ("半安定還元の判定: 『E[m] が有理的な m ≥ 3 があれば半安定』" ++
       "（Néron モデル／Néron–Ogg–Shafarevich）。" ++
       "★mathlib には EllipticCurve/Reduction.lean（IsMinimal・HasGoodReduction・" ++
       "HasMultiplicativeReduction・三分律）が**ある**が、判定そのものは無い" ++
       "（2026-08-29 に #check で確認）。**残る**") 8,
    .implicitStep
      ("★★★★★測定(2026-08-29): 原文 p.20 第 1 段の `d₀ = 23040` は" ++
       "`|GL₂(𝔽₃) × GL₂(𝔽₅)|` であり、mathlib の Matrix.card_GL_field から" ++
       "**decide で機械的に出る**——原文が『目視確認』した数値が確かめられた。" ++
       "★次数が 23040 を割ることは Lagrange だけである。" ++
       "★★残るのは**半安定還元の判定**だけである") 6 ]

end ABC3.Found.GenEll
