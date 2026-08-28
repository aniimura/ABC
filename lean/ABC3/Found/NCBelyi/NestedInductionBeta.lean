/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NCBelyi.MobiusConj
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★`β` を運ぶ入れ子帰納法（`Found`）

原典: S. Mochizuki, *Noncritical Belyi Maps* [NCBelyi]、物理 p.5。

原文 (NCBelyi p.5):
> hypothesis, and composing the resulting morphisms yields

## ★★★★★★★★★★★★★★★★★★★★★これは何か

`Lemma 2.4` の (b)（分離）を通すには、入れ子帰納法が **`β` も一緒に運ぶ**必要がある
——原文の

> In particular, replacing S by S′, β by f0(β), applying the induction hypothesis,
> and composing the resulting morphisms yields

の『**`β` を `f₀(β)` で置き換える**』がそれである。

★★★**本ファイルはその骨格を取る**（`nested_induction_descend_beta`）。
`Found/NCBelyi/DescendData.lean` の `nested_induction_descend'` は `Q S`／`P S` が
`β` に依らない形だったので、`β` に応じて降下先が変わる本件には使えなかった。

★機構は単純である——`P' S ≔ ∀ β, Q S β → P S β` と置けば
既存の `nested_induction`（測度 `(m(S), d(S))` の 2 重帰納）がそのまま回る。
★★帰納法の仮定が `∀ T` の形をしているので、`β` ごとに違う `T` を選んでよい。

## ★★★★★これで `Lemma 2.4` の (b) に残るのは 2 つ

| 部品 | 状態 |
|---|---|
| 数値の核（`f₀(β) ∉ S′`） | ★`SeparationStep.lean`（`§9-984`） |
| 測度が `ℚ`-Möbius で不変 | ★`MobiusRedDeg.lean`（`§9-983`） |
| `Gal`-安定性が `ℚ`-Möbius で不変 | ★`MobiusConj.lean`（`§9-985`） |
| 正規化（`|α| ≤ 1`、`|β| ≥ 4`） | ★`RatSeparation.lean`（第 409） |
| ★`β` を運ぶ帰納の骨格 | ★★**本ファイル** |
| ☆`Gal`-安定性が**多項式の像**で不変（`conjSet (p x) ⊆ p (conjSet x)`） | ☆**残る** |
| ☆合成を `ℙ¹` の有理写像として組み立てること | ☆**残る** |

★★1 つめの残りは `§9-985` の**多項式版**である
——`§9-985` は Möbius について同じことを変換多項式で示した。
多項式 `p` については `∏_{r ∈ conjSet x} (Y − p(r))` が `ℚ[Y]` に落ちること
（対称式）か、埋め込みの延長（`ℚ⟮p(x)⟯ ⊆ ℚ⟮x⟯` と `ℂ` が代数閉）が要る。
-/

namespace ABC3.Found.NCBelyi

open IntermediateField

/-! ## ★★★★★★★★★★★★★★★★★★★★★`β` を運ぶ帰納 -/

/-- ★★★★★★★★★★★★★★★★★★★★★**`β` を運ぶ入れ子帰納法**。

原文 (NCBelyi p.5):
> hypothesis, and composing the resulting morphisms yields

★原文の『replacing S by S′, **β by f₀(β)**, applying the induction hypothesis』
そのものである。

★★`P' S ≔ ∀ β, Q S β → P S β` と置けば既存の `nested_induction` がそのまま回る
——帰納法の仮定が `∀ T` の形だから、`β` ごとに違う降下先 `T` を選んでよい。 -/
theorem nested_induction_descend_beta {P Q : Finset ℂ → ℚ → Prop}
    (base : ∀ (S : Finset ℂ) (β : ℚ), maxRedDeg S = 0 → Q S β → P S β)
    (descend : ∀ (S : Finset ℂ) (β : ℚ), 0 < maxRedDeg S → Q S β →
      ∃ (T : Finset ℂ) (β' : ℚ),
        (maxRedDeg T < maxRedDeg S
          ∨ (maxRedDeg T = maxRedDeg S ∧ dSum T < dSum S))
        ∧ Q T β' ∧ (P T β' → P S β)) :
    ∀ (S : Finset ℂ) (β : ℚ), Q S β → P S β := by
  refine nested_induction (P := fun S => ∀ β : ℚ, Q S β → P S β) ?_ ?_
  · intro S h β hQ
    exact base S β h hQ
  · intro S hS ihm ihd β hQ
    obtain ⟨T, β', hmeas, hQ', himp⟩ := descend S β hS hQ
    refine himp ?_
    rcases hmeas with h | ⟨h1, h2⟩
    · exact ihm T h β' hQ'
    · exact ihd T h1 h2 β' hQ'

/-! ## ★★★★★★正規化の後も整である -/

/-- ★★★★★★**`ℚ`-Möbius の像も代数的**。

★正規化（`Lemma 2.3` の自己同型）を挟んだあとも `IsIntegral ℚ` が保たれる
——`σ x ∈ ℚ⟮x⟯` であり `ℚ⟮x⟯` は有限次だからである。 -/
theorem isIntegral_mobius (lam nu mu : ℚ) {x : ℂ} (hx : IsIntegral ℚ x) :
    IsIntegral ℚ ((nu : ℂ) * ((x - (lam : ℂ))⁻¹) + (mu : ℂ)) := by
  haveI : FiniteDimensional ℚ ℚ⟮x⟯ := IntermediateField.adjoin.finiteDimensional hx
  have hratmem : ∀ (K : IntermediateField ℚ ℂ) (q : ℚ), ((q : ℂ)) ∈ K := by
    intro K q
    have h := K.algebraMap_mem q
    simpa using h
  have hxm : x ∈ ℚ⟮x⟯ := mem_adjoin_simple_self ℚ x
  have hmem : (nu : ℂ) * ((x - (lam : ℂ))⁻¹) + (mu : ℂ) ∈ ℚ⟮x⟯ :=
    add_mem (mul_mem (hratmem _ nu) (inv_mem (sub_mem hxm (hratmem _ lam)))) (hratmem _ mu)
  have h1 : IsIntegral ℚ (⟨_, hmem⟩ : ℚ⟮x⟯) := IsIntegral.of_finite ℚ _
  simpa using h1.map (IsScalarTower.toAlgHom ℚ ℚ⟮x⟯ ℂ)

/-! ## ★出典の紐付け(`.src`)——★**条つきである。指標には数えない** -/

def nested_induction_descend_beta.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(β を運ぶ入れ子帰納法——replacing β by f₀(β))",
    sectionId := "ncbelyi-lemma-2-4" }

def isIntegral_mobius.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(ℚ-Möbius の像も代数的)",
    sectionId := "ncbelyi-lemma-2-4" }

def nested_induction_descend_beta.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "nested_induction(m(S)・d(S) の 2 重帰納、第 408)"
      (.inProject "ABC3" "ABC3.Found.NCBelyi.nested_induction") 3,
    .implicitStep
      ("★★★★★測定(2026-08-29): Lemma 2.4 の (b)(分離)を通すには帰納法が" ++
       "**β も一緒に運ぶ**必要がある(原文の『replacing S by S′, β by f₀(β)』)。" ++
       "★DescendData.lean の nested_induction_descend′ は Q S／P S が β に依らない形なので" ++
       "本件には使えなかった。★★P′ S ≔ ∀ β, Q S β → P S β と置けば既存の " ++
       "nested_induction がそのまま回る——帰納法の仮定が ∀ T の形だから、" ++
       "β ごとに違う降下先を選んでよい") 5,
    .implicitStep
      ("★★★これで Lemma 2.4 の (b) に残るのは**2 つ**である: " ++
       "(1) Gal-安定性が**多項式の像**で不変であること" ++
       "(conjSet (aeval x p) ⊆ (conjSet x).image (aeval · p))" ++
       "——§9-985 の多項式版。∏_{r ∈ conjSet x}(Y − p(r)) が ℚ[Y] に落ちること(対称式)か、" ++
       "埋め込みの延長(ℚ⟮p(x)⟯ ⊆ ℚ⟮x⟯ と ℂ が代数閉)が要る。" ++
       "(2) 合成を ℙ¹ の有理写像として組み立てること" ++
       "(原文の f(x) ∈ ℚ(x) は Möbius と多項式の交互合成)") 7 ]

end ABC3.Found.NCBelyi
