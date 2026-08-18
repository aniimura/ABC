/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop44Gl
import ABC3.Found.FrdI.Prop22Star

/-!
# [FrdI] Proposition 4.4, (iii) —— `Φ^birat` を `𝒟` 上の部分関手にする

★★★**2026-08-18 の測定で経路が変わったファイル**である。

原文 (FrdI p.85):
> Now assertion (iii) follows immediately from the existence of the functor

★(iii) は「(i) の関手の存在」＋ `Proposition 1.5, (ii)` から従う、と原典は言う。
★★そして (ii) の本文が `𝒟` 上の関手性の出所を明示している (FrdI p.83):

> functor "O×(−)" on D associated to the Frobenioid Cbirat [cf. Proposition 2.2, (ii),

★★★つまり **`𝒟` 上の関手性は `Proposition 2.2, (ii), (iii)` が供給する**。
「同じ底の対象は `𝒞^birat` で同型」を証明する筋ではない
(一度その筋を試して詰まった —— 記録は `frdi-decomposition.json` にある)。

## ★このファイルがすること

`Prop22Star.lean` の `otriStar`(= `Proposition 2.2, (ii)`)を
**`𝒞^birat` に適用する**。`𝒞^birat` が Frobenioid であることは
`birat_frobenioid_of_frobNormalized`(★(B) の追加仮定つき)から来る。
-/

universe v u v2 u2 w

namespace ABC3.Found.FrdI

open CategoryTheory Opposite

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (G : Frobenioid P)

/-! ## ★`𝒞^birat` を Frobenioid として取り出す -/

/-- ★★★**`𝒞^birat` の Frobenioid 構造**(★(B) の追加仮定つき)。

★仮定は「`𝒞^birat` の各対象が Frobenius-normalized」であり、
これは原文 `Definition 4.5, (i)` の **birationally Frobenius-normalized** そのもの。
★逸脱の分類は **(B)**(我々が仮定を追加した)。 -/
noncomputable def biratFrobenioid
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X) :
    Frobenioid (biratPre P G) :=
  birat_frobenioid_of_frobNormalized P G hfn

/-! ## ★★`Proposition 2.2, (ii)` を `𝒞^birat` に適用する -/

/-- ★★★★**`𝒪^▷(−)` は `𝒞^birat` の `𝒟*` 上の反変関手**。

★これが原文 `Proposition 4.4, (ii)` の
「the functor `𝒪^×(−)` **on 𝒟** associated to the Frobenioid `𝒞^birat`」の実体である。

★★**新しい証明は要らない** —— `Proposition 2.2, (ii)` の `otriStar` を
`𝒞^birat` に当てるだけ。`𝒟` 上の関手性はここから来る。 -/
noncomputable def otriStarBirat
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X) :
    (InducedCategory D (istrBase (biratPre P G)))ᵒᵖ ⥤ MonCat.{max (max u2 v) v2} :=
  otriStar (biratPre P G) (biratFrobenioid P G hfn).core (biratFrobenioid P G hfn)

/-! ## ★出典の紐付け(条つき) -/

/-- ★locator —— `Proposition 4.4, (iii)` の `𝒟` 上の関手性。 -/
def otriStarBirat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (ii) — 𝒟 上の関手 𝒪^×(−)",
    sectionId := "frdi-prop-4-4" }

end ABC3.Found.FrdI
