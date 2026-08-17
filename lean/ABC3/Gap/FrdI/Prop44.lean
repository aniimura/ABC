import ABC3.Found.FrdI.Prop44Core

/-!
# Gap — [FrdI] `Proposition 4.4, (ii)` の最後の 1 条

★**`Found/FrdI/Prop44Core.lean` で `𝒞^birat` の `FrobenioidCore` 21 条のうち
20 条と、`Definition 1.3, (i)(c)` の圏同値(`birat_plBkEquiv`)を証明した。**
★残るのは **`otriBase`(`Definition 1.3, (iii), (c)` の「全単射は `Base(φ)` にしか依らない」)
の 1 条だけ**である。

## ★原文

原文 (FrdI p.85):
> exercise to check that Cbirat is, in fact, a Frobenioid of group-like type. Moreover, it

★原文はここを **「routine exercise」** と書いて証明を置かない。

## ★★測定 —— 残る 1 条は何と同値か(2026-08-17)

`otriBase` は次の形である:

> `φ`・`φ' : A ⟶ B` が co-angular pre-step で `Base φ = Base φ'` のとき、
> `α ∈ 𝒪^▷(A)`・`β ∈ 𝒪^▷(B)` について
> `φ ≫ β = α ≫ φ` ⟹ `φ' ≫ β = α ≫ φ'`。

★★`𝒞^birat` では co-angular pre-step は**同型**なので、`p := φ⁻¹` と置くと
`β = p ≫ α ≫ φ` で、結論は `u ≫ α = α ≫ u`(`u := φ' ≫ p ∈ 𝒪^▷(A)`)に**同値**になる。
★★★すなわち **`otriBase` ⟺ `𝒪^▷(A^birat)`(= `𝒪^×(A^birat)`)が可換**。

★`Remark 1.3.1` は「`𝒪^▷(A)` は可換」を **(iii)(b),(c) から**導くが、
その (c) が `otriBase` 自身なので、**この経路は循環している**。

## ★★測定 —— どこまでは示せたか

★以下は `Prop44Core.lean` で**実際に証明した**もの:

| 事実 | 宣言 |
|---|---|
| `𝒞^birat` の自己射はすべて co-angular | `birat_isCoAngular_endo` |
| よって `𝒪^▷(A^birat)` の元はすべて**同型**(= 群になる) | (上と `birat_isIso_of_coaPre_birat`) |
| `𝒪^▷(A^birat)` の元は `[c]⁻¹ ≫ [s]`(`c`・`s : E ⟶ A` は `𝒞` の co-angular pre-step で `Base c = Base s`)の形 | `birat_hom_repr` |

★★さらに紙の上では次まで出ている(未実装):
**`𝒞` の `otriBase` を `c` と `s` に当てると `Θ_c = Θ_s` なので、
`u = [c]⁻¹ ≫ [s]` は `𝒪^▷(A)_𝒞` の像すべてと可換**である。
★したがって **`𝒪^▷(A)_𝒞` の像は `𝒪^▷(A^birat)` の中心に入る**。

★★★**足りないのは「`𝒪^▷(A)_𝒞` の像とその逆元が `𝒪^▷(A^birat)` を生成する」こと**
—— これが言えれば可換性は上の中心性から**直ちに**従う。
★モデル(数体上の直線束の Frobenioid)では `𝒪^▷(A)^gp = K^×` で成り立つが、
`Definition 1.3` の公理からこれを出す道は見つかっていない
(`preStepSpan` は `X ⟶ A`・`X ⟶ E` の span を与えるだけで、
必要な `A ⟶ E` 向きの co-angular pre-step は与えない)。

## ★影響範囲

★`Proposition 4.4, (ii)` が閉じないので、
`Prop 4.4` 全体の条なし `.src` は**まだ書けない**。
★`(i)`・`(iii)` の `Φ^birat`・`(iv)` の辞書は本条に依存しないので、
それらは独立に進められる。
-/

namespace ABC3.Gap.FrdI

/-- ★★**`Proposition 4.4, (ii)` に足りないもの**。

★分類は **② `missingMath`** —— 主張はモデルで真であり、
原典も「routine exercise」と書く。**足りているのは我々の構成**であって、
原典の飛躍だと名乗るには根拠が弱い。 -/
structure Gap_4_4_ii_otriBase : Prop where
  /-- 不足: `𝒪^▷(A)_𝒞` の像とその逆元が `𝒪^×(A^birat)` を生成すること
  (これがあれば中心性から可換性が出る)。 -/
  otriGpGenerates : True

def Gap_4_4_ii_otriBase.record : ABC3.Meta.GapRecord :=
  { source :=
      { paper := "FrdI", pdfPage := 85, item := "Proposition 4.4, (ii)",
        sectionId := "frdi-prop-4-4" },
    classification := ABC3.Meta.GapClass.missingMath,
    falsifier :=
      "★**これが埋まる条件**: `𝒪^▷(A)^gp -> 𝒪^×(A^birat)` が全射であることを " ++
      "`Definition 1.3` の公理から示すこと。★等価な形は " ++
      "「co-angular pre-step `c : E ⟶ A` が与えられたとき、" ++
      "`m ≫ c` と `m ≫ s` がともに `𝒪^▷(A)` に入るような `m : A ⟶ E` を取れる」。" ++
      "★あるいは `𝒪^×(A^birat)` の可換性を直接示す別経路でもよい。" ++
      "★★**③ に上がる条件**: 上のどちらも公理から出ないことを、" ++
      "`Definition 1.3` を満たすが `𝒪^×(A^birat)` が非可換になる Frobenioid の" ++
      "**反例**を作って示すこと。★このとき原文の「routine exercise」は誤りになる。" ++
      "★我々の測定(2026-08-17)では `𝒪^▷(A)_𝒞` の像が中心に入ることまでは出ており、" ++
      "反例があるなら中心拡大の交換子で見えるはずである。" }

end ABC3.Gap.FrdI
